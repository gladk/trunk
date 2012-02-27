/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "BeamRecorder.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
// to calculate strain of whole speciemen - first two bodies in subscribedBodies are Nodes. FIXME - make it clean!
#include <boost/lexical_cast.hpp>

/*
 *
	private :
		std::ofstream ofile; 
		std::vector<std::pair<std::list<unsigned int>, int> >	subscribedBodies; // [region]<list of bodies' ids, direction>
		bool			first;
		int			nodeGroupMask,beamGroupMask;

	public :
		std::string	 outputFile;
		unsigned int	 interval;
		std::list<std::pair<Vector3r,std::pair< Real, int > > > sections; // a list of section planes: <a midpoint, half length, direction >
		// int direction refers to: which coordinate to take from the Beam's.direction: x, y or z (0, 1 or 2)
 * */


BeamRecorder::BeamRecorder () : DataRecorder()
{
	outputFile = "";
	interval = 50;
	subscribedBodies.clear();
	first = true;
	// FIXME ...
        nodeGroupMask           = 1;
        beamGroupMask           = 2;
}


void BeamRecorder::preProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
	}
	else
	{
		sections_midpoints.clear();
		sections_halflength.clear();
		sections_directions.clear();
		
		std::list<std::pair<Vector3r,std::pair< Real, int > > >::iterator i=sections.begin();
		std::list<std::pair<Vector3r,std::pair< Real, int > > >::iterator e=sections.end();
		for( ; i!=e ; ++i)
		{
			sections_midpoints.push_back(i->first);
			sections_halflength.push_back(i->second.first);
			sections_directions.push_back(i->second.second);
		}
	}
}

void BeamRecorder::postProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
		sections.clear();
		std::vector<Vector3r>::iterator i1=sections_midpoints.begin();
		std::vector<Vector3r>::iterator e1=sections_midpoints.end();
		std::vector<Real>::iterator     i2=sections_halflength.begin();
		std::vector<int>::iterator      i3=sections_directions.begin();
		//std::vector<Vector3r>::iterator ea region_max.end();
		for( ; i1!=e1 ; ++i1 , ++i2 , ++i3 )
		{
			sections.push_back(std::make_pair(*i1,std::make_pair(*i2,*i3)));
		}

	//	outputFile += "_"+boost::lexical_cast<std::string>(Omega::instance().getTimeStep());
	//	std::cerr << "using dt for BeamRecorder output file: " << outputFile << "\n";

		ofile.open(outputFile.c_str());

		first = true;
		subscribedBodies.clear();
	}
}




bool BeamRecorder::isActivated(Scene*)
{
	return ((Omega::instance().getCurrentIteration() % interval == 0) && (ofile));
}


void BeamRecorder::action()
{
//	std::cerr << "BeamRecorder quits\n";
//	return;
	

	if(first)
	{
	// scan the sections and subscribe beams that cross each of them
	// FIXME - this is still problematic, because this should be done only
	// on the first run, regardless of saving/loading this file...
		first = false;
		subscribedBodies.resize( sections.size() );

		std::list<std::pair<Vector3r,std::pair<Real,int> > >::iterator i  =sections.begin();
		std::list<std::pair<Vector3r,std::pair<Real,int> > >::iterator end=sections.end();

		for(size_t section=0 ; i != end ; ++i, ++section )
		{
			Vector3r midpoint = i->first;
			Real     half     = i->second.first;
			int	 dir	  = i->second.second;

			if(!subscribedBodies[section].first.empty())
				std::cerr << "First, but not empty\n";

			subscribedBodies[section].second=dir; 

			BodyContainer::iterator bi    = ncb->bodies->begin();
			BodyContainer::iterator biEnd = ncb->bodies->end();
			for(  ; bi!=biEnd ; ++bi )
			{
				if( (*bi)->getGroupMask() & beamGroupMask )
				{
					Vector3r pos1 = (*(ncb->bodies))[static_cast<LatticeBeamParameters*>( (*bi)->physicalParameters.get() )->id1]->physicalParameters->se3.position; // beam first node
					Vector3r pos2 = (*(ncb->bodies))[static_cast<LatticeBeamParameters*>( (*bi)->physicalParameters.get() )->id2]->physicalParameters->se3.position; // beam second node
					
					if( pos1[dir] < pos2[dir] )
						std::swap(pos1,pos2); // make sure that pos1[dir] is bigger 
		
					if(        pos1[dir] > midpoint[dir]
						&& pos2[dir] < midpoint[dir] // beam crosses the plane!

						&& pos1[(dir+1)%3] > midpoint[(dir+1)%3]-half
						&& pos1[(dir+1)%3] < midpoint[(dir+1)%3]+half
						&& pos2[(dir+1)%3] > midpoint[(dir+1)%3]-half
						&& pos2[(dir+1)%3] < midpoint[(dir+1)%3]+half // both points are within halflength in direction (dir+1)%3

						&& pos1[(dir+2)%3] > midpoint[(dir+2)%3]-half
						&& pos1[(dir+2)%3] < midpoint[(dir+2)%3]+half
						&& pos2[(dir+2)%3] > midpoint[(dir+2)%3]-half
						&& pos2[(dir+2)%3] < midpoint[(dir+2)%3]+half // both points are within halflength in direction (dir+2)%3
						)
						subscribedBodies[section].first.push_back( (*bi)->getId() ); 
				}
			}
		}
	}


	for(size_t section=0 ; section<sections.size() ; ++section )
	{
		std::list<unsigned int>::iterator i   = subscribedBodies[section].first.begin();
		std::list<unsigned int>::iterator end = subscribedBodies[section].first.end();
		//int dir                               = subscribedBodies[section].second; // unused...
		
		//Real stress=0; // in direction 'dir'
		Real stress_x=0;
		Real stress_y=0;

		for( ; i != end ; ++i )
		{
			if( (*(ncb->bodies)).exists(*i) )
			{
				LatticeBeamParameters* beam 	= static_cast<LatticeBeamParameters*>( (*(ncb->bodies))[*i]->physicalParameters . get() );
				// direction is a unit vector, so direction[1] is cos(angle_to_dir)
				//Real s			= beam->strain() * std::abs(beam->direction[dir]) * beam->longitudalStiffness;
				Real s_x			= beam->strain() * std::abs(beam->direction[0]) * beam->longitudalStiffness;
				Real s_y			= beam->strain() * std::abs(beam->direction[1]) * beam->longitudalStiffness;
				stress_x += s_x;
				stress_y += s_y;

				(*(ncb->bodies))[*i]->shape->diffuseColor = Vector3r(1.0,0.0,((float)section)/3.0); // FIXME [1]
			}
		}

		ofile	<< lexical_cast<string>(stress_x) << " " << lexical_cast<string>(stress_y) << " \t ";
	}
	ofile << endl; 
		
	// [1]
	// mark colors what is monitored - FIXME - provide a way to do that for ALL DataRecorders in that parent class.
	// GLDrawSomething can just put a getClassName()
	 
}

YADE_PLUGIN((BeamRecorder));

YADE_REQUIRE_FEATURE(PHYSPAR);


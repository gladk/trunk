/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "NodeRecorder.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
// to calculate strain of whole speciemen - first two bodies in subscribedBodies are Nodes. FIXME - make it clean!
#include <boost/lexical_cast.hpp>


NodeRecorder::NodeRecorder () : DataRecorder()
{
	outputFile = "";
	interval = 50;
	subscribedBodies.clear();
	first = true;
	// FIXME ...
        nodeGroupMask           = 1;
        beamGroupMask           = 2;
	only_this_stiffness	= -1;
}


void NodeRecorder::preProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
	}
	else
	{
		regions_min.clear();
		regions_max.clear();
		std::list<std::pair<Vector3r,Vector3r> >::iterator i=regions.begin();
		std::list<std::pair<Vector3r,Vector3r> >::iterator e=regions.end();
		for( ; i!=e ; ++i)
		{
			regions_min.push_back(i->first);
			regions_max.push_back(i->second);
		}
	}
}

void NodeRecorder::postProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
		regions.clear();
		std::vector<Vector3r>::iterator ii=regions_min.begin();
		std::vector<Vector3r>::iterator ei=regions_min.end();
		std::vector<Vector3r>::iterator ia=regions_max.begin();
		//std::vector<Vector3r>::iterator ea region_max.end();
		for( ; ii!=ei ; ++ii , ++ia )
		{
			regions.push_back(std::make_pair(*ii,*ia));
		}

	//	outputFile += "_"+boost::lexical_cast<std::string>(Omega::instance().getTimeStep());
	//	std::cerr << "using dt for NodeRecorder output file: " << outputFile << "\n";

		ofile.open(outputFile.c_str());

		first = true;
		subscribedBodies.clear();
	}
}




bool NodeRecorder::isActivated(Scene*)
{
	return ((Omega::instance().getCurrentIteration() % interval == 0) && (ofile));
}


void NodeRecorder::action()
{
//	std::cerr << "NodeRecorder quits\n";
//	return;
	

	if(first)
	{
	// scan the regions and subscribe bodies for each of them
	// FIXME - this is still problematic, because this should be done only
	// on the first run, regardless of saving/loading this file...
		first = false;
		subscribedBodies.resize( regions.size() );

		std::list<std::pair<Vector3r,Vector3r> >::iterator i  =regions.begin();
		std::list<std::pair<Vector3r,Vector3r> >::iterator end=regions.end();

		for( int region=0 ; i != end ; ++i, ++region )
		{
			Vector3r min = i->first;
			Vector3r max = i->second;

			if(!subscribedBodies[region].empty())
				std::cerr << "First, but not empty\n";

			BodyContainer::iterator bi    = ncb->bodies->begin();
			BodyContainer::iterator biEnd = ncb->bodies->end();
			for(  ; bi!=biEnd ; ++bi )
			{
				if( (*bi)->getGroupMask() & nodeGroupMask )
				{
					Vector3r pos = (*bi)->physicalParameters->se3.position;
					if(        pos[0] > min[0] 
						&& pos[1] > min[1] 
						&& pos[2] > min[2] 
						&& pos[0] < max[0] 
						&& pos[1] < max[1] 
						&& pos[2] < max[2] 
						&& ((*bi)->getGroupMask() & nodeGroupMask)
					  )
					{
						if(only_this_stiffness > 0)
						{
							bool ok;ok=true;
							// check this node's beams, quite time consuming
							BodyContainer::iterator bi2    = ncb->bodies->begin();
							BodyContainer::iterator bi2End = ncb->bodies->end();
							for(  ; bi2!=bi2End ; ++bi2 )
							{
								if( (*bi2)->getGroupMask() & beamGroupMask )
								{
									LatticeBeamParameters* beam = dynamic_cast<LatticeBeamParameters*>((*bi2)->physicalParameters.get());
									if(beam)
									{
										if(beam->id1 == (*bi)->getId() || beam->id2 == (*bi)->getId())
											if(beam->longitudalStiffness != only_this_stiffness)
											{
												ok=false;
												break;
											}
									}
									else
									{
										std::cerr << "ERROR: that's not a beam!\n";
										exit(1);
									}
								}
							}
							if(ok)
								subscribedBodies[region].push_back((*bi)->getId());
						}
						else
						{
							subscribedBodies[region].push_back((*bi)->getId());
						}
					}
				}
			}
		}
	}

	for(size_t region=0 ; region<regions.size() ; ++region )
	{
		std::list<unsigned int>::iterator i   = subscribedBodies[region].begin();
		std::list<unsigned int>::iterator end = subscribedBodies[region].end();
		
		Real count=0;
		Vector3r sum=Vector3r(0,0,0);

		for( ; i != end ; ++i )
		{
			if( (*(ncb->bodies)).exists(*i) )
			{
				sum+=(*(ncb->bodies))[*i]->physicalParameters->se3.position;
				++count;
				
				//(*(ncb->bodies))[*i]->geometricalModel->diffuseColor = Vector3r(((float)((region+3)%5))/5.0,1.0,((float)region)/5.0); // FIXME [1]
				Vector3r col(0,0,0);
				switch(region%6)
				{       //                    0 0 0
					case 0 : col=Vector3r(0,0,1);break;
					case 1 : col=Vector3r(0,1,0);break;
					case 2 : col=Vector3r(0,1,1);break;
					case 3 : col=Vector3r(1,0,0);break;
					case 4 : col=Vector3r(1,0,1);break;
					case 5 : col=Vector3r(1,1,0);break;
					//                    1 1 1        
				}
				(*(ncb->bodies))[*i]->shape->diffuseColor = col; // FIXME [1]
			}
		}

		sum/=count;

		ofile	<< lexical_cast<string>(sum[0]) << " " 
			<< lexical_cast<string>(sum[1]) << " " 
			<< lexical_cast<string>(sum[2]) << " \t "; 
	}
	ofile << endl; 
		
	// [1]
	// mark colors what is monitored - FIXME - provide a way to do that for ALL DataRecorders in that parent class.
	// GLDrawSomething can just put a getClassName()
}

YADE_PLUGIN((NodeRecorder));

YADE_REQUIRE_FEATURE(PHYSPAR);


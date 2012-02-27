/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "NonLocalInitializer.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
#include<yade/pkg-lattice/LatticeSetParameters.hpp>
#include<yade/pkg-lattice/NonLocalDependency.hpp>


NonLocalInitializer::NonLocalInitializer () : GlobalEngine()
{
}

NonLocalInitializer::~NonLocalInitializer ()
{
}

bool NonLocalInitializer::calcNonLocal(Body* body1, Body* body2, BodyContainer* bodies, void* nonl)
//void NonLocalInitializer::calcNonLocal(Body* body1, Body* body2, BodyContainer* bodies, InteractionContainer* nonl)
//bool NonLocalInitializer::calcNonLocal(Body* body1, Body* body2, BodyContainer* bodies,//std::list<std::list<LatticeSetParameters::NonLocalInteraction> >
//std::vector<std::list<LatticeSetParameters::NonLocalInteraction , std::__malloc_alloc_template<sizeof(LatticeSetParameters::NonLocalInteraction)> > >& nonl)
{
	LatticeBeamParameters* beam1 	= static_cast<LatticeBeamParameters*>(body1->physicalParameters.get());
	LatticeBeamParameters* beam2 	= static_cast<LatticeBeamParameters*>(body2->physicalParameters.get());
	
//	static Real sqPi = std::sqrt(Mathr::PI);
	
	//bool exists = false;
	/*
	std::vector<LatticeSetParameters::NonLocalInteraction>::iterator nb = nonl.begin();
	std::vector<LatticeSetParameters::NonLocalInteraction>::iterator ne = nonl.end();
	for( ; nb != ne ; ++nb )
		if( 	   ( nb->id1 == body1->getId() && nb->id2 == body2->getId() )
			|| ( nb->id2 == body1->getId() && nb->id1 == body2->getId() ) )
			{
				exists = true;
				break;
			}
	*/
	
	//std::list<std::list<LatticeSetParameters::NonLocalInteraction> >
	/*
	std::vector<std::list<LatticeSetParameters::NonLocalInteraction , std::__malloc_alloc_template<sizeof(LatticeSetParameters::NonLocalInteraction)> > >::iterator nonlOuter = nonl.begin();
	
	static unsigned long int check =0;
	static unsigned long int idx = 0;
	++check;
	if(check > 4294967290ul )
	//if(check > 90000 )
	{
		check = 0;
		++idx;
	}
	if( idx >= nonl.size() )
	{
		cerr << "making next list, first is full.\n";
		nonl.resize( nonl.size()+1 );
	}
	for( int ii = 0 ; ii < idx ; ++ii )
		++nonlOuter;
	*/	
	/*
	std::list<std::list<LatticeSetParameters::NonLocalInteraction> >::iterator nonlEnd   = nonl.end();
	while( nonlOuter->size() == nonlOuter->max_size()-1 )
	{
		++nonlOuter;
		if( nonlOuter == nonlEnd )
		{
			cerr << "making next list, first is full.\n";
			nonl.resize( nonl.size()+1 );
			nonlOuter = nonl.end();
			--nonlOuter;
			nonlEnd = nonl.end();
		}
	}
	*/
	
	//if( ! nonl->find( body1->getId() , body2->getId() ) )
	//if( ! exists )
	{
		Real sqDist = (beam1->se3.position - beam2->se3.position).SquaredLength();
	//	Real dist = (beam1->se3.position - beam2->se3.position).length();
		Real len  = range;
		if( sqDist < std::pow(3 * len,2) )
	//	if( dist < 3 * len )
		// 3 * characteristic length give less than 0.1 % (2*  gives 1.83% and 1*  gives 37 % )
		{
			static unsigned long int total = 0;
			//Real dist = sqrt(sqDist);
			/*
			shared_ptr<Interaction> 		interaction(new Interaction( body1->getId() , body2->getId() ));
			shared_ptr<NonLocalDependency> 		nonLocal(new NonLocalDependency);
			
			nonLocal->gaussValue 			= std::exp( - std::pow( dist / len , 2) ) / ( len * sqPi );
			
			interaction->interactionPhysics 	= nonLocal;
			nonl->insert(interaction);
			*/
			LatticeSetParameters::NonLocalInteraction another;
			//another.gaussValue 			= std::exp( - std::pow( dist / len , 2) ) / ( len * sqPi );
	/////////		another.gaussValue 			= std::exp( - (sqDist/(len*len)) ) / ( len * sqPi );
			another.id1 				= body1->getId();
			another.id2 				= body2->getId();
			//nonlOuter->push_back(another);
			reinterpret_cast<LatticeSetParameters::NonLocalInteraction *>(nonl)[total] = another;
			++total;
			if( total > 121637000ul )
				cerr << "too much\n", exit(0);
			return true;
		}
		else
			return false;
	} 
}

void NonLocalInitializer::action(Scene * rootBody)
{
	return; // YADEview

	int beamGroupMask = YADE_CAST<LatticeSetParameters*>(rootBody->physicalParameters.get())->beamGroupMask;
	//std::list<std::list<LatticeSetParameters::NonLocalInteraction> >
	//std::vector<std::list<LatticeSetParameters::NonLocalInteraction , std::__malloc_alloc_template<sizeof(LatticeSetParameters::NonLocalInteraction)> > >& nonl = dynamic_cast<LatticeSetParameters*>(rootBody->physicalParameters.get())->nonl;
	void* nonl = YADE_CAST<LatticeSetParameters*>(rootBody->physicalParameters.get())->nonl;
	YADE_CAST<LatticeSetParameters*>(rootBody->physicalParameters.get())->range=range;

	BodyContainer::iterator bi    = rootBody->bodies->begin();
	BodyContainer::iterator bi2;
	BodyContainer::iterator biEnd = rootBody->bodies->end();
	
	//rootBody->interactions->clear();
	//nonl.clear();
	//nonl.resize(1);
	
	nonl = calloc(121637044ul,sizeof(LatticeSetParameters::NonLocalInteraction));
	if( nonl == 0 )
		cerr << "cannot allocate memory\n",exit(0);
	YADE_CAST<LatticeSetParameters*>(rootBody->physicalParameters.get())->nonl = nonl;
	
	unsigned long int& total = YADE_CAST<LatticeSetParameters*>(rootBody->physicalParameters.get())->total;
	total=0;

	int beam_counter = 0;
	float nodes_a=0;
	float nodes_all = rootBody->bodies->size();
	for(  ; bi!=biEnd ; ++bi )  // loop over all beams to create non-local
	{
		Body* body1 = (*bi).get(); // first_node
		
		if( ! ( body1->getGroupMask() & beamGroupMask ) )
			continue; // skip non-beams
			
		//total++;
		
	
		bi2 = bi;
		nodes_a+=1.0;
		
		for( ; bi2!=biEnd ; ++bi2 )
		{
			Body* body2 = (*bi2).get(); // all other beams (including this one)

			if( ! ( body2->getGroupMask() & beamGroupMask ) )
				continue; // skip non-beams
			
			//calcNonLocal(body1,body2,rootBody->bodies.get(),rootBody->interactions.get());
			if( calcNonLocal(body1,body2,rootBody->bodies.get(),nonl) )
				++total; 
		}
		if( ++beam_counter % 400 == 0 )
			cerr << "non-local dependency "<< total <<": " << ((nodes_a/nodes_all)*100.0)  << " %\n"; 
	}
//	cerr << "number of beams: " << total << "\n";
}

YADE_PLUGIN((NonLocalInitializer));

YADE_REQUIRE_FEATURE(PHYSPAR);


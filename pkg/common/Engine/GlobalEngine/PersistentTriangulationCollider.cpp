/*************************************************************************
*  Copyright (C) 2009 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>         *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include"PersistentTriangulationCollider.hpp"
#include<yade/core/Body.hpp>
#include<yade/core/Scene.hpp>
#include<yade/core/BodyContainer.hpp>
#include<limits>
#include<utility>
#include<vector>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/ElastMat.hpp>


YADE_REQUIRE_FEATURE(CGAL)
using namespace std;
		
// PersistentTriangulationCollider::PersistentTriangulationCollider() : Collider()
// {
// 	haveDistantTransient=false;
// 	isTriangulated = false;
// 	Tes = new ( TesselationWrapper );
// 
// 	nbObjects=0;
// 	xBounds.clear();
// 	yBounds.clear();
// 	zBounds.clear();
// 	minima.clear();
// 	maxima.clear();
// }


PersistentTriangulationCollider::~PersistentTriangulationCollider()
{
	delete Tes;
}

void PersistentTriangulationCollider::action ()
{
	shared_ptr<BodyContainer> bodies=scene->bodies;
	bool triangulationIteration = false;

//BEGIN VORONOI TESSELATION
	if ( !isTriangulated || Omega::instance().getCurrentIteration() % 50 == 0 )
	{
//TesselationWrapper Tes;
		Tes->clear();
		BodyContainer::iterator bi    = bodies->begin();
		BodyContainer::iterator biEnd = bodies->end();
		for ( ; bi!=biEnd ; ++bi )
		{
			if ( ( *bi )->isDynamic )
			{//means "is it a sphere (not a wall)"
				const Sphere* s = YADE_CAST<Sphere*> ( ( *bi )->shape.get() );
				Tes->insert ( (*bi)->state->pos[0],(*bi)->state->pos[1],(*bi)->state->pos[2], s->radius, ( *bi )->getId() );
			}
		}
		Tes->AddBoundingPlanes();
		isTriangulated = true;
		triangulationIteration = true;
		//}
//  else
//  {
//   //if (Omega::instance().getCurrentIteration() % 100 == 0) {
//   if (1) {
//   Tes->RemoveBoundingPlanes();
//   BodyContainer::iterator bi    = bodies->begin();
//   BodyContainer::iterator biEnd = bodies->end();
//   for ( ; bi!=biEnd ; ++bi )
//   {
//    if ( ( *bi )->isDynamic )
//    {//means "is it a sphere (not a wall)"
//     const Sphere* s = YADE_CAST<Sphere*> ( ( *bi )->shape.get() );
//     const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*> ( ( *bi )->physicalParameters.get() );
//     Tes->move ( p->se3.position[0],p->se3.position[1],p->se3.position[2], s->radius, ( *bi )->getId() );
// //    else {
// //    const Box* s = YADE_CAST<Box*>((*bi)->shape.get());
// //    const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*>((*bi)->physicalParameters.get());
// //    Tes.move(p->se3.position[0],p->se3.position[1],p->se3.position[2], s->radius, (*bi)->getId());
// //    }
//    }
//   }
//   Tes->AddBoundingPlanes();
//   }
	}
//ENDOF VORONOI TESSELATION


	interactions = scene->interactions;
	InteractionContainer::iterator I_end = interactions->end();
	for ( InteractionContainer::iterator I=interactions->begin(); I!=I_end; ++I )
	{
		// FIXME: eudoxos commented out as isReal and isNew is removed...
		//  if ( ( *I )->isReal ) ( *I )->isNew=false;
		//  if ( !haveDistantTransient ) ( *I )->isReal=false; // reset this flag, is used later... (??)
		( *I )->isNeighbor = false;// will be set again just below
	}

	if ( triangulationIteration )
	{
		std::pair<unsigned int,unsigned int> interaction_pair;
		unsigned int& id1 = interaction_pair.first;
		unsigned int& id2 = interaction_pair.second;

		unsigned int numberOfInteractions = Tes->NumberOfFacets(true);
		for ( unsigned int i=0; i<numberOfInteractions; ++i )
		{
			Tes->nextFacet ( interaction_pair );

			// look if the pair (id1,id2) already exists in the overleppingBB collection
			const shared_ptr<Interaction>& interaction=interactions->find ( body_id_t ( id1 ),body_id_t ( id2 ) );
			bool found= ( interaction!=0 );//Bruno's Hack

			// inserts the pair p=(id1,id2) if the two Aabb overlaps and if p does not exists in the overlappingBB
			if ( !found )
			{
				interactions->insert ( body_id_t ( id1 ),body_id_t ( id2 ) );
				//cerr << "inserted " << id1 << "-" << id2<<endl;
			}
			else interaction->isNeighbor = true;
			// removes the pair p=(id1,id2) if the two Aabb do not overlapp any more and if p already exists in the overlappingBB
		}

		vector< pair<unsigned int,unsigned int> > toErase;
		I_end = interactions->end();
		for ( InteractionContainer::iterator I=interactions->begin(); I!=I_end; ++I )
		{
			if (!(*I)->isNeighbor && !(haveDistantTransient && ( *I )->isReal()))//FIXME : isReal correct to check that the interaction has not been "requestErased"?
			{
				toErase.push_back ( pair<unsigned int,unsigned int> ( ( *I )->getId1() , ( *I )->getId2() ) );
				//cerr << "to delete " << ( *I )->getId1() << "-" << ( *I )->getId2() << "(isNeighbor=" << ( *I )->isNeighbor<< endl;
			}
			//interactions->erase ( ( *I )->getId1() , ( *I )->getId2() );
		}
		vector< pair<unsigned int,unsigned int> >::iterator it = toErase.begin();
		vector< pair<unsigned int,unsigned int> >::iterator it_end = toErase.end();
		for ( ;it!=it_end;++it )
		{
			interactions->erase ( it->first , it->second );
			//cerr << "deleted " << it->first << "-" << it->second<<endl;
		}
	}
}

YADE_PLUGIN((PersistentTriangulationCollider));



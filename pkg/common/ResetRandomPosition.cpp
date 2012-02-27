/*************************************************************************
*  Copyright (C) 2009 by Sergei Dorofeenko				 				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef pi
	#undef pi
#endif

#include<boost/random.hpp>
#include<yade/core/Body.hpp>
#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Facet.hpp>
#include<yade/pkg-common/InteractionLoop.hpp>
//#include<yade/pkg-dem/BodyMacroParameters.hpp>
#include"ResetRandomPosition.hpp"
#include<sstream>

YADE_PLUGIN((ResetRandomPosition));
CREATE_LOGGER(ResetRandomPosition);

boost::variate_generator<boost::mt19937,boost::uniform_real<> > 
	ResetRandomPosition::randomUnit(boost::mt19937(),boost::uniform_real<>(0,1));
boost::variate_generator<boost::mt19937,boost::uniform_real<> >
	ResetRandomPosition::randomSymmetricUnit(boost::mt19937(),boost::uniform_real<>(-1,1));

void ResetRandomPosition::action()
{
	if (first_run)
	{
		FOREACH(shared_ptr<Engine> eng, scene->engines)
		{
			bI=dynamic_cast<Collider*>(eng.get());
			if (bI) break;
		}
		if (!bI) 
		{
			LOG_FATAL("Can't find Collider." );
			return;
		}
		iGME=dynamic_cast<IGeomDispatcher*>(scene->engineByName("IGeomDispatcher").get());
		if (!iGME) 
		{
			InteractionLoop* iDsp=dynamic_cast<InteractionLoop*>(scene->engineByName("InteractionLoop").get());
			if (!iDsp)
			{
				LOG_FATAL("Can't find nor IGeomDispatcher nor InteractionLoop." );
				return;
			}
			iGME=dynamic_cast<IGeomDispatcher*>(iDsp->geomDispatcher.get());
			if (!iGME)
			{
				LOG_FATAL("Can't find IGeomDispatcher." );
				return;
			}
		}
		first_run=false;
		randomFacet= shared_ptr<RandomInt>(new RandomInt(boost::minstd_rand(),boost::uniform_int<>(0,factoryFacets.size()-1)));
	}

	shiftedBodies.clear();

	FOREACH(int id, subscribedBodies)
	{
		shared_ptr<Body> b = Body::byId(id);
		State* rb = b->state.get();
		Vector3r& position = rb->se3.position;

		if ( (position-point).dot(normal) < 0 )
		{
			Vector3r backup_pos = position;

			bool is_overlap;

			for (int attempt=0; attempt<maxAttempts; ++attempt)
			{
				position = (volumeSection) ? generatePositionInVolume():generatePositionOnSurface();

				const Real r = YADE_CAST<Sphere*>(b->shape.get())->radius;
				Bound bv;
				bv.min = Vector3r(position[0]-r, position[1]-r, position[2]-r);
				bv.max = Vector3r(position[0]+r, position[1]+r, position[2]+r);

				is_overlap=false;

				// Test overlap with already shifted bodies
				FOREACH(shared_ptr<Body> sb, shiftedBodies)
				{
					if (iGME->explicitAction(b,sb,/*force*/false)->geom)
					{
						is_overlap=true;
						break;
					}
				}
				if (is_overlap) continue; // new attempt

				// Test overlap with other bodies
				vector<Body::id_t> probedBodies=bI->probeBoundingVolume(bv);
				FOREACH(Body::id_t id, probedBodies){
					if (iGME->explicitAction(b,Body::byId(id),/*force*/false)->geom){
						is_overlap=true;
						break;
					}
				}
				if (is_overlap) continue; // new attempt
				break;
			}
			if (is_overlap)
			{
				LOG_WARN("Can't placing sphere during " << maxAttempts << " attempts.");
				position=backup_pos;
				return;
			}

			rb->vel		= Vector3r(//
					velocity[0]+velocityRange[0]*randomSymmetricUnit(),
					velocity[1]+velocityRange[1]*randomSymmetricUnit(),
					velocity[2]+velocityRange[2]*randomSymmetricUnit());
			rb->angVel = Vector3r(//
					angularVelocity[0]+angularVelocityRange[0]*randomSymmetricUnit(),
					angularVelocity[1]+angularVelocityRange[1]*randomSymmetricUnit(),
					angularVelocity[2]+angularVelocityRange[2]*randomSymmetricUnit());

			shiftedBodies.push_back(b);
		} 
	} // next sphere
}

Vector3r ResetRandomPosition::generatePositionOnSurface()
{
    Body::id_t facetId = factoryFacets[(*randomFacet)()];
    Real t1 = randomUnit();
    Real t2 = randomUnit()*(1-t1);

    shared_ptr<Body> facet = Body::byId(facetId);
	Facet* ifacet = static_cast<Facet*>(facet->shape.get());

	return t1*(ifacet->vertices[1]-ifacet->vertices[0])+t2*(ifacet->vertices[2]-ifacet->vertices[0])+ifacet->vertices[0]+facet->state->se3.position;
    
}

Vector3r ResetRandomPosition::generatePositionInVolume()
{
    Vector3r p1 = generatePositionOnSurface();
    Vector3r p2 = generatePositionOnSurface();
    Real t = randomUnit();
    return p1+t*(p2-p1);
    
}


//YADE_REQUIRE_FEATURE(PHYSPAR);


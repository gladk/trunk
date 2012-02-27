/*************************************************************************
*  Copyright (C) 2009 by Sergei Dorofeenko				 				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include <yade/pkg-common/PeriodicEngines.hpp>
#include <yade/pkg-common/Dispatching.hpp>
#include <yade/pkg-common/Collider.hpp>
#include <yade/core/Scene.hpp>
#include <vector>
#include <string>

using namespace std;

/// @brief Produces spheres over the course of a simulation. 
class ResetRandomPosition : public PeriodicEngine {
public:
	/// @brief Create one sphere per call.
	virtual void action();
private:
	/// @brief Pointer to Collider.
	/// It is necessary in order to probe the bounding volume for new sphere.
	Collider* bI;
	/// @brief Pointer to IGeomDispatcher.
	/// It is necessary in order to detect a real overlap with other bodies.
	IGeomDispatcher* iGME;

	std::vector<shared_ptr<Body> > shiftedBodies;

	bool first_run;
	//bool generateNewPosition(const shared_ptr<Body>& b, Vector3r& new_position);
   Vector3r generatePositionOnSurface();
   Vector3r generatePositionInVolume();

	typedef	boost::variate_generator<boost::minstd_rand,boost::uniform_int<> > RandomInt;
	shared_ptr<RandomInt> randomFacet;

	static boost::variate_generator<boost::mt19937,boost::uniform_real<> > 
		randomUnit;
	static boost::variate_generator<boost::mt19937,boost::uniform_real<> >
		randomSymmetricUnit;

	DECLARE_LOGGER;
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(ResetRandomPosition,GlobalEngine,"Creates spheres during simulation, placing them at random positions. Every time called, one new sphere will be created and inserted in the simulation.",
		((vector<Body::id_t>,factoryFacets,,,"The geometry of the section where spheres will be placed; they will be placed on facets or in volume between them depending on *volumeSection* flag."))
		((std::vector<int>,subscribedBodies,,,"Affected bodies."))
		((Vector3r,point,Vector3r::Zero(),,"??"))
		((Vector3r,normal,Vector3r(0,1,0),,"??"))
		((bool,volumeSection,((void)"define factory by facets.",false),,"Create new spheres inside factory volume rather than on its surface."))
		((int,maxAttempts,20,,"Max attempts to place sphere. If placing the sphere in certain random position would cause an overlap with any other physical body in the model, SpheresFactory will try to find another position."))
		((Vector3r,velocity,Vector3r::Zero(),,"Mean velocity of spheres."))
		((Vector3r,velocityRange,Vector3r::Zero(),,"Half size of a velocities distribution interval. New sphere will have random velocity within the range velocity±velocityRange."))
		((Vector3r,angularVelocity,Vector3r::Zero(),,"Mean angularVelocity of spheres."))
		((Vector3r,angularVelocityRange,Vector3r::Zero(),,"Half size of a angularVelocity distribution interval. New sphere will have random angularVelocity within the range angularVelocity±angularVelocityRange.")),
		first_run=true;
	);
};
REGISTER_SERIALIZABLE(ResetRandomPosition);



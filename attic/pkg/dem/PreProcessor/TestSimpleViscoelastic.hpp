/*************************************************************************
*  Copyright (C) 2006 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <yade/core/FileGenerator.hpp>

class TestSimpleViscoelastic : public FileGenerator {
    private:
	
	/// Collision time
	Real tc; 
	/// Normal restitution coefficient
	Real en;
	/// Shear restitution coefficient
       	Real es;
	/// Coulomb friction
	Real frictionAngle;
	/// Spheres radius
	Real R;
	/// Materials density
	Real density;
	/// Number of spheres
	Real nbSpheres;
	/// Height between spheres
	Real h;
	/// Initial velocity of spheres
	Vector3r velocity; 
	/// Initial angular velocity of spheres
	Vector3r angular_vel; 
	/// Gravity
	Vector3r gravity;
	/// Size of bottom box
	Vector3r groundSize;
	/// Base name of results file
	std::string outputBase;
	/// Saves interval
	int interval;
	
	bool rotationBlocked;

	void createActors(shared_ptr<Scene>& rootBody);
	void createSphere(shared_ptr<Body>& body, int i);
	void createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents);
	void positionRootBody(shared_ptr<Scene>& rootBody);

    public:
	TestSimpleViscoelastic();
	~TestSimpleViscoelastic();

    protected:
	bool generate();

	virtual void postProcessAttributes(bool deserializing);
	REGISTER_ATTRIBUTES(FileGenerator,(tc)(en)(es)(frictionAngle)(R)(density)(nbSpheres)(h)(velocity)(angular_vel)(rotationBlocked)(gravity)(groundSize)(outputBase)(interval));
	REGISTER_CLASS_NAME(TestSimpleViscoelastic);
	REGISTER_BASE_CLASS_NAME(FileGenerator);
};

REGISTER_SERIALIZABLE(TestSimpleViscoelastic);



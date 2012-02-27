/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/FileGenerator.hpp>

class SDECSpheresPlane : public FileGenerator
{
	private :
		Vector3r	 nbSpheres
				,disorder
				,groundSize
				,gravity;

		bool		 useSpheresAsGround
				,rotationBlocked;

		Real		 minRadius
				,density
				,maxRadius
				,dampingForce
				,dampingMomentum
				,spheresHeight
				,sphereYoungModulus
				,spherePoissonRatio
				,sphereFrictionDeg;

		int		 timeStepUpdateInterval;
	
		void createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents);
		void createSphere(shared_ptr<Body>& body, int i, int j, int k);
		void createGroundSphere(shared_ptr<Body>& body,Real radius, Real i, Real j, Real k);
		void createActors(shared_ptr<Scene>& rootBody);
		void positionRootBody(shared_ptr<Scene>& rootBody);

	// construction
	public :
		SDECSpheresPlane ();
		~SDECSpheresPlane ();
		bool generate();

	protected :
		virtual void postProcessAttributes(bool deserializing);
	REGISTER_ATTRIBUTES(FileGenerator,(nbSpheres)(minRadius)(maxRadius)(spheresHeight)(sphereYoungModulus)(spherePoissonRatio)(sphereFrictionDeg)(gravity)(density)(disorder)(groundSize)(useSpheresAsGround)(dampingForce)(dampingMomentum)(rotationBlocked)(timeStepUpdateInterval));
	REGISTER_CLASS_NAME(SDECSpheresPlane);
	REGISTER_BASE_CLASS_NAME(FileGenerator);
};

REGISTER_SERIALIZABLE(SDECSpheresPlane);




// Copyright (C) 2004 by Olivier Galizzi <olivier.galizzi@imag.fr>
// Copyright (C) 2004 by Janek Kozicki <cosurgi@berlios.de>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#ifndef YADE_PHYSPAR
	#error This file can be included only with the deprecated feature 'physpar' enabled!
#endif

#include<yade/lib-serialization/Serializable.hpp>
#include<yade/lib-multimethods/Indexable.hpp>
#include<yade/lib-base/yadeWm3Extra.hpp>
#include<vector>
#include<string>

class PhysicalParameters: public Serializable, public Indexable
{
	public:
			PhysicalParameters(): se3(Vector3r(0,0,0),Quaternionr(1,0,0,0)), refSe3(Vector3r(0,0,0),Quaternionr(1,0,0,0)), isDisplayed(true), blockedDOFs(0) {}
		Se3r se3; //! Se3 of the body in simulation space

		Se3r refSe3; //! Reference Se3 of body in display space
		// the following two are recomputed at every renderering cycle, no need to save/load them via registerAttributes
		Se3r dispSe3; //! Se3 of body in display space (!=se3 if scaling positions/orientations)
		bool isDisplayed; //! False if the body is not displayed in this cycle (clipped, for instance)

		unsigned blockedDOFs; //! Bitmask for degrees of freedom where velocity will be always zero, regardless of applied forces

		// mutex for updating the parameters from within the interaction loop
		boost::mutex updateMutex;

		enum {DOF_NONE=0,DOF_X=1,DOF_Y=2,DOF_Z=4,DOF_RX=8,DOF_RY=16,DOF_RZ=32};
		static const unsigned DOF_ALL=DOF_X|DOF_Y|DOF_Z|DOF_RX|DOF_RY|DOF_RZ; //! shorthand for all DOFs blocked
		static const unsigned DOF_XYZ=DOF_X|DOF_Y|DOF_Z; //! shorthand for all displacements blocked
		static const unsigned DOF_RXRYRZ=DOF_RX|DOF_RY|DOF_RZ; //! shorthand for all rotations blocked
		void setDOFfromVector3r(Vector3r disp,Vector3r rot=Vector3r(0,0,0))
		{
			//! set DOFs according to two Vector3r arguments (blocked is when disp[i]==1.0 or rot[i]==1.0)
			blockedDOFs = ((disp[0]==1.0)?DOF_X :0)|((disp[1]==1.0)?DOF_Y :0)|((disp[2]==1.0)?DOF_Z :0)|
			              ((rot [0]==1.0)?DOF_RX:0)|((rot [1]==1.0)?DOF_RY:0)|((rot [2]==1.0)?DOF_RZ:0);
			//Vector3<bool> would be actually better, I'm not checking now if it would work.
		};
		static unsigned axisDOF(int axis, bool rotationalDOF=false){return 1<<(axis+(rotationalDOF?3:0));} //! Return DOF_* constant for given axis∈{0,1,2} and rotationalDOF∈{false(default),true}; e.g. axisDOF(0,true)==DOF_RX
		std::vector<std::string> blockedDOFs_vec_get() const;
		void blockedDOFs_vec_set(const std::vector<std::string>& dofs);

	REGISTER_ATTRIBUTES(Serializable,(se3)(refSe3)(blockedDOFs));
	REGISTER_CLASS_AND_BASE(PhysicalParameters,Serializable Indexable);
	REGISTER_INDEX_COUNTER(PhysicalParameters);
};

REGISTER_SERIALIZABLE(PhysicalParameters);

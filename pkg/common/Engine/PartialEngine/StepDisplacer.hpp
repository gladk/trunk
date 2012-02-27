// 2008 © Václav Šmilauer <eudoxos@arcig.cz> 
#pragma once

#include<yade/core/PartialEngine.hpp>

class StepDisplacer: public PartialEngine {
	public:
		virtual void action();
	YADE_CLASS_BASE_DOC_ATTRS(StepDisplacer,PartialEngine,"Apply generalized displacement (displacement or rotation) stepwise on subscribed bodies.",
		((Se3r,deltaSe3,Se3r(Vector3r::Zero(),Quaternionr::Identity()),"Difference of position/orientation that will be added. Position is added to current :yref:`State::pos` using vector addition, orientation to :yref:`State::ori` using quaternion multiplication (rotation composition)."))
		((bool,setVelocities,false,"If true,  velocity and angularVelocity are modified in such a way that over one iteration (dt), the body will have prescribed se3 jump. In this case, se3 itself is not updated for :yref:`dynamic<Body::isDynamic>` bodies, since they would have the delta applied twice (here and in the :yref:`integrator<NewtonIntegrator>`). For non-dynamic bodies however, se3 *is* still updated."))
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(StepDisplacer);
	


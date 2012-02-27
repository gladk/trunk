// 2008 © Václav Šmilauer <eudoxos@arcig.cz> 
#pragma once

#include<yade/core/PartialEngine.hpp>
#include<yade/core/PhysicalParameters.hpp>
#include<yade/pkg-common/RigidBodyParameters.hpp>
#include<Python.h>

/* Engine interpolating between starting (current) and goal (given) se3, both position and orientation.
 *
 * The interpolation is linear, either in iteration space or (real/virtual) time space.
 *
 * A given python code can be executed when the ending se3 is reached and the engine deactivates itself.
 *
 * Note that eactly one body must be subscribed to this engine, since interpolating positions
 * on multiple bodies is meaningless.
 *
 * Attributes that should be assigned:
 * 	subscribedBodies (exactly one)
 * 	goal (7-tuple of position and orientation)
 * 	startVirt XOR startReal XOR startIter (at what time point to start)
 * 	goalVirt XOR goalReal XOR goalIter (at what time point to reach desired se3)
 * 	goalHook (optional: python command to be run when finished)
 * 	rotRelative (optional: rotation part of goal is relative to the rotation at the beginning)
 */
class Se3Interpolator: public PartialEngine {
	public:
		bool done,started;
		Se3r start;
		Se3r goal;
		bool rotRelative;
		Real goalVirt, goalReal; long goalIter;
		long startVirt, startReal; long startIter;
		string goalHook;
		virtual bool isActivated(Scene*){return !done;}
		Se3Interpolator(): done(false), started(false), rotRelative(false), goalVirt(0), goalReal(0), goalIter(0) {};
		virtual ~Se3Interpolator(){};
		virtual void action(Scene* mb);
	protected:
		REGISTER_ATTRIBUTES(PartialEngine,
			(done)
			(started)
			(start)
			(goal)
			(rotRelative)
			(goalVirt)
			(goalReal)
			(goalIter)
			(startVirt)
			(startReal)
			(startIter)
			(goalHook)
		);
		REGISTER_CLASS_AND_BASE(Se3Interpolator,PartialEngine);
		DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Se3Interpolator);
	



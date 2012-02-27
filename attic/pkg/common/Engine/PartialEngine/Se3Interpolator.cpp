
#include"Se3Interpolator.hpp"
#include<yade/pkg-common/PeriodicEngines.hpp>
#include<yade/lib-pyutil/gil.hpp>

YADE_PLUGIN((Se3Interpolator))
CREATE_LOGGER(Se3Interpolator);

void Se3Interpolator::action(Scene* mb){
	assert(subscribedBodies.size()==1);
	const shared_ptr<Body>&b=Body::byId(subscribedBodies[0],mb);
	Omega& O=Omega::instance();
	if(!started){
		started=true;
		start=b->physicalParameters->se3;
		if(rotRelative) goal.orientation=goal.orientation*start.orientation;
		int numGoalEndpoints=0;
		if(goalVirt>0){startVirt=O.getSimulationTime(); numGoalEndpoints++;}
		if(goalReal>0){startReal=PeriodicEngine::getClock(); numGoalEndpoints++;}
		if(goalIter>0){startIter=O.getCurrentIteration(); numGoalEndpoints++;}
		assert(numGoalEndpoints==1);
		LOG_DEBUG("Init'ed.");
	}
	Real t=-1; // t∈[0,1]
	if(goalVirt>0) t=(O.getSimulationTime()-startVirt)/(goalVirt-startVirt);
	if(goalReal>0) t=(PeriodicEngine::getClock()-startReal)/(goalReal-startReal);
	if(goalIter>0) t=(O.getCurrentIteration()-startIter)/(goalIter-startIter);
	assert(t>=0);
	LOG_DEBUG("t="<<t);
	t=min(t,1.);
	b->physicalParameters->se3.position=start.position+t*(goal.position-start.position);
	b->physicalParameters->se3.orientation=Quaternionr().Slerp(t,start.orientation,goal.orientation);
	if(t>=1.){
		done=true;
		LOG_DEBUG("Goal reached.");
		if(!goalHook.empty()) PyRunString(goalHook);
	}
}

YADE_REQUIRE_FEATURE(PHYSPAR);


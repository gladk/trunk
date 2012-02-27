// 2008 © Václav Šmilauer <eudoxos@arcig.cz> 
#include"UniaxialStrainer.hpp"
#include<boost/foreach.hpp>

#include<yade/core/Scene.hpp>
#include<yade/core/InteractionContainer.hpp>
#include<yade/pkg-common/Aabb.hpp>

YADE_PLUGIN((UniaxialStrainer));
/************************ UniaxialStrainer **********************/
CREATE_LOGGER(UniaxialStrainer);

void UniaxialStrainer::init(){
	needsInit=false;

	assert(posIds.size()>0);
	assert(negIds.size()>0);
	posCoords.clear(); negCoords.clear();
	FOREACH(body_id_t id,posIds){ const shared_ptr<Body>& b=Body::byId(id,scene); posCoords.push_back(b->state->pos[axis]);
		if(blockDisplacements && blockRotations) b->isDynamic=false;
		else{
			if(!blockDisplacements) b->state->blockedDOFs=State::axisDOF(axis); else b->state->blockedDOFs=State::DOF_XYZ;
			if(blockRotations) b->state->blockedDOFs|=State::DOF_RXRYRZ;
		}
	}
	FOREACH(body_id_t id,negIds){ const shared_ptr<Body>& b=Body::byId(id,scene); negCoords.push_back(b->state->pos[axis]);
		if(blockDisplacements && blockRotations) b->isDynamic=false;
		else{
			if(!blockDisplacements) b->state->blockedDOFs=State::axisDOF(axis); else b->state->blockedDOFs=State::DOF_XYZ;
			if(blockRotations) b->state->blockedDOFs|=State::DOF_RXRYRZ;
		}
	}

	assert(posIds.size()==posCoords.size() && negIds.size()==negCoords.size());

	originalLength=axisCoord(posIds[0])-axisCoord(negIds[0]);
	LOG_DEBUG("Reference particles: positive #"<<posIds[0]<<" at "<<axisCoord(posIds[0])<<"; negative #"<<negIds[0]<<" at "<<axisCoord(negIds[0]));
	LOG_INFO("Setting initial length to "<<originalLength<<" (between #"<<negIds[0]<<" and #"<<posIds[0]<<")");
	if(originalLength<=0) LOG_FATAL("Initial length is negative or zero (swapped reference particles?)! "<<originalLength);
	/* this happens is nan propagates from e.g. brefcom consitutive law in case 2 bodies have _exactly_ the same position
	 * (the the normal strain is 0./0.=nan). That is an user's error, however and should not happen. */
	if(isnan(originalLength)) LOG_FATAL("Initial length is NaN!");
	assert(originalLength>0 && !isnan(originalLength));

	assert(!isnan(strainRate) || !isnan(absSpeed));
	if(!isnan(std::numeric_limits<Real>::quiet_NaN())){ LOG_FATAL("NaN's are not properly supported (compiled, with -ffast-math?), which is required."); throw; }

	if(isnan(strainRate)){ strainRate=absSpeed/originalLength; LOG_INFO("Computed new strainRate "<<strainRate); }
	else {absSpeed=strainRate*originalLength;}

	if(!setSpeeds){
		initAccelTime_s=initAccelTime>=0 ? initAccelTime : scene->dt*(-initAccelTime);
		LOG_INFO("Strain speed will be "<<absSpeed<<", strain rate "<<strainRate<<", will be reached after "<<initAccelTime_s<<"s ("<<initAccelTime_s/scene->dt<<" steps).");
	} else {
		/* set speed such that it is linear on the strained axis; transversal speed is not set, which can perhaps create some problems.
			Note: all bodies in the simulation will have their speed set, since there is no way to tell which ones are part of the specimen
			and which are not.

			Speeds will be linearly interpolated beween axis positions p0,p1 and velocities v0,v1.
		*/
		initAccelTime_s=0;
		LOG_INFO("Strain speed will be "<<absSpeed<<", strain rate "<<strainRate<<"; velocities will be set directly at the beginning.");
		Real p0=axisCoord(negIds[0]), p1=axisCoord(posIds[0]); // limit positions
		Real v0,v1; // speeds at p0, p1
		switch(asymmetry){
			case -1: v0=-absSpeed; v1=0; break;
			case  0: v0=-absSpeed/2; v1=absSpeed/2; break;
			case  1: v0=0; v1=absSpeed; break;
			default: LOG_FATAL("Unknown asymmetry value "<<asymmetry<<" (should be -1,0,1)"); throw;
		}
		assert(p1>p0);
		// set speeds for particles on the boundary
		FOREACH(const shared_ptr<Body>& b, *scene->bodies){
			// skip bodies on the boundary, since those will have their positions updated directly
			if(std::find(posIds.begin(),posIds.end(),b->id)!=posIds.end() || std::find(negIds.begin(),negIds.end(),b->id)!=negIds.end()) { continue; }
			Real p=axisCoord(b->id);
			Real pNormalized=(p-p0)/(p1-p0);
			b->state->vel[axis]=pNormalized*(v1-v0)+v0;
		}
	}
	stressUpdateInterval=min(1000,max(1,(int)(1e-5/(abs(strainRate)*scene->dt))));
	LOG_INFO("Stress will be updated every "<<stressUpdateInterval<<" steps.");

	/* if we have default (<0) crossSectionArea, try to get it from root's Aabb;
	 * this will not work if there are foreign bodies in the simulation,
	 * in which case you must give the value yourself as engine attribute.
	 *
	 * A TODO option is to get crossSectionArea as average area of bounding boxes' of ABBBs
	 * of posIds and negIds perpendicular to axis. That might be better, except for cases where
	 * reference particles on either end do not coincide with the specimen cross-section.
	 *
	 * */
	if(crossSectionArea<=0){
		shared_ptr<Aabb> rbAABB;
		if (Omega::instance().getScene()->bound && (rbAABB=dynamic_pointer_cast<Aabb>(Omega::instance().getScene()->bound))){
			int axis2=(axis+1)%3, axis3=(axis+2)%3; // perpendicular axes indices
			Vector3r size=rbAABB->max-rbAABB->min;
			crossSectionArea=size[axis2]*size[axis3];
			LOG_INFO("Setting crossSectionArea="<<crossSectionArea<<", using axes #"<<axis2<<" and #"<<axis3<<".");
		} else {
			crossSectionArea=1.;
			LOG_WARN("No Axis Aligned Bounding BoxModel for scene, using garbage value ("<<crossSectionArea<<") for crossSectionArea!");
		}
	}
	assert(crossSectionArea>0);
}

void UniaxialStrainer::action(){
	if(needsInit) init();
	// postconditions for initParams
	assert(posIds.size()==posCoords.size() && negIds.size()==negCoords.size() && originalLength>0 && crossSectionArea>0);
	//nothing to do
	if(posIds.size()==0 || negIds.size()==0) return;
	// linearly increase strain to the desired value
	if(abs(currentStrainRate)<abs(strainRate)){
		Real t=Omega::instance().getSimulationTime();
		if(initAccelTime_s!=0) currentStrainRate=(t/initAccelTime_s)*strainRate;
		else currentStrainRate=strainRate;
	} else currentStrainRate=strainRate;
	// how much do we move (in total, symmetry handled below)
	Real dAX=currentStrainRate*originalLength*scene->dt;
	if(!isnan(stopStrain)){
		Real axialLength=axisCoord(posIds[0])-axisCoord(negIds[0]);
		Real newStrain=(axialLength+dAX)/originalLength-1;
		if((newStrain*stopStrain>0) && abs(newStrain)>=stopStrain){ // same sign of newStrain and stopStrain && over the limit from below in abs values
			dAX=originalLength*(stopStrain+1)-axialLength;
			LOG_INFO("Reached stopStrain "<<stopStrain<<", deactivating self and stopping in "<<idleIterations+1<<" iterations.");
			this->active=false;
			scene->stopAtIteration=Omega::instance().getCurrentIteration()+1+idleIterations;
		}
	}
	if(asymmetry==0) dAX*=.5; // apply half on both sides if straining symetrically
	if(asymmetry!=1){
		for(size_t i=0; i<negIds.size(); i++){
			negCoords[i]-=dAX;
			axisCoord(negIds[i])=negCoords[i]; // update current position
		}
	}
	if(asymmetry!=-1){
		for(size_t i=0; i<posIds.size(); i++){
			posCoords[i]+=dAX;
			axisCoord(posIds[i])=posCoords[i];
		}
	}

	Real axialLength=axisCoord(posIds[0])-axisCoord(negIds[0]);
	strain=axialLength/originalLength-1;

	// reverse if we're over the limit strain
	if(notYetReversed && limitStrain!=0 && ((currentStrainRate>0 && strain>limitStrain) || (currentStrainRate<0 && strain<limitStrain))) { currentStrainRate*=-1; notYetReversed=false; LOG_INFO("Reversed strain rate to "<<currentStrainRate); }

	// update forces and stresses
	if(Omega::instance().getCurrentIteration()%stressUpdateInterval==0) {
		computeAxialForce();
		avgStress=(sumPosForces+sumNegForces)/(2*crossSectionArea); // average nominal stress
	}
}

void UniaxialStrainer::computeAxialForce(){
	sumPosForces=sumNegForces=0;
	scene->forces.sync();
	FOREACH(body_id_t id, negIds) sumNegForces+=scene->forces.getForce(id)[axis];
	FOREACH(body_id_t id, posIds) sumPosForces-=scene->forces.getForce(id)[axis];
}

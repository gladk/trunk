// (c) 2007,2009 Vaclav Smilauer <eudoxos@arcig.cz> 

#include"Clump.hpp"
#include<algorithm>
#include<yade/core/Scene.hpp>
#include<yade/core/BodyContainer.hpp>
#include<yade/core/State.hpp>

YADE_PLUGIN((Clump)/*(ClumpMemberMover)*/);
CREATE_LOGGER(Clump);

#if 0
CREATE_LOGGER(ClumpMemberMover);
/**************************************************************************************
 ************************************* ClumpMemberMover ******************************
 **************************************************************************************/

/*! We only call clump's method, since it belongs there logically. It makes encapsulation of private members nicer, too.
 * _param pp passed by the dispatcher
 * _param clump passed by the dispatcher
 */
void ClumpMemberMover::action(){
	for(BodyContainer::iterator I=scene->bodies->begin(); I!=scene->bodies->end(); ++I){
		shared_ptr<Body> b = *I;
		if(b->isClump()){
			//LOG_TRACE("Applying movement to clump #"<<b->getId());
			static_cast<Clump*>(b.get())->moveMembers();
		}
	}
	//if(!clump->isDynamic) return; // perhaps clump that has been desactivated?!
}
#endif

/**************************************************************************************
 ******************************************** Clump ***********************************
 **************************************************************************************/

/*! @pre Body must be dynamic.
 * @pre Body must not be part or this clump already.
 * @pre Body must have valid (non-NULL) Body::physicalParameters
 * @pre Body must have id that is smaller than the clump's id (reason: processing order in NewtonIntegrator)
 */
void Clump::add(body_id_t subId){
	shared_ptr<Body> subBody=Body::byId(subId);

	// preconditions
	//assert(subBody->isDynamic);
	assert(state);
	assert(members.count(subId)==0);
	assert(subId<getId());

	// begin actual setup
	subBody->clumpId=getId();
	subBody->isDynamic=false;
	// for now, push just unitialized se3; will be updated by updateProperties
	members[subId]=Se3r();

	clumpId=getId(); // identifies a clump

	LOG_DEBUG("Added body #"<<subId<<" to clump #"<<getId());
}

/*! @pre Body with given id must be in the clump.
 */
void Clump::del(body_id_t subId){
	// erase the subBody; removing body that is not part of the clump is error
	assert(members.erase(subId)==1);
	// restore body's internal parameters;
	shared_ptr<Body> subBody=Body::byId(subId);
	subBody->clumpId=Body::ID_NONE;
	subBody->isDynamic=true;
	LOG_DEBUG("Removed body #"<<subId<<" from clump #"<<getId());
}

/*! @brief Calculate positions and orientations of members based on Clump's Se3; resets acceleration and angularAccelration to zero.
 *
 * This method is called by the ClumpMemberMover engine after each timestep.
 * @note Velocities of members are not updated, since members have isdynamic==false. It is possible, though, that someone needs to have a moving clump that is later broken apart and that liberated particle continue to move in the same way as they did within the clump. In that case, this will have to be completed.
 */
void Clump::moveMembers(){
	//const Se3r& mySe3(physicalParameters->se3);
	//const shared_ptr<RigidBodyParameters>& myRBP=static_pointer_cast<RigidBodyParameters>(physicalParameters);
	for(Clump::memberMap::iterator I=members.begin(); I!=members.end(); I++){
		// now, I->first is Body::id_t, I->second is Se3r of that body in the clump
		shared_ptr<Body> member=Body::byId(I->first);
		//const shared_ptr<RigidBodyParameters>& subRBP(YADE_PTR_CAST<RigidBodyParameters>(member->physicalParameters));
		State* subState=member->state.get();
		//LOG_TRACE("Old #"<<I->first<<"position: "<<subRBP->se3.position);
		subState->pos=state->pos+state->ori*I->second.position;
		subState->ori=state->ori*I->second.orientation;
		//LOG_TRACE("New #"<<I->first<<"position: "<<subRBP->se3.position);
		//LOG_TRACE("Clump #"<<getId()<<" moved #"<<I->first<<".");

		//! FIXME: we set velocity because of damping here; but since positions are integrated after all forces applied, these velocities will be used in the NEXT step for CundallNonViscousDamping. Does that matter?!
		//subState->vel=state->vel+state->angVel.Cross(I->second.position);
		subState->vel=state->vel+state->angVel.cross(subState->pos-state->pos);
		subState->angVel=state->angVel;
	}
	/* @bug Temporarily we reset acceleration and angularAcceleration of the clump here;
	 * should be a new negine that will take care of that?
	 */
	// const shared_ptr<RigidBodyParameters>& clumpRBP(YADE_PTR_CAST<RigidBodyParameters>(physicalParameters));
	#if 0
		if(Omega::instance().getCurrentIteration()%50==0){
			Real Erot=.5*clumpRBP->inertia[0]*pow(clumpRBP->angularVelocity[0],2)+.5*clumpRBP->inertia[1]*pow(clumpRBP->angularVelocity[1],2)+.5*clumpRBP->inertia[2]*pow(clumpRBP->angularVelocity[2],2);
			Real Etrans=.5*clumpRBP->mass*pow(clumpRBP->velocity.norm(),2);
			// (0,0,1) is gravity acceleration
			Real Epot=clumpRBP->se3.position.Dot(Vector3r(0,0,1))*clumpRBP->mass;
			LOG_TRACE("##"<<clumpId<<" energy "<<Erot+Etrans+Epot<<"\tv "<<Etrans<<"\tw "<<Erot<<"\tp "<<Epot);
		}
	#endif

	state->accel=state->angAccel=Vector3r::Zero();
}

/*! Clump's se3 will be updated (origin at centroid and axes coincident with principal inertia axes) and subSe3 modified in such a way that members positions in world coordinates will not change.

	The clump values that are changed are:
	-# Clump::members (holds position and orientation in clump's coordinate system)
	-# Clump::physicalParameters->mass (sum of masses of all members)
	-# Clump::physicalParameters->inertia (inertia of the aggregate - in clump coordinate system)
	-# Clump::physicalParameters->se3 (position and orientation of the clump; it is such that absolute positions and orientation of members will not chage)

	The algorithm is as follows:
	-# Clump::members values and Clump::physicalParameters::se3 are invalid from this point
	-# M=0; S=vector3r(0,0,0); I=zero tensor; (ALL calculations are in world coordinates!)
	-# loop over Clump::members (position x_i, mass m_i, inertia at subBody's centroid I_i) [this loop will be replaced by numerical integration (rasterization) for the intersecting case; the rest will be the same]
		- M+=m_i
		- S+=m_i*x_i (local static moments are zero (centroid)
		- get inertia tensor of subBody in world coordinates, by rotating the principal (local) tensor against subBody->se3->orientation; then translate it to world origin (parallel axes theorem), then I+=I_i_world
	-# clumpPos=S/M
	-# translate aggregate's inertia tensor; parallel axes on I (R=clumpPos): I^c_jk=I'_jk-M*(delta_jk R.R - R_j*R_k) [http://en.wikipedia.org/wiki/Moments_of_inertia#Parallel_axes_theorem]
	-# eigen decomposition of I, get principal inertia and rotation matrix of the clump
	-# se3->orientation=quaternion(rotation_matrix); se3->position=clumpPos
	-#	update subSe3s

	@todo I \em think that the order of transformation of inertia is:
		- from local to global: first rotate, then translate;
		- from global to local: first translate, then rotate,
	since rotation must be done with origin at the centroid... This needs to be verified, though.
	@todo all the rest of this routine needs to be verified!
	@todo implement the loop for intersecting bodies (may cut'n'paste from slum code, but that will work for spheres only!)

	@note User is responsible for calling this function when appropriate (after adding/removing bodies and before any subsequent simulation). This function can be rather slow by virtue of numerical integration.
	@note subBodie's velocities are not taken into account. This means that clump will be at still after being created, even if its composing particles did have some velocities. If this is concern for someone, it needs to be completed in the code below. See Clump::moveMembers for complementary issue.
	@todo Needs to be tested for physical correctness
	@param intersecting if true, evaluate mass and inertia numerically; otherwise, use analytical methods (parallel axes theorem) which disregard any intersections, but are much faster. */
void Clump::updateProperties(bool intersecting){
	LOG_DEBUG("Updating clump #"<<getId()<<" parameters");
	assert(members.size()>0);

	/* quantities suffixed by
		g: global (world) coordinates
		s: local subBody's coordinates
		c: local clump coordinates */
	double M=0; // mass
	Vector3r Sg(0,0,0); // static moment
	Matrix3r Ig(Matrix3r::Zero()), Ic(Matrix3r::Zero()); // tensors of inertia; is upper triangular, zeros instead of symmetric elements
	//Se3r& mySe3(physicalParameters->se3);
	//const shared_ptr<RigidBodyParameters>& clumpRBP(YADE_PTR_CAST<RigidBodyParameters>(physicalParameters));

	if(members.size()==1){
		LOG_DEBUG("Clump of size one will be treated specially.")
		memberMap::iterator I=members.begin();
		shared_ptr<Body> subBody=Body::byId(I->first);
		//const shared_ptr<RigidBodyParameters>& subRBP(YADE_PTR_CAST<RigidBodyParameters>(subBody->physicalParameters));
		State* subState=subBody->state.get();
		// se3 of the clump as whole is the same as the member's se3
		state->pos=subState->pos;
		state->ori=subState->ori;
		// relative member's se3 is identity
		I->second.position=Vector3r::Zero(); I->second.orientation=Quaternionr::Identity();
		state->inertia=subState->inertia;	
		state->mass=subState->mass;
		state->vel=Vector3r::Zero();
		state->angVel=Vector3r::Zero();
		return;
	}

	if(intersecting){
		LOG_WARN("Self-intersecting clumps not yet implemented, intersections will be ignored.");
		intersecting=false;}

	// begin non-intersecting loop here
	if(!intersecting){
		for(memberMap::iterator I=members.begin(); I!=members.end(); I++){
			// now, I->first is Body::id_t, I->second is Se3r of that body
			shared_ptr<Body> subBody=Body::byId(I->first);
			//const shared_ptr<RigidBodyParameters>& subRBP(YADE_PTR_CAST<RigidBodyParameters>(subBody->physicalParameters));
			State* subState=subBody->state.get();
			M+=subState->mass;
			Sg+=subState->mass*subState->pos;
			// transform from local to global coords
			// FIXME: verify this!
			Quaternionr subState_orientation_conjugate=subState->ori.conjugate();
			Matrix3r Imatrix=Matrix3r::Zero(); for(int i=0;i<3;i++) Imatrix(i,i)=subState->inertia[i]; // eigen: Imatrix.diagonal=subState->inertia;
			// TRWM3MAT(Imatrix); TRWM3QUAT(subRBP_orientation_conjugate);
			Ig+=Clump::inertiaTensorTranslate(Clump::inertiaTensorRotate(Imatrix,subState_orientation_conjugate),subState->mass,-1.*subState->pos);

			//TRWM3MAT(Clump::inertiaTensorRotate(Matrix3r(subRBP->inertia),subRBP_orientation_conjugate));
		}
	}
	TRVAR1(M);
	TRWM3MAT(Ig);
	TRWM3VEC(Sg);
	if(M==0){ state->mass=0; state->inertia=Vector3r(0,0,0); isDynamic=false; return; }

	state->pos=Sg/M; // clump's centroid
	// this will calculate translation only, since rotation is zero
	Matrix3r Ic_orientG=Clump::inertiaTensorTranslate(Ig, -M /* negative mass means towards centroid */, state->pos); // inertia at clump's centroid but with world orientation
	TRWM3MAT(Ic_orientG);

	Matrix3r R_g2c(Matrix3r::Zero()); //rotation matrix
	Ic_orientG(1,0)=Ic_orientG(0,1); Ic_orientG(2,0)=Ic_orientG(0,2); Ic_orientG(2,1)=Ic_orientG(1,2); // symmetrize
	//TRWM3MAT(Ic_orientG);
	matrixEigenDecomposition(Ic_orientG,R_g2c,Ic);
	/*! @bug eigendecomposition might be wrong. see http://article.gmane.org/gmane.science.physics.yade.devel/99 for message. It is worked around below, however.
	*/
	// has NaNs for identity matrix!
	TRWM3MAT(R_g2c);

	// these two should give the same result!
	//TRWM3MAT(Ic);
	//TRWM3MAT(Clump::inertiaTensorRotate(Ic_orientG,R_g2c));

	// set quaternion from rotation matrix
	state->ori=Quaternionr(R_g2c); state->ori.normalize();
	// now Ic is diagonal
	state->inertia=Vector3r(Ic(0,0),Ic(1,1),Ic(2,2));
	state->mass=M;


	// this block will be removed once EigenDecomposition works for diagonal matrices
	#if 1
		if(isnan(R_g2c(0,0))||isnan(R_g2c(0,1))||isnan(R_g2c(0,2))||isnan(R_g2c(1,0))||isnan(R_g2c(1,1))||isnan(R_g2c(1,2))||isnan(R_g2c(2,0))||isnan(R_g2c(2,1))||isnan(R_g2c(2,2))){
			LOG_FATAL("EigenDecomposition gave some NaNs, we will use imaginary values for clump inertia and orientation instead. I thought this may happen only for 1-member clumps which are now treated specially. Something is broken!");
			//FIXME: since EigenDecomposition is broken, use inertia of the first body instead;
			//!!!!! note that this is HIGHLY incorrect for all non-single clumps !!!!!
			memberMap::iterator I=members.begin();
			shared_ptr<Body> subBody=Body::byId(I->first);
			state->inertia=subBody->state->inertia*10.; // 10 is arbitrary; just to have inertia of clump bigger
			// orientation of the clump is broken as well, since is result of EigenDecomposition as well (rotation matrix)
			state->ori=Quaternionr::Identity();
		}
	#endif
	TRWM3VEC(state->inertia);

	// TODO: these might be calculated from members... but complicated... - someone needs that?!
	state->vel=state->angVel=Vector3r::Zero();

	// update subBodySe3s; subtract clump orientation (=apply its inverse first) to subBody's orientation
	// Conjugate is equivalent to Inverse for normalized quaternions
	for(memberMap::iterator I=members.begin(); I!=members.end(); I++){
		// now, I->first is Body::id_t, I->second is Se3r of that body
		shared_ptr<Body> subBody=Body::byId(I->first);
		//const shared_ptr<RigidBodyParameters>& subRBP(YADE_PTR_CAST<RigidBodyParameters>(subBody->physicalParameters));
		State* subState=subBody->state.get();
		I->second.orientation=state->ori.conjugate()*subState->ori;
		I->second.position=state->ori.conjugate()*(subState->pos-state->pos);
	}

}

/*! @brief Recalculates inertia tensor of a body after translation away from (default) or towards its centroid.
 *
 * @param I inertia tensor in the original coordinates; it is assumed to be upper-triangular (elements below the diagonal are ignored).
 * @param m mass of the body; if positive, translation is away from the centroid; if negative, towards centroid.
 * @param off offset of the new origin from the original origin
 * @return inertia tensor in the new coordinate system; the matrix is symmetric.
 */
Matrix3r Clump::inertiaTensorTranslate(const Matrix3r& I,const Real m, const Vector3r& off){
	Real ooff=off.dot(off);
	Matrix3r I2=I;
	//TRWM3VEC(off); TRVAR2(ooff,m); TRWM3MAT(I);
	// translation away from centroid
	/* I^c_jk=I'_jk-M*(delta_jk R.R - R_j*R_k) [http://en.wikipedia.org/wiki/Moments_of_inertia#Parallel_axes_theorem] */
	I2+=m*matrixFromElements(/* dIxx */ ooff-off[0]*off[0], /* dIxy */ -off[0]*off[1], /* dIxz */ -off[0]*off[2],
		/* sym */ 0., /* dIyy */ ooff-off[1]*off[1], /* dIyz */ -off[1]*off[2],
		/* sym */ 0., /* sym */ 0., /* dIzz */ ooff-off[2]*off[2]);
	I2(1,0)=I2(0,1); I2(2,0)=I2(0,2); I2(2,1)=I2(1,2);
	//TRWM3MAT(I2);
	return I2;
}

/*! @brief Recalculate body's inertia tensor in rotated coordinates.
 *
 * @param I inertia tensor in old coordinates
 * @param T rotation matrix from old to new coordinates
 * @return inertia tensor in new coordinates
 */
Matrix3r Clump::inertiaTensorRotate(const Matrix3r& I,const Matrix3r& T){
	/* [http://www.kwon3d.com/theory/moi/triten.html] */
	//TRWM3MAT(I); TRWM3MAT(T);
	return T.transpose()*I*T;
}

/*! @brief Recalculate body's inertia tensor in rotated coordinates.
 *
 * @param I inertia tensor in old coordinates
 * @param rot quaternion that describes rotation from old to new coordinates
 * @return inertia tensor in new coordinates
 */
Matrix3r Clump::inertiaTensorRotate(const Matrix3r& I, const Quaternionr& rot){
	Matrix3r T=rot.toRotationMatrix();
	return inertiaTensorRotate(I,T);
}



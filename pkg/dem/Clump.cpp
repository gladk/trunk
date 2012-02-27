// (c) 2007-2010 Vaclav Smilauer <eudoxos@arcig.cz> 

#include"Clump.hpp"
#include<algorithm>
#include<yade/core/Scene.hpp>
#include<yade/core/BodyContainer.hpp>
#include<yade/core/State.hpp>

YADE_PLUGIN((Clump));
CREATE_LOGGER(Clump);

void Clump::add(const shared_ptr<Body>& clumpBody, const shared_ptr<Body>& subBody){
	Body::id_t subId=subBody->getId();
	if(subBody->clumpId!=Body::ID_NONE) throw std::invalid_argument(("Body #"+lexical_cast<string>(subId)+" is already in clump #"+lexical_cast<string>(subBody->clumpId)).c_str());
	const shared_ptr<Clump> clump=YADE_PTR_CAST<Clump>(clumpBody->shape);
	if(clump->members.count(subId)!=0) throw std::invalid_argument(("Body #"+lexical_cast<string>(subId)+" is already part of this clump #"+lexical_cast<string>(clumpBody->id)).c_str());
	
	clump->members[subId]=Se3r(); // meaningful values will be put in by Clump::updateProperties
	subBody->clumpId=clumpBody->id;
	clumpBody->clumpId=clumpBody->id; // just to make sure
	clumpBody->setBounded(false); // disallow collisions with the clump itself
	LOG_DEBUG("Added body #"<<subBody->id<<" to clump #"<<clumpBody->id);
}

void Clump::del(const shared_ptr<Body>& clumpBody, const shared_ptr<Body>& subBody){
	// erase the subBody; removing body that is not part of the clump throws
	const shared_ptr<Clump> clump=YADE_PTR_CAST<Clump>(clumpBody->shape);
	if(clump->members.erase(subBody->id)!=1) throw std::invalid_argument(("Body #"+lexical_cast<string>(subBody->id)+" not part of clump #"+lexical_cast<string>(clumpBody->id)+"; not removing.").c_str());
	subBody->clumpId=Body::ID_NONE;
	LOG_DEBUG("Removed body #"<<subBody->id<<" from clump #"<<clumpBody->id);
}

void Clump::moveMembers(const shared_ptr<Body>& clumpBody, Scene* scene){
	const shared_ptr<Clump>& clump=YADE_PTR_CAST<Clump>(clumpBody->shape);
	const shared_ptr<State>& clumpState=clumpBody->state;

	FOREACH(MemberMap::value_type& B, clump->members){
		// B->first is Body::id_t, B->second is local Se3r of that body in the clump
		State* subState=Body::byId(B.first,scene)->state.get();
		Vector3r& subPos(B.second.position); Quaternionr& subOri(B.second.orientation);
		//LOG_TRACE("Old #"<<I->first<<"position: "<<subRBP->se3.position);

		// position update
		subState->pos=clumpState->pos+clumpState->ori*subPos;
		subState->ori=clumpState->ori*subOri;
		// velocity update
		subState->vel=clumpState->vel+clumpState->angVel.cross(subState->pos-clumpState->pos);
		subState->angVel=clumpState->angVel;

		//LOG_TRACE("New #"<<I->first<<"position: "<<subRBP->se3.position);
		//LOG_TRACE("Clump #"<<getId()<<" moved #"<<I->first<<".");
	}
}

/*! Clump's se3 will be updated (origin at centroid and axes coincident with principal inertia axes) and subSe3 modified in such a way that members positions in world coordinates will not change.

	TODO: numerical integration of inertia based on regular space sampling with relative tolerance WRT minimum sphere radius

	Note: velocities and angularVelocities of constituents are zeroed.

	OLD DOCS (will be cleaned up):

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

*/

void Clump::updateProperties(const shared_ptr<Body>& clumpBody, bool intersecting){
	LOG_DEBUG("Updating clump #"<<clumpBody->id<<" parameters");
	const shared_ptr<State> state(clumpBody->state);
	const shared_ptr<Clump> clump(YADE_PTR_CAST<Clump>(clumpBody->shape));
	
	if(clump->members.empty()){ throw std::runtime_error("Clump::updateProperties: clump has zero members."); }
	// trivial case
	if(clump->members.size()==1){
		LOG_DEBUG("Clump of size one will be treated specially.")
		MemberMap::iterator I=clump->members.begin();
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

	/* quantities suffixed by
		g: global (world) coordinates
		s: local subBody's coordinates
		c: local clump coordinates
	*/
	double M=0; // mass
	Vector3r Sg(0,0,0); // static moment, for getting clump's centroid
	Matrix3r Ig(Matrix3r::Zero()), Ic(Matrix3r::Zero()); // tensors of inertia; is upper triangular, zeros instead of symmetric elements

	if(intersecting){
		LOG_WARN("Self-intersecting clumps not yet implemented, intersections will be ignored.");
		intersecting=false;
	}

	// begin non-intersecting loop here
	if(!intersecting){
		FOREACH(MemberMap::value_type& I, clump->members){
			// I.first is Body::id_t, I.second is Se3r of that body
			shared_ptr<Body> subBody=Body::byId(I.first);
			State* subState=subBody->state.get();
			M+=subState->mass;
			Sg+=subState->mass*subState->pos;
			// transform from local to global coords
			Quaternionr subState_ori_conjugate=subState->ori.conjugate();
			Matrix3r Imatrix=Matrix3r::Zero(); Imatrix.diagonal()=subState->inertia; 
			// TRWM3MAT(Imatrix); TRWM3QUAT(subRBP_orientation_conjugate);
			Ig+=Clump::inertiaTensorTranslate(Clump::inertiaTensorRotate(Imatrix,subState_ori_conjugate),subState->mass,-1.*subState->pos);
			//TRWM3MAT(Clump::inertiaTensorRotate(Matrix3r(subRBP->inertia),subRBP_orientation_conjugate));
		}
	}
	TRVAR1(M); TRWM3MAT(Ig); TRWM3VEC(Sg);
	assert(M>0);

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

	// set quaternion from rotation matrix
	state->ori=Quaternionr(R_g2c); state->ori.normalize();
	// now Ic is diagonal
	state->inertia=Ic.diagonal();
	state->mass=M;


	// this block will be removed once EigenDecomposition works for diagonal matrices
	#if 1
		if(isnan(R_g2c(0,0))||isnan(R_g2c(0,1))||isnan(R_g2c(0,2))||isnan(R_g2c(1,0))||isnan(R_g2c(1,1))||isnan(R_g2c(1,2))||isnan(R_g2c(2,0))||isnan(R_g2c(2,1))||isnan(R_g2c(2,2))){
			throw std::logic_error("Clump::updateProperties: NaNs in eigen-decomposition of inertia matrix?!");
		}
	#endif
	TRWM3VEC(state->inertia);

	// TODO: these might be calculated from members... but complicated... - someone needs that?!
	state->vel=state->angVel=Vector3r::Zero();

	clumpBody->setAspherical(state->inertia[0]!=state->inertia[1] || state->inertia[0]!=state->inertia[2]);

	// update subBodySe3s; subtract clump orientation (=apply its inverse first) to subBody's orientation
	FOREACH(MemberMap::value_type& I, clump->members){
		// now, I->first is Body::id_t, I->second is Se3r of that body
		shared_ptr<Body> subBody=Body::byId(I.first);
		//const shared_ptr<RigidBodyParameters>& subRBP(YADE_PTR_CAST<RigidBodyParameters>(subBody->physicalParameters));
		State* subState=subBody->state.get();
		I.second.orientation=state->ori.conjugate()*subState->ori;
		I.second.position=state->ori.conjugate()*(subState->pos-state->pos);
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
	// TODO: replace by nicer eigen code

	// translation away from centroid
	/* I^c_jk=I'_jk-M*(delta_jk R.R - R_j*R_k) [http://en.wikipedia.org/wiki/Moments_of_inertia#Parallel_axes_theorem] */
	Matrix3r dI; dI<</* dIxx */ ooff-off[0]*off[0], /* dIxy */ -off[0]*off[1], /* dIxz */ -off[0]*off[2],
		/* sym */ 0., /* dIyy */ ooff-off[1]*off[1], /* dIyz */ -off[1]*off[2],
		/* sym */ 0., /* sym */ 0., /* dIzz */ ooff-off[2]*off[2];
	I2+=m*dI;
	I2(1,0)=I2(0,1); I2(2,0)=I2(0,2); I2(2,1)=I2(1,2); // symmetrize
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



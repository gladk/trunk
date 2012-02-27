/* LucScholtes2010  */

#include"CohesiveFrictionalPM.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/core/Omega.hpp>

YADE_PLUGIN((CFpmMat)(CFpmState)(CFpmPhys)(Ip2_CFpmMat_CFpmMat_CFpmPhys)(Law2_ScGeom_CFpmPhys_CohesiveFrictionalPM));

/********************** Law2_ScGeom_CFpmPhys_CohesiveFrictionalPM ****************************/
CREATE_LOGGER(Law2_ScGeom_CFpmPhys_CohesiveFrictionalPM);

void Law2_ScGeom_CFpmPhys_CohesiveFrictionalPM::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	ScGeom* geom = static_cast<ScGeom*>(ig.get()); 
	CFpmPhys* phys = static_cast<CFpmPhys*>(ip.get());
	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	Body* b1 = Body::byId(id1,scene).get();
	Body* b2 = Body::byId(id2,scene).get();
	
	Real displN = geom->penetrationDepth; // NOTE: the sign for penetrationdepth is different from ScGeom and Dem3DofGeom: geom->penetrationDepth>0 when spheres interpenetrate
	Real Dtensile=phys->FnMax/phys->kn;
	Real Dsoftening = phys->strengthSoftening*Dtensile; 
	
	/*to set the equilibrium distance between all cohesive elements when they first meet -> allows one to work with initial stress-free assembly*/
	if ( contact->isFresh(scene) ) { phys->initD = displN; phys->normalForce = Vector3r::Zero(); phys->shearForce = Vector3r::Zero();}
	Real D = displN - phys->initD; // interparticular distance is computed depending on the equilibrium distance

	/* Determination of interaction */
	if (D < 0){ //spheres do not touch 
	  if (!phys->isCohesive){ scene->interactions->requestErase(contact->getId1(),contact->getId2()); return; } // destroy the interaction before calculation
	  if ((phys->isCohesive) && (abs(D) > (Dtensile + Dsoftening))) { // spheres are bonded and the interacting distance is greater than the one allowed ny the defined cohesion
	    phys->isCohesive=false; 
	    // update body state with the number of broken bonds
	    CFpmState* st1=dynamic_cast<CFpmState*>(b1->state.get());
	    CFpmState* st2=dynamic_cast<CFpmState*>(b2->state.get());
	    st1->numBrokenCohesive+=1;
	    st2->numBrokenCohesive+=1;
	    //// the same thing but from ConcretePM
	    //const shared_ptr<Body>& body1=Body::byId(contact->getId1(),scene), body2=Body::byId(contact->getId2(),scene); assert(body1); assert(body2);
	    //const shared_ptr<CFpmState>& st1=YADE_PTR_CAST<CFpmState>(body1->state), st2=YADE_PTR_CAST<CFpmState>(body2->state);
	    //{ boost::mutex::scoped_lock lock(st1->updateMutex); st1->numBrokenCohesive+=1; }
	    //{ boost::mutex::scoped_lock lock(st2->updateMutex); st2->numBrokenCohesive+=1; }
	    // end of update
	    scene->interactions->requestErase(contact->getId1(),contact->getId2()); return;
	  }
	}	  
	
	/*NormalForce*/
	Real Fn=0, Dsoft=0;
		
	if ((D < 0) && (abs(D) > Dtensile)) { //to take into account strength softening
	  Dsoft = D+Dtensile; // Dsoft<0 for a negative value of Fn (attractive force)
	  Fn = -(phys->FnMax+(phys->kn/phys->strengthSoftening)*Dsoft); // computes FnMax - FnSoftening
	}
	else {
	  Fn = phys->kn*D;
	}
	phys->normalForce = Fn*geom->normal;  // NOTE normal is position2-position1 - It is directed from particle1 to particle2
	        
	/*ShearForce*/
	Vector3r& shearForce = phys->shearForce; 
		
	// using scGeom function rotateAndGetShear	
	State* st1 = Body::byId(id1,scene)->state.get();
	State* st2 = Body::byId(id2,scene)->state.get();

	geom->rotate(phys->shearForce);
	const Vector3r& dus = geom->shearIncrement();
	
	//Linear elasticity giving "trial" shear force
	shearForce -= phys->ks*dus;
	// needed for the next timestep
	phys->prevNormal = geom->normal;
		
	/* Morh-Coulomb criterion */
	Real maxFs = phys->FsMax + Fn*phys->tanFrictionAngle;

	if (shearForce.squaredNorm() > maxFs*maxFs){ 
	  shearForce*=maxFs/shearForce.norm(); // to fix the shear force to its yielding value
	}
	
	/* Apply forces */
	Vector3r f = phys->normalForce + shearForce;
	// these lines to adapt to periodic boundary conditions (NOTE applyForceAtContactPoint computes torque induced by normal and shear force too)
	if (!scene->isPeriodic)  
	applyForceAtContactPoint(f , geom->contactPoint , id2, st2->se3.position, id1, st1->se3.position);
	else { // in scg we do not wrap particles positions, hence "applyForceAtContactPoint" cannot be used when scene is periodic
		scene->forces.addForce(id1,-f);
		scene->forces.addForce(id2,f);
		scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
		scene->forces.addTorque(id2,(geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
	}
	
	/* Moment Rotation Law */
	// NOTE this part could probably be computed in ScGeom to avoid copy/paste multiplication !!!
	Quaternionr delta( b1->state->ori * phys->initialOrientation1.conjugate() *phys->initialOrientation2 * b2->state->ori.conjugate()); delta.normalize(); //relative orientation
	AngleAxisr aa(delta); // axis of rotation - this is the Moment direction UNIT vector; angle represents the power of resistant ELASTIC moment
	if(aa.angle() > Mathr::PI) aa.angle() -= Mathr::TWO_PI; // angle is between 0 and 2*pi, but should be between -pi and pi 
	  
	phys->cumulativeRotation = aa.angle();
	  
	//Find angle*axis. That's all.  But first find angle about contact normal. Result is scalar. Axis is contact normal.
	Real angle_twist(aa.angle() * aa.axis().dot(geom->normal) ); //rotation about normal
	Vector3r axis_twist(angle_twist * geom->normal);
	Vector3r moment_twist(axis_twist * phys->kr);
	  
	Vector3r axis_bending(aa.angle()*aa.axis() - axis_twist); //total rotation minus rotation about normal
	Vector3r moment_bending(axis_bending * phys->kr);
	Vector3r moment = moment_twist + moment_bending;
	
	Real MomentMax = phys->maxBend*std::fabs(phys->normalForce.norm());
	Real scalarMoment = moment.norm();

	/*Plastic moment */
	if(scalarMoment > MomentMax) 
	{
	  Real ratio = 0;
	  ratio *= MomentMax/scalarMoment; // to fix the moment to its yielding value
	  moment *= ratio;		
	  moment_twist *=  ratio;
	  moment_bending *= ratio;
	}

	phys->moment_twist = moment_twist;
	phys->moment_bending = moment_bending;
	  
	scene->forces.addTorque(id1,-moment);
	scene->forces.addTorque(id2, moment);

}

CREATE_LOGGER(Ip2_CFpmMat_CFpmMat_CFpmPhys);

void Ip2_CFpmMat_CFpmMat_CFpmPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction){
	
	/* avoid any updates if the interaction already exists */
	if(interaction->phys) return; 

	ScGeom* geom=dynamic_cast<ScGeom*>(interaction->geom.get());
	assert(geom);

	const shared_ptr<CFpmMat>& yade1 = YADE_PTR_CAST<CFpmMat>(b1);
	const shared_ptr<CFpmMat>& yade2 = YADE_PTR_CAST<CFpmMat>(b2);
			
	shared_ptr<CFpmPhys> contactPhysics(new CFpmPhys()); 
	
	/* From interaction physics */
	Real E1 	= yade1->young;
	Real E2 	= yade2->young;
	Real V1 	= yade1->poisson;
	Real V2 	= yade2->poisson;
	Real f1 	= yade1->frictionAngle;
	Real f2 	= yade2->frictionAngle;

	/* From interaction geometry */
	Real R1= geom->radius1;
	Real R2= geom->radius2;
	Real rMean = 0.5*(R1+R2); 
	Real crossSection = Mathr::PI*pow(min(R1,R2),2); 
	
	/* calculate stiffness */
	Real kNormal=0, kShear=0, kRotate=0;
	
	if(useAlphaBeta == true){
	  kNormal	= 2.*E1*R1*E2*R2/(E1*R1+E2*R2); // harmonic average of two stiffnesses
	  kShear 	= Alpha*kNormal;
	  kRotate	= Beta*kShear*rMean*rMean; 
	}
	else {
	  kNormal	= 2.*E1*R1*E2*R2/(E1*R1+E2*R2); // harmonic average of two stiffnesses
	  kShear	= 2.*E1*R1*V1*E2*R2*V2/(E1*R1*V1+E2*R2*V2); // harmonic average of two stiffnesses with ks=V*kn for each sphere
	  kRotate	= 0.;
	}

	/* Pass values calculated from above to CFpmPhys */
	contactPhysics->kn = kNormal;
	contactPhysics->ks = kShear;
	contactPhysics->kr = kRotate;
	contactPhysics->frictionAngle		= std::min(f1,f2); 
	contactPhysics->tanFrictionAngle	= std::tan(contactPhysics->frictionAngle); 
	contactPhysics->maxBend 		= eta*rMean;
	contactPhysics->prevNormal 		= geom->normal;
	contactPhysics->initialOrientation1	= Body::byId(interaction->getId1())->state->ori;
	contactPhysics->initialOrientation2	= Body::byId(interaction->getId2())->state->ori;
	
	///to set if the contact is cohesive or not
	if ( (scene->iter < cohesiveTresholdIteration) && ((tensileStrength>0) || (cohesion>0)) && (yade1->type == yade2->type)){ contactPhysics->isCohesive=true; }
	
	if ( contactPhysics->isCohesive ) {
	  contactPhysics->FnMax = tensileStrength*crossSection;
	  contactPhysics->strengthSoftening = strengthSoftening;
	  contactPhysics->FsMax = cohesion*crossSection;
	}
	
	interaction->phys = contactPhysics;
}

CFpmPhys::~CFpmPhys(){}

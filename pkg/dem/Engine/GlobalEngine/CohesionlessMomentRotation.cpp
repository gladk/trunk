/* CWBoon@2010  booncw@hotmail.com */


#include"CohesionlessMomentRotation.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/core/Omega.hpp>

/* Contact law have been verified with Plassiard et al. (2009) : A spherical discrete element model: calibration procedure and incremental response */


YADE_PLUGIN((Law2_SCG_MomentPhys_CohesionlessMomentRotation)(Ip2_MomentMat_MomentMat_MomentPhys)(MomentPhys));



/********************** Law2_SCG_MomentPhys_CohesionlessMomentRotation ****************************/
CREATE_LOGGER(Law2_SCG_MomentPhys_CohesionlessMomentRotation);

void Law2_SCG_MomentPhys_CohesionlessMomentRotation::go(shared_ptr<InteractionGeometry>& ig, shared_ptr<InteractionPhysics>& ip, Interaction* contact, Scene* rootBody){
	
	ScGeom* geom = static_cast<ScGeom*>(ig.get()); //InteractionGeometry
	MomentPhys* phys = static_cast<MomentPhys*>(ip.get()); //InteractionPhysics
	int id1 = contact->getId1(); //Id of body1
        int id2 = contact->getId2(); //Id of body2
	shared_ptr<BodyContainer>& bodies = rootBody->bodies;
	Body* b1 = ( *bodies ) [id1].get();
	Body* b2 = ( *bodies ) [id2].get();
	const Real& dt = scene->dt; //Get TimeStep


	/*NormalForce */
	Real displN=geom->penetrationDepth;  //get penetrationDepth from SCG
	if (displN<0){rootBody->interactions->requestErase(contact->getId1(),contact->getId2()); return;}  //! NOTE: the sign for penetrationdepth is different from SCG and Dem3DofGeom
	Real Fn =phys->kn*displN;   //scalar force
	phys->normalForce= Fn*geom->normal; //vector force NOTE: normal is position2-position1.  Therefore it is directed from particle1 to particle2
	

	Vector3r axis;
        Real angle;

	
	///////////////////////////////////////////////////////////////////////  SHEAR FORCE ///////////////////////////////////////////////
	/*ShearForce*/
	Vector3r& shearForce = phys->shearForce; //alias for phys->shearForce


	/// Here is the code with approximated rotations ///
	axis	 = phys->prevNormal.cross(geom->normal);  //axis of rotation between previous and current
	shearForce -= shearForce.cross(axis); //minus shearForce that is perpendicular (outward) to contact (not tangent)
	angle 	= dt*0.5*geom->normal.dot(Body::byId(id1)->state->angVel+Body::byId(id2)->state->angVel); 
	axis 	= angle*geom->normal;
	shearForce -= shearForce.cross(axis); //minus shearForce that is perpendicular outward to contact (not tangent)

	/// Here is the code with exact rotations 		 ///
	//Quaternionr q;
	//1st imagine we have to rotate shearForce ON the x-y plane as a result of contact point rotation
	//axis = phys->prevNormal.Cross(geom->normal); // axis of rotation of contact point
	//angle = Mathr::ACos(geom->normal.Dot(phys->prevNormal)); //angle of rotation of contact point
	//q=Quaternion(AngleAxisr(angleAxis)); //quaternion for rotation of contact point
	// shearForce = q*shearForce; //rotate shearForce from previous contact orientation to current contact orientation
	
	//Then imagine we have to find a direction by rotating  THROUGH the x-y plane, ABOUT the contact normal (i.e. axis).  The direction should be along the resultant angular velocity.
	//angle	= dt*0.5*geom->normal.Dot(Body::byId(id1)->state->angVel + Body::byId(id1)->state->angVel); //angle of rotation about the contact normal
	//axis	= geom->normal;
	//q.FromAxisAngle(axis,angle);
	//shearForce	= q*shearForce;
	/// 							 ///

	
	Vector3r x = geom->contactPoint;
	Vector3r c1x	= (x - b1->state->pos); //distance of contact point from centre1
	Vector3r c2x	= (x - b2 ->state->pos); //distance of contact ponit from centre2	
	/// The following definition of c1x and c2x is to avoid "granular ratcheting" 
	///  (see F. ALONSO-MARROQUIN, R. GARCIA-ROJO, H.J. HERRMANN, 
	///   "Micro-mechanical investigation of granular ratcheting, in Cyclic Behaviour of Soils and Liquefaction Phenomena",
	///   ed. T. Triantafyllidis (Balklema, London, 2004), p. 3-10 - and a lot more papers from the same authors)
	Vector3r _c1x_(0,0,0);
	Vector3r _c2x_(0,0,0);
	if (preventGranularRatcheting){
		_c1x_	= geom->radius1*geom->normal; //vector radius1
		_c2x_	= -geom->radius2*geom->normal; //vector radius2
	}else{
		_c1x_ = c1x;
		_c2x_ = c2x;
	}

	Vector3r relativeVelocity	= (b2->state->vel + b2->state->angVel.cross(_c2x_)) - (b1->state->vel + b1->state->angVel.cross(_c1x_));
	Vector3r shearVelocity=relativeVelocity- geom->normal.dot(relativeVelocity)*geom->normal; //1st find the scalar value of normal velocity.  Then change it to a vector.  Then minus.
	Vector3r shearDisplacement	= shearVelocity*dt; //in Dem3DofGeom displacementT() is cumulative
	shearForce -= phys->ks*shearDisplacement; //shearForce is cumulative here because shearDisplacement is not cumulative.  See CundallStrack.hpp, it is different!
      
	
	////////////////////////////////////////////////////////////   SHEAR FORCE END   ///////////////////////////////////////////////////////////////////////////////


	/* Coulomb's Law */
	Real maxFs = 0;
        Real Fs = shearForce.norm(); //scalar
	maxFs = std::max((Real) 0, Fn*phys->tanFrictionAngle); //Fn is scalar
      	if ( Fs  > maxFs )
	{
	shearForce *= maxFs/Fs;
	//if ( Fn<0 )  phys->normalForce = Vector3r::Zero();
	}

	phys->shearForce = shearForce; //for recording purposes

	/* Apply forces */
	Vector3r f = phys->normalForce + shearForce;
	rootBody->forces.addForce (id1,-f);
	rootBody->forces.addForce (id2, f);
	rootBody->forces.addTorque(id1,-(c1x).cross(f)); //about axis perpendicular to normal and force
	rootBody->forces.addTorque(id2, c2x.cross(f));

	/* Moment Rotation Law */
	Quaternionr delta( b1->state->ori * phys->initialOrientation1.conjugate() *phys->initialOrientation2 * b2->state->ori.conjugate()); //relative orientation
	AngleAxisr aa(angleAxisFromQuat(delta)); // axis of rotation - this is the Moment direction UNIT vector; angle represents the power of resistant ELASTIC moment
	if(aa.angle() > Mathr::PI) aa.angle() -= Mathr::TWO_PI; // angle is between 0 and 2*pi, but should be between -pi and pi 
	
	phys->cumulativeRotation = aa.angle();
	//Real elasticMoment = phys->kr * std::abs(aa.angle()); // positive value (*)	
	//Vector3r moment = aa.axis() * elasticMoment * (aa.angle()<0.0?-1.0:1.0); // restore sign. (*)
	

	//Find angle*axis. That's all.  But first find angle about contact normal. Result is scalar. Axis is contact normal.
	Real angle_twist(aa.angle() * aa.axis().dot(geom->normal) ); //rotation about normal
	Vector3r axis_twist(angle_twist * geom->normal);
	Vector3r moment_twist(axis_twist * phys->kr);
	
	Vector3r axis_bending(aa.angle()*aa.axis() - axis_twist); //total rotation minus rotation about normal
	Vector3r moment_bending(axis_bending * phys->kr);
	Vector3r moment = moment_twist + moment_bending;

	Real radiusm = 0.5*(geom->radius1 + geom->radius2); //radius is used here
	Real MomentMax = phys->Eta * radiusm* std::fabs(Fn);
	Real scalarMoment = moment.norm();


	if(scalarMoment > MomentMax) /*Plastic moment */
	{
		Real ratio=0;
		ratio *= MomentMax/scalarMoment;
		moment *= ratio;		
		moment_twist *=  ratio;
		moment_bending *= ratio;
	}

	phys->moment_twist = moment_twist;
	phys->moment_bending = moment_bending;
	
	rootBody->forces.addTorque(id1,-moment);
	rootBody->forces.addTorque(id2, moment);	
	/// Moment law	END ///

	phys->prevNormal = geom->normal;
	
}



CREATE_LOGGER(Ip2_MomentMat_MomentMat_MomentPhys);

void Ip2_MomentMat_MomentMat_MomentPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction){
	
	if(interaction->interactionPhysics) return; 

	ScGeom* scg=dynamic_cast<ScGeom*>(interaction->interactionGeometry.get());
			
	assert(scg);

	const shared_ptr<MomentMat>& sdec1 = YADE_PTR_CAST<MomentMat>(b1);
	const shared_ptr<MomentMat>& sdec2 = YADE_PTR_CAST<MomentMat>(b2);
			
	shared_ptr<MomentPhys> contactPhysics(new MomentPhys()); 

	/* From interaction physics */
	Real Ea 	= sdec1->young;
	Real Eb 	= sdec2->young;
	Real Va 	= sdec1->poisson;
	Real Vb 	= sdec2->poisson;
	Real fa 	= sdec1->frictionAngle;
	Real fb 	= sdec2->frictionAngle;
	Real etaA	= sdec1->eta;
	Real etaB 	= sdec2->eta;

	/* From interaction geometry */
	Real Da= scg->radius1;
	Real Db= scg->radius2;  			
	
	/* calculate stiffness */
	Real Kn=0;
	Real Ks=0;
	Real Kr=0;	
	
	Real radiusm = (Da + Db)/2; 

	if(userInputStiffness == true && useAlphaBeta == true){
	Kn = 0.5*Knormal*radiusm;
	Ks = Alpha*Knormal;
	Ks = 0.5*Ks*radiusm;
	Kr = Beta*Ks*radiusm*radiusm;  
	}
	else if (userInputStiffness == true && useAlphaBeta == false){//if userInputStiffness = true
	Kn = Knormal;
	Ks = Kshear;
	Kr = Krotate;
	}
	else{ 
	Kn = 2*Ea*Da*Eb*Db/(Ea*Da+Eb*Db);//harmonic average of two stiffnesses
	Ks = 2*Ea*Da*Va*Eb*Db*Vb/(Ea*Da*Va+Eb*Db*Va);//harmonic average of two stiffnesses with ks=V*kn for each sphere
	Kr=0;	
	}
	
	/* Pass values calculated from above to CSPhys */
	contactPhysics->kn = Kn;
	contactPhysics->ks = Ks;
	contactPhysics->kr = Kr;
	contactPhysics->Eta = std::min(etaA,etaB);
	contactPhysics->frictionAngle		= std::min(fa,fb); 
	contactPhysics->tanFrictionAngle	= std::tan(contactPhysics->frictionAngle); 
	contactPhysics->initialOrientation1	= Body::byId(interaction->getId1())->state->ori;
	contactPhysics->initialOrientation2	= Body::byId(interaction->getId2())->state->ori;
	contactPhysics->prevNormal 		= scg->normal; //This is also done in the Contact Law.  It is not redundant because this class is only called ONCE!

	interaction->interactionPhysics = contactPhysics;
}

MomentPhys::~MomentPhys(){}

Ip2_MomentMat_MomentMat_MomentPhys::~Ip2_MomentMat_MomentMat_MomentPhys(){};



/*************************************************************************
*  Copyright (C) 2006 by Bruno Chareyre				 	 *
*  bruno.chareyre@hmg.inpg.fr					  	 *
*									 *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
 
#include"TriaxialStressController.hpp"
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
#include<yade/core/State.hpp>

#include<yade/core/Scene.hpp>

CREATE_LOGGER(TriaxialStressController);
YADE_PLUGIN((TriaxialStressController));

TriaxialStressController::~TriaxialStressController()
{	
}

void TriaxialStressController::updateStiffness ()
{
	for (int i=0; i<6; ++i) stiffness[i] = 0;
	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for(  ; ii!=iiEnd ; ++ii ) if ((*ii)->isReal())
	{
		const shared_ptr<Interaction>& contact = *ii;
		Real fn = (static_cast<FrictPhys*>	(contact->interactionPhysics.get()))->normalForce.norm();
		if (fn!=0)
		{
			int id1 = contact->getId1(), id2 = contact->getId2();
			for (int index=0; index<6; ++index) if ( wall_id[index]==id1 || wall_id[index]==id2 )
			{
				FrictPhys* currentContactPhysics =
				static_cast<FrictPhys*> ( contact->interactionPhysics.get() );
				stiffness[index]  += currentContactPhysics->kn;
			}
		}
	}
}

void TriaxialStressController::controlExternalStress(int wall, Vector3r resultantForce, State* p, Real wall_max_vel)
{
	scene->forces.sync();
	Real translation=normal[wall].dot(getForce(scene,wall_id[wall])-resultantForce); 
	const bool log=false;
	if(log) LOG_DEBUG("wall="<<wall<<" actualForce="<<getForce(scene,wall_id[wall])<<", resultantForce="<<resultantForce<<", translation="<<translation);
	if (translation!=0)
	{
	   if (stiffness[wall]!=0)
	   {
			translation /= stiffness[wall];
			if(log) TRVAR2(translation,wall_max_vel*scene->dt)
			translation = std::min( std::abs(translation), wall_max_vel*scene->dt ) * Mathr::Sign(translation);
	   }
	   else translation = wall_max_vel * Mathr::Sign(translation)*scene->dt;
	}
	previousTranslation[wall] = (1-wallDamping)*translation*normal[wall] + 0.8*previousTranslation[wall];// formula for "steady-flow" evolution with fluctuations
	p->se3.position += previousTranslation[wall];
	// this is important is using VelocityBins. Otherwise the motion is never detected. Related to https://bugs.launchpad.net/yade/+bug/398089
	p->vel=previousTranslation[wall]/scene->dt;
	//if(log)TRVAR2(previousTranslation,p->se3.position);
}

void TriaxialStressController::action()
{
	// sync thread storage of ForceContainer
	scene->forces.sync();	
	if (first) {// sync boundaries ids in the table
		wall_id[0] = wall_bottom_id;
 		wall_id[1] = wall_top_id;
 		wall_id[2] = wall_left_id;
 		wall_id[3] = wall_right_id;
 		wall_id[4] = wall_front_id;
 		wall_id[5] = wall_back_id;}	

	if(thickness<0) thickness=2.0*YADE_PTR_CAST<Box>(Body::byId(wall_bottom_id,scene)->shape)->extents.y();
	State* p_bottom=Body::byId(wall_bottom_id,scene)->state.get();
	State* p_top=Body::byId(wall_top_id,scene)->state.get();
	State* p_left=Body::byId(wall_left_id,scene)->state.get();
	State* p_right=Body::byId(wall_right_id,scene)->state.get();
	State* p_front=Body::byId(wall_front_id,scene)->state.get();
	State* p_back=Body::byId(wall_back_id,scene)->state.get();
	height = p_top->se3.position.y() - p_bottom->se3.position.y() - thickness;
	width = p_right->se3.position.x() - p_left->se3.position.x() - thickness;
	depth = p_front->se3.position.z() - p_back->se3.position.z() - thickness;
	
	boxVolume = height * width * depth;
	if (first) {
		BodyContainer::iterator bi = scene->bodies->begin();
		BodyContainer::iterator biEnd = scene->bodies->end();
		spheresVolume = 0;
		for ( ; bi!=biEnd; ++bi )
		{
			if((*bi)->isClump()) continue;
			const shared_ptr<Body>& b = *bi;
			if ( b->isDynamic )
			{
				const shared_ptr<Sphere>& sphere =
						YADE_PTR_CAST<Sphere> ( b->shape );
				spheresVolume += 1.3333333*Mathr::PI*pow ( sphere->radius, 3 );
			}
		}
		max_vel1=3 * width /(height+width+depth)*max_vel;				
		max_vel2=3 * height /(height+width+depth)*max_vel;
		max_vel3 =3 * depth /(height+width+depth)*max_vel;
		first = false;
	}

	// NOT JUST at the first run, since sigma_iso may be changed
	// if the TriaxialCompressionEngine is used, sigma_iso is attributed to sigma1, sigma2 and sigma3
	if (isAxisymetric){
		sigma1=sigma2=sigma3=sigma_iso;
	}
	
	porosity = ( boxVolume - spheresVolume ) /boxVolume;
	position_top = p_top->se3.position.y();
	position_bottom = p_bottom->se3.position.y();
	position_right = p_right->se3.position.x();
	position_left = p_left->se3.position.x();
	position_front = p_front->se3.position.z();
	position_back = p_back->se3.position.z();
	
	// must be done _after_ height, width, depth have been calculated
	//Update stiffness only if it has been computed by StiffnessCounter (see "stiffnessUpdateInterval")
	if (Omega::instance().getCurrentIteration() % stiffnessUpdateInterval == 0 || Omega::instance().getCurrentIteration()<100) updateStiffness();
	bool isARadiusControlIteration = (Omega::instance().getCurrentIteration() % radiusControlInterval == 0);
	
	if (Omega::instance().getCurrentIteration() % computeStressStrainInterval == 0 ||
		 (internalCompaction && isARadiusControlIteration) )
		computeStressStrain();

	if (!internalCompaction) {
		Vector3r wallForce (0, sigma2*width*depth, 0);
		if (wall_bottom_activated) controlExternalStress(wall_bottom, -wallForce, p_bottom, max_vel2);
		if (wall_top_activated) controlExternalStress(wall_top, wallForce, p_top, max_vel2);
		
		wallForce = Vector3r(sigma1*height*depth, 0, 0);
		if (wall_left_activated) controlExternalStress(wall_left, -wallForce, p_left, max_vel1);
		if (wall_right_activated) controlExternalStress(wall_right, wallForce, p_right, max_vel1);
		
		wallForce = Vector3r(0, 0, sigma3*height*width);
		if (wall_back_activated) controlExternalStress(wall_back, -wallForce, p_back, max_vel3);
		if (wall_front_activated) controlExternalStress(wall_front, wallForce, p_front, max_vel3);
	}
	else //if internal compaction
	{
		if (isARadiusControlIteration) {
			//Real s = computeStressStrain(scene);
			if (sigma_iso<=meanStress) maxMultiplier = finalMaxMultiplier;
			if (meanStress==0) previousMultiplier = maxMultiplier;
			else {
				//     		previousMultiplier = 1+0.7*(sigma_iso-s)*(previousMultiplier-1.f)/(s-previousStress); // = (Dsigma/apparentModulus)*0.7
				//     		previousMultiplier = std::max(2-maxMultiplier, std::min(previousMultiplier, maxMultiplier));
				previousMultiplier = 1.+(sigma_iso-meanStress)/sigma_iso*(maxMultiplier-1.); // = (Dsigma/apparentModulus)*0.7
			}
			previousStress = meanStress;
			//Real apparentModulus = (s-previousStress)/(previousMultiplier-1.f);
			controlInternalStress(previousMultiplier);
		}
	}
}

void TriaxialStressController::computeStressStrain()
{
	scene->forces.sync();
	State* p_bottom=Body::byId(wall_bottom_id,scene)->state.get();
	State* p_top=Body::byId(wall_top_id,scene)->state.get();
	State* p_left=Body::byId(wall_left_id,scene)->state.get();
	State* p_right=Body::byId(wall_right_id,scene)->state.get();
	State* p_front=Body::byId(wall_front_id,scene)->state.get();
	State* p_back=Body::byId(wall_back_id,scene)->state.get();
	
 	height = p_top->se3.position.y() - p_bottom->se3.position.y() - thickness;
 	width = p_right->se3.position.x() - p_left->se3.position.x() - thickness;
 	depth = p_front->se3.position.z() - p_back->se3.position.z() - thickness;
	
	meanStress = 0;
	if (height0 == 0) height0 = height;
	if (width0 == 0) width0 = width;
	if (depth0 == 0) depth0 = depth;
	strain[0] = Mathr::Log(width0/width);
	strain[1] = Mathr::Log(height0/height);
	strain[2] = Mathr::Log(depth0/depth);
	volumetricStrain=strain[0]+strain[1]+strain[2];
	
	Real invXSurface = 1.f/(height*depth);
	Real invYSurface = 1.f/(width*depth);
	Real invZSurface = 1.f/(width*height);

 	force[wall_bottom]=getForce(scene,wall_id[wall_bottom]); stress[wall_bottom]=force[wall_bottom]*invYSurface;
	force[wall_top]=   getForce(scene,wall_id[wall_top]);    stress[wall_top]=force[wall_top]*invYSurface;
	force[wall_left]=  getForce(scene,wall_id[wall_left]);   stress[wall_left]=force[wall_left]*invXSurface;
	force[wall_right]= getForce(scene,wall_id[wall_right]);  stress[wall_right]= force[wall_right]*invXSurface;
	force[wall_front]= getForce(scene,wall_id[wall_front]);  stress[wall_front]=force[wall_front]*invZSurface;
        force[wall_back]=  getForce(scene,wall_id[wall_back]);   stress[wall_back]= force[wall_back]*invZSurface;

	for (int i=0; i<6; i++) meanStress-=stress[i].dot(normal[i]);
	meanStress/=6.;
}

void TriaxialStressController::controlInternalStress ( Real multiplier )
{
	spheresVolume *= pow ( multiplier,3 );
	BodyContainer::iterator bi    = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	for ( ; bi!=biEnd ; ++bi )
	{
		if ( ( *bi )->isDynamic )
		{
			( static_cast<Sphere*> ( ( *bi )->shape.get() ) )->radius *= multiplier;
				(*bi)->state->mass*=pow(multiplier,3);
				(*bi)->state->inertia*=pow(multiplier,5);

		}
	}
	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for (; ii!=iiEnd ; ++ii)
	{
		if ((*ii)->isReal()) {
			ScGeom* contact = static_cast<ScGeom*>((*ii)->interactionGeometry.get());
			if ((*(scene->bodies))[(*ii)->getId1()]->isDynamic)
				contact->radius1 = static_cast<Sphere*>((* (scene->bodies))[(*ii)->getId1()]->shape.get())->radius;
			if ((* (scene->bodies))[(*ii)->getId2()]->isDynamic)
				contact->radius2 = static_cast<Sphere*>((* (scene->bodies))[(*ii)->getId2()]->shape.get())->radius;
			const shared_ptr<FrictPhys>& contactPhysics = YADE_PTR_CAST<FrictPhys>((*ii)->interactionPhysics);
			contactPhysics->kn*=multiplier; contactPhysics->ks*=multiplier;
		}
	}
}

/*!
    \fn TriaxialStressController::ComputeUnbalancedForce( bool maxUnbalanced)
 */
Real TriaxialStressController::ComputeUnbalancedForce( bool maxUnbalanced)
{
	scene->forces.sync();
	//compute the mean contact force
	Real MeanForce = 0.f;
	long nForce = 0;
	shared_ptr<BodyContainer>& bodies = scene->bodies;

	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for(  ; ii!=iiEnd ; ++ii ) {
		if ((*ii)->isReal()) {
			const shared_ptr<Interaction>& contact = *ii;
			Real f = (static_cast<FrictPhys*> ((contact->interactionPhysics.get()))->normalForce+static_cast<FrictPhys*>(contact->interactionPhysics.get())->shearForce).squaredNorm();
			if (f!=0)
			{
			MeanForce += Mathr::Sqrt(f);
			++nForce;
			}
		}
	}
	if (nForce!=0) MeanForce /= nForce;


	if (!maxUnbalanced) {
		//compute mean Unbalanced Force
		Real MeanUnbalanced=0;
		long nBodies = 0;
		BodyContainer::iterator bi    = bodies->begin();
		BodyContainer::iterator biEnd = bodies->end();
		Real f;
		for(  ; bi!=biEnd ; ++bi ) {
			if ((*bi)->isDynamic) {
				f=getForce(scene,(*bi)->getId()).norm();
				MeanUnbalanced += f;
				if (f!=0) ++nBodies;
			}
		}
		if (nBodies != 0 && MeanForce != 0) MeanUnbalanced = MeanUnbalanced/nBodies/MeanForce;
		return  MeanUnbalanced;
	} else {
		//compute max Unbalanced Force
		Real MaxUnbalanced=0;
		BodyContainer::iterator bi    = bodies->begin();
		BodyContainer::iterator biEnd = bodies->end();
		for(  ; bi!=biEnd ; ++bi ) if ((*bi)->isDynamic)
				MaxUnbalanced = std::max(getForce(scene,(*bi)->getId()).norm(),MaxUnbalanced);
		if (MeanForce != 0) MaxUnbalanced = MaxUnbalanced/MeanForce;
		return MaxUnbalanced;
	}
}



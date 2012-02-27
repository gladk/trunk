/*************************************************************************
*  Copyright (C) 2009 by Luc Sibille                                     *
*  luc.sibille@univ-nantes.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "ThreeDTriaxialEngine.hpp"
#include<yade/core/Scene.hpp>
#include<yade/core/Omega.hpp>

#include<yade/lib-base/Math.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/lambda/lambda.hpp>
#include<yade/pkg-dem/Shop.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
#include<yade/pkg-common/ElastMat.hpp>

class Ip2_2xCohFrictMat_CohFrictPhys;

CREATE_LOGGER(ThreeDTriaxialEngine);
YADE_PLUGIN((ThreeDTriaxialEngine));


// ThreeDTriaxialEngine::ThreeDTriaxialEngine()
// {
// 	translationAxisy=Vector3r(0,1,0);
// 	translationAxisx=Vector3r(1,0,0);
// 	translationAxisz=Vector3r(0,0,1);
// 	strainRate1=0;
// 	currentStrainRate1=0;
// 	strainRate2=0;
// 	currentStrainRate2=0;
// 	strainRate2=0;
// 	currentStrainRate2=0;
// 	//StabilityCriterion=0.001;
// 	UnbalancedForce = 1;
// 	Key = "";
// 	//Iteration = 0;
// 	//testEquilibriumInterval = 20;
// 	firstRun=true;
// 	frictionAngleDegree = -1;
// 	updateFrictionAngle=false;
// 
// 	stressControl_1=false;
// 	stressControl_2=false;
// 	stressControl_3=false;
// 
//  	boxVolume=0;
// 
// }

ThreeDTriaxialEngine::~ThreeDTriaxialEngine()
{	
}


void ThreeDTriaxialEngine::action()
{
  
	if ( firstRun )
	{
		LOG_INFO ( "First run, will initialize!" );

		if (updateFrictionAngle) setContactProperties(frictionAngleDegree);
		
		height0 = height; depth0 = depth; width0 = width;

		if (stressControl_1){
		wall_right_activated=true; wall_left_activated=true; //are the right walls for direction 1?
		} else {
		wall_right_activated=false; wall_left_activated=false;
		}

		if (stressControl_2){
		wall_bottom_activated=true; wall_top_activated=true;
		} else {
		wall_bottom_activated=false; wall_top_activated=false;
		}

		if (stressControl_3){
		wall_front_activated=true; wall_back_activated=true; //are the right walls for direction 3?
		} else {
		wall_front_activated=false; wall_back_activated=false;
		}

		//internalCompaction=false;  //is needed to avoid a control for internal compaction by the TriaxialStressController engine
		
		isAxisymetric=false; //is needed to avoid a stress control according the parameter sigma_iso (but according to sigma1, sigma2 and sigma3)
	
		firstRun=false;

	}

	
	const Real& dt = scene->dt;

	if(!stressControl_1)  // control in strain if wanted
	{
		if ( currentStrainRate1 != strainRate1 ) currentStrainRate1 += ( strainRate1-currentStrainRate1 ) *0.0003;

		State* p_left=Body::byId(wall_left_id,scene)->state.get();
		p_left->pos += 0.5*currentStrainRate1*width*translationAxisx*dt;
		State* p_right=Body::byId(wall_right_id,scene)->state.get();
		p_right->pos -= 0.5*currentStrainRate1*width*translationAxisx*dt;

	} else {

		if ( currentStrainRate1 != strainRate1 ) currentStrainRate1 += ( strainRate1-currentStrainRate1 ) *0.0003;
		max_vel1 = 0.5*currentStrainRate1*width;
	}

	
	if(!stressControl_2)  // control in strain if wanted
	{
		if ( currentStrainRate2 != strainRate2 ) currentStrainRate2 += ( strainRate2-currentStrainRate2 ) *0.0003;

		State* p_bottom=Body::byId(wall_bottom_id,scene)->state.get();
		p_bottom->pos += 0.5*currentStrainRate2*height*translationAxisy*dt;
		State* p_top=Body::byId(wall_top_id,scene)->state.get();
		p_top->pos -= 0.5*currentStrainRate2*height*translationAxisy*dt;

	} else {

		if ( currentStrainRate2 != strainRate2 ) currentStrainRate2 += ( strainRate2-currentStrainRate2 ) *0.0003;
		max_vel2 = 0.5*currentStrainRate2*height;
	}


	if(!stressControl_3)  // control in strain if wanted
	{
		if ( currentStrainRate3 != strainRate3 ) currentStrainRate3 += ( strainRate3-currentStrainRate3 ) *0.0003;


		State* p_back=Body::byId(wall_back_id,scene)->state.get();
		p_back->pos += 0.5*currentStrainRate3*depth*translationAxisz*dt;
		State* p_front=Body::byId(wall_front_id,scene)->state.get();
		p_front->pos -= 0.5*currentStrainRate3*depth*translationAxisz*dt;

	} else {

		if ( currentStrainRate3 != strainRate3 ) currentStrainRate3 += ( strainRate3-currentStrainRate3 ) *0.0003;
		max_vel3 = 0.5*currentStrainRate3*depth;
	}
	
	TriaxialStressController::action(); // this function is called to perform the external stress control or the internal compaction

}

void ThreeDTriaxialEngine::setContactProperties(Real frictionDegree)
{
	scene = Omega::instance().getScene().get();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	FOREACH(const shared_ptr<Body>& b,*scene->bodies){
		if (b->isDynamic)
		YADE_PTR_CAST<FrictMat> (b->material)->frictionAngle = frictionDegree * Mathr::PI/180.0;
	}
				
	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions){
		if (!ii->isReal()) continue;
		const shared_ptr<FrictMat>& sdec1 = YADE_PTR_CAST<FrictMat>((*bodies)[(body_id_t) ((ii)->getId1())]->material);
		const shared_ptr<FrictMat>& sdec2 = YADE_PTR_CAST<FrictMat>((*bodies)[(body_id_t) ((ii)->getId2())]->material);
		//FIXME - why dynamic_cast fails here?
		FrictPhys* contactPhysics = YADE_CAST<FrictPhys*>((ii)->interactionPhysics.get());
		Real fa = sdec1->frictionAngle;
		Real fb = sdec2->frictionAngle;
		contactPhysics->frictionAngle			= std::min(fa,fb);
		contactPhysics->tangensOfFrictionAngle		= std::tan(contactPhysics->frictionAngle);
	}
  
} 



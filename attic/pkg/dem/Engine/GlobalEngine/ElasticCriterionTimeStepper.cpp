/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "ElasticCriterionTimeStepper.hpp"
#include<yade/pkg-dem/BodyMacroParameters.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/MacroMicroElasticRelationships.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/Sphere.hpp>


ElasticCriterionTimeStepper::ElasticCriterionTimeStepper() : TimeStepper() , sdecContactModel(new MacroMicroElasticRelationships)
{
	sdecGroupMask = 1;
}


ElasticCriterionTimeStepper::~ElasticCriterionTimeStepper()
{

}



void ElasticCriterionTimeStepper::findTimeStepFromBody(const shared_ptr<Body>& body)
{
	BodyMacroParameters * sdec	= dynamic_cast<BodyMacroParameters*>(body->physicalParameters.get());
	Sphere* sphere 		= dynamic_cast<Sphere*>(body->shape.get());

	if(! (sphere && sdec) )
		return; // not possible to compute!

	Real Dab  	= sphere->radius;
	Real rad3 	= std::pow(Dab,2); // radius to the power of 2, from sphere

	Real Eab 	= sdec->young;
	Real Vab 	= sdec->poisson;
	Real Dinit 	= 2*Dab; // assuming thet sphere is in contact with itself
	Real Sinit 	= Mathr::PI * std::pow( Dab , 2);

	Real alpha 	= sdecContactModel->alpha;
	Real beta 	= sdecContactModel->beta;
	Real gamma 	= sdecContactModel->gamma;

	Real Kn		= abs((Eab*Sinit/Dinit)*( (1+alpha)/(beta*(1+Vab) + gamma*(1-alpha*Dab) ) ));
	Real Ks		= abs(Kn*(1-alpha*Vab)/(1+Vab));

	Real dt = 0.1*min(
			  sqrt( sdec->mass       /  Kn       )
			, sqrt( sdec->inertia[0] / (Ks*rad3) )
		  );

	newDt = std::min(dt,newDt);
	computedSomething = true;
}


void ElasticCriterionTimeStepper::findTimeStepFromInteraction(const shared_ptr<Interaction>& interaction, shared_ptr<BodyContainer>& bodies)
{
	unsigned int id1 = interaction->getId1();
	unsigned int id2 = interaction->getId2();
		
	if( !( (*bodies)[id1]->getGroupMask() & (*bodies)[id2]->getGroupMask() & sdecGroupMask) )
		return; // skip other groups

	FrictPhys* sdecContact = dynamic_cast<FrictPhys*>(interaction->interactionPhysics.get());
	// if(!sdecContact) sdecContact=dynamic_cast<SDECLinkPhysics*>(interaction->interactionPhysics.get());
	ScGeom* interactionGeometry = dynamic_cast<ScGeom*>(interaction->interactionGeometry.get());
	BodyMacroParameters * body1	= dynamic_cast<BodyMacroParameters*>((*bodies)[id1]->physicalParameters.get());
	BodyMacroParameters * body2	= dynamic_cast<BodyMacroParameters*>((*bodies)[id2]->physicalParameters.get());

	if(! (sdecContact && interactionGeometry && body1 && body2))
		return;
	
	Real mass 	= std::min( body1->mass                   , body2->mass               );
//	if(mass == 0) 			// FIXME - remove that comment: zero mass and zero inertia are too stupid to waste time checking that.
//		mass 	= std::max( body1->mass                   , body2->mass               );
//	if(mass == 0)
//		return; // not possible to compute
	Real inertia 	= std::min( body1->inertia[0]             , body2->inertia[0]         );
//	if( inertia == 0)
//		inertia = std::max( body1->inertia[0]             , body2->inertia[0]         );
//	if( inertia == 0)
//		return;
	Real rad3 	= std::pow( std::max(interactionGeometry->radius1 , interactionGeometry->radius2 ) , 2); // radius to the power of 2, from sphere

	Real dt = 0.1*min(
			  sqrt( mass     / abs(sdecContact->initialKn)      )
			, sqrt( inertia  / abs(sdecContact->initialKs*rad3) )
		  );

	newDt = std::min(dt,newDt);
	computedSomething = true;
}


void ElasticCriterionTimeStepper::computeTimeStep()
{
	shared_ptr<BodyContainer>& bodies = ncb->bodies;
	shared_ptr<InteractionContainer>& interactions = ncb->interactions;

	newDt = Mathr::MAX_REAL;
	computedSomething = false; // this flag is to avoid setting timestep to MAX_REAL :)

	InteractionContainer::iterator ii    = interactions->begin();
	InteractionContainer::iterator iiEnd = interactions->end();
	for(  ; ii!=iiEnd ; ++ii )
		findTimeStepFromInteraction(*ii , bodies);

	if(! computedSomething)
	{
// no interactions at all? so let's try to estimate timestep by investigating bodies,
// simulating that a body in contact with itself. this happens only when there were not interactions at all.
		BodyContainer::iterator bi    = bodies->begin();
		BodyContainer::iterator biEnd = bodies->end();
		for(  ; bi!=biEnd ; ++bi )
		{
			shared_ptr<Body> b = *bi;
			if( b->getGroupMask() & sdecGroupMask)
				findTimeStepFromBody(b);
		}
	}	
	if(computedSomething)
	{
		Omega::instance().setTimeStep(newDt);		
		//cerr << "ElasticCriterionTimeStepper, timestep chosen is:" << Omega::instance().getTimeStep() << endl;
	}
}

YADE_PLUGIN((ElasticCriterionTimeStepper));

YADE_REQUIRE_FEATURE(PHYSPAR);


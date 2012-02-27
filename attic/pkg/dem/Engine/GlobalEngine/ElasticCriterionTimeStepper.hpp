/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/TimeStepper.hpp>

class Interaction;
class BodyContainer;
class MacroMicroElasticRelationships;

class ElasticCriterionTimeStepper : public TimeStepper
{
	private :
		Real		newDt;
		bool		computedSomething;
		shared_ptr<MacroMicroElasticRelationships> sdecContactModel;

		void findTimeStepFromBody(const shared_ptr<Body>&);
		void findTimeStepFromInteraction(const shared_ptr<Interaction>& , shared_ptr<BodyContainer>&);

	public :
		int sdecGroupMask; // FIXME - we should find a way to clean groupmask stuff

		ElasticCriterionTimeStepper();
		virtual ~ElasticCriterionTimeStepper();
	
		virtual void computeTimeStep(Scene* body);
	REGISTER_ATTRIBUTES(TimeStepper,(sdecGroupMask));
	REGISTER_CLASS_NAME(ElasticCriterionTimeStepper);
	REGISTER_BASE_CLASS_NAME(TimeStepper);
};

REGISTER_SERIALIZABLE(ElasticCriterionTimeStepper);



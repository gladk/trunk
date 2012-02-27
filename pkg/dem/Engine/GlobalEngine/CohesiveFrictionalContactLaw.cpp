/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre <bruno.chareyre@imag.fr>         *
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "CohesiveFrictionalContactLaw.hpp"
#include<yade/pkg-dem/CohFrictMat.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/CohFrictPhys.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>


YADE_PLUGIN((CohesiveFrictionalContactLaw)(Law2_ScGeom_CohFrictPhys_ElasticPlastic));
CREATE_LOGGER(Law2_ScGeom_CohFrictPhys_ElasticPlastic);

Vector3r translation_vect_ ( 0.10,0,0 );

void out ( Quaternionr q )
{
	AngleAxisr aa(angleAxisFromQuat(q));
	std::cout << " axis: " <<  aa.axis()[0] << " " << aa.axis()[1] << " " << aa.axis()[2] << ", angle: " << aa.angle() << " | ";
}

void outv ( Vector3r axis )
{
	std::cout << " axis: " <<  axis[0] << " " << axis[1] << " " << axis[2] << ", length: " << axis.norm() << " | ";
}



void CohesiveFrictionalContactLaw::action()
{
	if(!functor) functor=shared_ptr<Law2_ScGeom_CohFrictPhys_ElasticPlastic>(new Law2_ScGeom_CohFrictPhys_ElasticPlastic);
	functor->erosionActivated = erosionActivated;
	functor->detectBrokenBodies = detectBrokenBodies;
	functor->always_use_moment_law = always_use_moment_law;
	functor->shear_creep=shear_creep;
	functor->twist_creep=twist_creep;
	functor->creep_viscosity = creep_viscosity;
	functor->scene=scene;
	
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		functor->go(I->interactionGeometry, I->interactionPhysics, I.get(), scene);
	}
}


void Law2_ScGeom_CohFrictPhys_ElasticPlastic::go(shared_ptr<InteractionGeometry>& ig, shared_ptr<InteractionPhysics>& ip, Interaction* contact, Scene* ncb)
{
	const Real& dt = scene->dt;
// 		if (detectBrokenBodies  //Experimental, has no effect
// 		        && (*bodies)[contact->getId1()]->shape->getClassName() != "box"
// 		        && (*bodies)[contact->getId2()]->shape->getClassName() != "box") {
// 			YADE_CAST<CohFrictMat*> ((*bodies)[contact->getId1()]->material.get())->isBroken = false;
// 			YADE_CAST<CohFrictMat*> ((*bodies)[contact->getId2()]->material.get())->isBroken = false;}
	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	Body* b1 = Body::byId(id1,ncb).get();
	Body* b2 = Body::byId(id2,ncb).get();
	ScGeom* currentContactGeometry  = YADE_CAST<ScGeom*> (ig.get());
	CohFrictPhys* currentContactPhysics = YADE_CAST<CohFrictPhys*> (ip.get());

	Vector3r& shearForce    = currentContactPhysics->shearForce;

	if (contact->isFresh(ncb)) shearForce   = Vector3r::Zero();
	Real un     = currentContactGeometry->penetrationDepth;
	Real Fn    = currentContactPhysics->kn*un;
	currentContactPhysics->normalForce = Fn*currentContactGeometry->normal;
	if (un < 0 && (currentContactPhysics->normalForce.squaredNorm() > pow(currentContactPhysics->normalAdhesion,2)
	               || currentContactPhysics->normalAdhesion==0)) {
		// BREAK due to tension
		ncb->interactions->requestErase(contact->getId1(),contact->getId2());
		// contact->interactionPhysics was reset now; currentContactPhysics still hold the object, but is not associated with the interaction anymore
// 			currentContactPhysics->cohesionBroken = true;
// 			currentContactPhysics->normalForce = Vector3r::ZERO;
// 			currentContactPhysics->shearForce = Vector3r::ZERO;
	} else {
		State* de1 = Body::byId(id1,ncb)->state.get();
		State* de2 = Body::byId(id2,ncb)->state.get();
		///////////////////////// CREEP START ///////////
		if (shear_creep) shearForce -= currentContactPhysics->ks*(shearForce*dt/creep_viscosity);
		///////////////////////// CREEP END ////////////
		Vector3r dus = currentContactGeometry->rotateAndGetShear(shearForce,currentContactPhysics->prevNormal,de1,de2,dt);
		//Linear elasticity giving "trial" shear force
		shearForce -= currentContactPhysics->ks*dus;

		Real Fs = currentContactPhysics->shearForce.norm();
		Real maxFs = currentContactPhysics->shearAdhesion;
		if (!currentContactPhysics->cohesionDisablesFriction || maxFs==0)
			maxFs += Fn*currentContactPhysics->tangensOfFrictionAngle;
		maxFs = std::max((Real) 0, maxFs);
		if (Fs  > maxFs) {//Plasticity condition on shear force
			if (currentContactPhysics->fragile && !currentContactPhysics->cohesionBroken) {
				currentContactPhysics->SetBreakingState();
				maxFs = max((Real) 0, Fn*currentContactPhysics->tangensOfFrictionAngle);
			}
			maxFs = maxFs / Fs;
			if (maxFs>1) cerr << "maxFs>1!!" << endl;
			shearForce *= maxFs;
			if (Fn<0)  currentContactPhysics->normalForce = Vector3r::Zero();//Vector3r::Zero()
		}

		applyForceAtContactPoint(-currentContactPhysics->normalForce-shearForce, currentContactGeometry->contactPoint, id1, de1->se3.position, id2, de2->se3.position, ncb);

		/// Moment law        ///
		if (momentRotationLaw && (!currentContactPhysics->cohesionBroken || always_use_moment_law)) {
			// Not necessary. OK.
			//{// updates only orientation of contact (local coordinate system)
			// Vector3r axis = currentContactPhysics->prevNormal.UnitCross(currentContactGeometry->normal);
			// Real angle =  unitVectorsAngle(currentContactPhysics->prevNormal,currentContactGeometry->normal);
			// Quaternionr align(axis,angle);
			// currentContactPhysics->currentContactOrientation =  align * currentContactPhysics->currentContactOrientation;
			//}

			Quaternionr delta(b1->state->ori * currentContactPhysics->initialOrientation1.conjugate() *
			                  currentContactPhysics->initialOrientation2 * b2->state->ori.conjugate());
			if (twist_creep) {
				delta = delta * currentContactPhysics->twistCreep;
			}

			AngleAxisr aa(angleAxisFromQuat(delta)); // axis of rotation - this is the Moment direction UNIT vector; // angle represents the power of resistant ELASTIC moment
			if (aa.angle() > Mathr::PI) aa.angle() -= Mathr::TWO_PI;   // angle is between 0 and 2*pi, but should be between -pi and pi

			Real angle_twist(aa.angle() * aa.axis().dot(currentContactGeometry->normal));
			Vector3r axis_twist(angle_twist * currentContactGeometry->normal);

			if (twist_creep) {
				Real viscosity_twist = creep_viscosity * std::pow((2 * std::min(currentContactGeometry->radius1,currentContactGeometry->radius2)),2) / 16.0;
				Real angle_twist_creeped = angle_twist * (1 - dt/viscosity_twist);
				Quaternionr q_twist(AngleAxisr(angle_twist,currentContactGeometry->normal));
				//Quaternionr q_twist_creeped(currentContactGeometry->normal , angle_twist*0.996);
				Quaternionr q_twist_creeped(AngleAxisr(angle_twist_creeped,currentContactGeometry->normal));
				Quaternionr q_twist_delta(q_twist_creeped * q_twist.conjugate());
				currentContactPhysics->twistCreep = currentContactPhysics->twistCreep * q_twist_delta;
				// modify the initialRelativeOrientation to substract some twisting
				// currentContactPhysics->initialRelativeOrientation = currentContactPhysics->initialRelativeOrientation * q_twist_delta;
				//currentContactPhysics->initialOrientation1 = currentContactPhysics->initialOrientation1 * q_twist_delta;
				//currentContactPhysics->initialOrientation2 = currentContactPhysics->initialOrientation2 * q_twist_delta.Conjugate();
			}
			Vector3r moment_twist(axis_twist * currentContactPhysics->kr);
			Vector3r axis_bending(aa.angle()*aa.axis() - axis_twist);
			Vector3r moment_bending(axis_bending * currentContactPhysics->kr);
			Vector3r moment = moment_twist + moment_bending;
			currentContactPhysics->moment_twist = moment_twist;
			currentContactPhysics->moment_bending = moment_bending;
			scene->forces.addTorque(id1,-moment);
			scene->forces.addTorque(id2, moment);
		}
		/// Moment law END       ///
		currentContactPhysics->prevNormal = currentContactGeometry->normal;
	}
}

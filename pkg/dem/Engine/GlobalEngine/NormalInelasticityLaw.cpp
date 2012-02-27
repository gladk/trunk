/*************************************************************************
*  Copyright (C) 2008 by Jérôme DURIEZ                                   *
*  duriez@geo.hmg.inpg.fr                                                *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<yade/pkg-dem/NormalInelasticityLaw.hpp>

#include<yade/pkg-dem/NormalInelasticMat.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>

YADE_PLUGIN((Law2_ScGeom_NormalInelasticityPhys_NormalInelasticity));



void Law2_ScGeom_NormalInelasticityPhys_NormalInelasticity::go(shared_ptr<InteractionGeometry>& iG, shared_ptr<InteractionPhysics>& iP, Interaction* contact, Scene* scene)
{

	const Real& dt = scene->dt;

	int id1 = contact->getId1();
	int id2 = contact->getId2();
// 	cout << "contact entre " << id1 << " et " << id2;

	State* de1 = Body::byId(id1,scene)->state.get();
	State* de2 = Body::byId(id2,scene)->state.get();
	NormalInelasticMat* Mat1 = static_cast<NormalInelasticMat*>(Body::byId(id1,scene)->material.get());
	ScGeom* currentContactGeometry		= YADE_CAST<ScGeom*>(iG.get());
	NormalInelasticityPhys* currentContactPhysics = YADE_CAST<NormalInelasticityPhys*> (iP.get());

	Vector3r& shearForce 			= currentContactPhysics->shearForce;

	if (contact->isFresh(scene))
		{
		shearForce			= Vector3r::Zero();
		currentContactPhysics->previousun=0.0;
		currentContactPhysics->previousFn=0.0;
		currentContactPhysics->unMax=0.0;
		}


	// *** Computation of normal Force : depends of the history *** //
	if( currentContactGeometry->penetrationDepth < currentContactPhysics->unMax) 
		currentContactPhysics->kn *= Mat1->coeff_dech;

	Real Fn; // la valeur de Fn qui va etre calculee selon différentes manieres puis affectee
// 	cout << " Dans Law2 valeur de kn : " << currentContactPhysics->kn << endl;
	Real un = currentContactGeometry->penetrationDepth; // compte positivement lorsq vraie interpenetration
	
// 	cout << "un = " << un << " alors que unMax = "<< currentContactPhysics->unMax << " et previousun = " << currentContactPhysics->previousun << " et previousFn =" << currentContactPhysics->previousFn << endl;
	
	if(un >= currentContactPhysics->unMax)	// on est en charge, et au delà et sur la "droite principale"
	{
		Fn = currentContactPhysics->knLower*un;
		currentContactPhysics->unMax = std::abs(un);
// 		cout << "je suis dans le calcul normal " << endl;
	}
	else// a priori then we need a greater stifness. False in the case when we are already on the limit state, but a correction below will then do the jobthe law2 will do the job in the case with a correction
	{
		currentContactPhysics->kn = currentContactPhysics->knLower* Mat1->coeff_dech;
		Fn = currentContactPhysics->previousFn + currentContactPhysics->kn * (un-currentContactPhysics->previousun);	// Calcul incrémental de la nvlle force
// 		cout << "je suis dans l'autre calcul" << endl;
		if(std::abs(Fn) > std::abs(currentContactPhysics->knLower * un))		// verif qu'on ne depasse la courbe
		{
		    Fn = currentContactPhysics->knLower*un;
// 		    cout <<  "j'etais dans l'autre calcul mais j'ai corrige pour ne pas depasser la limite" << endl;
		}
		if(Fn < 0.0 )	// verif qu'on reste positif
		{
			Fn = 0;
// 			cout << "j'ai corrige pour ne pas etre negatif" << endl;
		}
	}
	currentContactPhysics->normalForce	= Fn*currentContactGeometry->normal;
// 	cout << "Fn appliquee " << Fn << endl << endl;
	// actualisation :
	currentContactPhysics->previousFn = Fn;
	currentContactPhysics->previousun = un;
// 	*** End of computation of normal force *** //

	
        if (   un < 0      )
        { // BREAK due to tension
		 scene->interactions->requestErase(contact->getId1(),contact->getId2());
		 // probably not useful anymore
                currentContactPhysics->normalForce = Vector3r::Zero();
                currentContactPhysics->shearForce = Vector3r::Zero();
        }
        else
        {
                Vector3r axis;
                Real angle;

//                 Here is the code with approximated rotations

		axis	 		= currentContactPhysics->prevNormal.cross(currentContactGeometry->normal);
		shearForce 	       -= shearForce.cross(axis);
		angle 			= dt*0.5*currentContactGeometry->normal.dot(de1->angVel+de2->angVel);
		axis 			= angle*currentContactGeometry->normal;
		shearForce 	       -= shearForce.cross(axis);

//                 Here is the code with exact rotationns

                // 		Quaternionr q;
                //
                // 		axis					= currentContactPhysics->prevNormal.cross(currentContactGeometry->normal);
                // 		angle					= acos(currentContactGeometry->normal.dot(currentContactPhysics->prevNormal));
                // 		q.fromAngleAxis(angle,axis);
                //
                // 		currentContactPhysics->shearForce	= currentContactPhysics->shearForce*q;
                //
                // 		angle					= dt*0.5*currentContactGeometry->normal.dot(de1->angVel+de2->angVel);
                // 		axis					= currentContactGeometry->normal;
                // 		q.fromAngleAxis(angle,axis);
                // 		currentContactPhysics->shearForce	= q*currentContactPhysics->shearForce;

                /// 							 ///

                Vector3r x				= currentContactGeometry->contactPoint;
                Vector3r c1x				= (x - de1->se3.position);
                Vector3r c2x				= (x - de2->se3.position);
                /// The following definition of c1x and c2x is to avoid "granular ratcheting" 
		///  (see F. ALONSO-MARROQUIN, R. GARCIA-ROJO, H.J. HERRMANN, 
		///   "Micro-mechanical investigation of granular ratcheting, in Cyclic Behaviour of Soils and Liquefaction Phenomena",
		///   ed. T. Triantafyllidis (Balklema, London, 2004), p. 3-10 - and a lot more papers from the same authors)
                Vector3r _c1x_	= currentContactGeometry->radius1*currentContactGeometry->normal;
                Vector3r _c2x_	= -currentContactGeometry->radius2*currentContactGeometry->normal;
                Vector3r relativeVelocity		= (de2->vel+de2->angVel.cross(_c2x_)) - (de1->vel+de1->angVel.cross(_c1x_));
                Vector3r shearVelocity			= relativeVelocity-currentContactGeometry->normal.dot(relativeVelocity)*currentContactGeometry->normal;
                Vector3r shearDisplacement		= shearVelocity*dt;


		shearForce -= currentContactPhysics->ks*shearDisplacement;

                Real maxFs = 0;
                Real Fs = currentContactPhysics->shearForce.norm();
                maxFs = std::max((Real) 0,Fn*currentContactPhysics->tangensOfFrictionAngle);
                
                if ( Fs  > maxFs )
                {

			maxFs = max((Real) 0, Fn * currentContactPhysics->tangensOfFrictionAngle);

			maxFs = maxFs / Fs;
			if (maxFs>1)
                        	cerr << "maxFs>1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                    	shearForce *= maxFs;
			if (Fn<0)  currentContactPhysics->normalForce = Vector3r::Zero();
                }


                //////// PFC3d SlipModel

                Vector3r f	= currentContactPhysics->normalForce + shearForce;
		scene->forces.addForce (id1,-f);
		scene->forces.addForce (id2, f);
		scene->forces.addTorque(id1,-c1x.cross(f));
		scene->forces.addTorque(id2, c2x.cross(f));


	// Moment law					 	 ///
		if(momentRotationLaw)
		{
			{// updates only orientation of contact (local coordinate system)
				Vector3r axis = currentContactPhysics->prevNormal.cross(currentContactGeometry->normal); axis.normalize();
				Real angle =  unitVectorsAngle(currentContactPhysics->prevNormal,currentContactGeometry->normal);
				Quaternionr align(AngleAxisr(angle,axis));
				currentContactPhysics->currentContactOrientation =  align * currentContactPhysics->currentContactOrientation;
			}

			Quaternionr delta( de1->se3.orientation * currentContactPhysics->initialOrientation1.conjugate() *
		                           currentContactPhysics->initialOrientation2 * de2->se3.orientation.conjugate());

			AngleAxisr aa(angleAxisFromQuat(delta)); // aa.axis() of rotation - this is the Moment direction UNIT vector; angle represents the power of resistant ELASTIC moment
			if(angle > Mathr::PI) angle -= Mathr::TWO_PI; // angle is between 0 and 2*pi, but should be between -pi and pi 

	//This indentation is a rewrite of original equations (the two commented lines), should work exactly the same.
//Real elasticMoment = currentContactPhysics->kr * std::abs(angle); // positive value (*)

			Real angle_twist(aa.angle() * aa.axis().dot(currentContactGeometry->normal) );
			Vector3r axis_twist(angle_twist * currentContactGeometry->normal);
			Vector3r moment_twist(axis_twist * currentContactPhysics->kr);

			Vector3r axis_bending(aa.angle()*aa.axis() - axis_twist);
			Vector3r moment_bending(axis_bending * currentContactPhysics->kr);

//Vector3r moment = aa.axis() * elasticMoment * (aa.angle()<0.0?-1.0:1.0); // restore sign. (*)

			Vector3r moment = moment_twist + moment_bending;

// 	Limitation par seuil plastique
			if (!momentAlwaysElastic)
			{
				Real normeMomentMax = currentContactPhysics->forMaxMoment * std::fabs(Fn);
				Real normeMoment = moment.norm();
				if(normeMoment>normeMomentMax)
					{
					moment *= normeMomentMax/normeMoment;
// 					nbreInteracMomPlastif++;
					}
			}
			scene->forces.addTorque(id1,-moment);
			scene->forces.addTorque(id2, moment);
		}
	// Moment law	END				 	 ///

                currentContactPhysics->prevNormal = currentContactGeometry->normal;
//             }
//         }
    }
}




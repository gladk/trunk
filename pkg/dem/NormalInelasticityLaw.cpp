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

YADE_PLUGIN((Law2_ScGeom6D_NormalInelasticityPhys_NormalInelasticity));



void Law2_ScGeom6D_NormalInelasticityPhys_NormalInelasticity::go(shared_ptr<IGeom>& iG, shared_ptr<IPhys>& iP, Interaction* contact)
{

	const Real& dt = scene->dt;

	int id1 = contact->getId1();
	int id2 = contact->getId2();
// 	cout << "contact entre " << id1 << " et " << id2;

	State* de1 = Body::byId(id1,scene)->state.get();
	State* de2 = Body::byId(id2,scene)->state.get();
	NormalInelasticMat* Mat1 = static_cast<NormalInelasticMat*>(Body::byId(id1,scene)->material.get());
	ScGeom6D* currentContactGeometry		= YADE_CAST<ScGeom6D*>(iG.get());
	NormalInelasticityPhys* currentContactPhysics = YADE_CAST<NormalInelasticityPhys*> (iP.get());

	Vector3r& shearForce 			= currentContactPhysics->shearForce;

	if (contact->isFresh(scene))
		{
		shearForce			= Vector3r::Zero();
		currentContactPhysics->previousun=0.0;
		currentContactPhysics->previousFn=0.0;
		currentContactPhysics->unMax=0.0;
		}


//	********	Computation of normal Force : depends of the history			*******	 //

	Real Fn; // la valeur de Fn qui va etre calculee selon différentes manieres puis affectee
// 	cout << " Dans Law2 valeur de kn : " << currentContactPhysics->kn << endl;
	Real un = currentContactGeometry->penetrationDepth; // >0 for real penetration
	
// 	cout << "un = " << un << " alors que unMax = "<< currentContactPhysics->unMax << " et previousun = " << currentContactPhysics->previousun << " et previousFn =" << currentContactPhysics->previousFn << endl;
	
	if(un >= currentContactPhysics->unMax)	// case of virgin load : on the "principal line" (limit state of the (un,Fn) space)
	{
		Fn = currentContactPhysics->knLower*un;
		currentContactPhysics->unMax = std::abs(un);
// 		cout << "je suis dans le calcul normal " << endl;
	}
	else// a priori then we need a greater stifness. False in the case when we are already on the limit state, but a correction below will then do the job
	{
		currentContactPhysics->kn = currentContactPhysics->knLower* Mat1->coeff_dech;
		// Incremental computation of the new normal force :
		Fn = currentContactPhysics->previousFn + currentContactPhysics->kn * (un-currentContactPhysics->previousun);
// 		cout << "je suis dans l'autre calcul" << endl;
		if(std::abs(Fn) > std::abs(currentContactPhysics->knLower * un))	// check if the limit state is not violated
		{
		    Fn = currentContactPhysics->knLower*un;
// 		    cout <<  "j'etais dans l'autre calcul mais j'ai corrige pour ne pas depasser la limite" << endl;
		}
		if(Fn < 0.0 )	// check to stay >=0
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



//	********	Tangential force				*******	 //
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

		
                Vector3r x				= currentContactGeometry->contactPoint;
                Vector3r c1x				= (x - de1->se3.position);
                Vector3r c2x				= (x - de2->se3.position);
                /// The following definition of c1x and c2x is to avoid "granular ratcheting" 
		///  (see F. ALONSO-MARROQUIN, R. GARCIA-ROJO, H.J. HERRMANN, 
		///   "Micro-mechanical investigation of granular ratcheting, in Cyclic Behaviour of Soils and Liquefaction Phenomena",
		///   ed. T. Triantafyllidis (Balklema, London, 2004), p. 3-10 - and a lot more papers from the same authors, or discussions on yade mailing lists)
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
                        	cerr << "maxFs>1!!!!!!!!!!!!!!!!!!!" << endl;
                    	shearForce *= maxFs;
			if (Fn<0)  currentContactPhysics->normalForce = Vector3r::Zero();
                }


                Vector3r f	= currentContactPhysics->normalForce + shearForce;
		scene->forces.addForce (id1,-f);
		scene->forces.addForce (id2, f);
		scene->forces.addTorque(id1,-c1x.cross(f));
		scene->forces.addTorque(id2, c2x.cross(f));
// 	*** End of computation of tangential force *** //


//	********	Moment law				*******	 //

		if(momentRotationLaw)
		{
			currentContactPhysics->moment_twist = (currentContactGeometry->getTwist()*currentContactPhysics->kr)*currentContactGeometry->normal ;
			currentContactPhysics->moment_bending = currentContactGeometry->getBending() * currentContactPhysics->kr;

			moment = currentContactPhysics->moment_twist + currentContactPhysics->moment_bending;

// 	Limitation by plastic threshold
			if (!momentAlwaysElastic)
			{
				Real normeMomentMax = currentContactPhysics->forMaxMoment * std::fabs(Fn);
				Real normeMoment = moment.norm();
				if(normeMoment>normeMomentMax)
					{
					moment *= normeMomentMax/normeMoment;
					}
			}
			scene->forces.addTorque(id1,-moment);
			scene->forces.addTorque(id2, moment);
		}
//	********	Moment law END				*******	 //

                currentContactPhysics->prevNormal = currentContactGeometry->normal;
    }
}




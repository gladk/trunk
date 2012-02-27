/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "ElasticCohesiveLaw.hpp"
#include<yade/pkg-dem/BodyMacroParameters.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
//#include<yade/pkg-dem/SDECLinkGeometry.hpp>
#include<yade/pkg-dem/SDECLinkPhysics.hpp>


#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>

#include<yade/lib-base/yadeWm3Extra.hpp>



ElasticCohesiveLaw::ElasticCohesiveLaw() : InteractionSolver()
{
	sdecGroupMask=1;
	first=true;
	momentRotationLaw = true;
}


void ElasticCohesiveLaw::action()
{
	shared_ptr<BodyContainer>& bodies = ncb->bodies;

	Real dt = Omega::instance().getTimeStep();

/// Permanents Links													///

	InteractionContainer::iterator ii    = ncb->interactions->begin();
	InteractionContainer::iterator iiEnd = ncb->interactions->end();
	for(  ; ii!=iiEnd ; ++ii )
	{
		if ((*ii)->isReal())
		{
			const shared_ptr<Interaction> contact2 = *ii;
	
			unsigned int id1 = contact2->getId1();
			unsigned int id2 = contact2->getId2();
			
			if( !( (*bodies)[id1]->getGroupMask() & (*bodies)[id2]->getGroupMask() & sdecGroupMask) ) continue;
	
			BodyMacroParameters* de1			= YADE_CAST<BodyMacroParameters*>((*bodies)[id1]->physicalParameters.get());
			BodyMacroParameters* de2			= YADE_CAST<BodyMacroParameters*>((*bodies)[id2]->physicalParameters.get());
			SDECLinkPhysics* currentContactPhysics		= YADE_CAST<SDECLinkPhysics*>(contact2->interactionPhysics.get());
			ScGeom* currentContactGeometry	= YADE_CAST<ScGeom*>(contact2->interactionGeometry.get());
			if((!currentContactPhysics)||(!currentContactGeometry)) continue;
	
			Real un 					= currentContactPhysics->equilibriumDistance-(de2->se3.position-de1->se3.position).Length();
			currentContactPhysics->normalForce		= currentContactPhysics->kn*un*currentContactGeometry->normal;
	
			if (first)
				currentContactPhysics->prevNormal 	= currentContactGeometry->normal;
			
			Vector3r axis;
			Real angle;
	
	/// Here is the code with approximated rotations 	 ///
	
			axis = currentContactPhysics->prevNormal.Cross(currentContactGeometry->normal);
			currentContactPhysics->shearForce	       -= currentContactPhysics->shearForce.Cross(axis);
			angle	 					= dt*0.5*currentContactGeometry->normal.Dot(de1->angularVelocity+de2->angularVelocity);
			axis 						= angle*currentContactGeometry->normal;
			currentContactPhysics->shearForce     	       -= currentContactPhysics->shearForce.Cross(axis);
	
	/// Here is the code without approximated rotations 	 ///
	
	// 		Quaternionr q;
	//
	// 		axis				= currentContactPhysics->prevNormal.cross(currentContactGeometry->normal);
	// 		angle				= acos(currentContactGeometry->normal.dot(currentContactPhysics->prevNormal));
	// 		q.fromAngleAxis(angle,axis);
	//
	// 		currentContactPhysics->shearForce	= q*currentContactPhysics->shearForce;
	//
	// 		angle				= dt*0.5*currentContactGeometry->normal.dot(de1->angularVelocity+de2->angularVelocity);
	// 		axis				= currentContactGeometry->normal;
	// 		q.fromAngleAxis(angle,axis);
	// 		currentContactPhysics->shearForce	= q*currentContactPhysics->shearForce;
	
	/// 							 ///
	
			Vector3r x				= de1->se3.position+(currentContactGeometry->radius1-0.5*un)*currentContactGeometry->normal;
			//Vector3r x				= (de1->se3.position+de2->se3.position)*0.5;
			//cout << currentContact->contactPoint << " || " << (de1->se3.position+de2->se3.position)*0.5 << endl;
			Vector3r c1x				= (x - de1->se3.position);
			Vector3r c2x				= (x - de2->se3.position);
	
			Vector3r relativeVelocity 		= (de2->velocity+de2->angularVelocity.Cross(c2x)) - (de1->velocity+de1->angularVelocity.Cross(c1x));
			Vector3r shearVelocity			= relativeVelocity-currentContactGeometry->normal.Dot(relativeVelocity)*currentContactGeometry->normal;
			Vector3r shearDisplacement		= shearVelocity*dt;
			currentContactPhysics->shearForce      -= currentContactPhysics->ks*shearDisplacement;
	
			Vector3r f = currentContactPhysics->normalForce + currentContactPhysics->shearForce;
			ncb->forces.addForce (id1,-f);
			ncb->forces.addForce (id2,+f);
			ncb->forces.addTorque(id1,-c1x.Cross(f));
			ncb->forces.addTorque(id2, c2x.Cross(f));
	
	
	
	/// Moment law					 	 ///
		if(momentRotationLaw)
		{
	
			if (first)
			{
				currentContactPhysics->prevRotation1 = de1->se3.orientation;
				currentContactPhysics->prevRotation2 = de2->se3.orientation;
				currentContactPhysics->averageRadius = (currentContactGeometry->radius1+currentContactGeometry->radius2)*0.5;
				currentContactPhysics->kr = currentContactPhysics->ks * currentContactPhysics->averageRadius * currentContactPhysics->averageRadius;
			}
	
			Vector3r n	= currentContactGeometry->normal;
			Vector3r prevN	= currentContactPhysics->prevNormal;
			Vector3r t1	= currentContactPhysics->shearForce;
			t1.Normalize();
			Vector3r t2	= n.UnitCross(t1);
	
	
			Quaternionr q_i_n,q_n_i;
	
			q_i_n=quaternionFromAxes(n,t1,t2);
			q_i_n=quaternionFromAxes(Vector3r(1,0,0),Vector3r(0,1,0),Vector3r(0,0,1)); // use identity matrix
			q_n_i = q_i_n.Inverse();
	
	/// Using Euler angle				 	 ///
	
	// 		Vector3r dBeta;
	// 		Vector3r orientation_Nc,orientation_Nc_old;
	// 		for(int i=0;i<3;i++)
	// 		{
	// 			int j = (i+1)%3;
	// 			int k = (i+2)%3;
	//
	// 			if (n[j]>=0)
	// 				orientation_Nc[k] = acos(n[i]); // what is Nc_new
	// 			else
	// 				orientation_Nc[k] = -acos(n[i]);
	//
	// 			if (prevN[j]>=0)
	// 				orientation_Nc_old[k] = acos(prevN[i]);
	// 			else
	// 				orientation_Nc_old[k] = -acos(prevN[i]);
	// 		}
	//
	// 		dBeta = orientation_Nc - orientation_Nc_old;
	//
	//
	// 		Vector3r dRotationA,dRotationB,da,db;
	// 		de1->se3.orientation.toEulerAngles(dRotationA);
	// 		de2->se3.orientation.toEulerAngles(dRotationB);
	//
	// 		currentContactPhysics->prevRotation1.toEulerAngles(da);
	// 		currentContactPhysics->prevRotation2.toEulerAngles(db);
	//
	// 		dRotationA -= da;
	// 		dRotationB -= db;
	//
	// 		Vector3r dUr = 	( currentContactGeometry->radius1*(  dRotationA  -  dBeta)
	// 				- currentContactGeometry->radius2*(  dRotationB  -  dBeta) ) * 0.5;
	
	/// Ending of use of eurler angle		 	 ///
	
	/// Using Quaternionr				 	 ///
	
			Quaternionr q,q1,q2;
			Vector3r dRotationAMinusDBeta,dRotationBMinusDBeta;
	
			q.Align(n,prevN);
			q1 = (de1->se3.orientation*currentContactPhysics->prevRotation1.Inverse())*q.Inverse();
			q2 = (de2->se3.orientation*currentContactPhysics->prevRotation2.Inverse())*q.Inverse();
			quaternionToEulerAngles(q1,dRotationAMinusDBeta);
			quaternionToEulerAngles(q2,dRotationBMinusDBeta);
			Vector3r dUr = ( currentContactGeometry->radius1*dRotationAMinusDBeta - currentContactGeometry->radius2*dRotationBMinusDBeta ) * 0.5;
	
	/// Ending of use of Quaternionr			 	 ///
	
			Vector3r dThetar = dUr/currentContactPhysics->averageRadius;
			currentContactPhysics->thetar += dThetar;
	
		//	Real fNormal = currentContactPhysics->normalForce.Length();
		//	Real normMPlastic = currentContactPhysics->heta*fNormal;
			Vector3r thetarn = q_i_n*currentContactPhysics->thetar; // rolling angle
			Vector3r mElastic = currentContactPhysics->kr * thetarn;
	
			//mElastic[0] = 0;  // No moment around normal direction
	
		//	Real normElastic = mElastic.Length();
	
	
			//if (normElastic<=normMPlastic)
			//{
			ncb->forces.addTorque(id1,-q_n_i*mElastic);
			ncb->forces.addTorque(id2, q_n_i*mElastic);
	
			//}  
			//else
			//{
			//	Vector3r mPlastic = mElastic;
			//	mPlastic.normalize();
			//	mPlastic *= normMPlastic;
			//	moments[id1]	-= q_n_i*mPlastic;
			//	moments[id2]	+= q_n_i*mPlastic;
			//	thetarn = mPlastic/currentContactPhysics->kr;
			//	currentContactPhysics->thetar = q_n_i*thetarn;
			//}
	
			currentContactPhysics->prevRotation1 = de1->se3.orientation;
			currentContactPhysics->prevRotation2 = de2->se3.orientation;
	
	/// Moment law	END				 	 ///
		}
	
			currentContactPhysics->prevNormal = currentContactGeometry->normal;
		}	
	}

	first = false; 
}


YADE_PLUGIN((ElasticCohesiveLaw));

YADE_REQUIRE_FEATURE(PHYSPAR);


/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<yade/pkg-dem/KinemCTDEngine.hpp>
#include<yade/core/State.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/core/Scene.hpp>

YADE_PLUGIN((KinemCTDEngine));

KinemCTDEngine::~KinemCTDEngine()
{
}

void KinemCTDEngine::action()
{

	leftbox = Body::byId(id_boxleft);
	rightbox = Body::byId(id_boxright);
	frontbox = Body::byId(id_boxfront);
	backbox = Body::byId(id_boxback);
	topbox = Body::byId(id_topbox);
	boxbas = Body::byId(id_boxbas);

// Computation of the current dimensions of the box : //
	Real Xleft = leftbox->state->pos.x() + (YADE_CAST<Box*>(leftbox->shape.get()))->extents.x();

	Real Xright = rightbox->state->pos.x() - (YADE_CAST<Box*>(rightbox->shape.get()))->extents.x();

	Real Zfront = frontbox->state->pos.z() - YADE_CAST<Box*>(frontbox->shape.get())->extents.z();
	Real Zback = backbox->state->pos.z() + (YADE_CAST<Box*>(backbox->shape.get()))->extents.z();

	Real Scontact = (Xright-Xleft)*(Zfront-Zback);	// that's so the value of section at the middle of the height of the box
// End of computation of the current dimensions of the box //

	scene->forces.sync(); Vector3r F_sup=scene->forces.getForce(id_topbox);
	Real current_NormalForce=(scene->forces.getForce(id_topbox)).y();
	current_sigma=current_NormalForce/(1000.0*Scontact);	// so we have the current value of sigma, in kPa

	if( ((compSpeed > 0) && (current_sigma < target_sigma)) || ((compSpeed < 0) && (current_sigma > target_sigma)) )
	{
		if(temoin!=0)
		{
// 			cout << "j'ai ici un temoin #0 visiblement. En effet temoin =" <<lexical_cast<string>(temoin) << endl;
			temoin=0;
// 			cout << "Maintenant (toujours dans le if temoin!=0), temoin =" <<lexical_cast<string>(temoin) << endl;
		}
		letMove();
	}
	else if (temoin==0)
	{
		stopMovement();
// 		cout << "Mouvement stoppe, temoin = " << lexical_cast<string>(temoin) << endl;
// 		cout << " Dans le if, temoin =" << lexical_cast<string>(temoin) << endl;
		string f;
		if (compSpeed > 0)
			f="sigmax_";
		else
			f="Sigmin_";

		Omega::instance().saveSimulation(Key + f +lexical_cast<string> (floor(target_sigma)) + "kPaReached.xml");
		temoin=1;
// 		cout << " Fin du if, temoin =" << lexical_cast<string>(temoin) << endl << endl;
	}

	
	for(unsigned int j=0;j<sigma_save.size();j++)
	{
		if( (  ( (compSpeed>0)&&(current_sigma > sigma_save[j]) ) || ((compSpeed<0)&&(current_sigma < sigma_save[j])) ) && (temoin_save[j]==0))
		{
			stopMovement();
			Omega::instance().saveSimulation(Key + "SigInt_" +lexical_cast<string> (floor(current_sigma)) + "kPareached.xml");
			temoin_save[j]=1;
		}
	}
// 	cout << "Fin de ApplyCondi, temoin = " << lexical_cast<string>(temoin) << endl;
		
}

void KinemCTDEngine::letMove()
{
	computeAlpha();
	const Real& dt = scene->dt;
	Real dh=-compSpeed*dt
		,dgamma=0
		;

	Real Ysup = topbox->state->pos.y();
	Real Ylat = leftbox->state->pos.y();

// 	Changes in vertical and horizontal position :

	topbox->state->pos += Vector3r(dgamma,dh,0);

	leftbox->state->pos += Vector3r(dgamma/2.0,dh/2.0,0);
	rightbox->state->pos += Vector3r(dgamma/2.0,dh/2.0,0);

	Real Ysup_mod = topbox->state->pos.y();
	Real Ylat_mod = leftbox->state->pos.y();

//	with the corresponding velocities :
	topbox->state->vel = Vector3r(dgamma/dt,dh/dt,0);
	leftbox->state->vel = Vector3r((dgamma/2.0)/dt,dh/(2.0*dt),0);
	rightbox->state->vel = Vector3r((dgamma/2.0)/dt,dh/(2.0*dt),0);

//	Then computation of the angle of the rotation to be done :
	if (alpha == Mathr::PI/2.0)	// Case of the very beginning
	{
		dalpha = - Mathr::ATan( dgamma / (Ysup_mod -Ylat_mod) );
	}
	else
	{
		Real A = (Ysup_mod - Ylat_mod) * 2.0*Mathr::Tan(alpha) / (2.0*(Ysup - Ylat) + dgamma*Mathr::Tan(alpha) );
		dalpha = Mathr::ATan( (A - Mathr::Tan(alpha))/(1.0 + A * Mathr::Tan(alpha)));
	}

	Quaternionr qcorr(AngleAxisr(dalpha,Vector3r::UnitZ()));

// On applique la rotation en changeant l'orientation des plaques et leurs vang
	leftbox->state->ori	= qcorr*leftbox->state->ori;
	leftbox->state->angVel	= Vector3r(0,0,1)*dalpha/dt;

	rightbox->state->ori	= qcorr*leftbox->state->ori;
	rightbox->state->angVel	= Vector3r(0,0,1)*dalpha/dt;

}

void KinemCTDEngine::computeAlpha()
{
	Quaternionr orientationLeftBox,orientationRightBox;
	orientationLeftBox = leftbox->state->ori;
	orientationRightBox = rightbox->state->ori;
	if(orientationLeftBox!=orientationRightBox)
	{
		cout << "WARNING !!! your lateral boxes have not the same orientation, you're not in the case of a box imagined for creating these engines" << endl;
	}
	AngleAxisr aa(angleAxisFromQuat(orientationLeftBox));
	alpha=Mathr::PI/2.0-aa.angle();		// right if the initial orientation of the body (on the beginning of the simulation) is q =(1,0,0,0) = FromAxisAngle((0,0,1),0)
}

void KinemCTDEngine::stopMovement()
{
	// annulation de la vitesse de la plaque du haut
	topbox->state->vel	=  Vector3r(0,0,0);

	// de la plaque gauche
	leftbox->state->vel	=  Vector3r(0,0,0);
	leftbox->state->angVel	=  Vector3r(0,0,0);

	// de la plaque droite
	rightbox->state->vel	=  Vector3r(0,0,0);
	rightbox->state->angVel	=  Vector3r(0,0,0);
}





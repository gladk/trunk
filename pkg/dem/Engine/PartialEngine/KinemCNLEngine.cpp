/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<yade/pkg-dem/KinemCNLEngine.hpp>

#include<yade/core/State.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
#include<yade/core/Scene.hpp>
#include <yade/lib-base/Math.hpp>

YADE_PLUGIN((KinemCNLEngine));

void KinemCNLEngine::action()
{
	if(LOG)	cout << "debut applyCondi du CNCEngine !!" << endl;
	leftbox = Body::byId(id_boxleft);
	rightbox = Body::byId(id_boxright);
	frontbox = Body::byId(id_boxfront);
	backbox = Body::byId(id_boxback);
	topbox = Body::byId(id_topbox);
	boxbas = Body::byId(id_boxbas);
	
	if(LOG)	cout << "gamma = " << lexical_cast<string>(gamma) << "  et gammalim = " << lexical_cast<string>(gammalim) << endl;
	if(gamma<=gammalim)
	{
		if(LOG)	cout << "Je suis bien dans la partie gamma < gammalim" << endl;
		if(temoin==0)

		{
			if(LOG) cout << "Je veux maintenir la Force a F_0 = : " << F_0 << endl; 
			temoin=1;
		}
		letMove();

	}
	else if (temoin<2)
	{
		stopMovement();		// INDISPENSABLE !
		it_stop=Omega::instance().getCurrentIteration();
		cout << "Cisaillement arrete : gammaLim atteint a l'iteration "<< it_stop << endl;
		temoin=2;
	}
	else if (temoin==2 && (Omega::instance().getCurrentIteration()==(it_stop+5000)) )
	{
		Omega::instance().saveSimulation(Key + "finCis" +lexical_cast<string> ( Omega::instance().getCurrentIteration() ) + ".xml");
		Omega::instance().stopSimulationLoop();
	}

	for(unsigned int j=0;j<gamma_save.size();j++)
	{
		if ((gamma > gamma_save[j]) && (temoin_save[j]==0))
		{
			stopMovement();		// reset of all the speeds before the save
			Omega::instance().saveSimulation(Key+"_"+lexical_cast<string> (floor(gamma*1000)) +"_" +lexical_cast<string> (floor(gamma*10000)-10*floor(gamma*1000))+ "mmsheared.xml");
			temoin_save[j]=1;
		}
	}

}

void KinemCNLEngine::letMove()
{
	shared_ptr<BodyContainer> bodies = scene->bodies;

	if(LOG)	cout << "It : " << Omega::instance().getCurrentIteration() << endl;
	computeDu();

	const Real& dt = scene->dt;
	Real dx = shearSpeed * dt;


	Real Ysup = topbox->state->pos.y();
	Real Ylat = leftbox->state->pos.y();

// 	Changes in vertical and horizontal position :

	
	topbox->state->pos += Vector3r(dx,deltaU,0);

	leftbox->state->pos += Vector3r(dx/2.0,deltaU/2.0,0);
	rightbox->state->pos += Vector3r(dx/2.0,deltaU/2.0,0);
	if(LOG)	cout << "deltaU reellemt applique :" << deltaU << endl;
	if(LOG)	cout << "qui nous a emmene en : y = " <<(topbox->state->pos).y() << endl;
	
	Real Ysup_mod = topbox->state->pos.y();
	Real Ylat_mod = leftbox->state->pos.y();

//	with the corresponding velocities :
	topbox->state->vel = Vector3r(shearSpeed,deltaU/dt,0);
	leftbox->state->vel = Vector3r(shearSpeed/2.0,deltaU/(2.0*dt),0);
	rightbox->state->vel = Vector3r(shearSpeed/2.0,deltaU/(2.0*dt),0);

//	Then computation of the angle of the rotation to be done :
	computeAlpha();
	if (alpha == Mathr::PI/2.0)	// Case of the very beginning
	{
		dalpha = - Mathr::ATan( dx / (Ysup_mod -Ylat_mod) );
	}
	else
	{
		Real A = (Ysup_mod - Ylat_mod) * 2.0*Mathr::Tan(alpha) / (2.0*(Ysup - Ylat) + dx*Mathr::Tan(alpha) );
		dalpha = Mathr::ATan( (A - Mathr::Tan(alpha))/(1.0 + A * Mathr::Tan(alpha)));
	}

	Quaternionr qcorr(AngleAxisr(dalpha,Vector3r::UnitZ()));

// On applique la rotation en changeant l'orientation des plaques, leurs vang et en affectant donc alpha
	leftbox->state->ori	= qcorr*leftbox->state->ori;
	leftbox->state->angVel	= Vector3r(0,0,1)*dalpha/dt;

	rightbox->state->ori	= qcorr*rightbox->state->ori;
	rightbox->state->angVel	= Vector3r(0,0,1)*dalpha/dt;

	gamma+=dx;
}
	

void KinemCNLEngine::computeAlpha()
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


void KinemCNLEngine::computeDu()
{
	scene->forces.sync(); Vector3r F_sup=scene->forces.getForce(id_topbox);
	
	if(firstRun)
	{
		if ( !myLdc )	// FIXME : reenable this feature since Law2 is a functor
		{
// 			vector<shared_ptr<Engine> >::iterator itFirst = scene->engines.begin();
// 			vector<shared_ptr<Engine> >::iterator itLast = scene->engines.end();
// 			for ( ;itFirst!=itLast; ++itFirst )
// 			{
// 				if ( ( *itFirst )->getClassName() == "Law2_ScGeom_NormalInelasticityPhys_NormalInelasticity" ) 
// 				{
// 					myLdc =  YADE_PTR_CAST<Law2_ScGeom_NormalInelasticityPhys_NormalInelasticity> ( *itFirst );
// 					coeff_dech = myLdc ->coeff_dech;
// 					if(LOG) cout << "Law2_ScGeom_NormalInelasticityPhys_NormalInelasticity engine found, with coeff_dech = " << coeff_dech << endl;
// 				}
// 			}
		}
		
		alpha=Mathr::PI/2.0;;
		F_0 = F_sup.y();
		cout << "F_0 initialise à : " << F_0 << endl;
		firstRun=false;
	}
	
	computeStiffness();
	if( (stiffness==0) )
	{
		deltaU=0;
		cerr << "stiffness(sample) = 0 => DNC en fait <=> à CNC..." << endl;
	}
	else
	{
		Real Ycourant = topbox->state->pos.y();
		deltaU = ( F_sup.y() - F_0 )/(stiffness);
		if(LOG) cout << "Lors du calcul de DU (utile pour deltaU) : F_0 = " << F_0 << "; Ycourant = " << Ycourant << endl;
	}


	if(LOG)	cout << "deltaU a permettre normalemt :" << deltaU << endl;

// 	Il va falloir prendre en compte la loi de contact qui induit une rigidite plus grande en decharge qu'en charge
	if(deltaU>0)
	{
		deltaU/=coeff_dech;
		if(LOG) cout << "Comme deltaU>0 => decharge, il est divise par " << coeff_dech << endl;
	}

	deltaU = (1-wallDamping)*deltaU;
	if(LOG)	cout << "deltaU apres amortissement :" << deltaU << endl;
	if(abs(deltaU) > max_vel*scene->dt)
	{
		if(LOG)	cout << "v induite pour cet it n° " <<Omega::instance().getCurrentIteration()<<" : " << deltaU/scene->dt << endl;
		deltaU=deltaU/abs(deltaU)*max_vel*scene->dt;
		
		if(LOG)	cout << "Correction appliquee pour ne pas depasser vmax(comp) = " << max_vel << endl;
	}

}

void KinemCNLEngine::stopMovement()
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

void KinemCNLEngine::computeStiffness()
{
	int nbre_contacts = 0;
	stiffness=0.0;
	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for(  ; ii!=iiEnd ; ++ii )
	{
		if ((*ii)->isReal())
		{
			const shared_ptr<Interaction>& contact = *ii;
			
			Real fn = (static_cast<FrictPhys*>	(contact->interactionPhysics.get()))->normalForce.norm();

			if (fn!=0)
			{
				int id1 = contact->getId1(), id2 = contact->getId2();
				if ( id_topbox==id1 || id_topbox==id2 )
					{
						FrictPhys* currentContactPhysics =
						static_cast<FrictPhys*> ( contact->interactionPhysics.get() );
						stiffness  += currentContactPhysics->kn;
						nbre_contacts += 1;
					}
			}
		}
	}
	if(LOG)	cout << "nbre billes en contacts : " << nbre_contacts << endl;
	if(LOG)	cout << "rigidite echantillon calculee : " << stiffness << endl;

}





/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include <yade/pkg-common/ParticleParameters.hpp>
 
#include<yade/pkg-dem/BodyMacroParameters.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
#include<yade/pkg-dem/SDECLinkPhysics.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>

#include <yade/pkg-common/Sphere.hpp>
#include "VolumicContactLaw.hpp"
//#include<yade/extra/TesselationWrapper.h>
#include<yade/pkg-dem/TesselationWrapper.hpp>
#include <time.h>

YADE_REQUIRE_FEATURE(NONSENSE__THIS_PLUGIN_WILL_NEVER_BE_BUILT_UNLESS_YOU_REMOVE_ABSOLUTE_PATHS_FROM_INCLUDES!!)


//  ../../YADE/trunk/examples/sphere_15.txt
//  ../../YADE/trunk/examples/sphere_33_mono.txt
//../../YADE/trunk/examples/sphere_100_poly.txt


VolumicContactLaw::VolumicContactLaw() : InteractionSolver()
//, T1(new TesselationWrapper)
{
	sdecGroupMask=1;
	momentRotationLaw = true;
}




// void VolumicContactLaw::speedTest()
// {
// //BEGIN SPEED TEST
// 	shared_ptr<BodyContainer>& bodies = ncb->bodies;
// 	//TesselationWrapper T;
// 	BodyContainer::iterator biBegin    = bodies->begin();
// 	BodyContainer::iterator biEnd = bodies->end();
// 	BodyContainer::iterator bi = biBegin;
// 	for(  ; bi!=biEnd ; ++bi )
// 	{
// 		if ((*bi)->isDynamic) {//means "is it a sphere (not a wall)"
// 			const Sphere* s = YADE_CAST<Sphere*>((*bi)->shape.get());
// 			const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*>((*bi)->physicalParameters.get());
// 			T.checkMinMax(p->se3.position[0],p->se3.position[1],p->se3.position[2], s->radius);
// 		}
// 	}
// 	
// 	clock_t	T1 = clock();	
// 	
// 	for (int j=0; j<30; j++)
// 	{
// 		T.clear2();
// 		T.bounded = false;
// 		for( bi = biBegin; bi!=biEnd ; ++bi )
// 		{
// 			if ((*bi)->isDynamic) {//means "is it a sphere (not a wall)"
// 				const Sphere* s = YADE_CAST<Sphere*>((*bi)->shape.get());
// 				const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*>((*bi)->physicalParameters.get());
// 				T.insert(p->se3.position[0],p->se3.position[1],p->se3.position[2], s->radius, (*bi)->getId());
// 			}
// 		}
// 		//T.AddBoundingPlanes();	
// 	}
// 	
// 	
// 	cerr << "Bouding planes apres : time = "<< difftime(clock(), T1)/ CLOCKS_PER_SEC << endl;
// 	T1 = clock();
// 	
// 	
// 	for (int j=0; j<30; j++)
// 	{	
// 		T.clear2();
// 		T.bounded = false;
// 		//T.AddBoundingPlanes();
// 		for( bi = biBegin; bi!=biEnd ; ++bi )
// 		{
// 			if ((*bi)->isDynamic) {//means "is it a sphere (not a wall)"
// 				const Sphere* s = YADE_CAST<Sphere*>((*bi)->shape.get());
// 				const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*>((*bi)->physicalParameters.get());
// 				T.insert(p->se3.position[0],p->se3.position[1],p->se3.position[2], s->radius, (*bi)->getId());
// 			}
// 		
// 		}
// 	
// 	}	
// 	
// 	cerr << "Bouding planes avant : time = "<< difftime(clock(), T1)/ CLOCKS_PER_SEC << endl;
// //END SPEED TEST
// 
// }


void VolumicContactLaw::action()
{
	shared_ptr<BodyContainer>& bodies = ncb->bodies;

	Real dt = Omega::instance().getTimeStep();

//	speedTest(ncb);	

	//BEGIN VORONOI TESSELATION
	TesselationWrapper T1;	
	T1.clear2();

	BodyContainer::iterator bi    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();

	for(  ; bi!=biEnd ; ++bi )
	{
		if ((*bi)->isDynamic) {//means "is it a sphere (not a wall)"
			const Sphere* s = YADE_CAST<Sphere*>((*bi)->shape.get());
			const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*>((*bi)->physicalParameters.get());
			T1.insert(p->se3.position[0],p->se3.position[1],p->se3.position[2], s->radius, (*bi)->getId());
	
		}
	}
// 	T1.ComputeVolumes();

 

// ====================================================================================================== JF
// 	Real bottom= 0, top=1, front=5 ,back=4 ,left=2 ,right=3 ;
		PhysicalParameters* p=static_cast<PhysicalParameters*> ( Body::byId ( 0 )->physicalParameters.get() );
	Real pminy = p->se3.position[1];

		PhysicalParameters* p1=static_cast<PhysicalParameters*> ( Body::byId ( 1 )->physicalParameters.get() );
	Real pmaxy = p1->se3.position[1];

		PhysicalParameters* p2=static_cast<PhysicalParameters*> ( Body::byId ( 4 )->physicalParameters.get() );
	Real pminz = p2->se3.position[2];

		PhysicalParameters* p3=static_cast<PhysicalParameters*> ( Body::byId ( 5 )->physicalParameters.get() );
	Real pmaxz = p3->se3.position[2];

		PhysicalParameters* p4=static_cast<PhysicalParameters*> ( Body::byId ( 2 )->physicalParameters.get() );
	Real pminx = p4->se3.position[0];

		PhysicalParameters* p5=static_cast<PhysicalParameters*> ( Body::byId ( 3 )->physicalParameters.get() );
	Real pmaxx = p5->se3.position[0];


 

	T1.ComputeVolumes(pminx,pmaxx,pminy,pmaxy,pminz,pmaxz,Omega::instance().getCurrentIteration());


	if (Omega::instance().getCurrentIteration() <= 0)
	{
 	BodyContainer::iterator bi1 = bodies->begin();
	int b2 = 5;
// 	les parois 
	compacite_init.push_back(0);//Le vecteur compacite_init intègre la compacité initiale dans les cellules de Voronoï
	compacite_init.push_back(0);// 	les parois 
	compacite_init.push_back(0);// 	les parois 
	compacite_init.push_back(0);// 	les parois 
	compacite_init.push_back(0);// 	les parois 
	compacite_init.push_back(0);// 	les parois 
	compacity_init = 0;
	compacity = 0;
	compacity_max = 0;
	interpenetration_max = 0;
	delta_interpenetration_max = 0;

// 		Real information =0;

	    for( ; bi1!=biEnd ; ++bi1 )
	    {
 	     if ((*bi1)->isDynamic) 
		{
		b2 = b2+1;
//		cout << "Itération = " << Omega::instance().getCurrentIteration() << endl;
		const Sphere* s = YADE_CAST<Sphere*>((*bi1)->shape.get());
 		Real vol_sphere = 4.18 * (pow(s->radius,3));//Volume de la sphère

	        compacite_init.push_back ((vol_sphere/T1.Volume((*bi1)->getId()))*0.999999);//La compacité dans la cellule de voronoï
// 		information = information + (vol_sphere/T1.Volume((*bi1)->getId()));

// 		cout << "compacite cellule id " << b2 << "   =  " << compacite_init[b2] << endl;
//  		cout << "numéro de sphere = " << (*bi1)->getId() << endl;
// 		cout << "volume de la cellule   =    " << T1.Volume((*bi1)->getId()) << endl;
// 		cout << "volume de la sphere   =    " << vol_sphere << endl;
// 		cout << " " <<endl;
 		}
	    }
// 	cout << "information   =          " << information <<endl;
	} 
	//ENDOF VORONOI TESSELATION
		

	
/// Non Permanents Links												///

	InteractionContainer::iterator ii    = ncb->interactions->begin();
	InteractionContainer::iterator iiEnd = ncb->interactions->end();
	for(  ; ii!=iiEnd ; ++ii )
	{
		if ((*ii)->isReal())
		{
			const shared_ptr<Interaction>& contact = *ii;
			int id1 = contact->getId1();
			int id2 = contact->getId2();
			
			if( !( (*bodies)[id1]->getGroupMask() & (*bodies)[id2]->getGroupMask() & sdecGroupMask) ) continue;

			const ScGeom*    currentContactGeometry= YADE_CAST<ScGeom*>(contact->interactionGeometry.get());
			FrictPhys* currentContactPhysics = YADE_CAST<FrictPhys*> (contact->interactionPhysics.get());
			if((!currentContactGeometry)||(!currentContactPhysics)) continue;
	
			BodyMacroParameters* de1 				= YADE_CAST<BodyMacroParameters*>((*bodies)[id1]->physicalParameters.get());
			BodyMacroParameters* de2 				= YADE_CAST<BodyMacroParameters*>((*bodies)[id2]->physicalParameters.get());
			
			Vector3r& shearForce 			= currentContactPhysics->shearForce;
	
// ============================================================================ 



		if (contact->isFresh(ncb)) // Si on a une nouvelle interpenetration detectée, on initialise 
			{
			shearForce			= Vector3r(0,0,0);
			currentContactPhysics->previousun=0.0;
			currentContactPhysics->previousFn=0.0;
			currentContactPhysics->unMax=0.0;
			}


		Real previousun = currentContactPhysics->previousun;
		Real previousFn = currentContactPhysics->previousFn;
//		Real kn = currentContactPhysics->kn;
		Real Fn; // la valeur de Fn qui va etre calculee selon différentes manieres puis affectee


		Real un=currentContactGeometry->penetrationDepth;
			
// ============================================================================ 

		if(un > currentContactPhysics->unMax)	// on est en charge, et au delà et sur la "droite principale"
			{
//			Fn = kn*un;
			currentContactPhysics->unMax = std::abs(un);
// 			cout << "je suis dans le calcul normal " << endl;
			}
//		else  On fait de la décharge
// 			{
// 			Fn = previousFn + coeff_dech * kn * (un-previousun);	// Calcul incrémental de la nvlle force
// // 			cout << "je suis dans l'autre calcul" << endl;
// 			if(std::abs(Fn) > std::abs(kn * un))		// verif qu'on ne depasse la courbe
// 				Fn = kn*un;
// 			if(Fn < 0.0 )	// verif qu'on reste positif FIXME ATTENTION ON S'EST FICHU DU NORMAL ADHESION !!!!
// 				{Fn = 0;
// // 				cout << "j'ai corrige pour ne pas etre negatif" << endl;
// 				}
// 			}
//   ==============================================================================
					

// 	 		cout << "Un_sphere_sphere= " << un << endl;
			
// 			currentContactPhysics->normalForce=currentContactPhysics->kn*std::max(un,(Real) 0) *currentContactGeometry->normal;

	
									
	//BEGIN USE VOLUME
	Real vol1 = ((*bodies)[id1]->isDynamic ? T1.Volume(id1) : 0);
	Real vol2 = ((*bodies)[id2]->isDynamic ? T1.Volume(id2) : 0);
	Real solid1 = ((*bodies)[id1]->isDynamic ? 4.18 * (pow(currentContactGeometry->radius1,3)) : 0);
	Real solid2 = ((*bodies)[id2]->isDynamic ? 4.18 * (pow(currentContactGeometry->radius2,3)) : 0);
	Real compacity = (solid1+solid2)/(vol1+vol2);

	Real delta_interp = (un-previousun); 


			
	Real compacity1 = (solid1/vol1);
	Real compacity2 = (solid2/vol2);
	Real rigidite = 20.5;
	Real ecrouissage = 0.24;
	Real beta2 = 2;
	Real r1 = currentContactGeometry->radius1;
	Real r2 = currentContactGeometry->radius2;
	Real rmoy= 2*(r1*r2)/(r1+r2);
	Real rigidite_max = (1*10000*rmoy*pow(3.14/6,1/3))/(3*(1-2*0.43));


  
	if ((*bodies)[id1]->isDynamic && (*bodies)[id2]->isDynamic)
	{


	Real compacite_init_mean=(compacite_init[id2]+compacite_init[id1])/2;


  	Real k12 = 15*(1-(ecrouissage/2))*(pow(compacity-compacite_init_mean,beta2))/(1-compacity);

 	Real k111 = (0.97 - 0.58 * ecrouissage) * exp (1.75 * (un/rmoy) *(1+(ecrouissage/2)));
 	Real k121 = (-4+15/(1+3*ecrouissage))*exp(-8*(1*un/rmoy));
 	Real k1111 = rmoy*rigidite*(k111 + k121 + k12);

	
		if (k1111>rigidite_max || k1111<0 ) k1111=rigidite_max;


 	 Real k22 = 15*(1-(ecrouissage/2))*(pow(compacity-compacite_init_mean,beta2))/(1-compacity);

 	 k111 = (0.97-0.58*ecrouissage)*exp(1.75*(1*un/rmoy)*(1+(ecrouissage/2)));
 	 k121 = (-4+15/(1+3*ecrouissage))*exp(-8*(1*un/rmoy));
 	 Real k2222 = rmoy*rigidite*(k111 + k121 +k22 );
// k111 + k121 + k22);


		if (k2222>rigidite_max || k2222<0) k2222=rigidite_max;


		if (k1111<=k2222) 
		{
 	 	Fn = (k2222 * delta_interp) + previousFn;

		currentContactPhysics->normalForce	= Fn*currentContactGeometry->normal;
		currentContactPhysics->previousFn = Fn;
		currentContactPhysics->previousun = un;

		}
		else if	(k1111>k2222) 
		{

 	 	Fn = (k1111 * delta_interp) + previousFn;

		currentContactPhysics->normalForce	= Fn*currentContactGeometry->normal;
		currentContactPhysics->previousFn = Fn;
		currentContactPhysics->previousun = un;
		}

	}
	else if ((*bodies)[id1]->isDynamic)
	{

	Real rmoy= 1*(r1*r2)/(r1+r2);	
  	rigidite_max = (1*10000*rmoy*pow(3.14/6,1/3))/(3*(1-2*0.43));

  	Real k12 = 15*(1-(ecrouissage/2))*(pow(compacity1-compacite_init[id1],beta2))/(1-compacity1);
	
//  Real k1 = k11 * k12;* (un-previousun)

 	Real k111 = (0.97-0.58*ecrouissage)*exp(1.75*(2*un/rmoy)*(1+(ecrouissage/2)));
 	Real k121 = (-4+15/(1+3*ecrouissage))*exp(-8*(2*un/rmoy));
	Real k1111 = rmoy*1*rigidite*(k111 + k121 + k12);
//k111 + k121 + k12);

	if (k1111>rigidite_max || k1111<0 ) k1111=rigidite_max;


// 	if (k1111>rigidite_max ){
// //	 cout << "SUPERIEUR A LA RIGIDITE DU MATERIAU SUPERIEUR A LA RIGIDITE DU MATERIAU" << endl;
// 	k1111=rigidite_max;
// 	}


 	Fn = (k1111 * delta_interp) + previousFn;

	currentContactPhysics->normalForce = Fn*currentContactGeometry->normal;
	currentContactPhysics->previousFn = Fn;
	currentContactPhysics->previousun = un;


	}
	else if ((*bodies)[id2]->isDynamic)
	{

	Real rmoy= 1*(r1*r2)/(r1+r2);
  	Real k22 = 15*(1-(ecrouissage/2))*(pow(compacity2-compacite_init[id2],beta2))/(1-compacity2);


	rigidite_max = (1*10000*rmoy*pow(3.14/6,1/3))/(3*(1-2*0.43));
		

 	Real k212 = (0.97-0.58*ecrouissage)*exp(1.75*(2*un/rmoy)*(1+(ecrouissage/2)));
	 Real k222 = (-4+15/(1+3*ecrouissage))*exp(-8*(2*un/rmoy));
 	Real k2222 = rmoy*1*rigidite*(k212+ k222 + k22);


	if (k2222>rigidite_max || k2222<0 ) k2222=rigidite_max;


	Fn = (k2222  *1*delta_interp) + previousFn;

	currentContactPhysics->normalForce	= Fn*currentContactGeometry->normal;
	currentContactPhysics->previousFn = Fn;
	currentContactPhysics->previousun = un;

	}
	else
	{

	cout << "WARNING" <<endl;
// 	cout << "numero id1 " << id1 << endl;
// 	cout << "numero id2 " << id2 << endl;
// 	currentContactPhysics->normalForce=currentContactPhysics->kn*std::max(un,(Real) 0) *currentContactGeometry->normal //loi normale
//  * max((Real) 1,pow(compacity*1.42, 3)); // + effet volumique

 	}
 

// 			currentContactPhysics->normalForce=currentContactPhysics->kn*std::max(un,(Real) 0) *currentContactGeometry->normal //loi normale
// 			* max((Real) 1,pow(compacity*1.42, 3)); // + effet volumique
			//END OF USE VOLUME		

			Vector3r axis;
			Real angle;
	
	/// Here is the code with approximated rotations 	 ///
			
			axis = currentContactPhysics->prevNormal.Cross(currentContactGeometry->normal);
			shearForce -= shearForce.Cross(axis);
			angle = dt*0.5*currentContactGeometry->normal.Dot(de1->angularVelocity+de2->angularVelocity);
			axis = angle*currentContactGeometry->normal;
			shearForce -= shearForce.Cross(axis);

		
	/// Here is the code with exact rotations 		 ///
	
	// 		Quaternionr q;
	//
	// 		axis					= currentContactPhysics->prevNormal.cross(currentContactGeometry->normal);
	// 		angle					= acos(currentContactGeometry->normal.dot(currentContactPhysics->prevNormal));
	// 		q.fromAngleAxis(angle,axis);
	//
	// 		currentContactPhysics->shearForce	= currentContactPhysics->shearForce*q;
	//
	// 		angle					= dt*0.5*currentContactGeometry->normal.dot(de1->angularVelocity+de2->angularVelocity);
	// 		axis					= currentContactGeometry->normal;
	// 		q.fromAngleAxis(angle,axis);
	// 		currentContactPhysics->shearForce	= q*currentContactPhysics->shearForce;
	
	/// 							 ///
	
			Vector3r x				= currentContactGeometry->contactPoint;
			Vector3r c1x				= (x - de1->se3.position);
			Vector3r c2x				= (x - de2->se3.position);
			Vector3r relativeVelocity		= (de2->velocity+de2->angularVelocity.Cross(c2x)) - (de1->velocity+de1->angularVelocity.Cross(c1x));
			Vector3r shearVelocity			= relativeVelocity-currentContactGeometry->normal.Dot(relativeVelocity)*currentContactGeometry->normal;
			Vector3r shearDisplacement		= shearVelocity*dt;
  			shearForce 			       -= currentContactPhysics->ks*shearDisplacement;


 
	
	// PFC3d SlipModel, is using friction angle. CoulombCriterion
			Real maxFs = currentContactPhysics->normalForce.SquaredLength() * std::pow(currentContactPhysics->tangensOfFrictionAngle,2);
			//cout << "TANGENTOfFrictionAngle     =    " << currentContactPhysics->tangensOfFrictionAngle << endl;

			if( shearForce.SquaredLength() > maxFs )
			{
				maxFs = Mathr::Sqrt(maxFs) / shearForce.Length();
				shearForce *= maxFs;
			}
	////////// PFC3d SlipModel
	
			Vector3r f				= currentContactPhysics->normalForce + shearForce;
			
			ncb->forces.addForce (id1,-f);
			ncb->forces.addForce (id2,+f);
			ncb->forces.addTorque(id1,-c1x.Cross(f));
			ncb->forces.addTorque(id2, c2x.Cross(f));
			
			currentContactPhysics->prevNormal = currentContactGeometry->normal;

		}
	}




// MGPOST
/*	

	if ( Omega::instance().getCurrentIteration() % 500 == 0 )//|| Omega::instance().getCurrentIteration() == 1 )
	{
//  		if (Omega::instance().getCurrentIteration() == 1) 
//  		{
		std::ofstream fichier ("monfichier.txt",ios::out);
//  		}

		fichier << "<?xml version='1.0'?>" << endl;
		fichier << "<mgpost mode='3D'>" <<endl;
		fichier << "<state time='0'>" <<endl;
//		std::ofstream fichier ("monfichier.txt",ios::app);
		
//		bool file_exists = std::ifstream (emma.txt); //if file does not exist, we will write colums titles
//		ofile.open(emma.txt, std::ios::app);
		//if (!file_exists) ofile<<"iteration s11 s22 s33 e11 e22 e33 unb_force porosity kineticE" << endl;
// 		if (!file_exists) ofile<<"iteration fn11 fn22 fn33 fn111 fn222 fn333 Position_right unb_f porosity Energie_Cinétique" << endl;

	BodyContainer::iterator bi    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();

	for(  ; bi!=biEnd ; ++bi )
	{

		if ((*bi)->isDynamic) {//means "is it a sphere (not a wall)"
 			const Sphere* s = YADE_CAST<Sphere*>((*bi)->shape.get());
 			const RigidBodyParameters* p = YADE_CAST<RigidBodyParameters*>((*bi)->physicalParameters.get());
// 			T1.insert(p->se3.position[0],p->se3.position[1],p->se3.position[2], s->radius, (*bi)->getId());

			Real vol = T1.Volume((*bi)->getId()) ;
			Real solid = 4.18 * (pow(s->radius,3));
	
		fichier << "<body>" <<endl;
		fichier << "<SPHER r='" <<s->radius << "'>" <<endl;
		fichier	<< "<position x='"<< p->se3.position[0] << "' y='" << p->se3.position[1] <<"' z='" << p->se3.position[2] <<"'/>" <<endl;
		fichier << "<velocity x='0.000000' y='0.000000' z='" << (solid/vol) <<"'/>" <<endl;
		fichier << "</SPHER>" <<endl;
		fichier << "</body>" <<endl;


		}

	}
// 	T1.ComputeVolumes();

		fichier << "<body>" <<endl;
		fichier << "<SPHER r='0.0'>" <<endl;
		fichier	<< "<position x='"<< p->se3.position[0] << "' y='" << p->se3.position[1] <<"' z='" << p->se3.position[2] <<"'/>" <<endl;
		fichier << "<velocity x='0.000000' y='0.000000' z='0.0'/>" <<endl;
		fichier << "</SPHER>" <<endl;
		fichier << "</body>" <<endl;

		fichier << "<body>" <<endl;
		fichier << "<SPHER r='0.0'>" <<endl;
		fichier	<< "<position x='"<< p->se3.position[0] << "' y='" << p->se3.position[1] <<"' z='" << p->se3.position[2] <<"'/>" <<endl;
		fichier << "<velocity x='0.000000' y='0.000000' z='1.0'/>" <<endl;
		fichier << "</SPHER>" <<endl;
		fichier << "</body>" <<endl;


// ====================================================================================================== JF
// 	Real bottom= 0, top=1, front=5 ,back=4 ,left=2 ,right=3 ;
		PhysicalParameters* p=static_cast<PhysicalParameters*> ( Body::byId ( 0 )->physicalParameters.get() );
	Real pminy = p->se3.position[1];

		fichier << "<body>" <<endl;
		fichier << "<WALLX pos='0.0'/>" <<endl;
		fichier << "<WALLX pos='1.0'/>" <<endl;
		fichier << "<WALLY pos='" << p->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLY pos='" << p->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLZ pos='0.0'/>" <<endl;
		fichier << "<WALLZ pos='1.0'/>" <<endl;
		fichier << "</body>" <<endl;


		PhysicalParameters* p1=static_cast<PhysicalParameters*> ( Body::byId ( 1 )->physicalParameters.get() );
	Real pmaxy = p1->se3.position[1];

		fichier << "<body>" <<endl;
		fichier << "<WALLX pos='0.0'/>" <<endl;
		fichier << "<WALLX pos='1.0'/>" <<endl;
		fichier << "<WALLY pos='" << p1->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLY pos='" << p1->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLZ pos='0.0'/>" <<endl;
		fichier << "<WALLZ pos='1.0'/>" <<endl;
		fichier << "</body>" <<endl;


		PhysicalParameters* p2=static_cast<PhysicalParameters*> ( Body::byId ( 4 )->physicalParameters.get() );
	Real pminz = p2->se3.position[2];

		fichier << "<body>" <<endl;
		fichier << "<WALLX pos='0.0'/>" <<endl;
		fichier << "<WALLX pos='1.0'/>" <<endl;
		fichier << "<WALLY pos='0.0'/>" <<endl;
		fichier << "<WALLY pos='1.0'/>" <<endl;
		fichier << "<WALLZ pos='" << p2->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLZ pos='" << p2->se3.position[1] << "'/>" <<endl;
		fichier << "</body>" <<endl;



		PhysicalParameters* p3=static_cast<PhysicalParameters*> ( Body::byId ( 5 )->physicalParameters.get() );
	Real pmaxz = p3->se3.position[2];

		fichier << "<body>" <<endl;
		fichier << "<WALLX pos='0.0'/>" <<endl;
		fichier << "<WALLX pos='1.0'/>" <<endl;
		fichier << "<WALLY pos='0.0'/>" <<endl;
		fichier << "<WALLY pos='1.0'/>" <<endl;
		fichier << "<WALLZ pos='" << p3->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLZ pos='" << p3->se3.position[1] << "'/>" <<endl;
		fichier << "</body>" <<endl;



		PhysicalParameters* p4=static_cast<PhysicalParameters*> ( Body::byId ( 2 )->physicalParameters.get() );
	Real pminx = p4->se3.position[0];


		fichier << "<body>" <<endl;
		fichier << "<WALLX pos='" << p4->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLX pos='" << p4->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLZ pos='0.0'/>" <<endl;
		fichier << "<WALLZ pos='1.0'/>" <<endl;
		fichier << "<WALLY pos='0.0'/>" <<endl;
		fichier << "<WALLY pos='1.0'/>" <<endl;
		fichier << "</body>" <<endl;



		PhysicalParameters* p5=static_cast<PhysicalParameters*> ( Body::byId ( 3 )->physicalParameters.get() );
	Real pmaxx = p5->se3.position[0];

		fichier << "<body>" <<endl;
		fichier << "<WALLX pos='" << p5->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLX pos='" << p5->se3.position[1] << "'/>" <<endl;
		fichier << "<WALLZ pos='0.0'/>" <<endl;
		fichier << "<WALLZ pos='1.0'/>" <<endl;
		fichier << "<WALLY pos='0.0'/>" <<endl;
		fichier << "<WALLY pos='1.0'/>" <<endl;
		fichier << "</body>" <<endl;

		fichier << "</state>" <<endl;
		fichier << "</mgpost>" <<endl;


	}*/
//T1.clear2();
}

 
YADE_PLUGIN((VolumicContactLaw));
YADE_REQUIRE_FEATURE(PHYSPAR);


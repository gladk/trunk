/*************************************************************************
*  Copyright (C) 2008 by Jérôme DURIEZ                                   *
*  duriez@geo.hmg.inpg.fr                                                *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Ip2_2xNormalInelasticMat_NormalInelasticityPhys.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/NormalInelasticityPhys.hpp>
#include<yade/pkg-dem/NormalInelasticMat.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>


//
//
//
/// Big WHAT THE ....? this code below is duplicated THREE times due to some weird IFs !
/// need to FIXME that.
/// but from all my testing it works currently. / janek
//
//

void Ip2_2xNormalInelasticMat_NormalInelasticityPhys::go(	  const shared_ptr<Material>& b1 // NormalInelasticMat
					, const shared_ptr<Material>& b2 // NormalInelasticMat
					, const shared_ptr<Interaction>& interaction)
{
	NormalInelasticMat* sdec1 = static_cast<NormalInelasticMat*>(b1.get());
	NormalInelasticMat* sdec2 = static_cast<NormalInelasticMat*>(b2.get());
	ScGeom* geom = YADE_CAST<ScGeom*>(interaction->geom.get());
	
	
	if(geom) // so it is ScGeom  - NON PERMANENT LINK
	{
		if(!interaction->phys)
		{
//std::cerr << " isNew, id1: " << interaction->getId1() << " id2: " << interaction->getId2()  << "\n";
			interaction->phys = shared_ptr<NormalInelasticityPhys>(new NormalInelasticityPhys());
			NormalInelasticityPhys* contactPhysics = YADE_CAST<NormalInelasticityPhys*>(interaction->phys.get());

			Real Ea 	= sdec1->young;
			Real Eb 	= sdec2->young;
			Real Va 	= sdec1->poisson;
			Real Vb 	= sdec2->poisson;
			Real Da 	= geom->radius1; // FIXME - multiply by factor of sphere interaction distance (so sphere interacts at bigger range that its geometrical size)
			Real Db 	= geom->radius2; // FIXME - as above
			Real fa 	= sdec1->frictionAngle;
			Real fb 	= sdec2->frictionAngle;

			Real Kn = 2.0*Ea*Da*Eb*Db/(Ea*Da+Eb*Db);//harmonic average of two stiffnesses
			
			Real Ks = 2.0*Ea*Da*Va*Eb*Db*Vb/(Ea*Da*Va+Eb*Db*Va);//harmonic average of two stiffnesses with ks=V*kn for each sphere

			// Jean-Patrick Plassiard, Noura Belhaine, Frederic
// 			Victor Donze, "Calibration procedure for spherical
			// discrete elements using a local moemnt law".
			Real Kr = betaR*std::pow((Da+Db)/2.0,2)*Ks;

			contactPhysics->frictionAngle			= std::min(fa,fb); // FIXME - this is actually a waste of memory space
			contactPhysics->tangensOfFrictionAngle		= std::tan(contactPhysics->frictionAngle);
			contactPhysics->forMaxMoment		= 1.0*(Da+Db)/2.0;	// 1.0 corresponding to ethaR which I don't know exactly where to define as a parameter...

			// Lot of suppress here around (>) r2276. Normally not bad but ? See Ip2_2xCohFrictMat_CohFrictPhys.cpp to re-find the initial §...
			
			contactPhysics->prevNormal 			= geom->normal;
			
			contactPhysics->knLower = Kn;
			contactPhysics->kn = Kn;			
			contactPhysics->ks = Ks;
			contactPhysics->kr = Kr;
		}
		
	}

};
YADE_PLUGIN((Ip2_2xNormalInelasticMat_NormalInelasticityPhys));




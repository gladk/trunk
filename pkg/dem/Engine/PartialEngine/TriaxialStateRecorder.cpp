/*************************************************************************
*  Copyright (C) 2006 by luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*  Copyright (C) 2008 by Bruno Chareyre                                  *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "TriaxialStateRecorder.hpp"
// #include <yade/pkg-dem/TriaxialCompressionEngine.hpp>
#include <yade/pkg-dem/TriaxialStressController.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include <yade/core/Omega.hpp>
#include <yade/core/Scene.hpp>
#include <boost/lexical_cast.hpp>
#include <yade/pkg-dem/ScGeom.hpp>
#include <yade/pkg-dem/FrictPhys.hpp>
#include <yade/pkg-dem/Shop.hpp>

CREATE_LOGGER(TriaxialStateRecorder);
TriaxialStateRecorder::~TriaxialStateRecorder() {};

void TriaxialStateRecorder::action ()
{
	// at the beginning of the file; write column titles
	if(out.tellp()==0)	out<<"iteration s11 s22 s33 e11 e22 e33 unb_force porosity kineticE"<<endl;
	
	if ( !triaxialStressController ){
		vector<shared_ptr<Engine> >::iterator itFirst = scene->engines.begin();
		vector<shared_ptr<Engine> >::iterator itLast = scene->engines.end();
		for ( ;itFirst!=itLast; ++itFirst ){
			if ( ( *itFirst )->getClassName() == "TriaxialCompressionEngine" || ( *itFirst )->getClassName() == "ThreeDTriaxialEngine" ){
				LOG_DEBUG ( "stress controller engine found" );
				triaxialStressController =  YADE_PTR_CAST<TriaxialStressController> ( *itFirst );
				//triaxialCompressionEngine = shared_ptr<TriaxialCompressionEngine> (static_cast<TriaxialCompressionEngine*> ( (*itFirst).get()));
			}
		}
		if ( !triaxialStressController ) LOG_DEBUG ( "stress controller engine NOT found" );
	}
	if ( ! ( Omega::instance().getCurrentIteration() % triaxialStressController->computeStressStrainInterval == 0 ) )
		triaxialStressController->computeStressStrain ();

	/// Compute porosity :
	Real Vs=0;
	Real V = ( triaxialStressController->height ) * ( triaxialStressController->width ) * ( triaxialStressController->depth );
	BodyContainer::iterator bi = scene->bodies->begin();
	BodyContainer::iterator biEnd = scene->bodies->end();
	for ( ; bi!=biEnd; ++bi ){
		if((*bi)->isClump()) continue;
		const shared_ptr<Body>& b = *bi;
		if ( b->isDynamic ){
			//Sorry, the next string was commented, because it gave a Warning "unused variable v". Anton Gladky
			//const Vector3r& v = b->state->vel;
			Vs += 1.3333333*Mathr::PI*pow ( YADE_PTR_CAST<Sphere>( b->shape)->radius, 3 );}
	}
	porosity = ( V - Vs ) /V;
	
	out << lexical_cast<string> ( Omega::instance().getCurrentIteration() ) << " "
 	<< lexical_cast<string> ( triaxialStressController->stress[triaxialStressController->wall_right][0] ) << " "
 	<< lexical_cast<string> ( triaxialStressController->stress[triaxialStressController->wall_top][1] ) << " "
 	<< lexical_cast<string> ( triaxialStressController->stress[triaxialStressController->wall_front][2] ) << " "
 	<< lexical_cast<string> ( triaxialStressController->strain[0] ) << " "
 	<< lexical_cast<string> ( triaxialStressController->strain[1] ) << " "
 	<< lexical_cast<string> ( triaxialStressController->strain[2] ) << " "
 	<< lexical_cast<string> ( triaxialStressController->ComputeUnbalancedForce () ) << " "
 	<< lexical_cast<string> ( porosity ) << " "
 	<< lexical_cast<string> ( Shop::kineticEnergy() )
 	<< endl;
}

YADE_PLUGIN((TriaxialStateRecorder));

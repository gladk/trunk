/*************************************************************************
*  Copyright (C) 2008 by Bruno Chareyre                                  *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#include<yade/pkg-common/ElastMat.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>
#include <yade/pkg-common/Sphere.hpp>
#include<yade/lib-triangulation/KinematicLocalisationAnalyser.hpp>
#include<yade/lib-triangulation/TriaxialState.h>
#include <yade/lib-triangulation/Tenseur3.h>
#include<boost/iostreams/filtering_stream.hpp>
#include<boost/iostreams/filter/bzip2.hpp>
#include<boost/iostreams/device/file.hpp>
#include<yade/pkg-dem/TriaxialCompressionEngine.hpp>
#include "MicroMacroAnalyser.hpp"

YADE_PLUGIN((MicroMacroAnalyser));
YADE_REQUIRE_FEATURE(CGAL)
CREATE_LOGGER(MicroMacroAnalyser);

MicroMacroAnalyser::~MicroMacroAnalyser(){ /*delete analyser;*/} //no need, its a shared_ptr now...

void MicroMacroAnalyser::postProcessAttributes(bool deserializing)
{
	if (deserializing) {
		bool file_exists = std::ifstream(outputFile.c_str());  //if file does not exist, we will write colums titles
		ofile.open(outputFile.c_str(), std::ios::app);
		if (!file_exists) ofile<<"iteration eps1w eps2w eps3w eps11g eps22g eps33g eps12g eps13g eps23g"<< endl;
	}
}

void MicroMacroAnalyser::action()
{
	//cerr << "MicroMacroAnalyser::action() (interval="<< interval <<", iteration="<< Omega::instance().getCurrentIteration()<<")" << endl;
	if (!triaxialCompressionEngine) {
		vector<shared_ptr<Engine> >::iterator itFirst = scene->engines.begin();
		vector<shared_ptr<Engine> >::iterator itLast = scene->engines.end();
		for (;itFirst!=itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "TriaxialCompressionEngine") {
				LOG_DEBUG("stress controller engine found");
				triaxialCompressionEngine =  YADE_PTR_CAST<TriaxialCompressionEngine> (*itFirst);}
		}
		if (!triaxialCompressionEngine) LOG_ERROR("stress controller engine not found");
	}
	if (triaxialCompressionEngine->strain[0]==0) return;// no deformation yet
	if (!initialized) {
		setState(1,true,false);
		//Check file here again, to make sure we write to the correct file when filename is modified after the scene is loaded
		bool file_exists = std::ifstream(outputFile.c_str());
		ofile.open(outputFile.c_str(), std::ios::app);
		if (!file_exists) ofile<<"iteration eps1w eps2w eps3w eps11g eps22g eps33g eps12g eps13g eps23g"<< endl;
		initialized=true;
	} else if (Omega::instance().getCurrentIteration() % interval == 0) {
		setState(2, true, compIncrt);
		if (compDeformation) {
			analyser->ComputeParticlesDeformation();
			//for (int i=0; i<analyser->ParticleDeformation.size();i++) cerr<< analyser->ParticleDeformation[i]<<endl;
			ostringstream oss;
			oss<<"deformation"<<incrtNumber++<<".vtk";
			analyser->DefToFile(oss.str().c_str());
		}
		CGT::Tenseur_sym3 epsg(analyser->grad_u_total);
		ofile << Omega::instance().getCurrentIteration() << analyser->Delta_epsilon(1,1)<<" "<<analyser->Delta_epsilon(2,2)<<" "<<analyser->Delta_epsilon(3,3)<<" "<<epsg(1,1)<<" "<<epsg(2,2)<< " "<<epsg(3,3)<<" "<<epsg(1,2)<<" "<<epsg(1,3)<<" "<<epsg(2,3)<<endl;
		analyser->SwitchStates();
	}
	//cerr << "ENDOF MicroMacro::action" << endl;
}

void MicroMacroAnalyser::setState(unsigned int state, bool save_states, bool computeIncrement)
{
	LOG_INFO("MicroMacroAnalyser::setState");
	CGT::TriaxialState& TS = makeState(state);
	if (state == 2) {
		analyser->Delta_epsilon(3,3) = analyser->TS1->eps3 - analyser->TS0->eps3;
		analyser->Delta_epsilon(1,1) = analyser->TS1->eps1 - analyser->TS0->eps1;
		analyser->Delta_epsilon(2,2) = analyser->TS1->eps2 - analyser->TS0->eps2;
		if (computeIncrement) {
			analyser->SetForceIncrements();
			analyser->SetDisplacementIncrements();
		}
	}
	if (save_states) {
		ostringstream oss;
		//oss<<stateFileName<<"_"<<Omega::instance().getCurrentIteration();
		oss<<stateFileName<<"_"<<stateNumber++;
		TS.to_file(oss.str().c_str(),/*use bz2?*/ true);
	}
	LOG_DEBUG("ENDOF MicroMacroAnalyser::setState");
}

//Copy simulation data in the triaxialState structure
CGT::TriaxialState& MicroMacroAnalyser::makeState(unsigned int state, const char* filename)
{
	Scene* scene=Omega::instance().getScene().get();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	CGT::TriaxialState* ts=0;
	if (state==1) ts = analyser->TS0;
	else if (state==2) ts = analyser->TS1;
	else LOG_ERROR("state must be 1 or 2, instead of " << state);
	CGT::TriaxialState& TS = *ts;
	
	TS.reset();
	long Ng = bodies->size();
	TS.mean_radius=0;
	TS.grains.resize(Ng);
	BodyContainer::iterator biBegin    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();
	BodyContainer::iterator bi = biBegin;
	Ng = 0;
	vector<body_id_t> fictiousVtx;
	for (; bi!=biEnd ; ++bi) {
		const body_id_t Idg = (*bi)->getId();
		TS.grains[Idg].id = Idg;
		if (!(*bi)->isDynamic) {
			TS.grains[Idg].isSphere = false; fictiousVtx.push_back(Idg);}
		else {//then it is a sphere (not a wall)
			++Ng;
			const Sphere* s = YADE_CAST<Sphere*> ((*bi)->shape.get());
			//const GranularMat* p = YADE_CAST<GranularMat*> ( ( *bi )->material.get() );
			const Vector3r& pos = (*bi)->state->pos;
			Real rad = s->radius;

			TS.grains[Idg].sphere = CGT::Sphere(CGT::Point(pos[0],pos[1],pos[2]),rad);
//    TS.grains[Idg].translation = trans;
			AngleAxisr aa(angleAxisFromQuat((*bi)->state->ori));
			Vector3r rotVec=aa.axis()*aa.angle();
			TS.grains[Idg].rotation = CGT::Vecteur(rotVec[0],rotVec[1],rotVec[2]);
			TS.box.base = CGT::Point(min(TS.box.base.x(), pos.x()-rad),
					min(TS.box.base.y(), pos.y()-rad),
					    min(TS.box.base.z(), pos.z()-rad));
			TS.box.sommet = CGT::Point(max(TS.box.sommet.x(), pos.x() +rad),
					max(TS.box.sommet.y(), pos.y() +rad),
					    max(TS.box.sommet.z(), pos.z() +rad));
			TS.mean_radius += TS.grains[Idg].sphere.weight();
		}
	}
	TS.mean_radius /= Ng;//rayon moyen
	LOG_INFO(" loaded : " << Ng << " grains with mean radius = " << TS.mean_radius);
	if (fictiousVtx.size()>=6){//boxes found, simulate them with big spheres
		CGT::Point& Pmin = TS.box.base; CGT::Point& Pmax = TS.box.sommet; 
		Real FAR = 1e4;
		TS.grains[fictiousVtx[0]].sphere = CGT::Sphere(CGT::Point(0.5*(Pmin.x()+Pmax.x()),Pmin.y()-FAR*(Pmax.x()-Pmin.x()),0.5*(Pmax.z()+Pmin.z())),FAR*(Pmax.x()-Pmin.x()));
		TS.grains[fictiousVtx[1]].sphere = CGT::Sphere(CGT::Point(0.5*(Pmin.x()+Pmax.x()),Pmax.y()+FAR*(Pmax.x()-Pmin.x()),0.5*(Pmax.z()+Pmin.z())),FAR*(Pmax.x()-Pmin.x()));
		TS.grains[fictiousVtx[2]].sphere =
		CGT::Sphere(CGT::Point(Pmin.x()-FAR*(Pmax.y()-Pmin.y()),0.5*(Pmax.y()+Pmin.y()),0.5*(Pmax.z()+Pmin.z())),FAR*(Pmax.y()-Pmin.y()));
		TS.grains[fictiousVtx[3]].sphere =
		CGT::Sphere(CGT::Point(Pmax.x()+FAR*(Pmax.y()-Pmin.y()),0.5*(Pmax.y()+Pmin.y()),0.5*(Pmax.z()+Pmin.z())),FAR*(Pmax.y()-Pmin.y()));
		TS.grains[fictiousVtx[4]].sphere =
		CGT::Sphere(CGT::Point(0.5*(Pmin.x()+Pmax.x()),0.5*(Pmax.y()+Pmin.y()),Pmin.z()-FAR*(Pmax.y()-Pmin.y())),FAR*(Pmax.y()-Pmin.y()));
		TS.grains[fictiousVtx[5]].sphere =
		CGT::Sphere(CGT::Point(0.5*(Pmin.x()+Pmax.x()),0.5*(Pmax.y()+Pmin.y()),Pmax.z()+FAR*(Pmax.y()-Pmin.y())),FAR*(Pmax.y()-Pmin.y()));
	}
	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for (; ii!=iiEnd ; ++ii) {
		if ((*ii)->isReal()) {
			CGT::TriaxialState::Contact *c = new CGT::TriaxialState::Contact;
			TS.contacts.push_back(c);
			CGT::TriaxialState::VectorGrain& grains = TS.grains;
			body_id_t id1 = (*ii)->getId1();
			body_id_t id2 = (*ii)->getId2();

			c->grain1 = & (TS.grains[id1]);
			c->grain2 = & (TS.grains[id2]);
			grains[id1].contacts.push_back(c);
			grains[id2].contacts.push_back(c);
			c->normal = CGT::Vecteur((YADE_CAST<ScGeom*> ((*ii)->interactionGeometry.get()))->normal.x(),
					(YADE_CAST<ScGeom*> ((*ii)->interactionGeometry.get()))->normal.y(),
					 (YADE_CAST<ScGeom*> ((*ii)->interactionGeometry.get()))->normal.z());
//    c->normal = ( grains[id2].sphere.point()-grains[id1].sphere.point() );
//    c->normal = c->normal/sqrt ( pow ( c->normal.x(),2 ) +pow ( c->normal.y(),2 ) +pow ( c->normal.z(),2 ) );
			c->position = CGT::Vecteur((YADE_CAST<ScGeom*> ((*ii)->interactionGeometry.get()))->contactPoint.x(),
					(YADE_CAST<ScGeom*> ((*ii)->interactionGeometry.get()))->contactPoint.y(),
					 (YADE_CAST<ScGeom*> ((*ii)->interactionGeometry.get()))->contactPoint.z());
//    c->position = 0.5* ( ( grains[id1].sphere.point()-CGAL::ORIGIN ) +
//          ( grains[id1].sphere.weight() *c->normal ) +
//          ( grains[id2].sphere.point()-CGAL::ORIGIN ) -
//          ( grains[id2].sphere.weight() *c->normal ) );
			c->fn = YADE_CAST<FrictPhys*> (((*ii)->interactionPhysics.get()))->normalForce.dot((YADE_CAST<ScGeom*> ((*ii)->interactionGeometry.get()))->normal);
			Vector3r fs = YADE_CAST<FrictPhys*> ((*ii)->interactionPhysics.get())->shearForce;
			c->fs = CGT::Vecteur(fs.x(),fs.y(),fs.z());
			c->old_fn = c->fn;
			c->old_fs = c->fs;
			c->frictional_work = 0;
		}
	}
	//Save various parameters if triaxialCompressionEngine is defined
	if (!triaxialCompressionEngine) {
		vector<shared_ptr<Engine> >::iterator itFirst = scene->engines.begin();
		vector<shared_ptr<Engine> >::iterator itLast = scene->engines.end();
		for (;itFirst!=itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "TriaxialCompressionEngine") {
				LOG_DEBUG("stress controller engine found");
				triaxialCompressionEngine =  YADE_PTR_CAST<TriaxialCompressionEngine> (*itFirst);}
		}
		if (!triaxialCompressionEngine) LOG_ERROR("stress controller engine not found");}
	
	if (triaxialCompressionEngine) {
		TS.wszzh = triaxialCompressionEngine->stress[triaxialCompressionEngine->wall_top][1];
		TS.wsxxd = triaxialCompressionEngine->stress[triaxialCompressionEngine->wall_right][0];
		TS.wsyyfa = triaxialCompressionEngine->stress[triaxialCompressionEngine->wall_front][2];
		TS.eps3 = triaxialCompressionEngine->strain[2];//find_parameter("eps3=", Statefile);
		TS.eps1 = triaxialCompressionEngine->strain[0];//find_parameter("eps1=", Statefile);
		TS.eps2 = triaxialCompressionEngine->strain[1];//find_parameter("eps2=", Statefile);
		TS.haut = triaxialCompressionEngine->height;//find_parameter("haut=", Statefile);
		TS.larg = triaxialCompressionEngine->width;//find_parameter("larg=", Statefile);
		TS.prof = triaxialCompressionEngine->depth;//find_parameter("prof=", Statefile);
		TS.porom = 0/*analyser->ComputeMacroPorosity() crasher?*/;//find_parameter("porom=", Statefile);
		TS.ratio_f = triaxialCompressionEngine-> ComputeUnbalancedForce(scene);  //find_parameter("ratio_f=", Statefile);
	}
	if (filename!=NULL) TS.to_file(filename);
	return TS;
}

// const vector<CGT::Tenseur3>& MicroMacroAnalyser::makeDeformationArray(const char* state_file1, const char* state_file0)
// {
// 	return analyser->ComputeParticlesDeformation(state_file1, state_file0);
// }

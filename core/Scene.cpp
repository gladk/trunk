/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Scene.hpp"
#include<yade/core/Engine.hpp>
#include<yade/core/Timing.hpp>
#include<yade/core/TimeStepper.hpp>

#include<yade/lib-base/Math.hpp>
#include<boost/foreach.hpp>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<boost/algorithm/string.hpp>

#include<yade/core/BodyContainer.hpp>
#include<yade/core/InteractionContainer.hpp>


// POSIX-only
#include<pwd.h>
#include<unistd.h>
#include<time.h>


YADE_PLUGIN((Scene));
CREATE_LOGGER(Scene);
// should be elsewhere, probably
bool TimingInfo::enabled=false;

void Scene::fillDefaultTags(){
	// fill default tags
	struct passwd* pw;
	char hostname[HOST_NAME_MAX];
	gethostname(hostname,HOST_NAME_MAX);
	pw=getpwuid(geteuid()); if(!pw) throw runtime_error("getpwuid(geteuid()) failed!");
	// a few default tags
	// real name: will have all non-ASCII characters replaced by ? since serialization doesn't handle that
	// the standard GECOS format is Real Name,,, - first comma and after will be discarded
	string gecos(pw->pw_gecos), gecos2; size_t p=gecos.find(","); if(p!=string::npos) boost::algorithm::erase_tail(gecos,gecos.size()-p); for(size_t i=0;i<gecos.size();i++){gecos2.push_back(((unsigned char)gecos[i])<128 ? gecos[i] : '?'); }
	tags.push_back(boost::algorithm::replace_all_copy(string("author=")+gecos2+" ("+string(pw->pw_name)+"@"+hostname+")"," ","~"));
	tags.push_back(string("isoTime="+boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time())));
	string id=boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time())+"p"+lexical_cast<string>(getpid());
	tags.push_back("id="+id);
	tags.push_back("d.id="+id);
	tags.push_back("id.d="+id);
}



void Scene::postLoad(Scene&){
	// FIXME: this should be no longer necessary with boost::serialization, but it must be checked
	/* since yade::serialization doesn't properly handle shared pointers, iterate over all bodies and make materials shared again, if id>=0 */
	FOREACH(const shared_ptr<Body>& b, *bodies){
		if(!b) continue; // erased body
		if(!b->material || b->material->id<0) continue; // not a shared material
		assert(b->material->id < (int)materials.size());
		b->material=materials[b->material->id];
	}
}



void Scene::moveToNextTimeStep(){
	if(needsInitializers){
		checkStateTypes();
		FOREACH(shared_ptr<Engine> e, initializers){ e->scene=this; if(e->dead || !e->isActivated()) continue; e->action(); } 
		needsInitializers=false;
	}
	if(!subStepping && subStep<0){
		// ** 1. ** prologue
		if(isPeriodic) cell->integrateAndUpdate(dt);
		//forces.reset(); // uncomment if ForceResetter is removed
		bool TimingInfo_enabled=TimingInfo::enabled; // cache the value, so that when it is changed inside the step, the engine that was just running doesn't get bogus values
		TimingInfo::delta last=TimingInfo::getNow(); // actually does something only if TimingInfo::enabled, no need to put the condition here
		// ** 2. ** engines
		FOREACH(const shared_ptr<Engine>& e, engines){
			e->scene=this;
			if(e->dead || !e->isActivated()) continue;
			e->action();
			if(TimingInfo_enabled) {TimingInfo::delta now=TimingInfo::getNow(); e->timingInfo.nsec+=now-last; e->timingInfo.nExec+=1; last=now;}
		}
		// ** 3. ** epilogue
		iter++;
		time+=dt;
	} else {
		/* IMPORTANT: take care to copy EXACTLY the same sequence as is in the block above !! */
		if(TimingInfo::enabled){ TimingInfo::enabled=false; LOG_INFO("O.timingEnabled disabled, since O.subStepping is used."); }
		if(subStep<-1 || subStep>(int)engines.size()){ LOG_ERROR("Invalid value of Scene::subStep ("<<subStep<<"), setting to -1 (prologue will be run)."); subStep=-1; }
		// if subStepping is disabled, it means we have not yet finished last step completely; in that case, do that here by running all remaining substeps at once
		// if subStepping is enabled, just run the step we need (the loop is traversed only once, with subs==subStep)
		int maxSubStep=subStep;
		if(!subStepping){ maxSubStep=engines.size(); LOG_INFO("Running remaining sub-steps ("<<subStep<<"…"<<maxSubStep<<") before disabling sub-stepping."); }
		for(int subs=subStep; subs<=maxSubStep; subs++){
			assert(subs>=-1 && subs<=(int)engines.size());
			// ** 1. ** prologue
			if(subs==-1){ if(isPeriodic) cell->integrateAndUpdate(dt); }
			// ** 2. ** engines
			else if(subs>=0 && subs<(int)engines.size()){ const shared_ptr<Engine>& e(engines[subs]); e->scene=this; if(!e->dead && e->isActivated()) e->action(); }
			// ** 3. ** epilogue
			else if(subs==(int)engines.size()){ iter++; time+=dt; /* gives -1 along with the increment afterwards */ subStep=-2; }
			// (?!)
			else { /* never reached */ assert(false); }
		}
		subStep++; // if not substepping, this will make subStep=-2+1=-1, which is what we want
	}
}



shared_ptr<Engine> Scene::engineByName(const string& s){
	FOREACH(shared_ptr<Engine> e, engines){
		if(e->getClassName()==s) return e;
	}
	return shared_ptr<Engine>();
}

bool Scene::timeStepperPresent(){
	int n=0;
	FOREACH(const shared_ptr<Engine>&e, engines){ if(dynamic_cast<TimeStepper*>(e.get())) n++; }
	if(n>1) throw std::runtime_error(string("Multiple ("+lexical_cast<string>(n)+") TimeSteppers in the simulation?!").c_str());
	return n>0;
}

bool Scene::timeStepperActive(){
	int n=0; bool ret=false;
	FOREACH(const shared_ptr<Engine>&e, engines){
		TimeStepper* ts=dynamic_cast<TimeStepper*>(e.get()); if(ts) { ret=ts->active; n++; }
	}
	if(n>1) throw std::runtime_error(string("Multiple ("+lexical_cast<string>(n)+") TimeSteppers in the simulation?!").c_str());
	return ret;
}

bool Scene::timeStepperActivate(bool a){
	int n=0;
	FOREACH(const shared_ptr<Engine> e, engines){
		TimeStepper* ts=dynamic_cast<TimeStepper*>(e.get());
		if(ts) { ts->setActive(a); n++; }
	}
	if(n>1) throw std::runtime_error(string("Multiple ("+lexical_cast<string>(n)+") TimeSteppers in the simulation?!").c_str());
	return n>0;
}



void Scene::checkStateTypes(){
	FOREACH(const shared_ptr<Body>& b, *bodies){
		if(!b || !b->material) continue;
		if(b->material && !b->state) throw std::runtime_error("Body #"+lexical_cast<string>(b->getId())+": has Body::material, but NULL Body::state.");
		if(!b->material->stateTypeOk(b->state.get())){
			throw std::runtime_error("Body #"+lexical_cast<string>(b->getId())+": Body::material type "+b->material->getClassName()+" doesn't correspond to Body::state type "+b->state->getClassName()+" (should be "+b->material->newAssocState()->getClassName()+" instead).");
		}
	}
}

void Scene::updateBound(){
	if(!bound) bound=shared_ptr<Bound>(new Bound);
	const Real& inf=std::numeric_limits<Real>::infinity();
	Vector3r mx(-inf,-inf,-inf);
	Vector3r mn(inf,inf,inf);
	FOREACH(const shared_ptr<Body>& b, *bodies){
		if(!b) continue;
		if(b->bound){
			for(int i=0; i<3; i++){
				if(!isinf(b->bound->max[i])) mx[i]=max(mx[i],b->bound->max[i]);
				if(!isinf(b->bound->min[i])) mn[i]=min(mn[i],b->bound->min[i]);
			}
		} else {
	 		mx=mx.cwise().max(b->state->pos);
 			mn=mn.cwise().min(b->state->pos);
		}
	}
	bound->min=mn; bound->max=mx;
}

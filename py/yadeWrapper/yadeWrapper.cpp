// 2007,2008 © Václav Šmilauer <eudoxos@arcig.cz> 

#include<sstream>
#include<map>
#include<vector>
#include<unistd.h>
#include<list>
#include<signal.h>

#include<boost/python.hpp>
#include<boost/python/raw_function.hpp>
#include<boost/python/suite/indexing/vector_indexing_suite.hpp>
#include<boost/bind.hpp>
#include<boost/lambda/bind.hpp>
#include<boost/thread/thread.hpp>
#include<boost/filesystem/operations.hpp>
#include<boost/date_time/posix_time/posix_time.hpp>
#include<boost/any.hpp>
#include<boost/python.hpp>
#include<boost/foreach.hpp>
#include<boost/algorithm/string.hpp>
#include<boost/version.hpp>

#include<yade/lib-base/Logging.hpp>
#include<yade/lib-serialization-xml/XMLFormatManager.hpp>
#include<yade/lib-pyutil/gil.hpp>
#include<yade/lib-pyutil/raw_constructor.hpp>
#include<yade/lib-pyutil/doc_opts.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/ThreadRunner.hpp>
#include<yade/core/FileGenerator.hpp>

#include<yade/pkg-dem/STLImporter.hpp>

#include<yade/core/Dispatcher.hpp>
#include<yade/core/GlobalEngine.hpp>
#include<yade/core/PartialEngine.hpp>
#include<yade/core/Functor.hpp>
#include<yade/pkg-common/ParallelEngine.hpp>

#include<yade/pkg-common/InteractionDispatchers.hpp>

// #include<yade/pkg-dem/Shop.hpp>
#include<yade/pkg-dem/Clump.hpp>

// local copy
#include<boost/math/nonfinite_num_facets.hpp>

#include<locale>
#include<boost/archive/codecvt_null.hpp>

using namespace boost;
using namespace std;

#include<yade/lib-serialization/ObjectIO.hpp>

#include<yade/extra/boost_python_len.hpp>

#ifdef YADE_LOG4CXX
	log4cxx::LoggerPtr logger=log4cxx::Logger::getLogger("yade.python");
#endif

/*
Python normally iterates over object it is has __getitem__ and __len__, which BodyContainer does.
However, it will not skip removed bodies automatically, hence this iterator which does just that.
*/
class pyBodyIterator{
	BodyContainer::iterator I, Iend;
	public:
	pyBodyIterator(const shared_ptr<BodyContainer>& bc){ I=bc->begin(); Iend=bc->end(); }
	pyBodyIterator pyIter(){return *this;}
	shared_ptr<Body> pyNext(){
		BodyContainer::iterator ret;
		while(I!=Iend){ ret=I; ++I; if(*ret) return *ret; }
		PyErr_SetNone(PyExc_StopIteration); python::throw_error_already_set(); /* never reached, but makes the compiler happier */ throw;
	}
};

class pyBodyContainer{
	public:
	const shared_ptr<BodyContainer> proxee;
	pyBodyIterator pyIter(){return pyBodyIterator(proxee);}
	pyBodyContainer(const shared_ptr<BodyContainer>& _proxee): proxee(_proxee){}
	shared_ptr<Body> pyGetitem(body_id_t _id){
		int id=(_id>=0 ? _id : proxee->size()+_id);
		if(id<0 || (size_t)id>=proxee->size()){ PyErr_SetString(PyExc_IndexError, "Body id out of range."); python::throw_error_already_set(); /* make compiler happy; never reached */ return shared_ptr<Body>(); }
		else return (*proxee)[id];
	}
	body_id_t append(shared_ptr<Body> b){
		// shoud be >=0, but Body is by default created with id 0... :-|
		if(b->getId()>=0){ PyErr_SetString(PyExc_IndexError,("Body already has id "+lexical_cast<string>(b->getId())+" set; appending such body (for the second time) is not allowed.").c_str()); python::throw_error_already_set(); }
		return proxee->insert(b);
	}
	vector<body_id_t> appendList(vector<shared_ptr<Body> > bb){
		/* prevent crash when adding lots of bodies (not clear why it happens exactly, bt is like this:

			#3  <signal handler called>
			#4  0x000000000052483f in boost::detail::atomic_increment (pw=0x8089) at /usr/include/boost/detail/sp_counted_base_gcc_x86.hpp:66
			#5  0x00000000005248b3 in boost::detail::sp_counted_base::add_ref_copy (this=0x8081) at /usr/include/boost/detail/sp_counted_base_gcc_x86.hpp:133
			#6  0x00000000005249ca in shared_count (this=0x7fff2e44db48, r=@0x7f08ffd692b8) at /usr/include/boost/detail/shared_count.hpp:227
			#7  0x00000000005258e3 in shared_ptr (this=0x7fff2e44db40) at /usr/include/boost/shared_ptr.hpp:165
			#8  0x0000000000505cff in BodyRedirectionVectorIterator::getValue (this=0x846f040) at /home/vaclav/yade/trunk/core/containers/BodyRedirectionVector.cpp:47
			#9  0x00007f0908af41ce in BodyContainerIteratorPointer::operator* (this=0x7fff2e44db60) at /home/vaclav/yade/build-trunk/include/yade-trunk/yade/core/BodyContainer.hpp:63
			#10 0x00007f0908af420a in boost::foreach_detail_::deref<BodyContainer, mpl_::bool_<false> > (cur=@0x7fff2e44db60) at /usr/include/boost/foreach.hpp:750
			#11 0x00007f0908adc5a9 in OpenGLRenderingEngine::renderGeometricalModel (this=0x77f1240, rootBody=@0x1f49220) at pkg/common/RenderingEngine/OpenGLRenderingEngine/OpenGLRenderingEngine.cpp:441
			#12 0x00007f0908adfb84 in OpenGLRenderingEngine::render (this=0x77f1240, rootBody=@0x1f49220, selection=-1) at pkg/common/RenderingEngine/OpenGLRenderingEngine/OpenGLRenderingEngine.cpp:232

		*/
		#if BOOST_VERSION<103500
			boost::try_mutex::scoped_try_lock lock(Omega::instance().renderMutex,true); // acquire lock on the mutex (true)
		#else
			boost::mutex::scoped_lock lock(Omega::instance().renderMutex);
		#endif
		vector<body_id_t> ret; FOREACH(shared_ptr<Body>& b, bb){ret.push_back(append(b));} return ret;
	}
	python::tuple appendClump(vector<shared_ptr<Body> > bb){
		// update clump members
		vector<body_id_t> ids(appendList(bb));
		// create and add clump itself
		shared_ptr<Clump> clump=shared_ptr<Clump>(new Clump());
		shared_ptr<Body> clumpAsBody=static_pointer_cast<Body>(clump);
		clump->isDynamic=true;
		proxee->insert(clumpAsBody);
		// add clump members to the clump
		FOREACH(body_id_t id, ids) clump->add(id);
		// update clump
		clump->updateProperties(false);
		return python::make_tuple(clump->getId(),ids);
	}
	vector<body_id_t> replace(vector<shared_ptr<Body> > bb){proxee->clear(); return appendList(bb);}
	long length(){return proxee->size();}
	void clear(){proxee->clear();}
	bool erase(body_id_t id){ return proxee->erase(id); }
};


class pyTags{
	public:
		pyTags(const shared_ptr<Scene> _mb): mb(_mb){}
		const shared_ptr<Scene> mb;
		bool hasKey(const string& key){ FOREACH(string val, mb->tags){ if(algorithm::starts_with(val,key+"=")){ return true;} } return false; }
		string getItem(const string& key){
			FOREACH(string& val, mb->tags){
				if(algorithm::starts_with(val,key+"=")){ string val1(val); algorithm::erase_head(val1,key.size()+1); return val1;}
			}
			PyErr_SetString(PyExc_KeyError,("Invalid key: "+key+".").c_str());
			python::throw_error_already_set(); /* make compiler happy; never reached */ return string();
		}
		void setItem(const string& key,const string& item){
			if(key.find("=")!=string::npos) {
				PyErr_SetString(PyExc_KeyError, "Key must not contain the '=' character (implementation limitation; sorry).");
				python::throw_error_already_set();
			}
			FOREACH(string& val, mb->tags){if(algorithm::starts_with(val,key+"=")){ val=key+"="+item; return; } }
			mb->tags.push_back(key+"="+item);
			}
		python::list keys(){
			python::list ret;
			FOREACH(string val, mb->tags){
				size_t i=val.find("=");
				if(i==string::npos) throw runtime_error("Tags must be in the key=value format (internal error?)");
				algorithm::erase_tail(val,val.size()-i); ret.append(val);
			}
			return ret;
		}
};


class pyInteractionIterator{
	InteractionContainer::iterator I, Iend;
	public:
	pyInteractionIterator(const shared_ptr<InteractionContainer>& ic){ I=ic->begin(); Iend=ic->end(); }
	pyInteractionIterator pyIter(){return *this;}
	shared_ptr<Interaction> pyNext(){
		InteractionContainer::iterator ret;
		while(I!=Iend){ ret=I; ++I; if((*ret)->isReal()) return *ret; }
		PyErr_SetNone(PyExc_StopIteration); python::throw_error_already_set();
		throw; // to avoid compiler warning; never reached
		//InteractionContainer::iterator ret=I; ++I; return *ret;
	}
};

class pyInteractionContainer{
	public:
		const shared_ptr<InteractionContainer> proxee;
		pyInteractionContainer(const shared_ptr<InteractionContainer>& _proxee): proxee(_proxee){}
		pyInteractionIterator pyIter(){return pyInteractionIterator(proxee);}
		shared_ptr<Interaction> pyGetitem(vector<body_id_t> id12){
			//if(!PySequence_Check(id12.ptr())) throw invalid_argument("Key must be a tuple");
			//if(python::len(id12)!=2) throw invalid_argument("Key must be a 2-tuple: id1,id2.");
			if(id12.size()==2){
				shared_ptr<Interaction> i=proxee->find(id12[0],id12[1]);
				if(i) return i; else { PyErr_SetString(PyExc_IndexError,"No such interaction"); python::throw_error_already_set(); /* make compiler happy; never reached */ return shared_ptr<Interaction>(); }
			}
			else if(id12.size()==1){ return (*proxee)[id12[0]];}
			else throw invalid_argument("2 integers (id1,id2) or 1 integer (nth) required.");
		}
		/* return nth _real_ iteration from the container (0-based index); this is to facilitate picking random interaction */
		shared_ptr<Interaction> pyNth(long n){
			long i=0; FOREACH(shared_ptr<Interaction> I, *proxee){ if(!I->isReal()) continue; if(i++==n) return I; }
			PyErr_SetString(PyExc_IndexError,(string("Interaction number out of range (")+lexical_cast<string>(n)+">="+lexical_cast<string>(i)+").").c_str());
			python::throw_error_already_set(); /* make compiler happy; never reached */ return shared_ptr<Interaction>();
		}
		long len(){return proxee->size();}
		void clear(){proxee->clear();}
		python::list withBody(long id){ python::list ret; FOREACH(const shared_ptr<Interaction>& I, *proxee){ if(I->isReal() && (I->getId1()==id || I->getId2()==id)) ret.append(I);} return ret;}
		python::list withBodyAll(long id){ python::list ret; FOREACH(const shared_ptr<Interaction>& I, *proxee){ if(I->getId1()==id || I->getId2()==id) ret.append(I);} return ret; }
		long countReal(){ long ret=0; FOREACH(const shared_ptr<Interaction>& I, *proxee){ if(I->isReal()) ret++; } return ret; }
		bool serializeSorted_get(){return proxee->serializeSorted;}
		void serializeSorted_set(bool ss){proxee->serializeSorted=ss;}
		void eraseNonReal(){ proxee->eraseNonReal(); }
};

class pyForceContainer{
		shared_ptr<Scene> scene;
	public:
		pyForceContainer(shared_ptr<Scene> _scene): scene(_scene) { }
		void checkId(long id){ if(id<0 || (size_t)id>=scene->bodies->size()){ PyErr_SetString(PyExc_IndexError, "Body id out of range."); python::throw_error_already_set(); /* never reached */ throw; } }
		Vector3r force_get(long id){  checkId(id); scene->forces.sync(); return scene->forces.getForce(id); }
		Vector3r torque_get(long id){ checkId(id); scene->forces.sync(); return scene->forces.getTorque(id); }
		Vector3r move_get(long id){   checkId(id); scene->forces.sync(); return scene->forces.getMove(id); }
		Vector3r rot_get(long id){    checkId(id); scene->forces.sync(); return scene->forces.getRot(id); }
		void force_add(long id, const Vector3r& f){  checkId(id); scene->forces.addForce (id,f); }
		void torque_add(long id, const Vector3r& t){ checkId(id); scene->forces.addTorque(id,t);}
		void move_add(long id, const Vector3r& t){   checkId(id); scene->forces.addMove(id,t);}
		void rot_add(long id, const Vector3r& t){    checkId(id); scene->forces.addRot(id,t);}
		long syncCount_get(){ return scene->forces.syncCount;}
		void syncCount_set(long count){ scene->forces.syncCount=count;}
};

class pyMaterialContainer{
		shared_ptr<Scene> scene;
	public:
		pyMaterialContainer(shared_ptr<Scene> _scene): scene(_scene) { }
		shared_ptr<Material> getitem_id(int _id){ int id=(_id>=0 ? _id : scene->materials.size()+_id); if(id<0 || (size_t)id>=scene->materials.size()){ PyErr_SetString(PyExc_IndexError, "Material id out of range."); python::throw_error_already_set(); /* never reached */ throw; } return Material::byId(id,scene); }
		shared_ptr<Material> getitem_label(string label){
			// translate runtime_error to KeyError (instead of RuntimeError) if the material doesn't exist
			try { return Material::byLabel(label,scene);	}
			catch (std::runtime_error& e){ PyErr_SetString(PyExc_KeyError,e.what()); python::throw_error_already_set(); /* never reached; avoids warning */ throw; }
		}
		int append(shared_ptr<Material> m){ scene->materials.push_back(m); m->id=scene->materials.size()-1; return m->id; }
		vector<int> appendList(vector<shared_ptr<Material> > mm){ vector<int> ret; FOREACH(shared_ptr<Material>& m, mm) ret.push_back(append(m)); return ret; }
		int len(){ return (int)scene->materials.size(); }
};

void termHandlerNormal(int sig){cerr<<"Yade: normal exit."<<endl; raise(SIGTERM);}
void termHandlerError(int sig){cerr<<"Yade: error exit."<<endl; raise(SIGTERM);}

class pyOmega{
	private:
		// can be safely removed now, since pyOmega makes an empty rootBody in the constructor, if there is none
		void assertScene(){if(!OMEGA.getScene()) throw std::runtime_error("No Scene instance?!"); }
		Omega& OMEGA;
	public:
	pyOmega(): OMEGA(Omega::instance()){
		shared_ptr<Scene> rb=OMEGA.getScene();
		if(!rb){
			OMEGA.init();
			rb=OMEGA.getScene();
		}
		assert(rb);
		if(!OMEGA.hasSimulationLoop()){
			OMEGA.createSimulationLoop();
		}
	};
	/* Create variables in python's __builtin__ namespace that correspond to labeled objects. At this moment, only engines and functors can be labeled (not bodies etc). */
	void mapLabeledEntitiesToVariables(){
		// not sure if we should map materials to variables by default...
		// a call to this functions would have to be added to pyMaterialContainer::append
		#if 0
			FOREACH(const shared_ptr<Material>& m, OMEGA.getScene()->materials){
				if(!m->label.empty()) { PyGILState_STATE gstate; gstate = PyGILState_Ensure(); PyRun_SimpleString(("__builtins__."+m->label+"=Omega().materials["+lexical_cast<string>(m->id)+"]").c_str()); PyGILState_Release(gstate); }
			}
		#endif
		FOREACH(const shared_ptr<Engine>& e, OMEGA.getScene()->engines){
			if(!e->label.empty()){
				pyRunString("__builtins__."+e->label+"=Omega().labeledEngine('"+e->label+"')");
			}
			if(isChildClassOf(e->getClassName(),"Dispatcher")){
				shared_ptr<Dispatcher> ee=dynamic_pointer_cast<Dispatcher>(e);
				FOREACH(const shared_ptr<Functor>& f, ee->functorArguments){
					if(!f->label.empty()){
						pyRunString("__builtins__."+f->label+"=Omega().labeledEngine('"+f->label+"')");
					}
				}
			}
			if(e->getClassName()=="InteractionDispatchers"){
				shared_ptr<InteractionDispatchers> ee=dynamic_pointer_cast<InteractionDispatchers>(e);
				list<shared_ptr<Functor> > eus;
				FOREACH(const shared_ptr<Functor>& eu,ee->geomDispatcher->functorArguments) eus.push_back(eu);
				FOREACH(const shared_ptr<Functor>& eu,ee->physDispatcher->functorArguments) eus.push_back(eu);
				FOREACH(const shared_ptr<Functor>& eu,ee->lawDispatcher->functorArguments) eus.push_back(eu);
				FOREACH(const shared_ptr<Functor>& eu,eus){
					if(!eu->label.empty()){
						pyRunString("__builtins__."+eu->label+"=Omega().labeledEngine('"+eu->label+"')");
					}
				}
			}
		}
	}

	long iter(){ return OMEGA.getScene()->currentIteration;}
	double simulationTime(){return OMEGA.getSimulationTime();}
	double realTime(){ return OMEGA.getComputationTime(); }
	double dt_get(){return OMEGA.getScene()->dt;}
	void dt_set(double dt){
		Scene* scene=OMEGA.getScene().get();
		// activate timestepper, if possible (throw exception if there is none)
		if(dt<=0){ if(!scene->timeStepperActivate(true)) /* not activated*/ throw runtime_error("No TimeStepper found in O.engines."); if(dt!=0) scene->dt=-dt; }
		else { scene->timeStepperActivate(false); scene->dt=dt; }
	}
	bool dynDt_get(){return OMEGA.getScene()->timeStepperActive();}
	long stopAtIter_get(){return OMEGA.getScene()->stopAtIteration; }
	void stopAtIter_set(long s){OMEGA.getScene()->stopAtIteration=s; }


	bool timingEnabled_get(){return TimingInfo::enabled;}
	void timingEnabled_set(bool enabled){TimingInfo::enabled=enabled;}
	// deprecated:
		unsigned long forceSyncCount_get(){ return OMEGA.getScene()->forces.syncCount;}
		void forceSyncCount_set(unsigned long count){ OMEGA.getScene()->forces.syncCount=count;}

	void run(long int numIter=-1,bool doWait=false){
		if(numIter>0) OMEGA.getScene()->stopAtIteration=OMEGA.getCurrentIteration()+numIter;
		OMEGA.startSimulationLoop();
		// timespec t1,t2; t1.tv_sec=0; t1.tv_nsec=40000000; /* 40 ms */
		// while(!OMEGA.isRunning()) nanosleep(&t1,&t2); // wait till we start, so that calling wait() immediately afterwards doesn't return immediately
		LOG_DEBUG("RUN"<<((OMEGA.getScene()->stopAtIteration-OMEGA.getCurrentIteration())>0?string(" ("+lexical_cast<string>(OMEGA.getScene()->stopAtIteration-OMEGA.getCurrentIteration())+" to go)"):string(""))<<"!");
		if(doWait) wait();
	}
	void pause(){Py_BEGIN_ALLOW_THREADS; OMEGA.stopSimulationLoop(); Py_END_ALLOW_THREADS; LOG_DEBUG("PAUSE!");}
	void step() { if(OMEGA.isRunning()) throw runtime_error("Called O.step() while simulation is running."); OMEGA.getScene()->moveToNextTimeStep(); /* LOG_DEBUG("STEP!"); run(1); wait(); */ }
	void wait(){
		if(OMEGA.isRunning()){LOG_DEBUG("WAIT!");} else return;
		timespec t1,t2; t1.tv_sec=0; t1.tv_nsec=40000000; /* 40 ms */ Py_BEGIN_ALLOW_THREADS; while(OMEGA.isRunning()) nanosleep(&t1,&t2); Py_END_ALLOW_THREADS;
		if(!OMEGA.simulationLoop->workerThrew) return;
		LOG_ERROR("Simulation error encountered."); OMEGA.simulationLoop->workerThrew=false; throw OMEGA.simulationLoop->workerException;
	}
	void load(std::string fileName) {
		Py_BEGIN_ALLOW_THREADS; OMEGA.joinSimulationLoop(); Py_END_ALLOW_THREADS; 
		OMEGA.setSimulationFileName(fileName);
		OMEGA.loadSimulation();
		OMEGA.createSimulationLoop();
		mapLabeledEntitiesToVariables();
		LOG_DEBUG("LOAD!");
	}
	void reload(){	load(OMEGA.getSimulationFileName());}
	void saveTmp(string mark=""){ save(":memory:"+mark);}
	void loadTmp(string mark=""){ load(":memory:"+mark);}
	void tmpToFile(string mark, string filename){
		if(OMEGA.memSavedSimulations.count(":memory:"+mark)==0) throw runtime_error("No memory-saved simulation named "+mark);
		iostreams::filtering_ostream out;
		if(algorithm::ends_with(filename,".bz2")) out.push(iostreams::bzip2_compressor());
		out.push(iostreams::file_sink(filename));
		if(!out.good()) throw runtime_error("Error while opening file `"+filename+"' for writing.");
		LOG_INFO("Saving :memory:"<<mark<<" to "<<filename);
		out<<OMEGA.memSavedSimulations[":memory:"+mark];
	}
	string tmpToString(string mark){
		if(OMEGA.memSavedSimulations.count(":memory:"+mark)==0) throw runtime_error("No memory-saved simulation named "+mark);
		return OMEGA.memSavedSimulations[":memory:"+mark];
	}

	void reset(){Py_BEGIN_ALLOW_THREADS; OMEGA.reset(); Py_END_ALLOW_THREADS; }
	void resetThisScene(){Py_BEGIN_ALLOW_THREADS; OMEGA.joinSimulationLoop(); Py_END_ALLOW_THREADS; OMEGA.resetScene(); OMEGA.createSimulationLoop();}
	void resetTime(){ OMEGA.getScene()->currentIteration=0; OMEGA.getScene()->simulationTime=0; OMEGA.timeInit(); }
	void switchScene(){ std::swap(OMEGA.scene,OMEGA.sceneAnother); }

	void save(std::string fileName){
		assertScene();
		OMEGA.saveSimulation(fileName);
		OMEGA.setSimulationFileName(fileName);
		LOG_DEBUG("SAVE!");
	}
	
	python::list miscParams_get(){
		python::list ret;
		FOREACH(shared_ptr<Serializable>& s, OMEGA.getScene()->miscParams){
			ret.append(s);
		}
		return ret;
	}

	void miscParams_set(vector<shared_ptr<Serializable> > ss){
		vector<shared_ptr<Serializable> >& miscParams=OMEGA.getScene()->miscParams;
		miscParams.clear();
		FOREACH(shared_ptr<Serializable> s, ss){
			miscParams.push_back(s);
		}
	}

	vector<shared_ptr<Engine> > engines_get(void){assertScene(); return OMEGA.getScene()->engines;}
	void engines_set(const vector<shared_ptr<Engine> >& egs){assertScene(); OMEGA.getScene()->engines.clear(); FOREACH(const shared_ptr<Engine>& e, egs) OMEGA.getScene()->engines.push_back(e); mapLabeledEntitiesToVariables(); }
	vector<shared_ptr<Engine> > initializers_get(void){assertScene(); return OMEGA.getScene()->initializers;}
	void initializers_set(const vector<shared_ptr<Engine> >& egs){assertScene(); OMEGA.getScene()->initializers.clear(); FOREACH(const shared_ptr<Engine>& e, egs) OMEGA.getScene()->initializers.push_back(e); mapLabeledEntitiesToVariables(); OMEGA.getScene()->needsInitializers=true; }

	python::object labeled_engine_get(string label){
		FOREACH(const shared_ptr<Engine>& eng, OMEGA.getScene()->engines){
			if(eng->label==label){ return python::object(eng); }
			shared_ptr<Dispatcher> me=dynamic_pointer_cast<Dispatcher>(eng);
			if(me){
				FOREACH(const shared_ptr<Functor>& eu, me->functorArguments){
					if(eu->label==label) return python::object(eu);
				}
			}
			shared_ptr<InteractionDispatchers> ee=dynamic_pointer_cast<InteractionDispatchers>(eng);
			if(ee){
				list<shared_ptr<Functor> > eus;
				FOREACH(const shared_ptr<Functor>& eu,ee->geomDispatcher->functorArguments) eus.push_back(eu);
				FOREACH(const shared_ptr<Functor>& eu,ee->physDispatcher->functorArguments) eus.push_back(eu);
				FOREACH(const shared_ptr<Functor>& eu,ee->lawDispatcher->functorArguments) eus.push_back(eu);
				FOREACH(const shared_ptr<Functor>& eu,eus){
					if(eu->label==label) return python::object(eu);
				}
			}
		}
		throw std::invalid_argument(string("No engine labeled `")+label+"'");
	}
	
	pyBodyContainer bodies_get(void){assertScene(); return pyBodyContainer(OMEGA.getScene()->bodies); }
	pyInteractionContainer interactions_get(void){assertScene(); return pyInteractionContainer(OMEGA.getScene()->interactions); }
	
	pyForceContainer forces_get(void){return pyForceContainer(OMEGA.getScene());}
	pyMaterialContainer materials_get(void){return pyMaterialContainer(OMEGA.getScene());}
	

	python::list listChildClassesNonrecursive(const string& base){
		python::list ret;
		for(map<string,DynlibDescriptor>::const_iterator di=Omega::instance().getDynlibsDescriptor().begin();di!=Omega::instance().getDynlibsDescriptor().end();++di) if (Omega::instance().isInheritingFrom((*di).first,base)) ret.append(di->first);
		return ret;
	}

	bool isChildClassOf(const string& child, const string& base){
		return (Omega::instance().isInheritingFrom_recursive(child,base));
	}

	python::list plugins_get(){
		const map<string,DynlibDescriptor>& plugins=Omega::instance().getDynlibsDescriptor();
		std::pair<string,DynlibDescriptor> p; python::list ret;
		FOREACH(p, plugins) ret.append(p.first);
		return ret;
	}

	pyTags tags_get(void){assertScene(); return pyTags(OMEGA.getScene());}

	void interactionContainer_set(string clss){
		Scene* rb=OMEGA.getScene().get();
		if(rb->interactions->size()>0) throw std::runtime_error("Interaction container not empty, will not change its class.");
		shared_ptr<InteractionContainer> ic=dynamic_pointer_cast<InteractionContainer>(ClassFactory::instance().createShared(clss));
		rb->interactions=ic;
	}
	string interactionContainer_get(string clss){ return OMEGA.getScene()->interactions->getClassName(); }

	void bodyContainer_set(string clss){
		Scene* rb=OMEGA.getScene().get();
		if(rb->bodies->size()>0) throw std::runtime_error("Body container not empty, will not change its class.");
		shared_ptr<BodyContainer> bc=dynamic_pointer_cast<BodyContainer>(ClassFactory::instance().createShared(clss));
		rb->bodies=bc;
	}
	string bodyContainer_get(string clss){ return OMEGA.getScene()->bodies->getClassName(); }
	#ifdef YADE_OPENMP
		int numThreads_get(){ return omp_get_max_threads();}
		void numThreads_set(int n){ int bcn=OMEGA.getScene()->forces.getNumAllocatedThreads(); if(bcn<n) LOG_WARN("ForceContainer has only "<<bcn<<" threads allocated. Changing thread number to on "<<bcn<<" instead of "<<n<<" requested."); omp_set_num_threads(min(n,bcn)); LOG_WARN("BUG: Omega().numThreads=n doesn't work as expected (number of threads is not changed globally). Set env var OMP_NUM_THREADS instead."); }
	#else
		int numThreads_get(){return 1;}
		void numThreads_set(int n){ LOG_WARN("Yade was compiled without openMP support, changing number of threads will have no effect."); }
	#endif
	#ifdef YADE_BOOST_SERIALIZATION
	void saveBoost(string filename){
		const shared_ptr<Scene>& scene=OMEGA.getScene();
		yade::ObjectIO::save(filename,"scene",scene);
	}
	void loadBoost(string filename){
		shared_ptr<Scene> scene(new Scene);
		yade::ObjectIO::load(filename,"scene",scene);
		OMEGA.setScene(scene);
	}
	#endif
	
	shared_ptr<Cell> cell_get(){ if(OMEGA.getScene()->isPeriodic) return OMEGA.getScene()->cell; return shared_ptr<Cell>(); }
	bool periodic_get(void){ return OMEGA.getScene()->isPeriodic; } 
	void periodic_set(bool v){ OMEGA.getScene()->isPeriodic=v; }

	void disableGdb(){
		signal(SIGSEGV,SIG_DFL);
		signal(SIGABRT,SIG_DFL);
	}
	void exitNoBacktrace(int status=0){
		if(status==0) signal(SIGSEGV,termHandlerNormal); /* unset the handler that runs gdb and prints backtrace */
		else signal(SIGSEGV,termHandlerError);
		// try to clean our mess
		Omega::instance().cleanupTemps();
		// flush all streams (so that in case we crash at exit, unflushed buffers are not lost)
		fflush(NULL);
		// attempt exit
		exit(status);
	}
	void runEngine(const shared_ptr<Engine>& e){ LOG_WARN("Omega().runEngine(): deprecated, use __call__ method of the engine instance directly instead; will be removed in the future."); e->scene=OMEGA.getScene().get(); e->action(); }
	std::string tmpFilename(){ return OMEGA.tmpFilename(); }
};

class pySTLImporter : public STLImporter {};

shared_ptr<Shape> Body_shape_deprec_get(const shared_ptr<Body>& b){ LOG_WARN("Body().mold and Body().geom attributes are deprecated, use 'shape' instead."); return b->shape; }
void Body_shape_deprec_set(const shared_ptr<Body>& b, shared_ptr<Shape> ig){ LOG_WARN("Body().mold and Body().geom attributes are deprecated, use 'shape' instead."); b->shape=ig; }

void FileGenerator_generate(const shared_ptr<FileGenerator>& fg, string outFile){ fg->setFileName(outFile); fg->setSerializationLibrary("XMLFormatManager"); bool ret=fg->generateAndSave(); LOG_INFO((ret?"SUCCESS:\n":"FAILURE:\n")<<fg->message); if(ret==false) throw runtime_error("Generator reported error: "+fg->message); };
void FileGenerator_load(const shared_ptr<FileGenerator>& fg){ string xml(Omega::instance().tmpFilename()+".xml.bz2"); LOG_DEBUG("Using temp file "<<xml); FileGenerator_generate(fg,xml); pyOmega().load(xml); }

BOOST_PYTHON_MODULE(wrapper)
{
	python::scope().attr("__doc__")="Wrapper for c++ internals of yade.";

	YADE_SET_DOCSTRING_OPTS;

	python::class_<pyOmega>("Omega")
		.add_property("iter",&pyOmega::iter,"Get current step number")
		.add_property("stopAtIter",&pyOmega::stopAtIter_get,&pyOmega::stopAtIter_set,"Get/set number of iteration after which the simulation will stop.")
		.add_property("time",&pyOmega::simulationTime,"Return virtual (model world) time of the simulation.")
		.add_property("realtime",&pyOmega::realTime,"Return clock (human world) time the simulation has been running.")
		.add_property("dt",&pyOmega::dt_get,&pyOmega::dt_set,"Current timestep (Δt) value.\n\n* assigning zero enables dynamic Δt control via a :yref:`TimeStepper` (raises an exception if there is no :yref:`TimeStepper` among :yref:`O.engines<Omega.engines>`)\n* assigning negative value enables dynamic Δt (as in the previous case) and sets positive timestep ``O.dt=|Δt|`` (will be used until the timestepper is run and updates it)\n* assigning positive value sets Δt to that value and disables dynamic Δt (via :yref:`TimeStepper`, if there is one).\n\n:yref:`dynDt<Omega.dynDt>` can be used to query whether dynamic Δt is in use.")
		.add_property("dynDt",&pyOmega::dynDt_get,"Whether a :yref:`TimeStepper` is used for dynamic Δt control. See :yref:`dt<Omega.dt>` on how to enable/disable :yref:`TimeStepper`.")
		.def("load",&pyOmega::load,"Load simulation from file.")
		.def("reload",&pyOmega::reload,"Reload current simulation")
		.def("save",&pyOmega::save,"Save current simulation to file (should be .xml or .xml.bz2)")
		.def("loadTmp",&pyOmega::loadTmp,(python::args("mark")=""),"Load simulation previously stored in memory by saveTmp. *mark* optionally distinguishes multiple saved simulations")
		.def("saveTmp",&pyOmega::saveTmp,(python::args("mark")=""),"Save simulation to memory (disappears at shutdown), can be loaded later with loadTmp. *mark* optionally distinguishes different memory-saved simulations.")
		.def("tmpToFile",&pyOmega::tmpToFile,(python::arg("fileName"),python::arg("mark")=""),"Save XML of :yref:`saveTmp<Omega.saveTmp>`'d simulation into *fileName*.")
		.def("tmpToString",&pyOmega::tmpToString,(python::arg("mark")=""),"Return XML of :yref:`saveTmp<Omega.saveTmp>`'d simulation as string.")
		.def("run",&pyOmega::run,(python::arg("nSteps")=-1,python::arg("wait")=false),"Run the simulation. *nSteps* how many steps to run, then stop (if positive); *wait* will cause not returning to python until simulation will have stopped.")
		.def("pause",&pyOmega::pause,"Stop simulation execution. (May be called from within the loop, and it will stop after the current step).")
		.def("step",&pyOmega::step,"Advance the simulation by one step. Returns after the step will have finished.")
		.def("wait",&pyOmega::wait,"Don't return until the simulation will have been paused. (Returns immediately if not running).")
		.def("reset",&pyOmega::reset,"Reset simulations completely (including another scene!).")
		.def("resetThisScene",&pyOmega::resetThisScene,"Reset current scene.")
		.def("switchScene",&pyOmega::switchScene,"Switch to alternative simulation (while keeping the old one). Calling the function again switches back to the first one. Note that most variables from the first simulation will still refer to the first simulation even after the switch\n(e.g. b=O.bodies[4]; O.switchScene(); [b still refers to the body in the first simulation here])")
		.def("labeledEngine",&pyOmega::labeled_engine_get,"Return instance of engine/functor with the given label. This function shouldn't be called by the user directly; every ehange in O.engines will assign respective global python variables according to labels.\n\nFor example::\n\tO.engines=[InsertionSortCollider(label='collider')]\n\tcollider.nBins=5 ## collider has become a variable after assignment to O.engines automatically)")
		.def("resetTime",&pyOmega::resetTime,"Reset simulation time: step number, virtual and real time. (Doesn't touch anything else, including timings).")
		.def("plugins",&pyOmega::plugins_get,"Return list of all plugins registered in the class factory.")
		.add_property("engines",&pyOmega::engines_get,&pyOmega::engines_set,"List of engines in the simulation (Scene::engines).")
		.add_property("miscParams",&pyOmega::miscParams_get,&pyOmega::miscParams_set,"MiscParams in the simulation (Scene::mistParams), usually used to save serializables that don't fit anywhere else, like GL functors")
		.add_property("initializers",&pyOmega::initializers_get,&pyOmega::initializers_set,"List of initializers (Scene::initializers).")
		.add_property("bodies",&pyOmega::bodies_get,"Bodies in the current simulation (container supporting index access by id and iteration)")
		.add_property("interactions",&pyOmega::interactions_get,"Interactions in the current simulation (container supporting index acces by either (id1,id2) or interactionNumber and iteration)")
		.add_property("materials",&pyOmega::materials_get,"Shared materials; they can be accessed by id or by label")
		.add_property("forces",&pyOmega::forces_get,"ForceContainer (forces, torques, displacements) in the current simulation.")
		.add_property("tags",&pyOmega::tags_get,"Tags (string=string dictionary) of the current simulation (container supporting string-index access/assignment)")
		.def("childClassesNonrecursive",&pyOmega::listChildClassesNonrecursive,"Return list of all classes deriving from given class, as registered in the class factory")
		.def("isChildClassOf",&pyOmega::isChildClassOf,"Tells whether the first class derives from the second one (both given as strings).")
		.add_property("timingEnabled",&pyOmega::timingEnabled_get,&pyOmega::timingEnabled_set,"Globally enable/disable timing services (see documentation of yade.timing).")
		.add_property("forceSyncCount",&pyOmega::forceSyncCount_get,&pyOmega::forceSyncCount_set,"Counter for number of syncs in ForceContainer, for profiling purposes.")
		.add_property("numThreads",&pyOmega::numThreads_get /* ,&pyOmega::numThreads_set*/ ,"Get maximum number of threads openMP can use.")
		.add_property("cell",&pyOmega::cell_get,"Periodic cell of the current scene (None if the scene is aperiodic).")
		.add_property("periodic",&pyOmega::periodic_get,&pyOmega::periodic_set,"Get/set whether the scene is periodic or not (True/False).")
		.def("exitNoBacktrace",&pyOmega::exitNoBacktrace,(python::arg("status")=0),"Disable SEGV handler and exit, optionally with given status number.")
		.def("disableGdb",&pyOmega::disableGdb,"Revert SEGV and ABRT handlers to system defaults.")
		#ifdef YADE_BOOST_SERIALIZATION
			.def("load2",&pyOmega::loadBoost,"[EXPERIMENTAL] load using boost::serialization (handles compression, XML/binary)")
			.def("save2",&pyOmega::saveBoost,"[EXPERIMENTAL] save using boost::serialization (handles compression, XML/binary)")
		#endif
		.def("runEngine",&pyOmega::runEngine,"Run given engine exactly once; simulation time, step number etc. will not be incremented (use only if you know what you do).")
		.def("tmpFilename",&pyOmega::tmpFilename,"Return unique name of file in temporary directory which will be deleted when yade exits.")
		;
	python::class_<pyTags>("TagsWrapper","Container emulating dictionary semantics for accessing tags associated with simulation. Tags are accesed by strings.",python::init<pyTags&>())
		.def("__getitem__",&pyTags::getItem)
		.def("__setitem__",&pyTags::setItem)
		.def("keys",&pyTags::keys)
		.def("has_key",&pyTags::hasKey);
	python::class_<pyBodyContainer>("BodyContainer",python::init<pyBodyContainer&>())
		.def("__getitem__",&pyBodyContainer::pyGetitem)
		.def("__len__",&pyBodyContainer::length)
		.def("__iter__",&pyBodyContainer::pyIter)
		.def("append",&pyBodyContainer::append,"Append one Body instance, return its id.")
		.def("append",&pyBodyContainer::appendList,"Append list of Body instance, return list of ids")
		.def("appendClumped",&pyBodyContainer::appendClump,"Append given list of bodies as a clump (rigid aggregate); return list of ids.")
		.def("clear", &pyBodyContainer::clear,"Remove all bodies (interactions not checked)")
		.def("erase", &pyBodyContainer::erase,"Erase body with the given id; all interaction will be deleted by InteractionDispatchers in the next step.")
		.def("replace",&pyBodyContainer::replace);
	python::class_<pyBodyIterator>("BodyIterator",python::init<pyBodyIterator&>())
		.def("__iter__",&pyBodyIterator::pyIter)
		.def("next",&pyBodyIterator::pyNext);
	python::class_<pyInteractionContainer>("InteractionContainer","Access to :yref:`interactions<Interaction>` of simulation, by using \n\n#. id's of both :yref:`Bodies<Body>` of the interactions, e.g. ``O.interactions[23,65]``\n#. iteraction over the whole container::\n\n\tfor i in O.interactions: print i.id1,i.id2\n\n.. note::\n\tIteration silently skips interactions that are not :yref:`real<Interaction.isReal>`.",python::init<pyInteractionContainer&>())
		.def("__iter__",&pyInteractionContainer::pyIter)
		.def("__getitem__",&pyInteractionContainer::pyGetitem)
		.def("__len__",&pyInteractionContainer::len)
		.def("countReal",&pyInteractionContainer::countReal,"Return number of interactions that are \"real\", i.e. they have phys and geom.")
		.def("nth",&pyInteractionContainer::pyNth,"Return n-th interaction from the container (usable for picking random interaction).")
		.def("withBody",&pyInteractionContainer::withBody,"Return list of real interactions of given body.")
		.def("withBodyAll",&pyInteractionContainer::withBodyAll,"Return list of all (real as well as non-real) interactions of given body.")
		.def("eraseNonReal",&pyInteractionContainer::eraseNonReal,"Erase all interactions that are not :yref:`real <InteractionContainer.isReal>`.")
		.add_property("serializeSorted",&pyInteractionContainer::serializeSorted_get,&pyInteractionContainer::serializeSorted_set)
		.def("clear",&pyInteractionContainer::clear,"Remove all interactions");
	python::class_<pyInteractionIterator>("InteractionIterator",python::init<pyInteractionIterator&>())
		.def("__iter__",&pyInteractionIterator::pyIter)
		.def("next",&pyInteractionIterator::pyNext);

	python::class_<pyForceContainer>("ForceContainer",python::init<pyForceContainer&>())
		.def("f",&pyForceContainer::force_get,(python::arg("id")),"Force applied on body.")
		.def("t",&pyForceContainer::torque_get,(python::arg("id")),"Torque applied on body.")
		.def("m",&pyForceContainer::torque_get,(python::arg("id")),"Deprecated alias for t (torque).")
		.def("move",&pyForceContainer::move_get,(python::arg("id")),"Displacement applied on body.")
		.def("rot",&pyForceContainer::rot_get,(python::arg("id")),"Rotation applied on body.")
		.def("addF",&pyForceContainer::force_add,(python::arg("id"),python::arg("f")),"Apply force on body (accumulates).")
		.def("addT",&pyForceContainer::torque_add,(python::arg("id"),python::arg("t")),"Apply torque on body (accumulates).")
		.def("addMove",&pyForceContainer::move_add,(python::arg("id"),python::arg("m")),"Apply displacement on body (accumulates).")
		.def("addRot",&pyForceContainer::rot_add,(python::arg("id"),python::arg("r")),"Apply rotation on body (accumulates).")
		.add_property("syncCount",&pyForceContainer::syncCount_get,&pyForceContainer::syncCount_set,"Number of synchronizations  of ForceContainer (cummulative); if significantly higher than number of steps, there might be unnecessary syncs hurting performance.")
		;

	python::class_<pyMaterialContainer>("MaterialContainer","Container for :yref:`Materials<Material>`. A material can be accessed using \n\n #. numerical index in range(0,len(cont)), like cont[2]; \n #. textual label that was given to the material, like cont['steel']. This etails traversing all materials and should not be used frequently.",python::init<pyMaterialContainer&>())
		.def("append",&pyMaterialContainer::append,"Add new shared material; changes its id and return it.")
		.def("append",&pyMaterialContainer::appendList,"Append list of Material instances, return list of ids.")
		.def("__getitem__",&pyMaterialContainer::getitem_id)
		.def("__getitem__",&pyMaterialContainer::getitem_label)
		.def("__len__",&pyMaterialContainer::len);

	python::class_<pySTLImporter>("STLImporter")
		.def("ymport",&pySTLImporter::import);

//////////////////////////////////////////////////////////////
///////////// proxyless wrappers 
	Serializable().pyRegisterClass(python::scope());

	python::class_<TimingDeltas, shared_ptr<TimingDeltas>, noncopyable >("TimingDeltas").add_property("data",&TimingDeltas::pyData,"Get timing data as list of tuples (label, execTime[nsec], execCount) (one tuple per checkpoint)").def("reset",&TimingDeltas::reset,"Reset timing information");

	python::scope().attr("O")=pyOmega();
}


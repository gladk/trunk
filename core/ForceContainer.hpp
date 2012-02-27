// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#include<string.h>
#include<vector>
#include<yade/lib-base/Math.hpp>
#include<yade/core/Body.hpp>

#include<boost/static_assert.hpp>
// make sure that (void*)&vec[0]==(void*)&vec
BOOST_STATIC_ASSERT(sizeof(Vector3r)==3*sizeof(Real));


#ifdef YADE_OPENMP

#include<omp.h>
/*! Container for Body External Variables (forces), typically forces and torques from interactions.
 * Values should be reset at every iteration by calling ForceContainer::reset();
 * If you want to add your own force type, you need to:
 *
 * 	1. Create storage vector
 * 	2. Create accessor function
 * 	3. Update the resize function
 * 	4. Update the reset function
 * 	5. update the sync function (for the multithreaded implementation)
 *
 * This class exists in two flavors: non-parallel and parallel. The parallel one stores
 * force increments separately for every thread and sums those when sync() is called.
 * The reason of this design is that the container is not truly random-access, but rather
 * is written to everywhere in one phase and read in the next one. Adding to force/torque
 * marks the container as dirty and sync() must be performed before reading the stored data.
 * Calling getForce/getTorque when the container is not synchronized throws an exception.
 *
 * It is intentional that sync() needs to be called exlicitly, since syncs are expensive and
 * the programmer should be aware of that. Sync is however performed only if the container
 * is dirty. Every full sync increments the syncCount variable, that should ideally equal
 * the number of steps (one per step).
 *
 * The number of threads (omp_get_max_threads) may not change once ForceContainer is constructed.
 *
 * The non-parallel flavor has the same interface, but sync() is no-op and synchronization
 * is not enforced at all.
 */

//! This is the parallel flavor of ForceContainer
class ForceContainer{
	private:
		typedef std::vector<Vector3r> vvector;
		std::vector<vvector> _forceData;
		std::vector<vvector> _torqueData;
		std::vector<vvector> _moveData;
		std::vector<vvector> _rotData;
		vvector _force, _torque, _move, _rot;
		std::vector<size_t> sizeOfThreads;
		size_t size;
		bool syncedSizes;
		int nThreads;
		bool synced,moveRotUsed;
		boost::mutex globalMutex;
		Vector3r _zero;

		inline void ensureSize(Body::id_t id, int threadN){
			assert(nThreads>omp_get_thread_num());
			if (sizeOfThreads[threadN]<=(size_t)id) resize(min((size_t)1.5*(id+100),(size_t)(id+2000)),threadN);
		}

		inline void ensureSynced(){ if(!synced) throw runtime_error("ForceContainer not thread-synchronized; call sync() first!"); }

		#if 0
			/*! Function to allow friend classes to get force even if not synced.
			* Dangerous! The caller must know what it is doing! (i.e. don't read after write
			* for a particular body id. */
			const Vector3r& getForceUnsynced (Body::id_t id){ensureSize(id); return _force[id];}
			const Vector3r& getTorqueUnsynced(Body::id_t id){ensureSize(id); return _force[id];}
		#endif
		// dummy function to avoid template resolution failure
		friend class boost::serialization::access; template<class ArchiveT> void serialize(ArchiveT & ar, unsigned int version){}
	public:
		ForceContainer(): size(0),syncedSizes(true),synced(true),moveRotUsed(false),_zero(Vector3r::Zero()),syncCount(0),lastReset(0){
			nThreads=omp_get_max_threads();
			for(int i=0; i<nThreads; i++){
				_forceData.push_back(vvector()); _torqueData.push_back(vvector());
				_moveData.push_back(vvector());  _rotData.push_back(vvector());
				sizeOfThreads.push_back(0);
			}
		}

		const Vector3r& getForce(Body::id_t id)         { ensureSynced(); return ((size_t)id<size)?_force[id]:_zero; }
		void  addForce(Body::id_t id, const Vector3r& f){ ensureSize(id,omp_get_thread_num()); synced=false;   _forceData[omp_get_thread_num()][id]+=f;}
		const Vector3r& getTorque(Body::id_t id)        { ensureSynced(); return ((size_t)id<size)?_torque[id]:_zero; }
		void addTorque(Body::id_t id, const Vector3r& t){ ensureSize(id,omp_get_thread_num()); synced=false;   _torqueData[omp_get_thread_num()][id]+=t;}
		const Vector3r& getMove(Body::id_t id)          { ensureSynced(); return ((size_t)id<size)?_move[id]:_zero; }
		void  addMove(Body::id_t id, const Vector3r& m) { ensureSize(id,omp_get_thread_num()); synced=false; moveRotUsed=true; _moveData[omp_get_thread_num()][id]+=m;}
		const Vector3r& getRot(Body::id_t id)           { ensureSynced(); return ((size_t)id<size)?_rot[id]:_zero; }
		void  addRot(Body::id_t id, const Vector3r& r)  { ensureSize(id,omp_get_thread_num()); synced=false; moveRotUsed=true; _rotData[omp_get_thread_num()][id]+=r;}
		/* To be benchmarked: sum thread data in getForce/getTorque upon request for each body individually instead of by the sync() function globally */
		// this function is used from python so that running simulation is not slowed down by sync'ing on occasions
		// since Vector3r writes are not atomic, it might (rarely) return wrong value, if the computation is running meanwhile
		Vector3r getForceSingle (Body::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_forceData [t][id]:_zero; } return ret; }
		Vector3r getTorqueSingle(Body::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_torqueData[t][id]:_zero; } return ret; }
		Vector3r getMoveSingle  (Body::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_moveData  [t][id]:_zero; } return ret; }
		Vector3r getRotSingle   (Body::id_t id){ Vector3r ret(Vector3r::Zero()); for(int t=0; t<nThreads; t++){ ret+=((size_t)id<sizeOfThreads[t])?_rotData   [t][id]:_zero; } return ret; }
		
		inline void syncSizesOfContainers() {
			if (syncedSizes) return;
			//check whether all containers have equal length, and if not resize it
			for(int i=0; i<nThreads; i++){
				if (sizeOfThreads[i]<size) resize(size,i);
			}
			_force.resize(size,Vector3r::Zero());
			_torque.resize(size,Vector3r::Zero());
			_move.resize(size,Vector3r::Zero());
			_rot.resize(size,Vector3r::Zero());
			syncedSizes=true;
		}
		/* Sum contributions from all threads, save to _force&_torque.
		 * Locks globalMutex, since one thread modifies common data (_force&_torque).
		 * Must be called before get* methods are used. Exception is thrown otherwise, since data are not consistent. */
		inline void sync(){
			if(synced) return;
			boost::mutex::scoped_lock lock(globalMutex);
			if(synced) return; // if synced meanwhile

			syncSizesOfContainers();

			for(long id=0; id<(long)size; id++){
				Vector3r sumF(Vector3r::Zero()), sumT(Vector3r::Zero());
				for(int thread=0; thread<nThreads; thread++){ sumF+=_forceData[thread][id]; sumT+=_torqueData[thread][id];}
				_force[id]=sumF; _torque[id]=sumT;
			}
			if(moveRotUsed){
				for(long id=0; id<(long)size; id++){
					Vector3r sumM(Vector3r::Zero()), sumR(Vector3r::Zero());
					for(int thread=0; thread<nThreads; thread++){ sumM+=_moveData[thread][id]; sumR+=_rotData[thread][id];}
					_move[id]=sumM; _rot[id]=sumR;
				}
			}
			synced=true; syncCount++;
		}
		unsigned long syncCount;
		long lastReset;

		void resize(size_t newSize, int threadN){
			_forceData [threadN].resize(newSize,Vector3r::Zero());
			_torqueData[threadN].resize(newSize,Vector3r::Zero());
			_moveData[threadN].resize(newSize,Vector3r::Zero());
			_rotData[threadN].resize(newSize,Vector3r::Zero());
			sizeOfThreads[threadN] = newSize;
			if (size<newSize) size=newSize;
			syncedSizes=false;
		}
		/*! Reset all data, also reset summary forces/torques and mark the container clean. */
		// perhaps should be private and friend Scene or whatever the only caller should be
		void reset(long iter){
			syncSizesOfContainers();
			for(int thread=0; thread<nThreads; thread++){
				memset(&_forceData [thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
				memset(&_torqueData[thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
				if(moveRotUsed){
					memset(&_moveData  [thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
					memset(&_rotData   [thread][0],0,sizeof(Vector3r)*sizeOfThreads[thread]);
				}
			}
			memset(&_force [0], 0,sizeof(Vector3r)*size);
			memset(&_torque[0], 0,sizeof(Vector3r)*size);
			if(moveRotUsed){
				memset(&_move  [0], 0,sizeof(Vector3r)*size);
				memset(&_rot   [0], 0,sizeof(Vector3r)*size);
			}
			synced=true; moveRotUsed=false;
			lastReset=iter;
		}
		//! say for how many threads we have allocated space
		const int& getNumAllocatedThreads() const {return nThreads;}
		const bool& getMoveRotUsed() const {return moveRotUsed;}
};

#else
//! This is the non-parallel flavor of ForceContainer
class ForceContainer {
	private:
		std::vector<Vector3r> _force;
		std::vector<Vector3r> _torque;
		std::vector<Vector3r> _move;
		std::vector<Vector3r> _rot;
		size_t size;
		inline void ensureSize(Body::id_t id){ if(size<=(size_t)id) resize(min((size_t)1.5*(id+100),(size_t)(id+2000)));}
		#if 0
			const Vector3r& getForceUnsynced (Body::id_t id){ return getForce(id);}
			const Vector3r& getTorqueUnsynced(Body::id_t id){ return getForce(id);}
		#endif
		bool moveRotUsed;
		// dummy function to avoid template resolution failure
		friend class boost::serialization::access; template<class ArchiveT> void serialize(ArchiveT & ar, unsigned int version){}
	public:
		ForceContainer(): size(0), moveRotUsed(false), syncCount(0), lastReset(0){}
		const Vector3r& getForce(Body::id_t id){ensureSize(id); return _force[id];}
		void  addForce(Body::id_t id,const Vector3r& f){ensureSize(id); _force[id]+=f;}
		const Vector3r& getTorque(Body::id_t id){ensureSize(id); return _torque[id];}
		void  addTorque(Body::id_t id,const Vector3r& t){ensureSize(id); _torque[id]+=t;}
		const Vector3r& getMove(Body::id_t id){ensureSize(id); return _move[id];}
		void  addMove(Body::id_t id,const Vector3r& f){ensureSize(id); moveRotUsed=true; _move[id]+=f;}
		const Vector3r& getRot(Body::id_t id){ensureSize(id); return _rot[id];}
		void  addRot(Body::id_t id,const Vector3r& f){ensureSize(id); moveRotUsed=true; _rot[id]+=f;}
		// single getters do the same as globally synced ones in the non-parallel flavor
		const Vector3r& getForceSingle (Body::id_t id){ ensureSize(id); return _force [id]; }
		const Vector3r& getTorqueSingle(Body::id_t id){ ensureSize(id); return _torque[id]; }
		const Vector3r& getMoveSingle  (Body::id_t id){ ensureSize(id); return _move  [id]; }
		const Vector3r& getRotSingle   (Body::id_t id){ ensureSize(id); return _rot   [id]; }

		//! Set all forces to zero
		void reset(long iter){
			memset(&_force [0],0,sizeof(Vector3r)*size);
			memset(&_torque[0],0,sizeof(Vector3r)*size);
			if(moveRotUsed){
				memset(&_move  [0],0,sizeof(Vector3r)*size);
				memset(&_rot   [0],0,sizeof(Vector3r)*size);
				moveRotUsed=false;
			}
			lastReset=iter;
		}
		//! No-op for API compatibility with the threaded version
		void sync(){return;}
		unsigned long syncCount;
		// interaction in which the container was last reset; used by NewtonIntegrator to detect whether ForceResetter was not forgotten
		long lastReset;
		/*! Resize the container; this happens automatically,
		 * but you may want to set the size beforehand to avoid resizes as the simulation grows. */
		void resize(size_t newSize){
			_force.resize(newSize,Vector3r::Zero());
			_torque.resize(newSize,Vector3r::Zero());
			_move.resize(newSize,Vector3r::Zero());
			_rot.resize(newSize,Vector3r::Zero());
			size=newSize;
		}
		const int getNumAllocatedThreads() const {return 1;}
		const bool& getMoveRotUsed() const {return moveRotUsed;}
};

#endif

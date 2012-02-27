// 2009 © Václav Šmilauer <eudoxos@arcig.cz>


#include"VelocityBins.hpp"
#include<yade/core/Scene.hpp>
#include<boost/foreach.hpp>
#ifndef FOREACH
#	define FOREACH BOOST_FOREACH
#endif
CREATE_LOGGER(VelocityBins);

bool VelocityBins::checkSize_incrementDists_shouldCollide(const Scene* scene){
	// number of particles increased, recollision necessary
	// smaller number of particles is handled in setBins
	if(bodyBins.size()<scene->bodies->size()) { bodyBins.resize(scene->bodies->size(),/* put new particles to the slowest bin*/ bins.size()-1); return true; }
	int i=0;
	const Real& dt=scene->dt;
	FOREACH(Bin& bin, bins){
		// NOTE: this mimics the integration scheme of NewtonIntegrator
		// if you use different integration method, it must be changed (or the infrastructure somehow modified to allow for that)
		bin.currDist+=dt*sqrt(bin.currMaxVelSq); i++;
		if(bin.currDist>bin.maxDist){
			LOG_TRACE("Collide: bin"<<i<<": max dist "<<bin.maxDist<<", current "<<bin.currDist);
			return true;
		}
	}
	return false;
}

void VelocityBins::setBins(Scene* scene, Real currMaxVelSq, Real refSweepLength){
	// initialization
		// sanity checks
		if(nBins<1 || nBins>100){ throw runtime_error("VelocityBins: Number of bins must be >=1 and <=100"); }
		if(binOverlap>=1 || binOverlap<=0){ LOG_ERROR("binOverlap set to 0.8 (was "<<binOverlap<<", not in range (0…1) )"); binOverlap=0.8;}
		// number of bins changed
		if(nBins!=bins.size()){
			LOG_INFO("New number of bins: "<<nBins);
			bins.resize(nBins);	
		}
		// number of bodies changed
		if(bodyBins.size()!=scene->bodies->size()) bodyBins.resize(scene->bodies->size(),-1);
		// set the new overall refMaxVelSq
		if(refMaxVelSq<0){ refMaxVelSq=currMaxVelSq; /* first time */}
		else {
			// there should be some maximum speed change parameter, so that bins do not change their limits (and therefore bodies, also!) too often, depending on 1 particle going crazy
			if(maxRefRelStep>0){
				Real limVelSq=pow(sqrt(refMaxVelSq)+(sqrt(currMaxVelSq)-sqrt(refMaxVelSq))*maxRefRelStep,2);
				// LOG_INFO("limVel="<<sqrt(limVelSq)<<",currMaxVel="<<sqrt(currMaxVelSq)<<",refMaxVel="<<sqrt(refMaxVelSq));
				if(currMaxVelSq>refMaxVelSq) refMaxVelSq=min(limVelSq,currMaxVelSq);
				else refMaxVelSq=max(limVelSq,currMaxVelSq);
				// LOG_INFO("new refMaxVelSq="<<sqrt(refMaxVelSq));
			}
			//}refMaxVelSq=min(max(refMaxVelSq/pow(1+maxRefRelStep,2),currMaxVelSq),refMaxVelSq*pow(1+maxRefRelStep,2));
			else refMaxVelSq=currMaxVelSq;
			if(refMaxVelSq==0) refMaxVelSq=currMaxVelSq;
		}
		LOG_TRACE("new refMaxVel: "<<sqrt(refMaxVelSq));
		// compute new minima/maxima for all bins, reset distance counters and real max velocity
		for(size_t i=0; i<nBins; i++){
			Bin& bin=bins[i];
			// 0th bin (fastest) has maximum the current maximum; slowest bin has minimum 0.
			bin.binMaxVelSq=(i==0       ? currMaxVelSq : refMaxVelSq/pow(binCoeff*binCoeff,(int)i));
			bin.binMinVelSq=(i==nBins-1 ? 0.           : refMaxVelSq/pow(binCoeff*binCoeff,(int)(i+1)));
			bin.maxDist=(i==0 ?
				(refMaxVelSq==0 ? 0: refSweepLength) :
				refSweepLength/pow(binCoeff,(int)i)
			);
			bin.currDist=0; bin.currMaxVelSq=0; bin.nBodies=0;
		}
	long moveFaster=0, moveSlower=0;
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue;
		Real velSq=VelocityBins::getBodyVelSq(b->state.get());
		binNo_t newBin=-1, oldBin=bodyBins[b->getId()];
		// we could compute logarithm and round it here, but perhaps this is faster
		for(size_t i=0; i<nBins; i++){
			// 1. do not allow to move up (slow down) more than one bin at time (reasonable?)
			// 2. for the current bin, put the lower margin lower by binOverlap^2, to avoid too many oscillations
			if((oldBin>=0 && oldBin==(binNo_t)(i-1)) || velSq>=bins[i].binMinVelSq*(oldBin==(binNo_t)i ? pow(binOverlap,2) : 1.) ){ newBin=(binNo_t)i; break;}
		}
		if(newBin<0){
			if(isnan(velSq)) throw runtime_error("Body #"+lexical_cast<string>(b->getId())+" has velocity==NaN!");
			throw logic_error("Body #"+lexical_cast<string>(b->getId())+", velSq="+lexical_cast<string>(velSq)+" was not put in any bin?!");
		}
		if(oldBin>=0) { if(newBin>oldBin) moveSlower++; else if(oldBin>newBin) moveFaster++; }
		bodyBins[b->getId()]=newBin;
		// LOG_TRACE("#"<<b->getId()<<": vel="<<sqrt(velSq)<<", bin "<<(int)newBin);
		bins[newBin].nBodies+=1;
	}
	#ifdef YADE_LOG4CXX
		// if debugging output
		if(logger->isDebugEnabled() && (scene->iter-histLast>=histInterval || histLast<0)){
			histLast=scene->iter;
			LOG_INFO(bodyBins.size()<<" bodies (moves: "<<moveFaster<<" faster, "<<moveSlower<<" slower), refMaxVel="<<sqrt(refMaxVelSq));
			for(size_t i=0; i<nBins; i++){
				int nChars=int(80*(bins[i].nBodies/Real(bodyBins.size())));
				string l; for(int j=0; j<nChars; j++) l+="#"; for(int j=nChars; j<80; j++) l+=" ";
				cerr<<"\t"<<i<<": |"<<l<<"| ("<<bins[i].nBodies<<") "<<sqrt(bins[i].binMinVelSq)<<"…"<<sqrt(bins[i].binMaxVelSq)<<endl;
			}
		}
	#endif
}

/* non-parallel implementations */
#ifdef YADE_OPENMP
	void VelocityBins::binVelSqInitialize(Real vSqInit){ FOREACH(Bin& bin, bins){ FOREACH(Real& vSq, bin.threadMaxVelSq) vSq=vSqInit; }}
	void VelocityBins::binVelSqUse(Body::id_t id, Real velSq){
		Real& maxVelSq(bins[bodyBins[id]].threadMaxVelSq[omp_get_thread_num()]);
		maxVelSq=max(maxVelSq,velSq);
	}
	void VelocityBins::binVelSqFinalize(){
		FOREACH(Bin& bin, bins){
			bin.currMaxVelSq=0;
			FOREACH(const Real& vSq, bin.threadMaxVelSq) bin.currMaxVelSq=max(bin.currMaxVelSq,vSq);
		}
	}
#else
	void VelocityBins::binVelSqInitialize(Real vSqInit){ FOREACH(Bin& bin, bins) bin.currMaxVelSq=vSqInit; }
	void VelocityBins::binVelSqUse(Body::id_t id, Real velSq){
		Real& maxVelSq(bins[bodyBins[id]].currMaxVelSq);
		maxVelSq=max(maxVelSq,velSq);
	}
	void VelocityBins::binVelSqFinalize(){}
#endif


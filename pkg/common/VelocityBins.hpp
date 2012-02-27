// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once 

#include<yade/core/Body.hpp> // for Body::id_t
#include<yade/core/State.hpp>
#include<yade/core/Interaction.hpp>
#include<vector>
#ifdef YADE_OPENMP
	#include<omp.h>
#endif

class Scene;

/* Class for putting bodies to velcoity bins, for optimization
of collision detection.
Each bin is characterized by its minimum/maximum velocity. 
See http://yade.wikia.com/wiki/Insertion_Sort_Collider_Stride#Enhancement_ideas:_velocity_bins for brief design overview.
*/
class VelocityBins{
	public:
	VelocityBins(int _nBins, Real _refMaxVelSq, Real _binCoeff=10, Real _binOverlap=0.8): refMaxVelSq(_refMaxVelSq), binCoeff(_binCoeff), binOverlap(_binOverlap), maxRefRelStep(-1), nBins(_nBins), histInterval(200), histLast(-1){}
	typedef signed char binNo_t;
	struct Bin{
		Bin(): binMinVelSq(-1), binMaxVelSq(-1), maxDist(0), currDist(0), currMaxVelSq(0), nBodies(0){
			#ifdef YADE_OPENMP
				threadMaxVelSq.resize(omp_get_max_threads(),0.);
			#endif
		};
		// limits for bin memebrship
		Real binMinVelSq, binMaxVelSq;
		// maximum distance that body in this bin can travel before it goes out of its swept bbox
		Real maxDist;
		// distance so far traveled by the fastest body in this bin (since last setBins)
		Real currDist;
		// maximum velSq over all bodies in this bin
		Real currMaxVelSq;
		// number of bodies in this bin (for informational purposes only)
		long nBodies;
		#ifdef YADE_OPENMP
			vector<Real> threadMaxVelSq;
		#endif
	};
	// bins themselves (their number is nBins)
	vector<Bin> bins;
	// each body has its bin number stored here
	vector<binNo_t> bodyBins;
	// reference overall maximum velocity (may be different from bins[0].maxVelSq)
	Real refMaxVelSq;
	// factor by which maximum velocity in each bin is smaller than in the higer one
	Real binCoeff;
	// relative overlap beween bins; body will not be moved from faster bin until its velocity is min*binOverlap; must be <=1
	Real binOverlap;
	// maximum relative change of reference max velocity per invocation (if <0, disabled; this is the default)
	Real maxRefRelStep;
	// number of bins; must be >=1 and <=100 (artificial upper limit)
	size_t nBins;
	// how often to show histogram, if LOG_DEBUG is enabled;
	long histInterval, histLast;
	// Assign bins to all bodies
	void setBins(Scene*, Real currMaxVelSq, Real refSweepLength);

	// Increment maximum per-bin distances and tell whether some bodies may be	already getting out of the swept bbox (in that case, we need to recompute bounding volumes and run the collider)
	// Also returns true if number of particles changed, in which case it adjusts internal storage accordingly
	bool checkSize_incrementDists_shouldCollide(const Scene*);
	
	/* NOTE: following 3 functions are separated because of multi-threaded operation of NewtonIntegrator
	in that case, every thread must have its own per-bin maximum and binVelSqFinalize will assign the
	overall maxima to the bins.

	The non-parallel implementation simply resets all Bin::realMaxVelSq, maximize if for each
	body in that bin and binVelSqFinalize() will not do nothing.
	*/
	// reset per-bin max velocities
	void binVelSqInitialize(Real velSqInit=0.);
	// use body max velocity -- called for every body at every step (from NewtonIntegrator, normally)
	void binVelSqUse(Body::id_t id, Real velSq);
	// actually assign max velocities to their respective bins
	void binVelSqFinalize();

	// get velSq for given body; this should be called from NewtonIntegrator as well,
	// to ensure that the same formulas are used (once we have angularVelocity + Aabb span,
	// for instance
	static Real getBodyVelSq(const State* state){
		return state->vel.squaredNorm();
	}

	DECLARE_LOGGER;
};

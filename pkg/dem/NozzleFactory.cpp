
#include<yade/pkg-dem/NozzleFactory.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


YADE_PLUGIN((NozzleFactory)(DragForceApplier));
CREATE_LOGGER(NozzleFactory);

void DragForceApplier::action(){
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue;
		Sphere* sphere=dynamic_cast<Sphere*>(b->shape.get());
		if(!sphere) continue;
		Vector3r dragForce=-b->state->vel.normalized()*(1/2.)*density*b->state->vel.squaredNorm()*.47*Mathr::PI*pow(sphere->radius,2);
		scene->forces.addForce(b->id,dragForce);
	}
}

void NozzleFactory::action(){
	// initialize random number generator with time seed
	static boost::minstd_rand randGen(TimingInfo::getNow(/* get the number even if timing is disabled globally */ true));
	static boost::variate_generator<boost::minstd_rand&, boost::uniform_real<Real> > randomUnit(randGen, boost::uniform_real<Real>(0,1));

	if(!collider){
		FOREACH(const shared_ptr<Engine>& e, scene->engines){ collider=dynamic_pointer_cast<Collider>(e); if(collider) break; }
		if(!collider) throw runtime_error("NozzleFactory: No Collider instance found in engines (needed for collision detection).");
	}
	goalMass+=massFlowRate*scene->dt; // totalMass that we want to attain in the current step

	normal.normalize();
	LOG_TRACE("totalMass="<<totalMass<<", goalMass="<<goalMass);

	while(totalMass<goalMass && (maxParticles<0 || numParticles<maxParticles)){
		// pick random radius (later: based on the particle size distribution)
		Real r=rMin+randomUnit()*(rMax-rMin);
		LOG_TRACE("Radius "<<r);
		Vector3r c=Vector3r::Zero();
		// until there is no overlap, pick a random position for the new particle
		int attempt;
		for(attempt=0; attempt<maxAttempt; attempt++){
			//TODO: update c (position on the circle in space)
			//Real angle=randomUnit()*2*Mathr::PI, rr=randomUnit()*(radius-r); // random polar coordinate inside the nozzle
			// this version places center in a box around the nozzle center (its size 1.15=2/sqrt(3) diagonal, which is not very exact
			c=center+Vector3r((randomUnit()-.5)*1.15*radius,(randomUnit()-.5)*1.15*radius,(randomUnit()-.5)*1.15*radius);
			
			LOG_TRACE("Center "<<c);
			Bound b; b.min=c-Vector3r(r,r,r); b.max=c+Vector3r(r,r,r);
			vector<Body::id_t> collidingParticles=collider->probeBoundingVolume(b);
			if(collidingParticles.size()==0) break;
			#ifdef YADE_DEBUG
				FOREACH(const Body::id_t& id, collidingParticles) LOG_TRACE(scene->iter<<":"<<attempt<<": collision with #" <<id);
			#endif
		}
		if(attempt==maxAttempt) {
			LOG_WARN("Unable to place new sphere after "<<maxAttempt<<" attempts, giving up.");
			return;
		}
		// pick random initial velocity (normal with some variation)
		// preliminary version that randomizes valocity magnitude but always makes initVel exactly aligned with normal
		Vector3r initVel=normal*(vMin+randomUnit()*(vMax-vMin)); // TODO: compute from vMin, vMax, vAngle, normal;
		//Vector3r initVel= normal*100;

		// create particle
		int mId=(materialId>=0 ? materialId : scene->materials.size()+materialId);
		if(mId<0 || (size_t) mId>=scene->materials.size()) throw std::invalid_argument(("NozzleFactory: invalid material id "+lexical_cast<string>(materialId)).c_str());
		const shared_ptr<Material>& material=scene->materials[mId];
		shared_ptr<Body> b(new Body);
		shared_ptr<Sphere> sphere(new Sphere); 
		shared_ptr<State> state(material->newAssocState());
		sphere->radius=r;  
		state->pos=state->refPos=c; 

		state->vel=initVel; 
		Real vol=(4/3.)*Mathr::PI*pow(r,3);
		state->mass=vol*material->density;
		state->inertia=(2./5.)*vol*r*r*material->density*Vector3r::Ones();

		b->shape=sphere; 
		b->state=state; 
		b->material=material;
		// insert particle in the simulation
		scene->bodies->insert(b);
		// increment total mass we've spit out
		totalMass+=state->mass;
		numParticles++;
		
	} 
	//std::cout<<"mass flow rate: "<<totalMass<<endl;totalMass =0.0;
};

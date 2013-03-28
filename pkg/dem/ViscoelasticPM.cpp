// 2009 © Sergei Dorofeenko <sega@users.berlios.de>
#include"ViscoelasticPM.hpp"
#include<yade/core/State.hpp>
#include<yade/pkg/dem/ScGeom.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg/common/Sphere.hpp>

YADE_PLUGIN((ViscElMat)(ViscElPhys)(Ip2_ViscElMat_ViscElMat_ViscElPhys)(Law2_ScGeom_ViscElPhys_Basic));

/* ViscElMat */
ViscElMat::~ViscElMat(){}

/* ViscElPhys */
ViscElPhys::~ViscElPhys(){}

/* Ip2_ViscElMat_ViscElMat_ViscElPhys */
void Ip2_ViscElMat_ViscElMat_ViscElPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) {
	// no updates of an existing contact 
	if(interaction->phys) return;
	ViscElMat* mat1 = static_cast<ViscElMat*>(b1.get());
	ViscElMat* mat2 = static_cast<ViscElMat*>(b2.get());
	const Real mass1 = Body::byId(interaction->getId1())->state->mass;
	const Real mass2 = Body::byId(interaction->getId2())->state->mass;
	const Real kn1 = mat1->kn*mass1; const Real cn1 = mat1->cn*mass1;
	const Real ks1 = mat1->ks*mass1; const Real cs1 = mat1->cs*mass1;
	const Real kn2 = mat2->kn*mass2; const Real cn2 = mat2->cn*mass2;
	const Real ks2 = mat2->ks*mass2; const Real cs2 = mat2->cs*mass2;
	
		
	ViscElPhys* phys = new ViscElPhys();
	
	if ((kn1>0) or (kn2>0)) {
		phys->kn = 1/( ((kn1>0)?1/kn1:0) + ((kn2>0)?1/kn2:0) );
	} else {
		phys->kn = 0;
	}
	if ((ks1>0) or (ks2>0)) {
		phys->ks = 1/( ((ks1>0)?1/ks1:0) + ((ks2>0)?1/ks2:0) );
	} else {
		phys->ks = 0;
	} 
	
	phys->cn = (cn1?1/cn1:0) + (cn2?1/cn2:0); phys->cn = phys->cn?1/phys->cn:0;
	phys->cs = (cs1?1/cs1:0) + (cs2?1/cs2:0); phys->cs = phys->cs?1/phys->cs:0;
	
	phys->tangensOfFrictionAngle = std::tan(std::min(mat1->frictionAngle, mat2->frictionAngle)); 
	phys->shearForce = Vector3r(0,0,0);
	
	if (mat1->Capillar and mat2->Capillar)  {
		if (mat1->Vb == mat2->Vb) {
			phys->Vb = mat1->Vb;
		} else {
			throw runtime_error("Vb should be equal for both particles!.");
		}
		
		if (mat1->gamma == mat2->gamma) {
			phys->gamma = mat1->gamma;
		} else {
			throw runtime_error("Gamma should be equal for both particles!.");
		}
	
		if (mat1->theta == mat2->theta) {
			phys->theta = (mat1->theta*M_PI/180.0);
		} else {
			throw runtime_error("Theta should be equal for both particles!.");
		}
		phys->Capillar=true;
	}
	
	interaction->phys = shared_ptr<ViscElPhys>(phys);
}

/* Law2_ScGeom_ViscElPhys_Basic */
void Law2_ScGeom_ViscElPhys_Basic::go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I){

	const ScGeom& geom=*static_cast<ScGeom*>(_geom.get());
	ViscElPhys& phys=*static_cast<ViscElPhys*>(_phys.get());

	const int id1 = I->getId1();
	const int id2 = I->getId2();
	
	if (geom.penetrationDepth<0) {
		if (phys.liqBridgeCreated and -geom.penetrationDepth<phys.sCrit and phys.Capillar) {
			//Capillar
			Real c0 = 0.96;
			Real c1 = 1.1;
			Real R = phys.R;
			Real s = -geom.penetrationDepth;
			
			Real beta = asin(pow(phys.Vb/((c0*R*R*R*(1+3*s/R)*(1+c1*sin(phys.theta)))), 1.0/4.0));
			Real r1 = (R*(1-cos(beta)) + s/2.0)/(cos(beta+phys.theta));
			Real r2 = R*sin(beta) + r1*(sin(beta+phys.theta)-1);
			Real Pc = phys.gamma*(1/r1 - 1/r2);

			Real fC = 2*M_PI*phys.gamma*R*sin(beta)*sin(beta+phys.theta) + M_PI*R*R*Pc*sin(beta)*sin(beta);
			
			/*
			std::cerr<<"R: "<<phys.R<<std::endl;
			std::cerr<<"s: "<<s<<std::endl;
			std::cerr<<"Vb: "<<phys.Vb<<std::endl;
			std::cerr<<"Theta: "<<phys.theta<<std::endl;
			std::cerr<<"Gamma: "<<phys.gamma<<std::endl;
			std::cerr<<"beta: "<<beta<<std::endl;
			std::cerr<<"r1: "<<r1<<std::endl;
			std::cerr<<"r2: "<<r2<<std::endl;
			std::cerr<<"Pc: "<<Pc<<std::endl;
			std::cerr<<"Scrit: "<<phys.sCrit<<std::endl;
			std::cerr<<"Fc: "<<fC<<std::endl<<std::endl;
			*/
			
			phys.normalForce = -fC*geom.normal;
		  if (I->isActive) {
				addForce (id1,-phys.normalForce,scene);
				addForce (id2, phys.normalForce,scene);
			};
			//std::cerr<<"Capillar: "<<phys.normalForce<<"; Pen Depth: "<< geom.penetrationDepth <<"; Schritt "<< phys.sCrit <<std::endl;
			return;
		} else {
			scene->interactions->requestErase(I);
			return;
		};
	};

	const BodyContainer& bodies = *scene->bodies;

	const State& de1 = *static_cast<State*>(bodies[id1]->state.get());
	const State& de2 = *static_cast<State*>(bodies[id2]->state.get());
	
	if (not(phys.liqBridgeCreated) and phys.Capillar) {
		phys.liqBridgeCreated = true;
		//std::cerr<<"First contact! Setting corresponding variables."<<std::endl;
		phys.sCrit = (1+0.5*phys.theta)*pow(phys.Vb,1/3.0);
		Sphere* s1=dynamic_cast<Sphere*>(bodies[id1]->shape.get());
		Sphere* s2=dynamic_cast<Sphere*>(bodies[id2]->shape.get());
		if (s1 and s2) {
			if (s1->radius == s2->radius) {
				phys.R = s1->radius;
			} else {
				throw runtime_error("We can calculate only monodisperse for the moment!.");
			}
		} else if (s1 and not(s2)) {
			phys.R = s1->radius;
		} else {
			phys.R = s2->radius;
		}
	}

	Vector3r& shearForce = phys.shearForce;
	if (I->isFresh(scene)) shearForce=Vector3r(0,0,0);
	const Real& dt = scene->dt;
	//Vector3r axis = phys.prevNormal.cross(geom.normal);
	//shearForce -= shearForce.cross(axis);
	//const Real angle = dt*0.5*geom.normal.dot(de1.angVel + de2.angVel);
	//axis = angle*geom.normal;
	//shearForce -= shearForce.cross(axis);
	shearForce = geom.rotate(shearForce);
	

	// Handle periodicity.
	const Vector3r shift2 = scene->isPeriodic ? scene->cell->intrShiftPos(I->cellDist): Vector3r::Zero(); 
	const Vector3r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(I->cellDist): Vector3r::Zero(); 

	const Vector3r c1x = (geom.contactPoint - de1.pos);
	const Vector3r c2x = (geom.contactPoint - de2.pos - shift2);
	
	const Vector3r relativeVelocity = (de1.vel+de1.angVel.cross(c1x)) - (de2.vel+de2.angVel.cross(c2x)) + shiftVel;
	const Real normalVelocity	= geom.normal.dot(relativeVelocity);
	const Vector3r shearVelocity	= relativeVelocity-normalVelocity*geom.normal;

	// As Chiara Modenese suggest, we store the elastic part 
	// and then add the viscous part if we pass the Mohr-Coulomb criterion.
	// See http://www.mail-archive.com/yade-users@lists.launchpad.net/msg01391.html
	shearForce += phys.ks*dt*shearVelocity; // the elastic shear force have a history, but
	Vector3r shearForceVisc = Vector3r::Zero(); // the viscous shear damping haven't a history because it is a function of the instant velocity 

	phys.normalForce = ( phys.kn * geom.penetrationDepth + phys.cn * normalVelocity ) * geom.normal;
	//phys.prevNormal = geom.normal;

	const Real maxFs = phys.normalForce.squaredNorm() * std::pow(phys.tangensOfFrictionAngle,2);
	if( shearForce.squaredNorm() > maxFs )
	{
		// Then Mohr-Coulomb is violated (so, we slip), 
		// we have the max value of the shear force, so 
		// we consider only friction damping.
		const Real ratio = sqrt(maxFs) / shearForce.norm();
		shearForce *= ratio;
	} 
	else 
	{
		// Then no slip occurs we consider friction damping + viscous damping.
		shearForceVisc = phys.cs*shearVelocity; 
	}
	
	if (I->isActive) {
		//std::cerr<<"Contact: "<<phys.normalForce<<std::endl;
		const Vector3r f = phys.normalForce + shearForce + shearForceVisc;
		addForce (id1,-f,scene);
		addForce (id2, f,scene);
		addTorque(id1,-c1x.cross(f),scene);
		addTorque(id2, c2x.cross(f),scene);
  }
}

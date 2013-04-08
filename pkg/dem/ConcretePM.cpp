// 2007,2008 © Václav Šmilauer <eudoxos@arcig.cz> 
#include"ConcretePM.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg/dem/DemXDofGeom.hpp>
#include<yade/pkg/dem/Shop.hpp>

YADE_PLUGIN((CpmState)(CpmMat)(Ip2_CpmMat_CpmMat_CpmPhys)(Ip2_FrictMat_CpmMat_FrictPhys)(CpmPhys)(Law2_SomeGeom_CpmPhys_Cpm)(Law2_Dem3DofGeom_CpmPhys_Cpm)(Law2_ScGeom_CpmPhys_Cpm)
	#ifdef YADE_OPENGL
		(Gl1_CpmPhys)
	#endif	
		(CpmStateUpdater));


/********************** Ip2_CpmMat_CpmMat_CpmPhys ****************************/


CREATE_LOGGER(Ip2_FrictMat_CpmMat_FrictPhys);
void Ip2_FrictMat_CpmMat_FrictPhys::go(const shared_ptr<Material>& pp1, const shared_ptr<Material>& pp2, const shared_ptr<Interaction>& interaction){
	const shared_ptr<FrictMat>& mat1 = YADE_PTR_CAST<FrictMat>(pp1);
	const shared_ptr<CpmMat>& mat2 = YADE_PTR_CAST<CpmMat>(pp2);
	Ip2_FrictMat_FrictMat_FrictPhys().go(mat1,mat2,interaction);
}




CREATE_LOGGER(Ip2_CpmMat_CpmMat_CpmPhys);
void Ip2_CpmMat_CpmMat_CpmPhys::go(const shared_ptr<Material>& pp1, const shared_ptr<Material>& pp2, const shared_ptr<Interaction>& interaction){
	// no updates of an already existing contact necessary
	if (interaction->phys) return;
	shared_ptr<CpmPhys> cpmPhys(new CpmPhys());
	interaction->phys = cpmPhys;
	CpmMat* mat1 = YADE_CAST<CpmMat*>(pp1.get());
	CpmMat* mat2 = YADE_CAST<CpmMat*>(pp2.get());

	// check unassigned values
	if (!mat1->neverDamage) {
		assert(!isnan(mat1->sigmaT));
		assert(!isnan(mat1->epsCrackOnset));
		assert(!isnan(mat1->crackOpening) || !isnan(mat1->relDuctility));
	}
	if (!mat2->neverDamage) {
		assert(!isnan(mat2->sigmaT));
		assert(!isnan(mat2->epsCrackOnset));
		assert(!isnan(mat2->crackOpening) || !isnan(mat2->relDuctility));
	}

	cpmPhys->damLaw = mat1->damLaw;
	// bodies sharing the same material; no averages necessary
	if (mat1->id>=0 && mat1->id == mat2->id) {
		cpmPhys->E = mat1->young;
		cpmPhys->G = mat1->young*mat1->poisson;
		cpmPhys->tanFrictionAngle = tan(mat1->frictionAngle);
		cpmPhys->undamagedCohesion = mat1->sigmaT;
		cpmPhys->isCohesive = (cohesiveThresholdIter < 0 || scene->iter < cohesiveThresholdIter);
		#define _CPATTR(a) cpmPhys->a=mat1->a
			_CPATTR(epsCrackOnset);
			_CPATTR(crackOpening);
			_CPATTR(relDuctility);
			_CPATTR(neverDamage);
			_CPATTR(dmgTau);
			_CPATTR(dmgRateExp);
			_CPATTR(plTau);
			_CPATTR(plRateExp);
			_CPATTR(isoPrestress);
		#undef _CPATTR
	} else {
		// averaging over both materials
		#define _AVGATTR(a) cpmPhys->a=.5*(mat1->a+mat2->a)
			cpmPhys->E = .5*(mat1->young + mat2->young);
			cpmPhys->G = .5*(mat1->poisson + mat2->poisson)*.5*(mat1->young + mat2->young);
			cpmPhys->tanFrictionAngle = tan(.5*(mat1->frictionAngle + mat2->frictionAngle));
			cpmPhys->undamagedCohesion = .5*(mat1->sigmaT + mat2->sigmaT);
			cpmPhys->isCohesive = (cohesiveThresholdIter < 0 || scene->iter < cohesiveThresholdIter);
			_AVGATTR(epsCrackOnset);
			_AVGATTR(crackOpening);
			_AVGATTR(relDuctility);
			cpmPhys->neverDamage = (mat1->neverDamage || mat2->neverDamage);
			_AVGATTR(dmgTau);
			_AVGATTR(dmgRateExp);
			_AVGATTR(plTau);
			_AVGATTR(plRateExp);
			_AVGATTR(isoPrestress);
		#undef _AVGATTR
	}

	// NOTE: some params are not assigned until in Law2_SomeGeom_CpmPhys_Cpm, since they need geometry as well; those are:
	// 	crossSection, kn, ks, refLength
}




/********************** CpmPhys ****************************/
CREATE_LOGGER(CpmPhys);

// !! at least one virtual function in the .cpp file
CpmPhys::~CpmPhys(){};

long CpmPhys::cummBetaIter=0, CpmPhys::cummBetaCount=0;

Real CpmPhys::solveBeta(const Real c, const Real N){
	#ifdef YADE_DEBUG
		cummBetaCount++;
	#endif
	const int maxIter = 20;
	const Real maxError = 1e-12;
	Real f, ret = 0.;
	for(int i = 0; i < maxIter; i++){
		#ifdef YADE_DEBUG
			cummBetaIter++;
		#endif
		Real aux = c*exp(N*ret)+exp(ret);
		f = log(aux);
		if (fabs(f) < maxError) return ret;
		Real df = (c*N*exp(N*ret)+exp(ret))/aux;
		ret -= f/df;
	}
	LOG_FATAL("No convergence after "<<maxIter<<" iters; c="<<c<<", N="<<N<<", ret="<<ret<<", f="<<f);
	throw runtime_error("CpmPhys::solveBeta failed to converge.");
}

Real CpmPhys::computeDmgOverstress(Real dt){
	if (dmgStrain >= epsN*omega) { // unloading, no viscous stress
		dmgStrain = epsN*omega;
		LOG_TRACE("Elastic/unloading, no viscous overstress");
		return 0.;
	}
	Real c = epsCrackOnset*(1-omega)*pow(dmgTau/dt,dmgRateExp)*pow(epsN*omega-dmgStrain,dmgRateExp-1.);
	Real beta = solveBeta(c,dmgRateExp);
	Real deltaDmgStrain = (epsN*omega-dmgStrain)*exp(beta);
	dmgStrain += deltaDmgStrain;
	LOG_TRACE("deltaDmgStrain="<<deltaDmgStrain<<", viscous overstress "<<(epsN*omega-dmgStrain)*E);
	/* σN=Kn(εN-εd); dmgOverstress=σN-(1-ω)*Kn*εN=…=Kn(ω*εN-εd) */
	return (epsN*omega-dmgStrain)*E;
}

Real CpmPhys::computeViscoplScalingFactor(Real sigmaTNorm, Real sigmaTYield,Real dt){
	if (sigmaTNorm<sigmaTYield) return 1.;
	Real c = undamagedCohesion*pow(plTau/(G*dt),plRateExp)*pow(sigmaTNorm-sigmaTYield,plRateExp-1.);
	Real beta = solveBeta(c,plRateExp);
	//LOG_DEBUG("scaling factor "<<1.-exp(beta)*(1-sigmaTYield/sigmaTNorm));
	return 1.-exp(beta)*(1-sigmaTYield/sigmaTNorm);
}

Real CpmPhys::funcG(const Real& kappaD, const Real& epsCrackOnset, const Real& epsFracture, const bool& neverDamage, const int& damLaw) {
	if (kappaD<epsCrackOnset || neverDamage) return 0;
	switch (damLaw) {
		case 0: // linear
			return (1.-epsCrackOnset/kappaD)/(1.-epsCrackOnset/epsFracture);
		case 1: // exponential
			return 1.-(epsCrackOnset/kappaD)*exp(-(kappaD-epsCrackOnset)/epsFracture);
	}
	throw runtime_error("CpmPhys::funcG: wrong damLaw\n");
}

Real CpmPhys::funcGDKappa(const Real& kappaD, const Real& epsCrackOnset, const Real& epsFracture, const bool& neverDamage, const int& damLaw) {
	switch (damLaw) {
		case 0: // linear
			return epsCrackOnset / ((1.-epsCrackOnset/epsFracture)*kappaD*kappaD);
		case 1: // exponential
			return epsCrackOnset/kappaD * (1./kappaD + 1./epsFracture) * exp(-(kappaD-epsCrackOnset)/epsFracture);
	}
	throw runtime_error("CpmPhys::funcGDKappa: wrong damLaw\n");
}

Real CpmPhys::funcGInv(const Real& omega, const Real& epsCrackOnset, const Real& epsFracture, const bool& neverDamage, const int& damLaw) {
	if (omega==0. || neverDamage) return 0;
	switch (damLaw) {
		case 0: // linear
			return epsCrackOnset / (1. - omega*(1. - epsCrackOnset/epsFracture));
		case 1: // exponential
			// Newton's iterations
			Real fg,dfg,decr,ret=epsCrackOnset,tol=1e-3;
			int maxIter = 100;
			for (int i=0; i<maxIter; i++) {
				fg = - omega + 1. - epsCrackOnset/ret * exp(-(ret-epsCrackOnset)/epsFracture);
				//dfg = (epsCrackOnset/ret/ret - epsCrackOnset*(ret-epsCrackOnset)/ret/epsFracture/epsFracture) * exp(-(ret-epsCrackOnset)/epsFracture);
				dfg = CpmPhys::funcGDKappa(ret,epsCrackOnset,epsFracture,neverDamage,damLaw);
				decr = fg/dfg;
				ret -= decr;
				//printf("i %d fg %e dfg %e decr %e ret %e\n",i,fg,dfg,decr,ret);
				if (fabs(decr/epsCrackOnset) < tol) {
					return ret;
				}
			}
			throw runtime_error("CpmPhys::funcGInv: no convergence\n");
	}
	throw runtime_error("CpmPhys::funcGInv: wrong damLaw\n");
}

void CpmPhys::setDamage(Real dmg) {
	if (neverDamage) { return; }
	omega = dmg;
	kappaD = CpmPhys::funcGInv(dmg,epsCrackOnset,epsFracture,neverDamage,damLaw);
}

void CpmPhys::setRelResidualStrength(Real r) {
	if (neverDamage) { return; }
	if (r == 1.) {
		relResidualStrength = r;
		kappaD = omega = 0.;
		return;
	}
	Real k = epsFracture;
	Real g,dg,f,df,tol=1e-3,e0i=1./epsCrackOnset,decr;
	int maxIter = 100;
	int i;
	for (i=0; i<maxIter; i++) {
		g = CpmPhys::funcG(k,epsCrackOnset,epsFracture,neverDamage,damLaw);
		dg = CpmPhys::funcGDKappa(k,epsCrackOnset,epsFracture,neverDamage,damLaw);
		f = -r + (1-g)*k*e0i;
		df = e0i*(1-g-k*dg);
		decr = f/df;
		k -= decr;
		if (fabs(decr) < tol) {
			kappaD = k;
			omega = CpmPhys::funcG(k,epsCrackOnset,epsFracture,neverDamage,damLaw);
			relResidualStrength = r;
			return;
		}
	}
	throw runtime_error("CpmPhys::setRelResidualStrength: no convergence\n");
}








/********************** Law2_SomeGeom_CpmPhys_Cpm ****************************/
CREATE_LOGGER(Law2_Dem3DofGeom_CpmPhys_Cpm);
CREATE_LOGGER(Law2_ScGeom_CpmPhys_Cpm);

#ifdef YADE_CPM_FULL_MODEL_AVAILABLE
	#include"../../../brefcom-mm.hh"
#endif

// #undef CPM_MATERIAL_MODEL (force trunk version of the model)


#ifdef CPM_MATERIAL_MODEL
	#define CPM_MATERIAL_MODEL_A \
		Real& epsNPl(phys->epsNPl);\
		const Real& dt = scene->dt;\
		const Real& dmgTau(phys->dmgTau);\
		const Real& plTau(phys->plTau);\
		const Real& yieldLogSpeed(this->yieldLogSpeed);\
		const int& yieldSurfType(this->yieldSurfType);\
		const Real& yieldEllipseShift(this->yieldEllipseShift);\
		const Real& epsSoft(this->epsSoft);\
		const Real& relKnSoft(this->relKnSoft);
#else
	#define CPM_MATERIAL_MODEL_A
#endif


// constitutive law 
#ifdef CPM_MATERIAL_MODEL
	#define CPM_MATERIAL_MODEL_B \
		/* complicated version */ \
		if (epsSoft >= 0) epsN += phys->isoPrestress/E; \
		else { /* take softening into account for the prestress */ \
			Real sigmaSoft=E*epsSoft;\
			if (phys->isoPrestress >= sigmaSoft) epsN += phys -> isoPrestress/E; /* on the non-softened branch yet */ \
			/* otherwise take the regular and softened branches separately (different moduli) */ \
			else epsN += sigmaSoft/E+(phys->isoPrestress-sigmaSoft)/(E*relKnSoft);\
		} \
		CPM_MATERIAL_MODEL
#else
	#define CPM_MATERIAL_MODEL_B \
		/* simplified public model */ \
		epsN += phys->isoPrestress/E; \
		/* very simplified version of the constitutive law */ \
		kappaD = max(max(0.,epsN),kappaD); /* internal variable, max positive strain (non-decreasing) */ \
		/* Real epsFracture = crackOpening/geom->refLength; */ \
		omega = isCohesive? phys->funcG(kappaD,epsCrackOnset,epsFracture,neverDamage,damLaw) : 1.; /* damage variable (non-decreasing, as funcG is also non-decreasing) */ \
		sigmaN = (1-(epsN>0?omega:0))*E*epsN; /* damage taken in account in tension only */ \
		sigmaT = G*epsT; /* trial stress */ \
		Real yieldSigmaT = max((Real)0.,undamagedCohesion*(1-omega)-sigmaN*tanFrictionAngle); /* Mohr-Coulomb law with damage */ \
		if (sigmaT.squaredNorm() > yieldSigmaT*yieldSigmaT) { \
			Real scale = yieldSigmaT/sigmaT.norm(); \
			sigmaT *= scale; /* stress return */ \
			epsT *= scale; \
			/* epsPlSum += yieldSigmaT*geom->slipToStrainTMax(yieldSigmaT/G);*/ /* adjust strain */ \
		} \
		relResidualStrength = isCohesive? (kappaD<epsCrackOnset? 1. : (1-omega)*(kappaD)/epsCrackOnset) : 0;
#endif


#ifdef YADE_DEBUG
	#define CPM_YADE_DEBUG_A \
		if(isnan(epsN)){\
			/*LOG_FATAL("refLength="<<geom->refLength<<"; pos1="<<geom->se31.position<<"; pos2="<<geom->se32.position<<"; displacementN="<<geom->displacementN());*/ \
			throw runtime_error("!! epsN==NaN !!");\
		}
#else
	#define CPM_YADE_DEBUG_A
#endif


#define YADE_VERIFY(condition) if(!(condition)){LOG_FATAL("Verification `"<<#condition<<"' failed!"); LOG_FATAL("in interaction #"<<I->getId1()<<"+#"<<I->getId2()); Omega::instance().saveSimulation("/tmp/verificationFailed.xml"); throw;}
#define NNAN(a) YADE_VERIFY(!isnan(a));
#define NNANV(v) YADE_VERIFY(!isnan(v[0])); assert(!isnan(v[1])); assert(!isnan(v[2]));



#define CPM_GO_A \
	/* just the first time */ \
	if (I->isFresh(scene)) { \
		Vector3r pos1 = scene->bodies->operator[](I->id1)->state->pos; \
		Vector3r pos2 = scene->bodies->operator[](I->id2)->state->pos; \
		Real minRad = (geom->refR1 <= 0? geom->refR2 : (geom->refR2 <=0? geom->refR1 : min(geom->refR1,geom->refR2))); \
		Vector3r shift2 = scene->isPeriodic? Vector3r(scene->cell->hSize*I->cellDist.cast<Real>()) : Vector3r::Zero();\
		phys->refLength = (pos2 - pos1 + shift2).norm(); \
		phys->crossSection = Mathr::PI*pow(minRad,2); \
		phys->refPD = geom->refR1 + geom->refR2 - phys->refLength; \
		phys->kn = phys->crossSection*phys->E/phys->refLength; \
		phys->ks = phys->crossSection*phys->G/phys->refLength; \
		phys->epsFracture = isnan(phys->crackOpening)? phys->epsCrackOnset*phys->relDuctility : phys->crackOpening/(2*minRad); /* *geom->refLength */; \
	} \
	\
	/* shorthands */ \
	Real& epsN(phys->epsN); \
	Vector3r& epsT(phys->epsT); \
	Real& kappaD(phys->kappaD); \
	/* Real& epsPlSum(phys->epsPlSum); */ \
	const Real& E(phys->E); \
	const Real& undamagedCohesion(phys->undamagedCohesion); \
	const Real& tanFrictionAngle(phys->tanFrictionAngle); \
	const Real& G(phys->G); \
	const Real& crossSection(phys->crossSection); \
	const Real& omegaThreshold(this->omegaThreshold); \
	const Real& epsCrackOnset(phys->epsCrackOnset); \
	Real& relResidualStrength(phys->relResidualStrength); \
	/*const Real& crackOpening(phys->crackOpening); */ \
	/*const Real& relDuctility(phys->relDuctility); */ \
	const Real& epsFracture(phys->epsFracture); \
	const int& damLaw(phys->damLaw); \
	const bool& neverDamage(phys->neverDamage); \
	Real& omega(phys->omega); \
	Real& sigmaN(phys->sigmaN); \
	Vector3r& sigmaT(phys->sigmaT); \
	Real& Fn(phys->Fn); \
	Vector3r& Fs(phys->Fs); /* for python access */ \
	const bool& isCohesive(phys->isCohesive); \
	\
	CPM_MATERIAL_MODEL_A

#define CPM_GO_B \
	/* debugging */ \
	CPM_YADE_DEBUG_A \
	\
	NNAN(epsN); NNANV(epsT); \
	\
	/* constitutive law */ \
	CPM_MATERIAL_MODEL_B \
	\
	sigmaN -= phys->isoPrestress; \
	\
	NNAN(kappaD); NNAN(epsFracture); NNAN(omega); \
	NNAN(sigmaN); NNANV(sigmaT); NNAN(crossSection); \
	\
	/* handle broken contacts */ \
	if (epsN>0. && ((isCohesive && omega>omegaThreshold) || !isCohesive)) { \
		/* if (isCohesive) { */ \
			 const shared_ptr<Body>& body1 = Body::byId(I->getId1(),scene), body2 = Body::byId(I->getId2(),scene); assert(body1); assert(body2); \
			 const shared_ptr<CpmState>& st1 = YADE_PTR_CAST<CpmState>(body1->state), st2 = YADE_PTR_CAST<CpmState>(body2->state); \
			/* nice article about openMP::critical vs. scoped locks: http://www.thinkingparallel.com/2006/08/21/scoped-locking-vs-critical-in-openmp-a-personal-shootout/ */ \
			{ boost::mutex::scoped_lock lock(st1->updateMutex); st1->numBrokenCohesive += 1; /* st1->epsPlBroken += epsPlSum; */ } \
			{ boost::mutex::scoped_lock lock(st2->updateMutex); st2->numBrokenCohesive += 1; /* st2->epsPlBroken += epsPlSum; */ } \
		/* } */ \
		scene->interactions->requestErase(I); \
		return; \
	} \
	\
	Fn = sigmaN*crossSection; phys->normalForce = -Fn*geom->normal; \
	Fs = sigmaT*crossSection; phys->shearForce = -Fs;






void Law2_Dem3DofGeom_CpmPhys_Cpm::go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I){
	Dem3DofGeom* geom = static_cast<Dem3DofGeom*>(_geom.get());
	CpmPhys* phys = static_cast<CpmPhys*>(_phys.get());
	
	CPM_GO_A

	epsN = geom->strainN();
	epsT = geom->strainT();

	CPM_GO_B

	applyForceAtContactPoint(-phys->normalForce - phys->shearForce, geom->contactPoint, I->getId1(), geom->se31.position, I->getId2(), geom->se32.position);
}





void Law2_ScGeom_CpmPhys_Cpm::go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I){
	ScGeom* geom=static_cast<ScGeom*>(_geom.get());
	CpmPhys* phys=static_cast<CpmPhys*>(_phys.get());

	CPM_GO_A
	
	epsN = - (-phys->refPD + geom->penetrationDepth) / phys->refLength;
	epsT = geom->rotate(epsT);
	epsT += geom->shearIncrement() / (phys->refLength + phys->refPD) ; 

	CPM_GO_B

	Body::id_t id1 = I->getId1();
 	Body::id_t id2 = I->getId2();

	State* b1 = Body::byId(id1,scene)->state.get();
	State* b2 = Body::byId(id2,scene)->state.get();	

	Vector3r f = -phys->normalForce - phys->shearForce;
	if (!scene->isPeriodic) {
		applyForceAtContactPoint(f, geom->contactPoint , id1, b1->se3.position, id2, b2->se3.position);
	} else {
		scene->forces.addForce(id1,f);
		scene->forces.addForce(id2,-f);
		scene->forces.addTorque(id1,(geom->radius1+.5*(phys->refPD-geom->penetrationDepth))*geom->normal.cross(f));
		scene->forces.addTorque(id2,(geom->radius2+.5*(phys->refPD-geom->penetrationDepth))*geom->normal.cross(f));
	}
}








#ifdef YADE_OPENGL
	/********************** Gl1_CpmPhys ****************************/
	#include<yade/lib/opengl/OpenGLWrapper.hpp>
	#include<yade/lib/opengl/GLUtils.hpp>

	CREATE_LOGGER(Gl1_CpmPhys);

	bool Gl1_CpmPhys::contactLine=true;
	bool Gl1_CpmPhys::dmgLabel=true;
	bool Gl1_CpmPhys::dmgPlane=false;
	bool Gl1_CpmPhys::epsNLabel=true;
	bool Gl1_CpmPhys::epsT=false;
	bool Gl1_CpmPhys::epsTAxes=false;
	bool Gl1_CpmPhys::normal=false;
	Real Gl1_CpmPhys::colorStrainRatio=-1;


	void Gl1_CpmPhys::go(const shared_ptr<IPhys>& ip, const shared_ptr<Interaction>& i, const shared_ptr<Body>& b1, const shared_ptr<Body>& b2, bool wireFrame){
		const shared_ptr<CpmPhys>& phys = static_pointer_cast<CpmPhys>(ip);
		const shared_ptr<GenericSpheresContact>& geom = YADE_PTR_CAST<GenericSpheresContact>(i->geom);
		// FIXME: get the scene for periodicity; ugly!
		Scene* scene=Omega::instance().getScene().get();

		//Vector3r lineColor(phys->omega,1-phys->omega,0.0); /* damaged links red, undamaged green */
		Vector3r lineColor = Shop::scalarOnColorScale(1.-phys->relResidualStrength);

		if(colorStrainRatio>0) lineColor = Shop::scalarOnColorScale(phys->epsN/(phys->epsCrackOnset*colorStrainRatio));

		// FIXME: should be computed by the renderer; for now, use the real values
		//Vector3r pos1=geom->se31.position, pos2=geom->se32.position;
		Vector3r pos1 = scene->bodies->operator[](i->id1)->state->pos, pos2 = scene->bodies->operator[](i->id2)->state->pos;
		if (scene->isPeriodic) {
			Vector3r dPos = pos2 - pos1;
			pos1=scene->cell->wrapShearedPt(pos1);
			Vector3r shift2 = scene->isPeriodic? Vector3r(scene->cell->hSize*i->cellDist.cast<Real>()) : Vector3r::Zero();\
			pos2 = pos1 + dPos + shift2;
			//phys->refLength = (pos2 - pos1 + shift2).norm();
			//pos2=pos1+(geom->se32.position-geom->se31.position);
		}
		/*
		if (scene->isPeriodic) {
			Vector3r temp = pos2 - pos1;
			pos1 = scene->cell->wrapShearedPt(pos1);
			pos2 = pos1 + temp;
		}
		*/

		if (contactLine) GLUtils::GLDrawLine(pos1,pos2,lineColor);
		if (dmgLabel) { GLUtils::GLDrawNum(phys->omega,0.5*(pos1+pos2),lineColor); }
		else if (epsNLabel) { GLUtils::GLDrawNum(phys->epsN,0.5*(pos1+pos2),lineColor); }
		if (phys->omega>0 && dmgPlane) {
			Real halfSize = sqrt(1-phys->relResidualStrength)*.5*.705*sqrt(phys->crossSection);
			Vector3r midPt = .5*Vector3r(pos1+pos2);
			glDisable(GL_CULL_FACE);
			glPushMatrix();
				glTranslatev(midPt);
				Quaternionr q; q.setFromTwoVectors(Vector3r::UnitZ(),geom->normal);
				AngleAxisr aa(q);
				glRotatef(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
				glBegin(GL_POLYGON);
					glColor3v(lineColor); 
					glVertex3d(halfSize,0.,0.);
					glVertex3d(.5*halfSize,.866*halfSize,0.);
					glVertex3d(-.5*halfSize,.866*halfSize,0.);
					glVertex3d(-halfSize,0.,0.);
					glVertex3d(-.5*halfSize,-.866*halfSize,0.);
					glVertex3d(.5*halfSize,-.866*halfSize,0.);
				glEnd();
			glPopMatrix();
		}

		Vector3r cp = static_pointer_cast<GenericSpheresContact>(i->geom)->contactPoint;
		if (scene->isPeriodic) {cp = scene->cell->wrapShearedPt(cp);}
		if (epsT) {
			Real maxShear = (phys->undamagedCohesion-phys->sigmaN*phys->tanFrictionAngle)/phys->G;
			Real relShear = phys->epsT.norm()/maxShear;
			Real scale = phys->refLength;
			Vector3r dirShear = phys->epsT; dirShear.normalize();
			if(epsTAxes){
				GLUtils::GLDrawLine(cp-Vector3r(scale,0,0),cp+Vector3r(scale,0,0));
				GLUtils::GLDrawLine(cp-Vector3r(0,scale,0),cp+Vector3r(0,scale,0));
				GLUtils::GLDrawLine(cp-Vector3r(0,0,scale),cp+Vector3r(0,0,scale));
			}
			GLUtils::GLDrawArrow(cp,cp+dirShear*relShear*scale,Vector3r(1.,0.,0.));
			GLUtils::GLDrawLine(cp+dirShear*relShear*scale,cp+dirShear*scale,Vector3r(.3,.3,.3));

			/* normal strain */ GLUtils::GLDrawArrow(cp,cp+geom->normal*(phys->epsN/maxShear),Vector3r(0.,1.,0.));
		}
		//if(normal) GLUtils::GLDrawArrow(cp,cp+geom->normal*.5*phys->equilibriumDist,Vector3r(0.,1.,0.));
	}
#endif










/********************** CpmStateUpdater ****************************/
CREATE_LOGGER(CpmStateUpdater);
//Real CpmStateUpdater::maxOmega=0.;
//Real CpmStateUpdater::avgRelResidual=0.;

void CpmStateUpdater::update(Scene* _scene){
	Scene *scene = _scene? _scene : Omega::instance().getScene().get();
	vector<BodyStats> bodyStats; bodyStats.resize(scene->bodies->size());
	assert(bodyStats[0].nCohLinks == 0); // should be initialized by dfault ctor
	avgRelResidual = 0; Real nAvgRelResidual = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if (!I) continue;
		if (!I->isReal()) continue;
		shared_ptr<CpmPhys> phys = dynamic_pointer_cast<CpmPhys>(I->phys);
		if (!phys) continue;
		const Body::id_t id1 = I->getId1(), id2 = I->getId2();
		GenericSpheresContact* geom = YADE_CAST<GenericSpheresContact*>(I->geom.get());

		const Vector3r& n = geom->normal;
		const Real& Fn = phys->Fn;
		const Vector3r& Fs = phys->Fs;
		//stress[i,j] += geom->refLength*(Fn*n[i]*n[j]+0.5*(Fs[i]*n[j]+Fs[j]*n[i]));
		//stress += geom->refLength*(Fn*outer(n,n)+.5*(outer(Fs,n)+outer(n,Fs)));
		Matrix3r stress = phys->refLength*(Fn*n*n.transpose()+.5*(Fs*n.transpose()+n*Fs.transpose()));
		
		bodyStats[id1].stress += stress;
		bodyStats[id2].stress += stress;
		bodyStats[id1].nLinks++; bodyStats[id2].nLinks++;
		
		if (!phys->isCohesive) continue;
		bodyStats[id1].nCohLinks++; bodyStats[id1].dmgSum += (1-phys->relResidualStrength); // bodyStats[id1].epsPlSum += phys->epsPlSum;
		bodyStats[id2].nCohLinks++; bodyStats[id2].dmgSum += (1-phys->relResidualStrength); // bodyStats[id2].epsPlSum += phys->epsPlSum;
		maxOmega = max(maxOmega,phys->omega);
		avgRelResidual += phys->relResidualStrength;
		nAvgRelResidual += 1;
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				bodyStats[id1].dmgRhs += n*n.transpose()*(1-phys->relResidualStrength);
				bodyStats[id2].dmgRhs += n*n.transpose()*(1-phys->relResidualStrength);
			}
		}
	}
	Matrix3r identity = Matrix3r::Identity();
	Real tr;
	FOREACH(shared_ptr<Body> B, *scene->bodies){
		if (!B) continue;
		const Body::id_t& id = B->getId();
		// add damaged contacts that have already been deleted
		CpmState* state = dynamic_cast<CpmState*>(B->state.get());
		if (!state) continue;
		state->stress = bodyStats[id].stress;
		int cohLinksWhenever = bodyStats[id].nCohLinks+state->numBrokenCohesive;
		if (cohLinksWhenever>0) {
			state->normDmg = (bodyStats[id].dmgSum+state->numBrokenCohesive)/cohLinksWhenever;
			// state->normEpsPl = (bodyStats[id].epsPlSum+state->epsPlBroken)/cohLinksWhenever;
			if (state->normDmg>1) {
				LOG_WARN("#"<<id<<" normDmg="<<state->normDmg<<" nCohLinks="<<bodyStats[id].nCohLinks<<", numBrokenCohesive="<<state->numBrokenCohesive<<", dmgSum="<<bodyStats[id].dmgSum<<", numAllCohLinks"<<cohLinksWhenever);
			}
			Matrix3r& dmgRhs = bodyStats[id].dmgRhs;
			dmgRhs *= 15./cohLinksWhenever;
			tr = 3*state->normDmg;
			state->damageTensor = .5 * (dmgRhs - tr*identity);
		}
		else { state->normDmg = 0; /*state->normEpsPl=0;*/ state->damageTensor = Matrix3r::Zero(); }
		B->shape->color = Vector3r(state->normDmg,1-state->normDmg,B->state->blockedDOFs==State::DOF_ALL?0:1);
		nAvgRelResidual += 0.5*state->numBrokenCohesive; // add half or broken interactions, other body has the other half
		Sphere* sphere = dynamic_cast<Sphere*>(B->shape.get());
		if(!sphere) continue;
		Real& r = sphere->radius;
		state->stress = bodyStats[id].stress/(4/3.*Mathr::PI*r*r*r/.62)*.5;
	}
	avgRelResidual /= nAvgRelResidual;
}





#undef YADE_VERIFY
#undef NNAN
#undef NNANV

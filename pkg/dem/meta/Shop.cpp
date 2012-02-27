// 2007 © Václav Šmilauer <eudoxos@arcig.cz>
#include"Shop.hpp"

#include<limits>

#include<boost/filesystem/convenience.hpp>
#include<boost/tokenizer.hpp>
#include<boost/tuple/tuple.hpp>

#include<yade/core/Scene.hpp>
#include<yade/core/Body.hpp>

#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>

#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/ElastMat.hpp>

#include<yade/pkg-common/Bo1_Sphere_Aabb.hpp>
#include<yade/pkg-common/Bo1_Box_Aabb.hpp>
#include<yade/pkg-dem/NewtonIntegrator.hpp>
#include<yade/pkg-dem/Ig2_Sphere_Sphere_ScGeom.hpp>
#include<yade/pkg-dem/Ig2_Box_Sphere_ScGeom.hpp>
#include<yade/pkg-dem/Ip2_FrictMat_FrictMat_FrictPhys.hpp>
//#include<yade/pkg-dem/SimpleViscoelasticBodyParameters.hpp>
#include<yade/pkg-dem/ViscoelasticPM.hpp>
/*class Bo1_Sphere_Aabb;
class Bo1_Box_Aabb;
class SceneShape; */


#include<yade/pkg-common/ForceResetter.hpp>

#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/pkg-common/InteractionDispatchers.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/GravityEngines.hpp>

#include<yade/pkg-dem/GlobalStiffnessTimeStepper.hpp>
#include<yade/pkg-dem/ElasticContactLaw.hpp>

#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>


#include<yade/pkg-dem/Tetra.hpp>

#include<boost/foreach.hpp>
#ifndef FOREACH
	#define FOREACH BOOST_FOREACH
#endif


#define _SPEC_CAST(orig,cast) template<> void Shop::setDefault<orig>(string key, orig val){setDefault(key,cast(val));}
_SPEC_CAST(const char*,string);
_SPEC_CAST(char*,string);
#undef _SPEC_CAST

CREATE_LOGGER(Shop);
// YADE_PLUGIN((Shop));

map<string,boost::any> Shop::defaults;

/*! Flip periodic cell by given number of cells.

Still broken, some interactions are missed. Should be checked.
*/

Matrix3r Shop::flipCell(const Matrix3r& _flip){
	Scene* scene=Omega::instance().getScene().get(); const shared_ptr<Cell>& cell(scene->cell); const Matrix3r& trsf(cell->trsf);
	Vector3r size=cell->getSize();
	Matrix3r flip;
	if(_flip==Matrix3r::Zero()){
		bool hasNonzero=false;
		for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
			if(i==j){ flip(i,j)=0; continue; }
			flip(i,j)=-floor(.5+trsf(i,j)/(size[j]/size[i]));
			if(flip(i,j)!=0) hasNonzero=true;
		}
		if(!hasNonzero) {LOG_TRACE("No flip necessary."); return Matrix3r::Zero();}
		LOG_DEBUG("Computed flip matrix: upper "<<flip(0,1)<<","<<flip(0,2)<<","<<flip(1,2)<<"; lower "<<flip(1,0)<<","<<flip(2,0)<<","<<flip(2,1));
	} else {
		flip=_flip;
	}

	// current cell coords of bodies
	vector<Vector3i > oldCells; oldCells.resize(scene->bodies->size());
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue; cell->wrapShearedPt(b->state->pos,oldCells[b->getId()]);
	}

	// change cell trsf here
	Matrix3r trsfInc;
	for(int i=0; i<3; i++) for(int j=0; j<3; j++){
		if(i==j) { if(flip(i,j)!=0) LOG_WARN("Non-zero diagonal term at ["<<i<<","<<j<<"] is meaningless and will be ignored."); trsfInc(i,j)=0; continue; }
		// make sure non-diagonal entries are "integers"
		if(flip(i,j)!=double(int(flip(i,j)))) LOG_WARN("Flip matrix entry "<<flip(i,j)<<" at ["<<i<<","<<j<<"] not integer?! (will be rounded)");
		trsfInc(i,j)=int(flip(i,j))*size[j]/size[i];
	}
	TRWM3MAT(cell->trsf);
	TRWM3MAT(trsfInc);
	cell->trsf+=trsfInc;
	cell->postProcessAttributes(/*deserializing*/true);

	// new cell coords of bodies
	vector<Vector3i > newCells; newCells.resize(scene->bodies->size());
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue;
		cell->wrapShearedPt(b->state->pos,newCells[b->getId()]);
	}

	// remove all potential interactions
	scene->interactions->eraseNonReal();
	// adjust Interaction::cellDist for real interactions; 
	FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
		body_id_t id1=i->getId1(),id2=i->getId2();
		// this must be the same for both old and new interaction: cell2-cell1+cellDist
		// c₂-c₁+c₁₂=k; c'₂+c₁'+c₁₂'=k   (new cell coords have ')
		// c₁₂'=(c₂-c₁+c₁₂)-(c₂'-c₁')
		i->cellDist=(oldCells[id2]-oldCells[id1]+i->cellDist)-(newCells[id2]-newCells[id1]);
	}


	// force reinitialization of the collider
	bool colliderFound=false;
	FOREACH(const shared_ptr<Engine>& e, scene->engines){
		Collider* c=dynamic_cast<Collider*>(e.get());
		if(c){ colliderFound=true; c->invalidatePersistentData(); }
	}
	if(!colliderFound) LOG_WARN("No collider found while flipping cell; continuing simulation might give garbage results.");
	return flip;
}

/* Apply force on contact point to 2 bodies; the force is oriented as it applies on the first body and is reversed on the second.
 */
void Shop::applyForceAtContactPoint(const Vector3r& force, const Vector3r& contPt, body_id_t id1, const Vector3r& pos1, body_id_t id2, const Vector3r& pos2, Scene* rootBody){
	rootBody->forces.addForce(id1,force);
	rootBody->forces.addForce(id2,-force);
	rootBody->forces.addTorque(id1,(contPt-pos1).cross(force));
	rootBody->forces.addTorque(id2,-(contPt-pos2).cross(force));
}


/*! Compute sum of forces in the whole simulation and averages stiffness.

Designed for being used with periodic cell, where diving the resulting components by
areas of the cell will give average stress in that direction.

Requires all .isReal() interaction to have interactionPhysics deriving from NormShearPhys.
*/
Vector3r Shop::totalForceInVolume(Real& avgIsoStiffness, Scene* _rb){
	Scene* rb=_rb ? _rb : Omega::instance().getScene().get();
	Vector3r force(Vector3r::Zero()); Real stiff=0; long n=0;
	FOREACH(const shared_ptr<Interaction>&I, *rb->interactions){
		if(!I->isReal()) continue;
		NormShearPhys* nsi=YADE_CAST<NormShearPhys*>(I->interactionPhysics.get());
		force+=Vector3r(abs(nsi->normalForce[0]+nsi->shearForce[0]),abs(nsi->normalForce[1]+nsi->shearForce[1]),abs(nsi->normalForce[2]+nsi->shearForce[2]));
		stiff+=(1/3.)*nsi->kn+(2/3.)*nsi->ks; // count kn in one direction and ks in the other two
		n++;
	}
	avgIsoStiffness= n>0 ? (1./n)*stiff : -1;
	return force;
}

Real Shop::unbalancedForce(bool useMaxForce, Scene* _rb){
	Scene* rb=_rb ? _rb : Omega::instance().getScene().get();
	rb->forces.sync();
	// get maximum force on a body and sum of all forces (for averaging)
	Real sumF=0,maxF=0,currF; int nb=0;
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b->isDynamic) continue;
		currF=rb->forces.getForce(b->id).norm(); maxF=max(currF,maxF); sumF+=currF; nb++;
	}
	Real meanF=sumF/nb;
	// get max force on contacts
	// Real maxContactF=0; 										//It gave a warning. Anton Gladky.
	sumF=0; nb=0;
	FOREACH(const shared_ptr<Interaction>& I, *rb->interactions){
		if(!I->isReal()) continue;
		shared_ptr<NormShearPhys> nsi=YADE_PTR_CAST<NormShearPhys>(I->interactionPhysics); assert(nsi);
		sumF+=(nsi->normalForce+nsi->shearForce).norm(); nb++;
	}
	sumF/=nb;
	return (useMaxForce?maxF:meanF)/(sumF);
}

Real Shop::kineticEnergy(Scene* _rb){
	Scene* rb=_rb ? _rb : Omega::instance().getScene().get();
	Real ret=0.;
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b->isDynamic) continue;
		// ½(mv²+ωIω)
		ret+=.5*(b->state->mass*b->state->vel.squaredNorm()+b->state->angVel.dot(diagMult(b->state->inertia,b->state->angVel)));
	}
	return ret;
}





template <typename valType> valType Shop::getDefault(const string& key) {
	ensureInit();
	try{return boost::any_cast<valType>(defaults[key]);}
	catch(boost::bad_any_cast& e){
		LOG_FATAL("Cast error while getting key `"<<key<<"' of type `"<<typeid(valType).name()<<"' ("<<e.what()<<").");
		if(!hasDefault(key) || (defaults[key].type()==typeid(void))){LOG_FATAL("Key `"<<key<<"' not defined in the map.");}
		else{LOG_INFO("(key `"<<key<<"' exists and if of type `"<<defaults[key].type().name()<<"').");}
		throw;
	}
}

void Shop::init(){
	//LOG_INFO("Container length is "<<defaults.size()<<endl);
	/*for(map<string,boost::any>::iterator I=defaults.begin(); I!=defaults.end(); I++){
		LOG_INFO("Key `"<<I->first<<", type `"<<I->second.type().name()<<"'.");
	}*/

	defaults["container_is_not_empty"]=boost::any(0); // prevent loops from ensureInit();

	setDefault("phys_density",2e3);
	setDefault("phys_young",30e9);
	setDefault("phys_poisson",.3);
	setDefault("phys_frictionAngle",0.5236); //30˚
	setDefault("phys_se3_orientation",Quaternionr::Identity());

	setDefault("aabb_randomColor",false);
	setDefault("aabb_color",Vector3r(0,1,0));

	setDefault("mold_randomColor",true);
	setDefault("mold_color",Vector3r(1,0,0));

	setDefault("shape_randomColor",true);
	setDefault("shape_color",Vector3r(0,0,1));

	setDefault("shape_wire",false);
	setDefault("shape_shadowCaster",true);

	setDefault("param_damping",.2);
	setDefault("param_gravity",Vector3r(0,0,-10));
	setDefault<int>("param_timeStepUpdateInterval",300);

	setDefault("param_pythonInitExpr",string("print 'Hello world!'"));
	setDefault("param_pythonRunExpr",string(""));

}

/*! Create root body. */
shared_ptr<Scene> Shop::rootBody(){
	return shared_ptr<Scene>(new Scene);
}


/*! Assign default set of actors (initializers and engines) to an existing Scene.
 */
void Shop::rootBodyActors(shared_ptr<Scene> rootBody){
	// initializers	
	rootBody->initializers.clear();

	shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
	boundDispatcher->add(new Bo1_Sphere_Aabb);
	boundDispatcher->add(new Bo1_Box_Aabb);
	boundDispatcher->add(new Bo1_Tetra_Aabb);
	rootBody->initializers.push_back(boundDispatcher);

	//engines
	rootBody->engines.clear();

	if(getDefault<int>("param_timeStepUpdateInterval")>0){
		shared_ptr<GlobalStiffnessTimeStepper> sdecTimeStepper(new GlobalStiffnessTimeStepper);
		sdecTimeStepper->timeStepUpdateInterval=getDefault<int>("param_timeStepUpdateInterval");
		sdecTimeStepper->timeStepUpdateInterval=300;
		rootBody->engines.push_back(sdecTimeStepper);
	}

	rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));

	rootBody->engines.push_back(boundDispatcher);

	rootBody->engines.push_back(shared_ptr<Engine>(new InsertionSortCollider));

	shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);
	interactionGeometryDispatcher->add(new Ig2_Sphere_Sphere_ScGeom);
	interactionGeometryDispatcher->add(new Ig2_Box_Sphere_ScGeom);
	interactionGeometryDispatcher->add(new Ig2_Tetra_Tetra_TTetraGeom);
	rootBody->engines.push_back(interactionGeometryDispatcher);

	shared_ptr<InteractionPhysicsDispatcher> interactionPhysicsDispatcher(new InteractionPhysicsDispatcher);
	interactionPhysicsDispatcher->add(new Ip2_FrictMat_FrictMat_FrictPhys);
	rootBody->engines.push_back(interactionPhysicsDispatcher);
		
	shared_ptr<ElasticContactLaw> constitutiveLaw(new ElasticContactLaw);
	rootBody->engines.push_back(constitutiveLaw);

	if(getDefault<Vector3r>("param_gravity")!=Vector3r(0,0,0)){
		shared_ptr<GravityEngine> gravityCondition(new GravityEngine);
		gravityCondition->gravity=getDefault<Vector3r>("param_gravity");
		rootBody->engines.push_back(gravityCondition);
	}
	
	shared_ptr<NewtonIntegrator> newton(new NewtonIntegrator);
	newton->damping=max(getDefault<double>("param_damping"),0.);
	rootBody->engines.push_back(newton);
}


shared_ptr<FrictMat> Shop::defaultGranularMat(){
	shared_ptr<FrictMat> mat(new FrictMat);
	mat->density=getDefault<double>("phys_density");
	mat->young=getDefault<double>("phys_young");
	mat->poisson=getDefault<double>("phys_poisson");
	mat->frictionAngle=getDefault<double>("phys_frictionAngle");
	return mat;
}

/*! Create body - sphere. */
shared_ptr<Body> Shop::sphere(Vector3r center, Real radius, shared_ptr<Material> mat){
	shared_ptr<Body> body(new Body);
	body->isDynamic=true;
	body->material=mat ? mat : static_pointer_cast<Material>(defaultGranularMat());
	body->state->pos=center;
	body->state->mass=4.0/3.0*Mathr::PI*radius*radius*radius*body->material->density;
	body->state->inertia=Vector3r(2.0/5.0*body->state->mass*radius*radius,2.0/5.0*body->state->mass*radius*radius,2.0/5.0*body->state->mass*radius*radius);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->shape=shared_ptr<Sphere>(new Sphere(radius));
	return body;
}

/*! Create body - box. */
shared_ptr<Body> Shop::box(Vector3r center, Vector3r extents, shared_ptr<Material> mat){
	shared_ptr<Body> body(new Body);
	body->isDynamic=true;
	body->material=mat ? mat : static_pointer_cast<Material>(defaultGranularMat());
	body->state->pos=center;
	Real mass=8.0*extents[0]*extents[1]*extents[2]*body->material->density;
	body->state->mass=mass;
	body->state->inertia=Vector3r(mass*(4*extents[1]*extents[1]+4*extents[2]*extents[2])/12.,mass*(4*extents[0]*extents[0]+4*extents[2]*extents[2])/12.,mass*(4*extents[0]*extents[0]+4*extents[1]*extents[1])/12.);
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->shape=shared_ptr<Box>(new Box(extents));
	return body;
}

/*! Create body - tetrahedron. */
shared_ptr<Body> Shop::tetra(Vector3r v_global[4], shared_ptr<Material> mat){
	shared_ptr<Body> body(new Body);
	body->isDynamic=true;
	body->material=mat ? mat : static_pointer_cast<Material>(defaultGranularMat());
	Vector3r centroid=(v_global[0]+v_global[1]+v_global[2]+v_global[3])*.25;
	Vector3r v[4]; for(int i=0; i<4; i++) v[i]=v_global[i]-centroid;
	body->state->pos=centroid;
	body->state->mass=body->material->density*TetrahedronVolume(v);
	// inertia will be calculated below, by TetrahedronWithLocalAxesPrincipal
	body->bound=shared_ptr<Aabb>(new Aabb);
	body->shape=shared_ptr<Tetra>(new Tetra(v[0],v[1],v[2],v[3]));
	// make local axes coincident with principal axes
	TetrahedronWithLocalAxesPrincipal(body);
	return body;
}


void Shop::saveSpheresToFile(string fname){
	const shared_ptr<Scene>& rootBody=Omega::instance().getScene();
	ofstream f(fname.c_str());
	if(!f.good()) throw runtime_error("Unable to open file `"+fname+"'");

	FOREACH(shared_ptr<Body> b, *rootBody->bodies){
		if (!b->isDynamic) continue;
		shared_ptr<Sphere>	intSph=dynamic_pointer_cast<Sphere>(b->shape);
		if(!intSph) continue;
		const Vector3r& pos=b->state->pos;
		f<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" "<<intSph->radius<<endl; // <<" "<<1<<" "<<1<<endl;
	}
	f.close();
}

Real Shop::getSpheresVolume(){
	const shared_ptr<Scene>& rootBody=Omega::instance().getScene();
	Real vol=0;
	FOREACH(shared_ptr<Body> b, *rootBody->bodies){
		if (!b->isDynamic) continue;
		shared_ptr<Sphere> intSph=YADE_PTR_CAST<Sphere>(b->shape);
		vol += 4.18879020*pow(intSph->radius,3);
	}
	return vol;
}


vector<pair<Vector3r,Real> > Shop::loadSpheresFromFile(string fname, Vector3r& minXYZ, Vector3r& maxXYZ){
	if(!boost::filesystem::exists(fname)) {
		throw std::invalid_argument(string("File with spheres `")+fname+"' doesn't exist.");
	}
	vector<pair<Vector3r,Real> > spheres;
	ifstream sphereFile(fname.c_str());
	if(!sphereFile.good()) throw std::runtime_error("File with spheres `"+fname+"' couldn't be opened.");
	Vector3r C;
	Real r=0;
	string line;
	size_t lineNo=0;
	while(std::getline(sphereFile, line, '\n')){
		lineNo++;
		boost::tokenizer<boost::char_separator<char> > toks(line,boost::char_separator<char>(" \t"));
		int i=0;
		FOREACH(const string& s, toks){
			if(i<3) C[i]=lexical_cast<Real>(s);
			if(i==3) r=lexical_cast<Real>(s);
			i++;
		}
		if(i==0) continue; // empty line, skipped (can be the last one)
		// Wenjie's format: 5 per line are boxes, which should be skipped
		if(i==5) continue;
		if((i!=4) and (i!=6)) {
			LOG_ERROR("Line "+lexical_cast<string>(lineNo)+" in the spheres file "+fname+" has "+lexical_cast<string>(i)+" columns instead of 0,4,5 or 6.");
			LOG_ERROR("The result may be garbage!");
		}
		for(int j=0; j<3; j++) { minXYZ[j]=(spheres.size()>0?min(C[j]-r,minXYZ[j]):C[j]-r); maxXYZ[j]=(spheres.size()>0?max(C[j]+r,maxXYZ[j]):C[j]+r);}
		spheres.push_back(pair<Vector3r,Real>(C,r));
	}
	return spheres;
}


vector<pair<Vector3r,Real> > Shop::loadSpheresSmallSdecXyz(Vector3r& minXYZ, Vector3r& maxXYZ){
	int i=-1;
	vector<pair<Vector3r,Real> > spheres;
	while(Shop::smallSdecXyzData[++i].r>0){
		Vector3r C(smallSdecXyzData[i].C0,smallSdecXyzData[i].C1,smallSdecXyzData[i].C2); Real r(smallSdecXyzData[i].r);
		for(int j=0; j<3; j++) { minXYZ[j]=(spheres.size()>0?min(C[j]-r,minXYZ[j]):C[j]-r); maxXYZ[j]=(spheres.size()>0?max(C[j]+r,maxXYZ[j]):C[j]+r);}
		spheres.push_back(pair<Vector3r,Real>(C,r));
		//LOG_DEBUG("Loaded sphere with radius "<<r<<" at "<<C);
	}
	return spheres;
}

Real Shop::PWaveTimeStep(shared_ptr<Scene> _rb){
	shared_ptr<Scene> rb=_rb;
	if(!rb)rb=Omega::instance().getScene();
	Real dt=std::numeric_limits<Real>::infinity();
	FOREACH(const shared_ptr<Body>& b, *rb->bodies){
		if(!b->material || !b->shape) continue;
		shared_ptr<ElastMat> ebp=dynamic_pointer_cast<ElastMat>(b->material);
		shared_ptr<Sphere> s=dynamic_pointer_cast<Sphere>(b->shape);
		if(!ebp || !s) continue;
		Real density=b->state->mass/((4/3.)*Mathr::PI*pow(s->radius,3));
		dt=min(dt,s->radius/sqrt(ebp->young/density));
	}
	return dt;
}
/* Project 3d point into 2d using spiral projection along given axis;
 * the returned tuple is
 * 	
 *  (height relative to the spiral, distance from axis, theta )
 *
 * dH_dTheta is the inclination of the spiral (height increase per radian),
 * theta0 is the angle for zero height (by given axis).
 */
boost::tuple<Real,Real,Real> Shop::spiralProject(const Vector3r& pt, Real dH_dTheta, int axis, Real periodStart, Real theta0){
	int ax1=(axis+1)%3,ax2=(axis+2)%3;
	Real r=sqrt(pow(pt[ax1],2)+pow(pt[ax2],2));
	Real theta;
	if(r>Mathr::ZERO_TOLERANCE){
		theta=acos(pt[ax1]/r);
		if(pt[ax2]<0) theta=Mathr::TWO_PI-theta;
	}
	else theta=0;
	Real hRef=dH_dTheta*(theta-theta0);
	long period;
	if(isnan(periodStart)){
		Real h=Shop::periodicWrap(pt[axis]-hRef,hRef-Mathr::PI*dH_dTheta,hRef+Mathr::PI*dH_dTheta,&period);
		return boost::make_tuple(r,h,theta);
	}
	else{
		// Real hPeriodStart=(periodStart-theta0)*dH_dTheta;
		//TRVAR4(hPeriodStart,periodStart,theta0,theta);
		//Real h=Shop::periodicWrap(pt[axis]-hRef,hPeriodStart,hPeriodStart+2*Mathr::PI*dH_dTheta,&period);
		theta=Shop::periodicWrap(theta,periodStart,periodStart+2*Mathr::PI,&period);
		Real h=pt[axis]-hRef+period*2*Mathr::PI*dH_dTheta;
		//TRVAR3(pt[axis],pt[axis]-hRef,period);
		//TRVAR2(h,theta);
		return boost::make_tuple(r,h,theta);
	}
}

shared_ptr<Interaction> Shop::createExplicitInteraction(body_id_t id1, body_id_t id2, bool force){
	InteractionGeometryDispatcher* geomMeta=NULL;
	InteractionPhysicsDispatcher* physMeta=NULL;
	shared_ptr<Scene> rb=Omega::instance().getScene();
	if(rb->interactions->find(body_id_t(id1),body_id_t(id2))!=0) throw runtime_error(string("Interaction #")+lexical_cast<string>(id1)+"+#"+lexical_cast<string>(id2)+" already exists.");
	FOREACH(const shared_ptr<Engine>& e, rb->engines){
		if(!geomMeta) { geomMeta=dynamic_cast<InteractionGeometryDispatcher*>(e.get()); if(geomMeta) continue; }
		if(!physMeta) { physMeta=dynamic_cast<InteractionPhysicsDispatcher*>(e.get()); if(physMeta) continue; }
		InteractionDispatchers* id(dynamic_cast<InteractionDispatchers*>(e.get()));
		if(id){ geomMeta=id->geomDispatcher.get(); physMeta=id->physDispatcher.get(); }
		if(geomMeta&&physMeta){break;}
	}
	if(!geomMeta) throw runtime_error("No InteractionGeometryDispatcher in engines or inside InteractionDispatchers.");
	if(!physMeta) throw runtime_error("No InteractionPhysicsDispatcher in engines or inside InteractionDispatchers.");
	shared_ptr<Body> b1=Body::byId(id1,rb), b2=Body::byId(id2,rb);
	if(!b1) throw runtime_error(("No body #"+lexical_cast<string>(id1)).c_str());
	if(!b2) throw runtime_error(("No body #"+lexical_cast<string>(id2)).c_str());
	shared_ptr<Interaction> i=geomMeta->explicitAction(b1,b2,/*force*/force);
	assert(force && i);
	if(!i) return i;
	physMeta->explicitAction(b1->material,b2->material,i);
	rb->interactions->insert(i);
	return i;
}

Shop::sphereGeomStruct Shop::smallSdecXyzData[]={
	{0.027814,0.028311,0.053055,0.006921},
	{0.040115,0.039488,0.051634,0.009072},
	{0.014935,0.042338,0.065370,0.014248},
	{0.040157,0.020724,0.049777,0.007933},
	{0.054315,0.030277,0.034105,0.015247},
	{0.084782,0.035042,0.040466,0.006506},
	{0.063445,0.049730,0.067288,0.011911},
	{0.068606,0.025283,0.065262,0.012022},
	{0.088031,0.014291,0.039471,0.014521},
	{0.116754,0.037915,0.051292,0.008274},
	{0.100791,0.035164,0.067461,0.013108},
	{0.114599,0.012268,0.055589,0.015359},
	{0.138767,0.018895,0.036148,0.009280},
	{0.145595,0.044496,0.031125,0.015939},
	{0.132246,0.042180,0.051209,0.007763},
	{0.143392,0.023461,0.060496,0.015920},
	{0.171105,0.022212,0.032836,0.012108},
	{0.181104,0.037549,0.046727,0.010873},
	{0.166356,0.047869,0.055943,0.008384},
	{0.182075,0.018453,0.050915,0.008701},
	{0.199011,0.025231,0.029222,0.016017},
	{0.223922,0.034850,0.046495,0.009033},
	{0.201425,0.039383,0.074087,0.010160},
	{0.208859,0.020393,0.059165,0.015109},
	{0.247185,0.029358,0.040632,0.009137},
	{0.256247,0.054160,0.035898,0.007448},
	{0.241821,0.048789,0.051785,0.013900},
	{0.262536,0.027938,0.052724,0.007258},
	{0.266457,0.015729,0.044737,0.007849},
	{0.300361,0.040203,0.038619,0.010238},
	{0.269716,0.038512,0.048249,0.006284},
	{0.284778,0.015856,0.057268,0.010625},
	{0.308568,0.015046,0.037486,0.015678},
	{0.324288,0.040734,0.028459,0.015762},
	{0.312379,0.037454,0.070796,0.013332},
	{0.325630,0.017879,0.071709,0.009211},
	{0.337472,0.021731,0.049983,0.015831},
	{0.363910,0.040206,0.054198,0.013145},
	{0.347960,0.050663,0.050631,0.006258},
	{0.353319,0.027257,0.066410,0.007567},
	{0.382606,0.035870,0.041879,0.009000},
	{0.403064,0.040715,0.052848,0.013380},
	{0.378583,0.049140,0.071144,0.009526},
	{0.392096,0.013326,0.048449,0.016326},
	{0.419694,0.020328,0.050734,0.012238},
	{0.435554,0.040547,0.055823,0.013440},
	{0.416138,0.032889,0.071377,0.008786},
	{0.438445,0.022231,0.064865,0.007190},
	{0.449096,0.034012,0.035415,0.006513},
	{0.465516,0.041105,0.049831,0.016459},
	{0.451849,0.034788,0.071794,0.010092},
	{0.470798,0.014474,0.061308,0.013016},
	{0.022520,0.051083,0.034275,0.015762},
	{0.051991,0.077988,0.036367,0.009543},
	{0.020986,0.077522,0.050688,0.015394},
	{0.045377,0.061146,0.050244,0.013259},
	{0.070165,0.058855,0.037378,0.014762},
	{0.079072,0.088030,0.037625,0.009745},
	{0.070175,0.083486,0.059729,0.008680},
	{0.080809,0.052527,0.057405,0.008263},
	{0.102876,0.053069,0.047621,0.012600},
	{0.129932,0.071275,0.041493,0.015762},
	{0.100905,0.082761,0.049471,0.015577},
	{0.114885,0.069784,0.065561,0.008883},
	{0.124029,0.048642,0.037433,0.006910},
	{0.151400,0.072980,0.032672,0.007510},
	{0.125754,0.087963,0.056533,0.007065},
	{0.151978,0.057727,0.065740,0.009624},
	{0.170038,0.050106,0.034273,0.009151},
	{0.180330,0.067662,0.044758,0.013741},
	{0.163195,0.080867,0.050099,0.008542},
	{0.193257,0.052669,0.058994,0.009924},
	{0.213500,0.053456,0.046639,0.012293},
	{0.214894,0.078031,0.032394,0.009344},
	{0.204112,0.073752,0.062427,0.013778},
	{0.209893,0.055436,0.071207,0.007226},
	{0.229179,0.063602,0.037779,0.008378},
	{0.253466,0.070318,0.046778,0.010045},
	{0.243394,0.085523,0.061896,0.013644},
	{0.268095,0.055682,0.049778,0.010865},
	{0.272490,0.056346,0.032166,0.007300},
	{0.293444,0.075937,0.050821,0.007124},
	{0.268497,0.087864,0.061002,0.010071},
	{0.285287,0.053447,0.056283,0.007652},
	{0.303399,0.065465,0.032905,0.015762},
	{0.331186,0.075873,0.030308,0.014024},
	{0.305661,0.077322,0.067721,0.013776},
	{0.321014,0.057758,0.055215,0.012434},
	{0.340583,0.055745,0.032542,0.006766},
	{0.363513,0.081185,0.053753,0.009367},
	{0.342074,0.072811,0.056793,0.013501},
	{0.362170,0.061143,0.054690,0.007397},
	{0.370272,0.074686,0.039972,0.007302},
	{0.402953,0.071810,0.034403,0.012804},
	{0.382151,0.086219,0.067255,0.009054},
	{0.382421,0.059723,0.048574,0.013805},
	{0.418492,0.059752,0.050291,0.011257},
	{0.431649,0.077456,0.029368,0.012777},
	{0.423405,0.082143,0.054675,0.011680},
	{0.438239,0.060029,0.044749,0.009130},
	{0.451129,0.062181,0.031179,0.009710},
	{0.479146,0.076782,0.041030,0.013576},
	{0.450502,0.080910,0.050514,0.015762},
	{0.469945,0.064905,0.054090,0.006331},
	{0.013563,0.097184,0.042932,0.007007},
	{0.045430,0.112786,0.026860,0.013686},
	{0.026238,0.108238,0.071820,0.008664},
	{0.045567,0.084792,0.054663,0.010547},
	{0.055851,0.098322,0.045255,0.008878},
	{0.092159,0.113570,0.039595,0.014076},
	{0.042022,0.114607,0.055621,0.014833},
	{0.076167,0.098666,0.069346,0.010212},
	{0.091668,0.095217,0.029729,0.006767},
	{0.114684,0.116416,0.048755,0.009589},
	{0.101121,0.106787,0.066807,0.014051},
	{0.117768,0.092385,0.068330,0.006541},
	{0.136812,0.102557,0.043292,0.008163},
	{0.140060,0.121797,0.037341,0.012236},
	{0.129633,0.101076,0.058761,0.006790},
	{0.148062,0.103611,0.061993,0.009901},
	{0.154918,0.097447,0.035701,0.011902},
	{0.173991,0.107853,0.041188,0.010506},
	{0.163383,0.113509,0.070590,0.007969},
	{0.183313,0.096713,0.054555,0.008896},
	{0.193900,0.091794,0.032836,0.015762},
	{0.221431,0.105734,0.041347,0.016142},
	{0.210572,0.112748,0.064643,0.009085},
	{0.220128,0.094827,0.066506,0.011310},
	{0.244598,0.102360,0.048021,0.008201},
	{0.260796,0.107346,0.041782,0.007826},
	{0.247778,0.114675,0.056955,0.007342},
	{0.263854,0.105240,0.057300,0.008131},
	{0.275329,0.080384,0.035952,0.016171},
	{0.282596,0.108994,0.045510,0.014327},
	{0.273227,0.123662,0.063561,0.010354},
	{0.286850,0.088067,0.054854,0.007260},
	{0.305680,0.088521,0.040830,0.008724},
	{0.323658,0.122036,0.045958,0.009372},
	{0.305577,0.108660,0.055999,0.010069},
	{0.326432,0.092403,0.071328,0.010172},
	{0.332822,0.094105,0.046980,0.010736},
	{0.360057,0.114979,0.050602,0.008675},
	{0.345196,0.113199,0.049286,0.006350},
	{0.359002,0.099534,0.049942,0.006821},
	{0.379041,0.093289,0.031702,0.012829},
	{0.380592,0.110656,0.050707,0.012044},
	{0.372838,0.103784,0.069441,0.009079},
	{0.398837,0.089939,0.069654,0.008117},
	{0.404951,0.099259,0.042662,0.015929},
	{0.430865,0.122250,0.034910,0.014087},
	{0.408865,0.117363,0.058603,0.007815},
	{0.433797,0.100725,0.068087,0.010539},
	{0.451913,0.097993,0.031985,0.009479},
	{0.468554,0.117154,0.051922,0.008497},
	{0.450291,0.108384,0.052179,0.011764},
	{0.468896,0.094189,0.071862,0.008497},
	{0.017769,0.123023,0.046288,0.012482},
	{0.039600,0.152973,0.031490,0.013071},
	{0.022486,0.159910,0.058277,0.006730},
	{0.053035,0.135514,0.059822,0.006873},
	{0.060820,0.135597,0.032391,0.009756},
	{0.075793,0.140311,0.044451,0.007184},
	{0.059415,0.151698,0.048706,0.013209},
	{0.075397,0.146441,0.067826,0.009066},
	{0.096374,0.134144,0.029818,0.009090},
	{0.114668,0.150443,0.049586,0.009896},
	{0.098456,0.146495,0.056580,0.008196},
	{0.115068,0.135950,0.063642,0.006871},
	{0.120480,0.129505,0.039205,0.007619},
	{0.156971,0.148691,0.034191,0.009631},
	{0.137556,0.146645,0.052913,0.006461},
	{0.141950,0.131573,0.059911,0.010688},
	{0.168096,0.130776,0.044977,0.013464},
	{0.181562,0.156313,0.051610,0.015347},
	{0.157127,0.150002,0.051947,0.008174},
	{0.193242,0.120660,0.050909,0.012412},
	{0.201771,0.131287,0.031430,0.011359},
	{0.221674,0.149773,0.053281,0.006817},
	{0.205379,0.141341,0.053543,0.011532},
	{0.224583,0.129576,0.068306,0.012332},
	{0.238902,0.136528,0.043920,0.011950},
	{0.258910,0.151108,0.049705,0.013475},
	{0.232642,0.144977,0.060278,0.006814},
	{0.257772,0.129747,0.049843,0.007917},
	{0.272859,0.130994,0.031434,0.013547},
	{0.298637,0.152991,0.041282,0.015762},
	{0.278487,0.158034,0.062419,0.008408},
	{0.294879,0.136446,0.061848,0.007854},
	{0.305700,0.128154,0.049784,0.009982},
	{0.317822,0.154906,0.054591,0.007666},
	{0.313082,0.139892,0.064913,0.007314},
	{0.318412,0.124214,0.060868,0.006584},
	{0.345403,0.128047,0.055609,0.007438},
	{0.364150,0.139196,0.044629,0.013056},
	{0.335338,0.144636,0.042778,0.015825},
	{0.370092,0.124629,0.059323,0.007443},
	{0.381964,0.126156,0.033148,0.011118},
	{0.405718,0.145022,0.044830,0.007344},
	{0.384403,0.150897,0.061572,0.015826},
	{0.394762,0.130323,0.049303,0.009909},
	{0.418237,0.126941,0.051554,0.007327},
	{0.427303,0.148241,0.032287,0.012279},
	{0.408184,0.153492,0.057158,0.007815},
	{0.434164,0.125257,0.064771,0.012118},
	{0.445472,0.139238,0.049416,0.009855},
	{0.469080,0.147719,0.054616,0.015762},
	{0.446519,0.149106,0.070188,0.007289},
	{0.463066,0.128741,0.072069,0.010712},
	{0.030541,0.171854,0.034568,0.007470},
	{0.036934,0.187996,0.038951,0.010436},
	{0.017296,0.187394,0.069932,0.015775},
	{0.045470,0.171678,0.064724,0.015762},
	{0.056514,0.169207,0.032678,0.010404},
	{0.079003,0.186246,0.029530,0.016079},
	{0.059706,0.192913,0.050156,0.012943},
	{0.081987,0.176494,0.057278,0.010079},
	{0.098980,0.166199,0.053356,0.010173},
	{0.117318,0.178623,0.040798,0.015290},
	{0.102820,0.182572,0.068939,0.007422},
	{0.119261,0.172887,0.062372,0.007118},
	{0.134901,0.165613,0.030499,0.008886},
	{0.147868,0.177545,0.033698,0.009023},
	{0.136653,0.175577,0.056404,0.008907},
	{0.141868,0.158596,0.054271,0.006316},
	{0.165708,0.172844,0.039254,0.010243},
	{0.182639,0.183677,0.050062,0.011541},
	{0.165084,0.185889,0.065932,0.012227},
	{0.180948,0.173447,0.066522,0.007376},
	{0.196132,0.172251,0.035240,0.011530},
	{0.227201,0.182137,0.045706,0.008997},
	{0.204299,0.184844,0.052950,0.010342},
	{0.225988,0.175305,0.062124,0.008828},
	{0.237177,0.160547,0.051402,0.009195},
	{0.261161,0.186504,0.033019,0.015762},
	{0.241624,0.187436,0.069292,0.007769},
	{0.257173,0.165912,0.069549,0.008701},
	{0.275377,0.158582,0.034952,0.008984},
	{0.288726,0.191148,0.034408,0.009808},
	{0.271403,0.185700,0.067431,0.011799},
	{0.284218,0.172191,0.051535,0.010347},
	{0.306049,0.174027,0.032682,0.007870},
	{0.316942,0.179608,0.046528,0.010610},
	{0.313026,0.194025,0.059273,0.009027},
	{0.316419,0.168708,0.067588,0.006369},
	{0.344253,0.166378,0.030810,0.010546},
	{0.367887,0.186357,0.033459,0.007112},
	{0.358319,0.179340,0.050728,0.013841},
	{0.366596,0.160155,0.045909,0.007602},
	{0.383817,0.160452,0.034727,0.012399},
	{0.396469,0.185253,0.029237,0.015979},
	{0.377815,0.185629,0.073245,0.013364},
	{0.396205,0.166460,0.049105,0.007346},
	{0.413309,0.167440,0.053605,0.007275},
	{0.434009,0.184506,0.031699,0.013620},
	{0.411121,0.187057,0.053929,0.011981},
	{0.435082,0.168647,0.049551,0.010284},
	{0.451053,0.158269,0.035414,0.012610},
	{0.470872,0.191256,0.051955,0.013847},
	{0.447897,0.178522,0.059981,0.008314},
	{0.465203,0.170594,0.069944,0.011281},
	{0.018220,0.199513,0.032480,0.012472},
	{0.043929,0.213054,0.038912,0.014059},
	{0.022453,0.225231,0.059937,0.012276},
	{0.039438,0.194934,0.058858,0.008970},
	{0.062291,0.203776,0.030515,0.008161},
	{0.067233,0.212277,0.050888,0.007657},
	{0.063555,0.218942,0.074236,0.011822},
	{0.086676,0.195765,0.056391,0.009775},
	{0.100062,0.203963,0.032083,0.010760},
	{0.123439,0.223656,0.039869,0.007878},
	{0.087363,0.218666,0.047798,0.013687},
	{0.109855,0.201133,0.048087,0.008215},
	{0.129021,0.205306,0.043706,0.011683},
	{0.157516,0.221570,0.036912,0.010458},
	{0.126303,0.227505,0.059162,0.010490},
	{0.143861,0.194185,0.066210,0.010561},
	{0.172910,0.201043,0.033978,0.013203},
	{0.185035,0.221440,0.045513,0.013181},
	{0.168070,0.220355,0.072675,0.009360},
	{0.180093,0.198384,0.061991,0.007046},
	{0.208314,0.209539,0.044961,0.009495},
	{0.220276,0.226073,0.034315,0.013522},
	{0.202082,0.218715,0.068875,0.006812},
	{0.227664,0.201788,0.056524,0.012609},
	{0.232983,0.200019,0.028712,0.015762},
	{0.256565,0.219723,0.031699,0.007469},
	{0.227652,0.224321,0.057241,0.006644},
	{0.254299,0.203686,0.071224,0.012930},
	{0.262164,0.203438,0.047374,0.006460},
	{0.291834,0.224737,0.049349,0.006948},
	{0.268860,0.219381,0.069067,0.007715},
	{0.289943,0.200233,0.060624,0.012722},
	{0.301505,0.208536,0.039247,0.012306},
	{0.322382,0.217100,0.051740,0.012086},
	{0.314574,0.219134,0.070936,0.007251},
	{0.325878,0.206084,0.069072,0.007512},
	{0.344014,0.201827,0.047775,0.012892},
	{0.351810,0.217320,0.036368,0.007867},
	{0.338587,0.227804,0.063659,0.010701},
	{0.356571,0.197024,0.069754,0.008434},
	{0.376453,0.200295,0.047115,0.014105},
	{0.394956,0.226437,0.029396,0.010725},
	{0.378323,0.226833,0.061314,0.014536},
	{0.389287,0.205171,0.066511,0.007194},
	{0.416665,0.206194,0.024533,0.013492},
	{0.438684,0.208898,0.032405,0.009490},
	{0.404012,0.213687,0.050773,0.015762},
	{0.439529,0.206631,0.054037,0.011117},
	{0.452197,0.194638,0.042409,0.007397},
	{0.478578,0.218454,0.044278,0.008032},
	{0.461476,0.216922,0.055979,0.012027},
	{0.472499,0.206672,0.070146,0.008643},
	{0.031034,0.238557,0.044665,0.009415},
	{0.046212,0.247434,0.048494,0.008580},
	{0.029656,0.251629,0.061475,0.011572},
	{0.046678,0.229877,0.059527,0.012161},
	{0.060543,0.228584,0.035366,0.008957},
	{0.081148,0.255496,0.031889,0.011593},
	{0.066991,0.258165,0.055249,0.015762},
	{0.078630,0.234901,0.067187,0.009147},
	{0.093096,0.239315,0.031956,0.008521},
	{0.123349,0.260861,0.042759,0.012372},
	{0.097630,0.259258,0.071577,0.008194},
	{0.104871,0.242539,0.055778,0.015789},
	{0.130936,0.239003,0.044952,0.008368},
	{0.153591,0.241938,0.041655,0.008205},
	{0.143130,0.253503,0.053056,0.011112},
	{0.158308,0.235817,0.066131,0.006732},
	{0.169649,0.237061,0.031901,0.009847},
	{0.177988,0.267639,0.051141,0.008419},
	{0.168590,0.251903,0.050804,0.009913},
	{0.187601,0.240303,0.056261,0.008679},
	{0.202598,0.238631,0.030252,0.008355},
	{0.221493,0.253090,0.049019,0.008616},
	{0.198234,0.260140,0.040252,0.015762},
	{0.221897,0.228125,0.071226,0.008950},
	{0.236460,0.237179,0.031506,0.006307},
	{0.246367,0.254381,0.050248,0.011434},
	{0.230061,0.245247,0.060087,0.007429},
	{0.243784,0.239481,0.061145,0.007205},
	{0.270130,0.238534,0.030357,0.015762},
	{0.292685,0.259379,0.051551,0.008561},
	{0.270568,0.256201,0.053461,0.012554},
	{0.288088,0.246768,0.065776,0.009704},
	{0.307205,0.228427,0.050927,0.006869},
	{0.321366,0.246094,0.045009,0.011972},
	{0.314218,0.248123,0.067437,0.006555},
	{0.317918,0.234420,0.056996,0.005112},
	{0.345996,0.237475,0.032189,0.012907},
	{0.358078,0.266408,0.046713,0.010906},
	{0.341845,0.249863,0.050455,0.009550},
	{0.352246,0.244796,0.066207,0.009874},
	{0.370773,0.231259,0.033624,0.012678},
	{0.389736,0.252726,0.029076,0.015191},
	{0.374368,0.254883,0.060882,0.013568},
	{0.399580,0.236433,0.066670,0.009396},
	{0.419086,0.229820,0.030486,0.008709},
	{0.431898,0.253013,0.026908,0.014166},
	{0.424434,0.259566,0.058997,0.006353},
	{0.425774,0.238512,0.052362,0.015762},
	{0.449984,0.229928,0.034152,0.014447},
	{0.468335,0.247681,0.035305,0.009986},
	{0.446284,0.265609,0.058199,0.014044},
	{0.479957,0.238820,0.056505,0.015762},
	{0.016058,0.272296,0.031907,0.015903},
	{0.043188,0.296346,0.035622,0.013645},
	{0.020265,0.293955,0.054775,0.015873},
	{0.041510,0.268037,0.070523,0.007849},
	{0.062225,0.280373,0.040949,0.011078},
	{0.086682,0.281592,0.035069,0.013702},
	{0.069169,0.284177,0.059909,0.009020},
	{0.085718,0.266880,0.066273,0.006910},
	{0.097403,0.266171,0.048894,0.009477},
	{0.107654,0.281428,0.038660,0.007576},
	{0.088316,0.284050,0.066490,0.010354},
	{0.113503,0.280552,0.060462,0.015014},
	{0.122859,0.278804,0.036141,0.006759},
	{0.142878,0.297979,0.037101,0.011706},
	{0.133777,0.283064,0.049774,0.007955},
	{0.154649,0.263474,0.065692,0.007507},
	{0.169175,0.273460,0.031813,0.013607},
	{0.179325,0.294941,0.051789,0.010091},
	{0.165815,0.277188,0.064486,0.010220},
	{0.192842,0.268948,0.065506,0.008519},
	{0.203664,0.282860,0.030067,0.008961},
	{0.221520,0.290705,0.031569,0.010600},
	{0.204075,0.289620,0.057034,0.015762},
	{0.217062,0.273218,0.067345,0.007561},
	{0.236623,0.274054,0.029708,0.008988},
	{0.259386,0.293049,0.032016,0.011916},
	{0.240332,0.286268,0.050913,0.015762},
	{0.252585,0.269128,0.068044,0.010690},
	{0.277228,0.263086,0.036137,0.007242},
	{0.297808,0.298823,0.038507,0.011507},
	{0.273109,0.288314,0.070360,0.007871},
	{0.283880,0.280089,0.053191,0.014002},
	{0.311353,0.272674,0.042400,0.016116},
	{0.318185,0.295940,0.050900,0.009579},
	{0.300031,0.294218,0.056245,0.006682},
	{0.324041,0.283074,0.063804,0.007747},
	{0.338338,0.274337,0.045421,0.010405},
	{0.366658,0.294673,0.036475,0.013011},
	{0.345468,0.300165,0.061263,0.006730},
	{0.371053,0.276011,0.062724,0.007746},
	{0.382159,0.266197,0.045056,0.007040},
	{0.385302,0.283479,0.044022,0.010008},
	{0.382756,0.291084,0.070614,0.012573},
	{0.403108,0.264928,0.054184,0.015762},
	{0.419100,0.274729,0.031563,0.011467},
	{0.437640,0.289426,0.049922,0.012611},
	{0.416213,0.289680,0.065434,0.013463},
	{0.425331,0.272106,0.052713,0.007637},
	{0.454008,0.282624,0.038875,0.008275},
	{0.483070,0.298139,0.038400,0.016342},
	{0.446217,0.292641,0.069129,0.007245},
	{0.473937,0.269688,0.049149,0.015337},
	{0.027818,0.310268,0.030259,0.007775},
	{0.037764,0.320592,0.051731,0.008781},
	{0.016776,0.324398,0.064621,0.015762},
	{0.043065,0.304671,0.073634,0.008404},
	{0.060248,0.307914,0.031150,0.007448},
	{0.077271,0.327681,0.031663,0.013927},
	{0.055259,0.322950,0.050938,0.008327},
	{0.079365,0.305257,0.051792,0.015762},
	{0.095660,0.308012,0.031429,0.010463},
	{0.103681,0.332264,0.031625,0.012290},
	{0.089797,0.330875,0.053651,0.011579},
	{0.102216,0.307408,0.054272,0.006891},
	{0.135063,0.305375,0.051629,0.006373},
	{0.146687,0.332542,0.040487,0.015166},
	{0.119733,0.321591,0.052005,0.015762},
	{0.139154,0.311001,0.062753,0.006747},
	{0.166611,0.304212,0.025702,0.013467},
	{0.179450,0.331099,0.030605,0.013540},
	{0.168140,0.316792,0.050046,0.012886},
	{0.181651,0.308754,0.061301,0.006449},
	{0.196784,0.307080,0.030792,0.016081},
	{0.215900,0.329030,0.028037,0.012976},
	{0.194417,0.327273,0.054083,0.014566},
	{0.222938,0.309720,0.049824,0.007765},
	{0.236225,0.315755,0.036277,0.011069},
	{0.259154,0.323023,0.045475,0.012979},
	{0.229226,0.318793,0.058722,0.006414},
	{0.263876,0.304803,0.055229,0.008220},
	{0.280036,0.308037,0.046953,0.010221},
	{0.289039,0.323201,0.031727,0.013078},
	{0.277826,0.320270,0.068751,0.009886},
	{0.281755,0.301785,0.064415,0.008406},
	{0.311687,0.313114,0.034665,0.007177},
	{0.321013,0.324686,0.029636,0.008514},
	{0.307138,0.327658,0.072817,0.007383},
	{0.319851,0.313892,0.064538,0.009449},
	{0.333424,0.313519,0.041662,0.007831},
	{0.360179,0.317042,0.032555,0.009090},
	{0.345335,0.325054,0.051853,0.011632},
	{0.356818,0.307255,0.054092,0.008452},
	{0.376257,0.314516,0.044839,0.010566},
	{0.404942,0.322076,0.039531,0.009398},
	{0.377646,0.330732,0.066428,0.015762},
	{0.393248,0.303006,0.052180,0.011230},
	{0.414496,0.303557,0.046583,0.009348},
	{0.429851,0.323348,0.043622,0.015876},
	{0.411263,0.319949,0.059797,0.008998},
	{0.432707,0.308279,0.070597,0.011926},
	{0.453339,0.305155,0.030681,0.015166},
	{0.467869,0.322571,0.052192,0.010042},
	{0.449922,0.317745,0.059766,0.010027},
	{0.461253,0.304033,0.060778,0.007790},
	{0.021541,0.332873,0.040513,0.010232},
	{0.038000,0.356554,0.031870,0.007826},
	{0.032187,0.360165,0.053724,0.015074},
	{0.051751,0.340670,0.059341,0.011595},
	{0.050035,0.338842,0.039691,0.008214},
	{0.074340,0.355322,0.033149,0.012791},
	{0.064951,0.354774,0.051051,0.007431},
	{0.073704,0.344766,0.063534,0.010807},
	{0.090705,0.345601,0.029750,0.006412},
	{0.113418,0.360863,0.031428,0.008821},
	{0.096799,0.362099,0.049082,0.008369},
	{0.107463,0.341185,0.072071,0.014400},
	{0.121217,0.345472,0.047617,0.008563},
	{0.151164,0.358997,0.045262,0.011416},
	{0.135273,0.360729,0.057086,0.008467},
	{0.150044,0.351765,0.063675,0.007622},
	{0.167267,0.337703,0.045162,0.006559},
	{0.185291,0.357837,0.029177,0.013202},
	{0.174685,0.353137,0.051201,0.010915},
	{0.182603,0.341175,0.071944,0.009153},
	{0.199074,0.341536,0.033569,0.008593},
	{0.223648,0.364790,0.032933,0.013250},
	{0.206142,0.359909,0.070232,0.011880},
	{0.218376,0.338725,0.063741,0.006775},
	{0.235978,0.341262,0.044598,0.015762},
	{0.248271,0.357981,0.028623,0.010326},
	{0.239381,0.361261,0.054185,0.006676},
	{0.262678,0.342288,0.073448,0.011112},
	{0.270573,0.339357,0.051890,0.007719},
	{0.291651,0.365260,0.034029,0.012619},
	{0.278330,0.357467,0.061963,0.013218},
	{0.287271,0.340846,0.047013,0.009740},
	{0.311721,0.350588,0.047719,0.015762},
	{0.329835,0.361562,0.031124,0.010858},
	{0.308858,0.362142,0.074170,0.007884},
	{0.318657,0.331854,0.062466,0.008662},
	{0.334080,0.339898,0.036091,0.011614},
	{0.361053,0.363146,0.028462,0.015405},
	{0.342858,0.361137,0.049372,0.011564},
	{0.357992,0.345258,0.056926,0.010458},
	{0.376633,0.335585,0.028650,0.016006},
	{0.388485,0.356145,0.042465,0.011050},
	{0.367467,0.367712,0.056170,0.013400},
	{0.398534,0.339188,0.050584,0.007905},
	{0.410905,0.345867,0.031138,0.016090},
	{0.433557,0.364385,0.032172,0.009789},
	{0.414823,0.362718,0.050261,0.009697},
	{0.425989,0.342978,0.064208,0.012469},
	{0.446622,0.347487,0.046030,0.008548},
	{0.463905,0.355130,0.053326,0.007996},
	{0.447393,0.368651,0.049700,0.012945},
	{0.474977,0.341544,0.046688,0.010746},
	{0.019773,0.386775,0.049869,0.014542},
	{0.049900,0.390242,0.042917,0.008807},
	{0.022183,0.402617,0.069178,0.007405},
	{0.034110,0.379357,0.068905,0.008232},
	{0.067211,0.383096,0.027225,0.015626},
	{0.078560,0.396236,0.052905,0.006480},
	{0.068770,0.402250,0.062282,0.008351},
	{0.075963,0.374658,0.055577,0.014744},
	{0.104738,0.376287,0.048997,0.007889},
	{0.119282,0.396959,0.048600,0.010216},
	{0.096405,0.402327,0.056212,0.008008},
	{0.107422,0.388063,0.061078,0.009161},
	{0.131436,0.376425,0.038024,0.015884},
	{0.148514,0.406019,0.048436,0.013951},
	{0.132288,0.405409,0.061120,0.006654},
	{0.151652,0.388678,0.066685,0.006268},
	{0.169202,0.379342,0.050424,0.015870},
	{0.188174,0.396136,0.022829,0.007254},
	{0.170203,0.400366,0.065913,0.009355},
	{0.194966,0.372889,0.055640,0.010498},
	{0.212711,0.382905,0.044145,0.006430},
	{0.220250,0.394315,0.029698,0.010635},
	{0.195739,0.397199,0.044541,0.015762},
	{0.211654,0.377701,0.056428,0.006888},
	{0.234965,0.382146,0.029519,0.007750},
	{0.263901,0.397112,0.046806,0.006980},
	{0.234022,0.391910,0.051822,0.015536},
	{0.260195,0.374135,0.045281,0.015762},
	{0.280399,0.376799,0.056025,0.006714},
	{0.283192,0.402145,0.026391,0.009786},
	{0.279561,0.389831,0.070096,0.006789},
	{0.299585,0.373604,0.057431,0.011953},
	{0.312740,0.372343,0.033714,0.009518},
	{0.335479,0.387674,0.036866,0.008966},
	{0.313317,0.390191,0.058539,0.009610},
	{0.327009,0.381826,0.052052,0.007696},
	{0.349111,0.381525,0.050691,0.009803},
	{0.370847,0.389607,0.040148,0.007174},
	{0.336900,0.404187,0.051731,0.013298},
	{0.364431,0.386995,0.051672,0.006272},
	{0.380486,0.372616,0.037016,0.007843},
	{0.387432,0.394596,0.035861,0.010668},
	{0.383440,0.392839,0.070329,0.010113},
	{0.397354,0.378289,0.069987,0.008162},
	{0.412381,0.379269,0.031555,0.014769},
	{0.431738,0.399853,0.031097,0.013491},
	{0.416807,0.388592,0.068089,0.010234},
	{0.434246,0.385754,0.051562,0.008510},
	{0.459761,0.377812,0.034232,0.008546},
	{0.469310,0.399879,0.030659,0.015762},
	{0.448851,0.389468,0.070969,0.016061},
	{0.472731,0.373287,0.052274,0.012192},
	{0.025845,0.411569,0.031706,0.012222},
	{0.053452,0.432332,0.034575,0.012284},
	{0.013101,0.425034,0.052342,0.015517},
	{0.035661,0.409240,0.053396,0.011700},
	{0.069697,0.410872,0.045637,0.008799},
	{0.085224,0.435507,0.032725,0.015801},
	{0.056265,0.422981,0.058037,0.013129},
	{0.076985,0.417420,0.069616,0.007327},
	{0.086433,0.412353,0.047787,0.008139},
	{0.113150,0.426126,0.051781,0.007572},
	{0.095559,0.431375,0.071602,0.012141},
	{0.110469,0.407765,0.062669,0.008395},
	{0.131849,0.404271,0.035325,0.007326},
	{0.154057,0.433609,0.051576,0.006804},
	{0.130968,0.431234,0.049111,0.011155},
	{0.146668,0.411976,0.072037,0.006514},
	{0.173190,0.416526,0.046090,0.012972},
	{0.190589,0.427391,0.031914,0.009879},
	{0.165638,0.431320,0.062630,0.009368},
	{0.181578,0.418447,0.070458,0.012485},
	{0.201938,0.420097,0.050532,0.007909},
	{0.210786,0.430524,0.043687,0.007384},
	{0.199643,0.440108,0.058556,0.013086},
	{0.219101,0.413160,0.043518,0.011656},
	{0.235474,0.419530,0.052019,0.007861},
	{0.254554,0.428819,0.029869,0.009518},
	{0.251638,0.431825,0.063212,0.015327},
	{0.253052,0.411818,0.042596,0.011772},
	{0.268275,0.411836,0.031385,0.007133},
	{0.299229,0.425848,0.044101,0.011809},
	{0.277241,0.422348,0.053397,0.012319},
	{0.285524,0.396817,0.050771,0.014608},
	{0.311540,0.402649,0.035156,0.015936},
	{0.331451,0.434548,0.035449,0.006888},
	{0.310673,0.434227,0.066608,0.012917},
	{0.322590,0.420054,0.049005,0.008242},
	{0.340380,0.414409,0.031798,0.009373},
	{0.362861,0.429363,0.031830,0.012967},
	{0.340153,0.437578,0.060896,0.007843},
	{0.370512,0.415453,0.055970,0.010546},
	{-1.,-1.,-1.,-1. } /* sentinel: non-positive radius */
};

Vector3r Shop::inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2)
{
	return v0+((v2-v0)*(v1-v0).norm()+(v1-v0)*(v2-v0).norm())/((v1-v0).norm()+(v2-v1).norm()+(v0-v2).norm());
}

void Shop::getViscoelasticFromSpheresInteraction( Real m, Real tc, Real en, Real es, shared_ptr<ViscElMat> b)
{
    b->kn = m/tc/tc * ( Mathr::PI*Mathr::PI + Mathr::Pow(Mathr::Log(en),2) );
    b->cn = -2.0*m/tc * Mathr::Log(en);
    b->ks = 2.0/7.0 * m/tc/tc * ( Mathr::PI*Mathr::PI + Mathr::Pow(Mathr::Log(es),2) );
    b->cs = -2.0/7.0 * m/tc * Mathr::Log(es);

    if (Mathr::FAbs(b->cn) <= Mathr::ZERO_TOLERANCE ) b->cn=0;
    if (Mathr::FAbs(b->cs) <= Mathr::ZERO_TOLERANCE ) b->cs=0;
}

/* This function is copied almost verbatim from scientific python, module Visualization, class ColorScale
 *
 */
Vector3r Shop::scalarOnColorScale(Real x, Real xmin, Real xmax){
	Real xnorm=min((Real)1.,max((x-xmin)/(xmax-xmin),(Real)0.));
	if(xnorm<.25) return Vector3r(0,.4*xnorm,1);
	if(xnorm<.5)  return Vector3r(0,1,1.-4.*(xnorm-.25));
	if(xnorm<.75) return Vector3r(4*(xnorm-.5),1.,0);
	return Vector3r(1,1-4*(xnorm-.75),0);
}

/* Wrap floating point number into interval (x0,x1〉such that it is shifted
 * by integral number of the interval range. If given, *period will hold
 * this number. The wrapped value is returned.
 */
Real Shop::periodicWrap(Real x, Real x0, Real x1, long* period){
	Real xNorm=(x-x0)/(x1-x0);
	Real xxNorm=xNorm-floor(xNorm);
	if(period) *period=(long)floor(xNorm);
	return x0+xxNorm*(x1-x0);
}

void Shop::getStressForEachBody(vector<Shop::bodyState>& bodyStates){
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	bodyStates.resize(scene->bodies->size());
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		const NormShearPhys* phys = YADE_CAST<NormShearPhys*>(I->interactionPhysics.get());
		Dem3DofGeom* geom=YADE_CAST<Dem3DofGeom*>(I->interactionGeometry.get());	//For the moment only for Dem3DofGeom!!!
		if(!phys) continue;
		const body_id_t id1=I->getId1(), id2=I->getId2();
		
		Real minRad=(geom->refR1<=0?geom->refR2:(geom->refR2<=0?geom->refR1:min(geom->refR1,geom->refR2)));
		Real crossSection=Mathr::PI*pow(minRad,2);
		
		Vector3r normalStress=((1./crossSection)*geom->normal.dot(phys->normalForce))*geom->normal;
		Vector3r shearStress;
		for(int i=0; i<3; i++){
			int ix1=(i+1)%3,ix2=(i+2)%3;
			shearStress[i]=geom->normal[ix1]*phys->shearForce[ix1]+geom->normal[ix2]*phys->shearForce[ix2];
			shearStress[i]/=crossSection;
		}
		
		bodyStates[id1].normStress+=normalStress;
		bodyStates[id2].normStress+=normalStress;
		
		bodyStates[id1].shearStress+=shearStress;
		bodyStates[id2].shearStress+=shearStress;
	}
}

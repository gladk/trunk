#include<yade/pkg-dem/DomainLimiter.hpp>
#include<yade/pkg-dem/Shop.hpp>

YADE_PLUGIN((DomainLimiter)(LawTester)
	#ifdef YADE_OPENGL
		(GlExtra_LawTester)
	#endif
);

void DomainLimiter::action(){
	std::list<Body::id_t> out;
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue;
		const Vector3r& p(b->state->pos);
		if(p[0]<lo[0] || p[0]>hi[0] || p[1]<lo[1] || p[1]>hi[1] || p[2]<lo[2] || p[2]>hi[2]) out.push_back(b->id);
	}
	FOREACH(Body::id_t id, out){
		scene->bodies->erase(id);
		nDeleted++;
	}
}

#include<yade/pkg-dem/DemXDofGeom.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/L3Geom.hpp>
#include<yade/pkg-common/NormShearPhys.hpp>
#include<yade/pkg-common/LinearInterpolate.hpp>
#include<yade/lib-pyutil/gil.hpp>

CREATE_LOGGER(LawTester);

void LawTester::postLoad(LawTester&){
	if(ids.size()==0) return; // uninitialized object, don't do nothing at all
	if(ids.size()!=2) throw std::invalid_argument("LawTester.ids: exactly two values must be given.");
	if(path.empty() && rotPath.empty()) throw invalid_argument("LawTester.{path,rotPath}: at least one point must be given.");
	if(pathSteps.empty()) throw invalid_argument("LawTester.pathSteps: at least one value must be given.");
	size_t pathSize=max(path.size(),rotPath.size());
	// update path points
	_pathU.clear(); _pathU.push_back(Vector3r::Zero());
	_pathR.clear(); _pathR.push_back(Vector3r::Zero());
	for(size_t i=0; i<pathSize; i++) {
		_pathU.push_back(i<path.size()?path[i]:(path.empty()?Vector3r::Zero():*(path.rbegin())));
		_pathR.push_back(i<rotPath.size()?rotPath[i]:(rotPath.empty()?Vector3r::Zero():*(rotPath.rbegin())));
	}
	// update time points from distances, repeat last distance if shorter than path
	_pathT.clear(); _pathT.push_back(0);
	for(size_t i=0; i<pathSteps.size(); i++) _pathT.push_back(_pathT[i]+pathSteps[i]);
	int lastDist=pathSteps[pathSteps.size()-1];
	for(size_t i=pathSteps.size(); i<pathSize; i++) _pathT.push_back(*(_pathT.rbegin())+lastDist);
}

void LawTester::action(){
	Vector2r foo; // avoid undefined ~Vector2r with clang?
	if(ids.size()!=2) throw std::invalid_argument("LawTester.ids: exactly two values must be given.");
	LOG_DEBUG("=================== LawTester step "<<step<<" ========================");
	const shared_ptr<Interaction> Inew=scene->interactions->find(ids[0],ids[1]);
	string strIds("##"+lexical_cast<string>(ids[0])+"+"+lexical_cast<string>(ids[1]));
	// interaction not found at initialization
	if(!I && (!Inew || !Inew->isReal())) throw std::runtime_error("LawTester: interaction "+strIds+" does not exist"+(Inew?" (to be honest, it does exist, but it is not real).":"."));
	// interaction was deleted meanwhile
	if(I && (!Inew || !Inew->isReal())) throw std::runtime_error("LawTester: interaction "+strIds+" was deleted"+(Inew?" (is not real anymore).":"."));
	// different interaction object
	if(I && Inew && I!=Inew) throw std::logic_error("LawTester: interacion "+strIds+" is a different object now?!");
	assert(Inew);
	bool doInit=(!I);
	if(doInit) I=Inew;

	id1=I->getId1(); id2=I->getId2();
	// test object types
	GenericSpheresContact* gsc=dynamic_cast<GenericSpheresContact*>(I->geom.get());
	ScGeom* scGeom=dynamic_cast<ScGeom*>(I->geom.get());
	Dem3DofGeom* d3dGeom=dynamic_cast<Dem3DofGeom*>(I->geom.get());
	L3Geom* l3Geom=dynamic_cast<L3Geom*>(I->geom.get());
	L6Geom* l6Geom=dynamic_cast<L6Geom*>(I->geom.get());
	//NormShearPhys* phys=dynamic_cast<NormShearPhys*>(I->phys.get());			//Disabled because of warning
	if(!gsc) throw std::invalid_argument("LawTester: IGeom of "+strIds+" not a GenericSpheresContact.");
	if(!scGeom && !d3dGeom && !l3Geom) throw std::invalid_argument("LawTester: IGeom of "+strIds+" is neither ScGeom, nor Dem3DofGeom, nor L3Geom (or L6Geom).");
	assert(!((bool)scGeom && (bool)d3dGeom && (bool)l3Geom)); // nonsense
	// get body objects
	State *state1=Body::byId(id1,scene)->state.get(), *state2=Body::byId(id2,scene)->state.get();
	scene->forces.sync();
	//if(state1->blockedDOFs!=State::DOF_RXRYRZ) { LOG_INFO("Blocking rotational DOFs for #"<<id1); state1->blockedDOFs=State::DOF_RXRYRZ;}
	//if(state2->blockedDOFs!=State::DOF_RXRYRZ) { LOG_INFO("Blocking rotational DOFs for #"<<id2); state2->blockedDOFs=State::DOF_RXRYRZ;}
	if(state1->blockedDOFs!=State::DOF_ALL) { LOG_INFO("Blocking all DOFs for #"<<id1); state1->blockedDOFs=State::DOF_ALL;}
	if(state2->blockedDOFs!=State::DOF_ALL) { LOG_INFO("Blocking all DOFs for #"<<id2); state2->blockedDOFs=State::DOF_ALL;}


	if(step>*(_pathT.rbegin())){
		LOG_INFO("Last step done, setting zero velocities on #"<<id1<<", #"<<id2<<".");
		state1->vel=state1->angVel=state2->vel=state2->angVel=Vector3r::Zero();
		if(doneHook.empty()){ LOG_INFO("No doneHook set, dying."); dead=true; }
		else{ LOG_INFO("Running doneHook: "<<doneHook);	pyRunString(doneHook);}
		return;
	}
	/* initialize or update local axes and trsf */
	rotGeom=Vector3r(NaN,NaN,NaN); // this is meaningful only with l6geom
	if(!l3Geom){
		axX=gsc->normal; /* just in case */ axX.normalize();
		if(doInit){ // initialization of the new interaction
			// take vector in the y or z direction, depending on its length; arbitrary, but one of them is sure to be non-zero
			axY=axX.cross(axX[1]>axX[2]?Vector3r::UnitY():Vector3r::UnitZ());
			axY.normalize();
			axZ=axX.cross(axY);
			LOG_DEBUG("Initial axes x="<<axX<<", y="<<axY<<", z="<<axZ);
		} else { // udpate of an existing interaction
			if(scGeom){
				scGeom->rotate(axY); scGeom->rotate(axZ);
				scGeom->rotate(shearTot);
				shearTot+=scGeom->shearIncrement();
				ptGeom=Vector3r(-scGeom->penetrationDepth,shearTot.dot(axY),shearTot.dot(axZ));
			}
			else{ // d3dGeom
				throw runtime_error("LawTester: Dem3DofGeom not yet supported.");
				// essentially copies code from ScGeom, which is not very nice indeed; oh well…
				Vector3r vRel=(state2->vel+state2->angVel.cross(-gsc->refR2*gsc->normal))-(state1->vel+state1->angVel.cross(gsc->refR1*gsc->normal));
			}
		}
		// update the transformation
		// the matrix is orthonormal, since axX, axY are normalized and and axZ is their corss-product
		trsf.col(0)=axX; trsf.col(1)=axY; trsf.col(2)=axZ;
	} else {
		trsf=l3Geom->trsf;
		axX=trsf.col(0); axY=trsf.col(1); axZ=trsf.col(2);
		ptGeom=l3Geom->u;
		if(l6Geom) rotGeom=l6Geom->phi;
	}
	trsfQ=Quaternionr(trsf);
	contPt=gsc->contactPoint;
	refLength=gsc->refR1+gsc->refR2;
	renderLength=.5*refLength;
	
	// here we go ahead, finally
	Vector3r u=linearInterpolate<Vector3r,int>(step,_pathT,_pathU,_interpPos);
	Vector3r phi=linearInterpolate<Vector3r,int>(step,_pathT,_pathR,_interpPosRot);
	Vector3r dU=u-uPrev; uPrev=u;
	Vector3r dPhi=phi-phiPrev; phiPrev=phi;
	if(displIsRel){
		LOG_DEBUG("Relative displacement diff is "<<dU<<" (will be normalized by "<<gsc->refR1+gsc->refR2<<")");
		dU*=refLength;
	}
	LOG_DEBUG("Absolute diff is: displacement "<<dU<<", rotation "<<dPhi);
	ptOurs+=dU;
	rotOurs+=dPhi;

	// reset velocities where displacement is controlled
	//for(int i=0; i<3; i++){ if(forceControl[i]==0){ state1.vel[i]=0; state2.vel[i]=0; }
	
	// shear is applied as rotation of id2: dε=r₁dθ → dθ=dε/r₁;
	Vector3r vel[2],angVel[2];
	//State* states[]={state1,state2};
	for(int i=0; i<2; i++){
		int sign=(i==0?-1:1);
		Real weight=(i==0?1-idWeight:idWeight);
		Real radius=(i==0?gsc->refR1:gsc->refR2);
		// signed and weighted displacement/rotation to be applied on this sphere (reversed for #0)
		// some rotations must cancel the sign, by multiplying by sign again
		Vector3r ddU=sign*dU*weight;
		Vector3r ddPhi=sign*dPhi*weight;
		vel[i]=angVel[i]=Vector3r::Zero();

		// normal displacement

		vel[i]+=axX*ddU[0]/scene->dt;

		// shear rotation

		//   multiplication by sign cancels sign in ddU, since rotation is non-symmetric (to increase shear, both spheres have the same rotation)
		//   (unlike shear displacement, which is symmetric)
		// rotation around Z (which gives y-shear) must be inverted: +ry gives +εzm while -rz gives +εy
		Real rotZ=-sign*rotWeight*ddU[1]/radius, rotY=sign*rotWeight*ddU[2]/radius;
		angVel[i]+=(rotY*axY+rotZ*axZ)/scene->dt;

		// shear displacement

		// angle that is traversed by a sphere in order to give desired ddU when displaced on the branch of r1+r2
		// FIXME: is the branch value correct here?!
		Real arcAngleY=atan((1-rotWeight)*ddU[1]/radius), arcAngleZ=atan((1-rotWeight)*ddU[2]/radius);
		// same, but without the atan, which can be disregarded for small increments:
		//    Real arcAngleY=(1-rotWeight)*ddU[1]/radius, arcAngleZ=(1-rotWeight)*ddU[2]/radius; 
		vel[i]+=axY*radius*sin(arcAngleY)/scene->dt; vel[i]+=axZ*radius*sin(arcAngleZ)/scene->dt;

		// compensate distance increase caused by motion along the perpendicular axis
		// cos(argAngle*) is always positive, regardless of the orientation
		// and the compensation is always in the -εx sense (-sign → +1 for #0, -1 for #1)
		vel[i]+=-sign*axX*radius*((1-cos(arcAngleY))+(1-cos(arcAngleZ)))/scene->dt;

		// rotation
		angVel[i]+=trsf*ddPhi;

		LOG_DEBUG("vel="<<vel[i]<<", angVel="<<angVel[i]<<", rotY,rotZ="<<rotY<<","<<rotZ<<", arcAngle="<<arcAngleY<<","<<arcAngleZ<<", sign="<<sign<<", weight="<<weight);

	}
	state1->vel=vel[0]; state1->angVel=angVel[0];
	state2->vel=vel[1]; state2->angVel=angVel[1];
	LOG_DEBUG("Body #"<<id1<<", setting vel="<<vel[0]<<", angVel="<<angVel[0]);
	LOG_DEBUG("Body #"<<id2<<", setting vel="<<vel[1]<<", angVel="<<angVel[1]);

	/* find out where are we at in the path, run hooks if approriate */
	// _pathT has the first (zero) value added by us, so we skip it
	int nPathT=_pathT.size();
	for(int i=1; i<nPathT; i++){
		// i-th point on _pathT is (i-1)th on path; run corresponding hook, if it exists
		if(step==_pathT[i] && hooks.size()>(i-1) && !hooks[i-1].empty()) pyRunString(hooks[i-1]);
	}
	step++;
}

#ifdef YADE_OPENGL
#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/lib-opengl/GLUtils.hpp>
#include<yade/pkg-common/GLDrawFunctors.hpp>
#include<yade/pkg-common/OpenGLRenderer.hpp>
#include<GL/glu.h>

CREATE_LOGGER(GlExtra_LawTester);

void GlExtra_LawTester::render(){
	// scene object changed (after reload, for instance), for re-initialization
	if(tester && tester->scene!=scene) tester=shared_ptr<LawTester>();

	if(!tester){ FOREACH(shared_ptr<Engine> e, scene->engines){ tester=dynamic_pointer_cast<LawTester>(e); if(tester) break; } }
	if(!tester){ LOG_ERROR("No LawTester in O.engines, killing myself."); dead=true; return; }

	//if(tester->renderLength<=0) return;
	glColor3v(Vector3r(1,0,1));

	// switch to local coordinates
	glTranslatev(tester->contPt);
	#if EIGEN_WORLD_VERSION==2
		glMultMatrixd(Eigen::Transform3d(tester->trsf).data());
	#elif EIGEN_WORLD_VERSION==3
		glMultMatrixd(Eigen::Affine3d(tester->trsf).data());
	#endif


	glDisable(GL_LIGHTING); 
	//glColor3v(Vector3r(1,0,1));
	//glBegin(GL_LINES); glVertex3v(Vector3r::Zero()); glVertex3v(.1*Vector3r::Ones()); glEnd(); 
	//GLUtils::GLDrawText(string("This is the contact point!"),Vector3r::Zero(),Vector3r(1,0,1));

	// local axes
	glLineWidth(2.);
	for(int i=0; i<3; i++){
		Vector3r pt=Vector3r::Zero(); pt[i]=.5*tester->renderLength; Vector3r color=.3*Vector3r::Ones(); color[i]=1;
		GLUtils::GLDrawLine(Vector3r::Zero(),pt,color);
		GLUtils::GLDrawText(string(i==0?"x":(i==1?"y":"z")),pt,color);
	}

	// put the origin to the initial (no-shear) point, so that the current point appears at the contact point
	glTranslatev(Vector3r(0,tester->ptOurs[1],tester->ptOurs[2]));


	const int t(tester->step); const vector<int>& TT(tester->_pathT); const vector<Vector3r>& VV(tester->_pathU);
	size_t numSegments=TT.size();
	const Vector3r colorBefore=Vector3r(.7,1,.7), colorAfter=Vector3r(1,.7,.7);

	// scale displacement, if they have the strain meaning
	Real scale=1;
	if(tester->displIsRel) scale=tester->refLength;

	// find maximum displacement, draw axes in the shear plane
	Real displMax=0;
	FOREACH(const Vector3r& v, VV) displMax=max(v.squaredNorm(),displMax);
	displMax=1.2*scale*sqrt(displMax);

	glLineWidth(1.);
	GLUtils::GLDrawLine(Vector3r(0,-displMax,0),Vector3r(0,displMax,0),Vector3r(.5,0,0));
	GLUtils::GLDrawLine(Vector3r(0,0,-displMax),Vector3r(0,0,displMax),Vector3r(.5,0,0));

	// draw path
	glLineWidth(4.);
	for(size_t segment=0; segment<numSegments-1; segment++){
		// different colors before and after the current point
		Real t0=TT[segment],t1=TT[segment+1];
		const Vector3r &from=-VV[segment]*scale, &to=-VV[segment+1]*scale;
		// current segment
		if(t>t0 && t<t1){
			Real norm=(t-t0)/(t1-t0);
			GLUtils::GLDrawLine(from,from+(to-from)*norm,colorBefore);
			GLUtils::GLDrawLine(from+(to-from)*norm,to,colorAfter);
		} else {	// other segment
			GLUtils::GLDrawLine(from,to,t<t0?colorAfter:colorBefore);
		}
	}

	glLineWidth(1.);
}

#endif /* YADE_OPENGL */

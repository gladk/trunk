// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>

#include"OpenGLRenderingEngine.hpp"
#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/lib-opengl/GLUtils.hpp>
#include<yade/core/Timing.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/Aabb.hpp>

#ifdef __APPLE__
#  include <OpenGL/glu.h>
#  include <OpenGL/gl.h>
#  include <GLUT/glut.h>
#else
#  include <GL/glu.h>
#  include <GL/gl.h>
#  include <GL/glut.h>
#endif

YADE_PLUGIN((OpenGLRenderingEngine));
YADE_REQUIRE_FEATURE(OPENGL)
CREATE_LOGGER(OpenGLRenderingEngine);

bool OpenGLRenderingEngine::initDone=false;
const int OpenGLRenderingEngine::numClipPlanes;

void OpenGLRenderingEngine::init(){

	map<string,DynlibDescriptor>::const_iterator di = Omega::instance().getDynlibsDescriptor().begin();
	map<string,DynlibDescriptor>::const_iterator diEnd = Omega::instance().getDynlibsDescriptor().end();
	for(;di!=diEnd;++di){
		if (Omega::instance().isInheritingFrom((*di).first,"GlStateFunctor")) addStateFunctor((*di).first);
		if (Omega::instance().isInheritingFrom((*di).first,"GlBoundFunctor")) addBoundingVolumeFunctor((*di).first);
		if (Omega::instance().isInheritingFrom((*di).first,"GlShapeFunctor")) addInteractingGeometryFunctor((*di).first);
		if (Omega::instance().isInheritingFrom((*di).first,"GlInteractionGeometryFunctor")) addInteractionGeometryFunctor((*di).first);
		if (Omega::instance().isInheritingFrom((*di).first,"GlInteractionPhysicsFunctor")) addInteractionPhysicsFunctor((*di).first);
	}
	postProcessAttributes(true);

	clipPlaneNormals.resize(numClipPlanes);

	initgl();

	static bool glutInitDone=false;
	if(!glutInitDone){
		glutInit(&Omega::instance().origArgc,Omega::instance().origArgv);
		/* transparent spheres (still not working): glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE | GLUT_ALPHA); glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE); */
		glutInitDone=true;
	}

	initDone=true;
	// glGetError crashes at some machines?! Was never really useful, anyway.
	// reported http://www.mail-archive.com/yade-users@lists.launchpad.net/msg01482.html
	#if 0
		int e=glGetError();
		if(e!=GL_NO_ERROR) throw runtime_error((string("OpenGLRenderingEngine::init returned GL error ")+lexical_cast<string>(e)).c_str());
	#endif
}

void OpenGLRenderingEngine::setBodiesRefSe3(){
	LOG_DEBUG("(re)initializing reference positions and orientations.");
	FOREACH(const shared_ptr<Body>& b, *scene->bodies) if(b && b->state) { b->state->refPos=b->state->pos; b->state->refOri=b->state->ori; }
	scene->cell->refSize=scene->cell->getSize();
	// scene->cell->refShear=scene->cell->shear;
}


void OpenGLRenderingEngine::initgl(){
	LOG_DEBUG("(re)initializing GL for gldraw methods.\n");
	BOOST_FOREACH(vector<string>& s,stateFunctorNames) (static_pointer_cast<GlStateFunctor>(ClassFactory::instance().createShared(s[1])))->initgl();
	BOOST_FOREACH(vector<string>& s,boundFunctorNames)	(static_pointer_cast<GlBoundFunctor>(ClassFactory::instance().createShared(s[1])))->initgl();
	BOOST_FOREACH(vector<string>& s,shapeFunctorNames)	(static_pointer_cast<GlShapeFunctor>(ClassFactory::instance().createShared(s[1])))->initgl();
	BOOST_FOREACH(vector<string>& s,interactionGeometryFunctorNames) (static_pointer_cast<GlInteractionGeometryFunctor>(ClassFactory::instance().createShared(s[1])))->initgl();
	BOOST_FOREACH(vector<string>& s,interactionPhysicsFunctorNames) (static_pointer_cast<GlInteractionPhysicsFunctor>(ClassFactory::instance().createShared(s[1])))->initgl();
}

void OpenGLRenderingEngine::renderWithNames(const shared_ptr<Scene>& _scene){
	scene=_scene;
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b || !b->shape) continue;
		glPushMatrix();
		const Se3r& se3=b->state->se3;
		AngleAxisr aa(angleAxisFromQuat(se3.orientation));
		glTranslatef(se3.position[0],se3.position[1],se3.position[2]);
		glRotatef(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
		//if(b->shape->getClassName() != "LineSegment"){ // FIXME: a body needs to say: I am selectable ?!?!
			glPushName(b->getId());
			shapeDispatcher(b->shape,b->state,wire || b->shape->wire,viewInfo);
			glPopName();
		//}
		glPopMatrix();
	}
};

bool OpenGLRenderingEngine::pointClipped(const Vector3r& p){
	if(numClipPlanes<1) return false;
	for(int i=0;i<numClipPlanes;i++) if(clipPlaneActive[i]&&(p-clipPlaneSe3[i].position).dot(clipPlaneNormals[i])<0) return true;
	return false;
}


void OpenGLRenderingEngine::setBodiesDispInfo(){
	if(scene->bodies->size()!=bodyDisp.size()) bodyDisp.resize(scene->bodies->size());
	bool scaleRotations=(rotScale!=1.0);
	bool scaleDisplacements=(dispScale!=Vector3r::Ones());
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b || !b->state) continue;
		size_t id=b->getId();
		const Vector3r& pos=b->state->pos; const Vector3r& refPos=b->state->refPos;
		const Quaternionr& ori=b->state->ori; const Quaternionr& refOri=b->state->refOri;
		Vector3r cellPos=(!scene->isPeriodic ? pos : scene->cell->wrapShearedPt(pos)); // inside the cell if periodic, same as pos otherwise
		bodyDisp[id].isDisplayed=!pointClipped(cellPos);	
		// if no scaling and no periodic, return quickly
		if(!(scaleDisplacements||scaleRotations||scene->isPeriodic)){ bodyDisp[id].pos=pos; bodyDisp[id].ori=ori; continue; }
		// apply scaling
		bodyDisp[id].pos=cellPos; // point of reference (inside the cell for periodic)
		if(scaleDisplacements) bodyDisp[id].pos+=diagMult(dispScale,Vector3r(pos-refPos)); // add scaled translation to the point of reference
		if(!scaleRotations) bodyDisp[id].ori=ori;
		else{
			Quaternionr relRot=refOri.conjugate()*ori;
			AngleAxisr aa(angleAxisFromQuat(relRot));
			aa.angle()*=rotScale;
			bodyDisp[id].ori=refOri*Quaternionr(aa);
		}
	}
}

// draw periodic cell, if active
void OpenGLRenderingEngine::drawPeriodicCell(){
	if(!scene->isPeriodic) return;
	glColor3v(Vector3r(1,1,0));
	glPushMatrix();
		Vector3r size=scene->cell->getSize();
		if(dispScale!=Vector3r::Ones()) size+=diagMult(dispScale,Vector3r(size-scene->cell->refSize));
		glTranslatev(scene->cell->shearPt(.5*size)); // shear center (moves when sheared)
		glMultMatrixd(scene->cell->getGlShearTrsfMatrix());
		glScalev(size);
		glutWireCube(1);
	glPopMatrix();
}

void OpenGLRenderingEngine::resetSpecularEmission(){
	const GLfloat specular[4]={.3,.3,.3,1};
	const GLfloat emission[4]={0,0,0,0};
	glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
	glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,emission);
}

void OpenGLRenderingEngine::render(const shared_ptr<Scene>& _scene,body_id_t selection /* not sure. maybe a list of selections, or maybe bodies themselves should remember if they are selected? */) {

	if(!initDone) init();
	assert(initDone);
	current_selection = selection;

	scene=_scene;

	// assign scene inside functors

	// just to make sure, since it is not initialized by default
	if(!scene->bound) scene->bound=shared_ptr<Aabb>(new Aabb);

	// recompute emissive light colors for highlighted bodies
	Real now=TimingInfo::getNow(/*even if timing is disabled*/true)*1e-9;
	highlightEmission0[0]=highlightEmission0[1]=highlightEmission0[2]=.8*normSquare(now,1);
	highlightEmission1[0]=highlightEmission1[1]=highlightEmission0[2]=.5*normSaw(now,2);
		
	// clipping
	assert(clipPlaneNormals.size()==(size_t)numClipPlanes);
	for(size_t i=0;i<(size_t)numClipPlanes; i++){
		// someone could have modified those from python and truncate the vectors; fill those here in that case
		if(i==clipPlaneSe3.size()) clipPlaneSe3.push_back(Se3r(Vector3r::Zero(),Quaternionr::Identity()));
		if(i==clipPlaneActive.size()) clipPlaneActive.push_back(false);
		if(i==clipPlaneNormals.size()) clipPlaneNormals.push_back(Vector3r::UnitX());
		// end filling stuff modified from python
		if(clipPlaneActive[i]) clipPlaneNormals[i]=clipPlaneSe3[i].orientation*Vector3r(0,0,1);
		/* glBegin(GL_LINES);glVertex3v(clipPlaneSe3[i].position);glVertex3v(clipPlaneSe3[i].position+clipPlaneNormals[i]);glEnd(); */
	}

	// set displayed Se3 of body (scaling) and isDisplayed (clipping)
	setBodiesDispInfo();

	// set light source
	const GLfloat pos[4]	= {lightPos[0],lightPos[1],lightPos[2],1.0};
	const GLfloat ambientColor[4]={0.5,0.5,0.5,1.0};	
	//const GLfloat specularColor[4]={0.5,0.5,0.5,1.0};	
	glClearColor(bgColor[0],bgColor[1],bgColor[2],1.0);
	glLightfv(GL_LIGHT0, GL_POSITION,pos);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColor);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1); // important: do lighting calculations on both sides of polygons
	//glLightfv(GL_LIGHT0, GL_SPECULAR, specularColor);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glEnable(GL_CULL_FACE);
	// http://www.sjbaker.org/steve/omniv/opengl_lighting.html
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

	resetSpecularEmission();

	drawPeriodicCell();

	
	if (dof || id) renderDOF_ID();
	if (bound) renderBoundingVolume();
	if (shape){
		renderShape();
	}
	if (intrAllWire) renderAllInteractionsWire();
	if (intrGeom) renderInteractionGeometry();
	if (intrPhys) renderInteractionPhysics();
}

void OpenGLRenderingEngine::renderAllInteractionsWire(){
	FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
		glColor3v(i->isReal()? Vector3r(0,1,0) : Vector3r(.5,0,1));
		Vector3r p1=Body::byId(i->getId1(),scene)->state->pos;
		const Vector3r& size=scene->cell->getSize();
		Vector3r shift2(i->cellDist[0]*size[0],i->cellDist[1]*size[1],i->cellDist[2]*size[2]);
		// in sheared cell, apply shear on the mutual position as well
		shift2=scene->cell->shearPt(shift2);
		Vector3r rel=Body::byId(i->getId2(),scene)->state->pos+shift2-p1;
		if(scene->isPeriodic) p1=scene->cell->wrapShearedPt(p1);
		glBegin(GL_LINES); glVertex3v(p1);glVertex3v(Vector3r(p1+rel));glEnd();
	}
}

void OpenGLRenderingEngine::renderDOF_ID(){	
	const GLfloat ambientColorSelected[4]={10.0,0.0,0.0,1.0};	
	const GLfloat ambientColorUnselected[4]={0.5,0.5,0.5,1.0};	
	FOREACH(const shared_ptr<Body> b, *scene->bodies){
		if(!b) continue;
		if(b->shape && ((b->getGroupMask() & mask) || b->getGroupMask()==0)){
			if(!id && b->state->blockedDOFs==0) continue;
			if(current_selection==b->getId()){glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColorSelected);}
			{ // write text
				glColor3f(1.0-bgColor[0],1.0-bgColor[1],1.0-bgColor[2]);
				unsigned DOF = b->state->blockedDOFs;
				std::string sDof = std::string("") 
										+ (((DOF & State::DOF_X )!=0)?"X":" ")
										+ (((DOF & State::DOF_Y )!=0)?"Y":" ")
										+ (((DOF & State::DOF_Z )!=0)?"Z":" ")
										+ (((DOF & State::DOF_RX)!=0)?"RX":"  ")
										+ (((DOF & State::DOF_RY)!=0)?"RY":"  ")
										+ (((DOF & State::DOF_RZ)!=0)?"RZ":"  ");
				std::string sId = boost::lexical_cast<std::string>(b->getId());
				std::string str;
				if(dof && id) sId += " ";
				if(id) str += sId;
				if(dof) str += sDof;
				const Vector3r& h(current_selection==b->getId() ? highlightEmission0 : Vector3r(1,1,1));
				glColor3v(h);
				GLUtils::GLDrawText(str,bodyDisp[b->id].pos,h);
			}
			if(current_selection == b->getId()){glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColorUnselected);}
		}
	}
}

void OpenGLRenderingEngine::renderInteractionGeometry(){	
	{
		boost::mutex::scoped_lock lock(scene->interactions->drawloopmutex);
		FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
			if(!I->interactionGeometry) continue;
			const shared_ptr<Body>& b1=Body::byId(I->getId1(),scene), b2=Body::byId(I->getId2(),scene);
			if(!(bodyDisp[I->getId1()].isDisplayed||bodyDisp[I->getId2()].isDisplayed)) continue;
			glPushMatrix(); interactionGeometryDispatcher(I->interactionGeometry,I,b1,b2,intrWire); glPopMatrix();
		}
	}
}


void OpenGLRenderingEngine::renderInteractionPhysics(){	
	{
		boost::mutex::scoped_lock lock(scene->interactions->drawloopmutex);
		FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
			if(!I->interactionPhysics) continue;
			const shared_ptr<Body>& b1=Body::byId(I->getId1(),scene), b2=Body::byId(I->getId2(),scene);
			body_id_t id1=I->getId1(), id2=I->getId2();
			if(!(bodyDisp[id1].isDisplayed||bodyDisp[id2].isDisplayed)) continue;
			glPushMatrix(); interactionPhysicsDispatcher(I->interactionPhysics,I,b1,b2,intrWire); glPopMatrix();
		}
	}
}

void OpenGLRenderingEngine::renderBoundingVolume(){	
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b || !b->bound) continue;
		if(!bodyDisp[b->getId()].isDisplayed) continue;
		if(b->bound && ((b->getGroupMask()&mask) || b->getGroupMask()==0)){
			glPushMatrix(); boundDispatcher(b->bound,scene.get()); glPopMatrix();
		}
	}
	// since we remove the functor as Scene doesn't inherit from Body anymore, hardcore the rendering routine here
	// for periodic scene, renderPeriodicCell is called separately
	if(!scene->isPeriodic){
		if(!scene->bound) scene->updateBound();
		glColor3v(Vector3r(0,1,0));
		Vector3r size=scene->bound->max-scene->bound->min;
		Vector3r center=.5*(scene->bound->min+scene->bound->max);
		glPushMatrix();
			glTranslatev(center);
			glScalev(size);
			glutWireCube(1);
		glPopMatrix();
	}
}


void OpenGLRenderingEngine::renderShape()
{
	// Additional clipping planes: http://fly.srk.fer.hr/~unreal/theredbook/chapter03.html
	#if 0
		GLdouble clip0[4]={.0,1.,0.,0.}; // y<0
		glClipPlane(GL_CLIP_PLANE0,clip0);
		glEnable(GL_CLIP_PLANE0);
	#endif

	//const GLfloat ambientColorSelected[4]={10.0,0.0,0.0,1.0};	
	//const GLfloat ambientColorUnselected[4]={0.5,0.5,0.5,1.0};

	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b || !b->shape) continue;
		if(!bodyDisp[b->getId()].isDisplayed) continue;
		Vector3r pos=bodyDisp[b->getId()].pos;
		Quaternionr ori=bodyDisp[b->getId()].ori;
		if(!b->shape || !((b->getGroupMask()&mask) || b->getGroupMask()==0)) continue;
		glPushMatrix();
			AngleAxisr aa(angleAxisFromQuat(ori));	
			glTranslatef(pos[0],pos[1],pos[2]);
			glRotatef(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
			if(current_selection==b->getId() || b->shape->highlight){
				// set hightlight
				const Vector3r& h(current_selection==b->getId() ? highlightEmission0 : highlightEmission1);
				glMaterialv(GL_FRONT_AND_BACK,GL_EMISSION,h);
				glMaterialv(GL_FRONT_AND_BACK,GL_SPECULAR,h);
				//
				shapeDispatcher(b->shape,b->state,wire || b->shape->wire,viewInfo);
				// reset highlight
				resetSpecularEmission();
			} else {
				// no highlight; in case previous functor fiddled with glMaterial
				resetSpecularEmission();
				shapeDispatcher(b->shape,b->state,wire || b->shape->wire,viewInfo);
			}
		glPopMatrix();
		if(current_selection==b->getId() || b->shape->highlight){
			if(!b->bound || wire || b->shape->wire) GLUtils::GLDrawInt(b->getId(),pos);
			else {
				// move the label towards the camera by the bounding box so that it is not hidden inside the body
				const Vector3r& mn=b->bound->min; const Vector3r& mx=b->bound->max; const Vector3r& p=pos;
				Vector3r ext(viewDirection[0]>0?p[0]-mn[0]:p[0]-mx[0],viewDirection[1]>0?p[1]-mn[1]:p[1]-mx[1],viewDirection[2]>0?p[2]-mn[2]:p[2]-mx[2]); // signed extents towards the camera
				Vector3r dr=-1.01*(viewDirection.dot(ext)*viewDirection);
				GLUtils::GLDrawInt(b->getId(),pos+dr,Vector3r::Ones());
			}
		}
		// if the body goes over the cell margin, draw it in positions where the bbox overlaps with the cell in wire
		// precondition: pos is inside the cell.
		if(b->bound && scene->isPeriodic){
			const Vector3r& cellSize(scene->cell->getSize());
			pos=scene->cell->unshearPt(pos); // remove the shear component
			// traverse all periodic cells around the body, to see if any of them touches
			Vector3r halfSize=b->bound->max-b->bound->min; halfSize*=.5;
			Vector3r pmin,pmax;
			Vector3i i;
			for(i[0]=-1; i[0]<=1; i[0]++) for(i[1]=-1;i[1]<=1; i[1]++) for(i[2]=-1; i[2]<=1; i[2]++){
				if(i[0]==0 && i[1]==0 && i[2]==0) continue; // middle; already rendered above
				Vector3r pos2=pos+Vector3r(cellSize[0]*i[0],cellSize[1]*i[1],cellSize[2]*i[2]); // shift, but without shear!
				pmin=pos2-halfSize; pmax=pos2+halfSize;
				if(pmin[0]<=cellSize[0] && pmax[0]>=0 &&
					pmin[1]<=cellSize[1] && pmax[1]>=0 &&
					pmin[2]<=cellSize[2] && pmax[2]>=0) {
					Vector3r pt=scene->cell->shearPt(pos2);
					if(pointClipped(pt)) continue;
					glPushMatrix();
						glTranslatev(pt);
						glRotatef(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
						shapeDispatcher(b->shape,b->state,/*wire*/ true, viewInfo);
					glPopMatrix();
				}
			}
		}
	}
}



void OpenGLRenderingEngine::postProcessAttributes(bool deserializing){
	if(!deserializing) return;
	for(unsigned int i=0;i<boundFunctorNames.size();i++) boundDispatcher.add1DEntry(boundFunctorNames[i][0],boundFunctorNames[i][1]);
	for(unsigned int i=0;i<shapeFunctorNames.size();i++) shapeDispatcher.add1DEntry(shapeFunctorNames[i][0],shapeFunctorNames[i][1]);
	for(unsigned int i=0;i<interactionGeometryFunctorNames.size();i++) interactionGeometryDispatcher.add1DEntry(interactionGeometryFunctorNames[i][0],interactionGeometryFunctorNames[i][1]);
	for(unsigned int i=0;i<interactionPhysicsFunctorNames.size();i++) interactionPhysicsDispatcher.add1DEntry(interactionPhysicsFunctorNames[i][0],interactionPhysicsFunctorNames[i][1]);	
}
void OpenGLRenderingEngine::addInteractionGeometryFunctor(const string& str2){
	string str1 = (static_pointer_cast<GlInteractionGeometryFunctor>(ClassFactory::instance().createShared(str2)))->renders();
	vector<string> v; v.push_back(str1); v.push_back(str2); interactionGeometryFunctorNames.push_back(v);
}
void OpenGLRenderingEngine::addInteractionPhysicsFunctor(const string& str2){
	string str1 = (static_pointer_cast<GlInteractionPhysicsFunctor>(ClassFactory::instance().createShared(str2)))->renders();
	vector<string> v; v.push_back(str1); v.push_back(str2); interactionPhysicsFunctorNames.push_back(v);
}
void OpenGLRenderingEngine::addStateFunctor(const string& str2){
	string str1 = (static_pointer_cast<GlStateFunctor>(ClassFactory::instance().createShared(str2)))->renders();
	vector<string> v; v.push_back(str1); v.push_back(str2); stateFunctorNames.push_back(v);
}
void OpenGLRenderingEngine::addBoundingVolumeFunctor(const string& str2){
	string str1 = (static_pointer_cast<GlBoundFunctor>(ClassFactory::instance().createShared(str2)))->renders();
	vector<string> v; v.push_back(str1); v.push_back(str2); boundFunctorNames.push_back(v);
}
void OpenGLRenderingEngine::addInteractingGeometryFunctor(const string& str2){
	string str1 = (static_pointer_cast<GlShapeFunctor>(ClassFactory::instance().createShared(str2)))->renders();
	vector<string> v; v.push_back(str1); v.push_back(str2); shapeFunctorNames.push_back(v);
}

/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2005 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"GLViewer.hpp"
#include"QtGUI.hpp"
#include"YadeCamera.hpp"
#include"YadeQtMainWindow.hpp"

#include<yade/lib-opengl/FpsTracker.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/core/Body.hpp>
#include<yade/core/Scene.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/core/DisplayParameters.hpp>
#include<boost/filesystem/operations.hpp>
#include<boost/algorithm/string.hpp>
#include<boost/version.hpp>
#include<boost/python.hpp>
#include<yade/lib-serialization/ObjectIO.hpp>
using namespace boost;

#ifdef YADE_GL2PS
	#include<gl2ps.h>
#endif

CREATE_LOGGER(GLViewer);

GLLock::GLLock(GLViewer* _glv): boost::try_mutex::scoped_lock(Omega::instance().renderMutex), glv(_glv){
	glv->makeCurrent();
}
GLLock::~GLLock(){ glv->doneCurrent(); }


void GLViewer::updateGL(void){/*GLLock lock(this); */QGLViewer::updateGL();}

/* additionally try: doneCurrent; glFlush; glSwapBuffers after paintGL */
void GLViewer::paintGL(void){
	/* paintGL encapsulated preDraw, draw and postDraw within QGLViewer. If the mutex cannot be locked,
	 * we just return without repainting */
	boost::try_mutex::scoped_try_lock lock(Omega::instance().renderMutex);
	#if BOOST_VERSION<103500
		if(lock.locked())
	#else
		if(lock.owns_lock())
	#endif
		{
			this->makeCurrent();
			QGLViewer::paintGL();
		}
	this->doneCurrent();
}

GLViewer::~GLViewer(){ /* get the GL mutex when closing */ GLLock lock(this); }

GLViewer::GLViewer(int id, shared_ptr<OpenGLRenderingEngine> _renderer, QWidget * parent, QGLWidget * shareWidget) : QGLViewer(parent,"glview",shareWidget)//, qglThread(this,rendererInit)
{
	isMoving=false;
	renderer=_renderer; renderer->init();
	drawGridXYZ[0]=drawGridXYZ[1]=drawGridXYZ[2]=false;
	drawScale=true;
	timeDispMask=TIME_REAL|TIME_VIRT|TIME_ITER;
	viewId = id;
	cut_plane = 0;
	cut_plane_delta = -2;
	grid_subdivision = false;
	resize(550,550);

	if (id==0) setCaption("Primary view");
	else if(id==-1) { setCaption("Player view"); viewId=0; }
	else setCaption("Secondary view number "+lexical_cast<string>(id));
	show();
	
	mouseMovesCamera();
	manipulatedClipPlane=-1;

	if(manipulatedFrame()==0) setManipulatedFrame(new qglviewer::ManipulatedFrame());

	xyPlaneConstraint=shared_ptr<qglviewer::LocalConstraint>(new qglviewer::LocalConstraint());
	//xyPlaneConstraint->setTranslationConstraint(qglviewer::AxisPlaneConstraint::AXIS,qglviewer::Vec(0,0,1));
	//xyPlaneConstraint->setRotationConstraint(qglviewer::AxisPlaneConstraint::FORBIDDEN,qglviewer::Vec(0,0,1));
	manipulatedFrame()->setConstraint(NULL);

	setKeyDescription(Qt::Key_C,"Set scene center so that all bodies are visible; if a body is selected, center around this body.");
	setKeyDescription(Qt::Key_C & Qt::ALT,"Set scene center to median body position");
	setKeyDescription(Qt::Key_D,"Toggle time display mask");
	setKeyDescription(Qt::Key_G,"Toggle grid");
	setKeyDescription(Qt::Key_X,"Toggle YZ grid (or: align manipulated clip plane normal with +X)");
	setKeyDescription(Qt::Key_Y,"Toggle XZ grid (or: align manipulated clip plane normal with +Y)");
	setKeyDescription(Qt::Key_Z,"Toggle XY grid (or: align manipulated clip plane normal with +Z)");
	setKeyDescription(Qt::Key_Period,"Toggle grid subdivision by 10");
	setKeyDescription(Qt::Key_S & Qt::ALT,   "Save QGLViewer state to /tmp/qglviewerState.xml");
	setKeyDescription(Qt::Key_T,"Switch orthographic / perspective camera");
	setKeyDescription(Qt::Key_O,"Set narrower field of view");
	setKeyDescription(Qt::Key_P,"Set wider field of view");
	setKeyDescription(Qt::Key_R,"Revolve around scene center");
	setKeyDescription(Qt::Key_V,"Save PDF of the current view to /tmp/yade-snapshot-0001.pdf (whichever number is available first). (Must be compiled with the gl2ps feature.)");
	setKeyDescription(Qt::Key_Plus,    "Cut plane increase");
	setKeyDescription(Qt::Key_Minus,   "Cut plane decrease");
	setKeyDescription(Qt::Key_Slash,   "Cut plane step decrease");
	setKeyDescription(Qt::Key_Asterisk,"Cut plane step increase");
	// unset default shortcuts
 	setPathKey(-Qt::Key_F1);
 	setPathKey(-Qt::Key_F2);
	setKeyDescription(Qt::Key_Escape,"Manipulate scene (default)");
	setKeyDescription(Qt::Key_F1,"Manipulate clipping plane #1 (2,...)");
	setKeyDescription(Qt::Key_1,"Make the manipulated clipping plane parallel with plane #1 (2,...)");
	setKeyDescription(Qt::Key_1 & Qt::ALT,"Add/remove plane #1 (2,...) to/from the bound group");
	setKeyDescription(Qt::Key_0,"Clear the bound group");
	setKeyDescription(Qt::Key_Space,"Activate/deactivate the manipulated clipping plane");

	centerScene();
}

void GLViewer::mouseMovesCamera(){
	camera()->frame()->setWheelSensitivity(-1.0f);

	setMouseBinding(Qt::SHIFT + Qt::LeftButton, SELECT);
	//setMouseBinding(Qt::RightButton, NO_CLICK_ACTION);
	setMouseBinding(Qt::SHIFT + Qt::LeftButton + Qt::RightButton, FRAME, ZOOM);
	setMouseBinding(Qt::SHIFT + Qt::MidButton, FRAME, TRANSLATE);
	setMouseBinding(Qt::SHIFT + Qt::RightButton, FRAME, ROTATE);
	setWheelBinding(Qt::ShiftButton , FRAME, ZOOM);

	setMouseBinding(Qt::LeftButton + Qt::RightButton, CAMERA, ZOOM);
	setMouseBinding(Qt::MidButton, CAMERA, ZOOM);
	setMouseBinding(Qt::LeftButton, CAMERA, ROTATE);
	setMouseBinding(Qt::RightButton, CAMERA, TRANSLATE);
	setWheelBinding(Qt::NoButton, CAMERA, ZOOM);
};

void GLViewer::mouseMovesManipulatedFrame(qglviewer::Constraint* c){
	setMouseBinding(Qt::LeftButton + Qt::RightButton, FRAME, ZOOM);
	setMouseBinding(Qt::MidButton, FRAME, ZOOM);
	setMouseBinding(Qt::LeftButton, FRAME, ROTATE);
	setMouseBinding(Qt::RightButton, FRAME, TRANSLATE);
	setWheelBinding(Qt::NoButton , FRAME, ZOOM);
	manipulatedFrame()->setConstraint(c);
}

bool GLViewer::isManipulating(){
	return isMoving || manipulatedClipPlane>=0;
}

void GLViewer::resetManipulation(){
	mouseMovesCamera();
	setSelectedName(-1);
	isMoving=false;
	manipulatedClipPlane=-1;
}

void GLViewer::startClipPlaneManipulation(int planeNo){
	assert(planeNo<renderer->numClipPlanes);
	resetManipulation();
	mouseMovesManipulatedFrame(xyPlaneConstraint.get());
	manipulatedClipPlane=planeNo;
	const Se3r se3(renderer->clipPlaneSe3[planeNo]);
	manipulatedFrame()->setPositionAndOrientation(qglviewer::Vec(se3.position[0],se3.position[1],se3.position[2]),qglviewer::Quaternion(se3.orientation.x(),se3.orientation.y(),se3.orientation.z(),se3.orientation.w()));
	string grp=strBoundGroup();
	displayMessage("Manipulating clip plane #"+lexical_cast<string>(planeNo+1)+(grp.empty()?grp:" (bound planes:"+grp+")"));
}

void GLViewer::useDisplayParameters(size_t n){
	LOG_DEBUG("Loading display parameters from #"<<n);
	vector<shared_ptr<DisplayParameters> >& dispParams=Omega::instance().getScene()->dispParams;
	if(dispParams.size()<=(size_t)n){LOG_ERROR("Display parameters #"<<n<<" don't exist (number of entries "<<dispParams.size()<<")"); return;}
	const shared_ptr<DisplayParameters>& dp=dispParams[n];
	string val;
	if(dp->getValue("OpenGLRenderingEngine",val)){ istringstream oglre(val);
		#ifdef YADE_SERIALIZE_USING_BOOST 
			yade::ObjectIO::load<typeof(renderer),boost::archive::xml_iarchive>(oglre,"renderer",renderer);
		#else
			IOFormatManager::loadFromStream("XMLFormatManager",oglre,"renderer",renderer);
		#endif
	}
	else { LOG_WARN("OpenGLRenderingEngine configuration not found in display parameters, skipped.");}
	if(dp->getValue("GLViewer",val)){ GLViewer::setState(val);}
	else { LOG_WARN("GLViewer configuration not found in display parameters, skipped."); }
}

void GLViewer::saveDisplayParameters(size_t n){
	LOG_DEBUG("Saving display parameters to #"<<n);
	vector<shared_ptr<DisplayParameters> >& dispParams=Omega::instance().getScene()->dispParams;
	if(dispParams.size()<=n){while(dispParams.size()<=n) dispParams.push_back(shared_ptr<DisplayParameters>(new DisplayParameters));} assert(n<dispParams.size());
	shared_ptr<DisplayParameters>& dp=dispParams[n];
	ostringstream oglre;
	#ifdef YADE_SERIALIZE_USING_BOOST
		yade::ObjectIO::save<typeof(renderer),boost::archive::xml_oarchive>(oglre,"renderer",renderer);
	#else
		IOFormatManager::saveToStream("XMLFormatManager",oglre,"renderer",renderer);
	#endif
	dp->setValue("OpenGLRenderingEngine",oglre.str());
	dp->setValue("GLViewer",GLViewer::getState());
}

string GLViewer::getState(){
	QString origStateFileName=stateFileName();
	string tmpFile=Omega::instance().tmpFilename();
	setStateFileName(tmpFile); saveStateToFile(); setStateFileName(origStateFileName);
	LOG_DEBUG("State saved to temp file "<<tmpFile);
	// read tmp file contents and return it as string
	// this will replace all whitespace by space (nowlines will disappear, which is what we want)
	ifstream in(tmpFile.c_str()); string ret; while(!in.eof()){string ss; in>>ss; ret+=" "+ss;}; in.close();
	boost::filesystem::remove(boost::filesystem::path(tmpFile));
	return ret;
}

void GLViewer::setState(string state){
	string tmpFile=Omega::instance().tmpFilename();
	std::ofstream out(tmpFile.c_str());
	if(!out.good()){ LOG_ERROR("Error opening temp file `"<<tmpFile<<"', loading aborted."); return; }
	out<<state; out.close();
	LOG_DEBUG("Will load state from temp file "<<tmpFile);
	QString origStateFileName=stateFileName(); setStateFileName(tmpFile); restoreStateFromFile(); setStateFileName(origStateFileName);
	boost::filesystem::remove(boost::filesystem::path(tmpFile));
}

void GLViewer::keyPressEvent(QKeyEvent *e)
{
	last_user_event = boost::posix_time::second_clock::local_time();

	if(false){}
	/* special keys: Escape and Space */
	//else if(e->key()==Qt::Key_F9 || e->key()==Qt::Key_F10 || e->key()==Qt::Key_F11 || e->key()==Qt::Key_F12){ YadeQtMainWindow::self->closeView(this); }
	else if(e->key()==Qt::Key_Escape){
		if(!isManipulating()){ /* CRASH, deleting this: YadeQtMainWindow::self->closeView(this); */ return; }
		else { resetManipulation(); displayMessage("Manipulating scene."); }
	}
	else if(e->key()==Qt::Key_Space){
		if(manipulatedClipPlane>=0) {displayMessage("Clip plane #"+lexical_cast<string>(manipulatedClipPlane+1)+(renderer->clipPlaneActive[manipulatedClipPlane]?" de":" ")+"activated"); renderer->clipPlaneActive[manipulatedClipPlane]=!renderer->clipPlaneActive[manipulatedClipPlane]; }
	}
	/* function keys */
	else if(e->key()==Qt::Key_F1 || e->key()==Qt::Key_F2 || e->key()==Qt::Key_F3 /* || ... */ ){
		int n=0; if(e->key()==Qt::Key_F1) n=1; else if(e->key()==Qt::Key_F2) n=2; else if(e->key()==Qt::Key_F3) n=3; assert(n>0); int planeId=n-1;
		if(planeId>=renderer->numClipPlanes) return;
		if(planeId!=manipulatedClipPlane) startClipPlaneManipulation(planeId);
	}
	/* numbers */
	else if(e->key()==Qt::Key_0 &&(e->state() & AltButton)) { boundClipPlanes.clear(); displayMessage("Cleared bound planes group.");}
	else if(e->key()==Qt::Key_1 || e->key()==Qt::Key_2 || e->key()==Qt::Key_3 /* || ... */ ){
		int n=0; if(e->key()==Qt::Key_1) n=1; else if(e->key()==Qt::Key_2) n=2; else if(e->key()==Qt::Key_3) n=3; assert(n>0); int planeId=n-1;
		if(planeId>=renderer->numClipPlanes) return; // no such clipping plane
		if(e->state() & AltButton){
			if(boundClipPlanes.count(planeId)==0) {boundClipPlanes.insert(planeId); displayMessage("Added plane #"+lexical_cast<string>(planeId+1)+" to the bound group: "+strBoundGroup());}
			else {boundClipPlanes.erase(planeId); displayMessage("Removed plane #"+lexical_cast<string>(planeId+1)+" from the bound group: "+strBoundGroup());}
		}
		else if(manipulatedClipPlane>=0 && manipulatedClipPlane!=planeId) {
			const Quaternionr& o=renderer->clipPlaneSe3[planeId].orientation;
			manipulatedFrame()->setOrientation(qglviewer::Quaternion(o.x(),o.y(),o.z(),o.w()));
			displayMessage("Copied orientation from plane #1");
		}
	}
	else if(e->key()==Qt::Key_7 || e->key()==Qt::Key_8 || e->key()==Qt::Key_9){
		int nn=-1; if(e->key()==Qt::Key_7)nn=0; else if(e->key()==Qt::Key_8)nn=1; else if(e->key()==Qt::Key_9)nn=2; assert(nn>=0); size_t n=(size_t)nn;
		if(e->state() & AltButton) saveDisplayParameters(n);
		else useDisplayParameters(n);
	}
	/* letters alphabetically */
	else if(e->key()==Qt::Key_C && (e->state() & AltButton)){ displayMessage("Median centering"); centerMedianQuartile(); }
	else if(e->key()==Qt::Key_C){
		// center around selected body
		if(selectedName() >= 0 && (*(Omega::instance().getScene()->bodies)).exists(selectedName())) setSceneCenter(manipulatedFrame()->position());
		// make all bodies visible
		else centerScene();
	}
	else if(e->key()==Qt::Key_D &&(e->state() & AltButton)){ body_id_t id; if((id=Omega::instance().getScene()->selectedBody)>=0){ const shared_ptr<Body>& b=Body::byId(id); b->isDynamic=!b->isDynamic; LOG_INFO("Body #"<<id<<" now "<<(b->isDynamic?"":"NOT")<<" dynamic"); } }
	else if(e->key()==Qt::Key_D) {timeDispMask+=1; if(timeDispMask>(TIME_REAL|TIME_VIRT|TIME_ITER))timeDispMask=0; }
	else if(e->key()==Qt::Key_G) {bool anyDrawn=drawGridXYZ[0]||drawGridXYZ[1]||drawGridXYZ[2]; for(int i=0; i<3; i++)drawGridXYZ[i]=!anyDrawn; }
	else if (e->key()==Qt::Key_M && selectedName() >= 0){
		if(!(isMoving=!isMoving)){displayMessage("Moving done."); mouseMovesCamera();}
		else{ displayMessage("Moving selected object"); mouseMovesManipulatedFrame();}
	}
	else if (e->key() == Qt::Key_T) camera()->setType(camera()->type()==qglviewer::Camera::ORTHOGRAPHIC ? qglviewer::Camera::PERSPECTIVE : qglviewer::Camera::ORTHOGRAPHIC);
	else if(e->key()==Qt::Key_O) camera()->setFieldOfView(camera()->fieldOfView()*0.9);
	else if(e->key()==Qt::Key_P) camera()->setFieldOfView(camera()->fieldOfView()*1.1);
	else if(e->key()==Qt::Key_R){ // reverse the clipping plane; revolve around scene center if no clipping plane selected
		if(manipulatedClipPlane>=0 && manipulatedClipPlane<renderer->numClipPlanes){
			/* here, we must update both manipulatedFrame orientation and renderer->clipPlaneSe3 orientation in the same way */
			Quaternionr& ori=renderer->clipPlaneSe3[manipulatedClipPlane].orientation;
			ori=Quaternionr(AngleAxisr(Mathr::PI,Vector3r(0,1,0)))*ori; 
			manipulatedFrame()->setOrientation(qglviewer::Quaternion(qglviewer::Vec(0,1,0),Mathr::PI)*manipulatedFrame()->orientation());
			displayMessage("Plane #"+lexical_cast<string>(manipulatedClipPlane+1)+" reversed.");
		}
		else {
			camera()->setRevolveAroundPoint(sceneCenter());
		}
	}
	else if(e->key()==Qt::Key_S){
		LOG_INFO("Saving QGLViewer state to /tmp/qglviewerState.xml");
		setStateFileName("/tmp/qglviewerState.xml"); saveStateToFile(); setStateFileName(QString::null);
	}
	else if(e->key()==Qt::Key_X || e->key()==Qt::Key_Y || e->key()==Qt::Key_Z){
		int axisIdx=(e->key()==Qt::Key_X?0:(e->key()==Qt::Key_Y?1:2));
		if(manipulatedClipPlane<0){ drawGridXYZ[axisIdx]=!drawGridXYZ[axisIdx]; }
		else{ // align clipping plane with world axis
			// x: (0,1,0),pi/2; y: (0,0,1),pi/2; z: (1,0,0),0
			qglviewer::Vec axis(0,0,0); axis[(axisIdx+1)%3]=1;
			manipulatedFrame()->setOrientation(qglviewer::Quaternion(axis,axisIdx==2?0:Mathr::PI/2));
		}
	}
	else if(e->key()==Qt::Key_Period) grid_subdivision = !grid_subdivision;
#ifdef YADE_GL2PS
	else if(e->key()==Qt::Key_V){
		for(int i=0; ;i++){
			std::ostringstream fss; fss<<"/tmp/yade-snapshot-"<<setw(4)<<setfill('0')<<i<<".pdf";
			if(!boost::filesystem::exists(fss.str())){ nextFrameSnapshotFilename=fss.str(); break; }
		}
		LOG_INFO("Will save snapshot to "<<nextFrameSnapshotFilename);
	}
#endif

//////////////////////////////////////////////
// FIXME that all should be in some nice GUI
// 
// Cutting plane
	else if( e->key()==Qt::Key_Plus ){
			cut_plane = std::min(1.0, cut_plane + std::pow(10.0,(double)cut_plane_delta));
			static_cast<YadeCamera*>(camera())->setCuttingDistance(cut_plane);
			displayMessage("Cut plane: "+lexical_cast<std::string>(cut_plane));
	}else if( e->key()==Qt::Key_Minus ){
			cut_plane = std::max(0.0, cut_plane - std::pow(10.0,(double)cut_plane_delta));
			static_cast<YadeCamera*>(camera())->setCuttingDistance(cut_plane);
			displayMessage("Cut plane: "+lexical_cast<std::string>(cut_plane));
	}else if( e->key()==Qt::Key_Slash ){
			cut_plane_delta -= 1;
			displayMessage("Cut plane increment: 1e"+(cut_plane_delta>0?std::string("+"):std::string(""))+lexical_cast<std::string>(cut_plane_delta));
	}else if( e->key()==Qt::Key_Asterisk ){
			cut_plane_delta = std::min(1+cut_plane_delta,-1);
			displayMessage("Cut plane increment: 1e"+(cut_plane_delta>0?std::string("+"):std::string(""))+lexical_cast<std::string>(cut_plane_delta));
	}

// FIXME END
//////////////////////////////////////////////
//
	else if(e->key()!=Qt::Key_Escape && e->key()!=Qt::Key_Space) QGLViewer::keyPressEvent(e);
	updateGL();
}
/* Center the scene such that periodic cell is contained in the view */
void GLViewer::centerPeriodic(){
	Scene* scene=Omega::instance().getScene().get();
	assert(scene->isPeriodic);
	Vector3r center=.5*scene->cell->getSize();
	Vector3r halfSize=.5*scene->cell->getSize();
	float radius=std::max(halfSize[0],std::max(halfSize[1],halfSize[2]));
	LOG_DEBUG("Periodic scene center="<<center<<", halfSize="<<halfSize<<", radius="<<radius);
	setSceneCenter(qglviewer::Vec(center[0],center[1],center[2]));
	setSceneRadius(radius*1.5);
	showEntireScene();
}

/* Calculate medians for x, y and z coordinates of all bodies;
 *then set scene center to median position and scene radius to 2*inter-quartile distance.
 *
 * This function eliminates the effect of lonely bodies that went nuts and enlarge
 * the scene's Aabb in such a way that fitting the scene to see the Aabb makes the
 * "central" (where most bodies is) part very small or even invisible.
 */
void GLViewer::centerMedianQuartile(){
	Scene* scene=Omega::instance().getScene().get();
	if(scene->isPeriodic){ centerPeriodic(); return; }
	long nBodies=scene->bodies->size();
	if(nBodies<4) {
		LOG_INFO("Less than 4 bodies, median makes no sense; calling centerScene() instead.");
		return centerScene();
	}
	std::vector<Real> coords[3];
	for(int i=0;i<3;i++)coords[i].reserve(nBodies);
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		for(int i=0; i<3; i++) coords[i].push_back(b->state->pos[i]);
	}
	Vector3r median,interQuart;
	for(int i=0;i<3;i++){
		sort(coords[i].begin(),coords[i].end());
		median[i]=*(coords[i].begin()+nBodies/2);
		interQuart[i]=*(coords[i].begin()+3*nBodies/4)-*(coords[i].begin()+nBodies/4);
	}
	LOG_DEBUG("Median position is"<<median<<", inter-quartile distance is "<<interQuart);
	setSceneCenter(qglviewer::Vec(median[0],median[1],median[2]));
	setSceneRadius(2*(interQuart[0]+interQuart[1]+interQuart[2])/3.);
	showEntireScene();
}

void GLViewer::centerScene(){
	Scene* rb=Omega::instance().getScene().get();
	if (!rb) return;
	if(rb->isPeriodic){ centerPeriodic(); return; }

	if(rb->bodies->size()<renderer->selectBodyLimit){LOG_INFO("Less than "+lexical_cast<string>(renderer->selectBodyLimit)+" bodies, moving possible. Select with shift, press 'm' to move.");}
	else{LOG_INFO("More than "+lexical_cast<string>(renderer->selectBodyLimit)+" (OpenGLRenderingEngine::selectBodyLimit) bodies. Moving not possible.");}
	Vector3r min,max;	
	if(rb->bound){
		min=rb->bound->min; max=rb->bound->max;
		if(std::max(max[0]-min[0],std::max(max[1]-min[1],max[2]-min[2]))<=0){
			// Aabb is not yet calculated...
			LOG_DEBUG("rootBody's bound not yet calculated or has one dimension zero, attempt get that from bodies' positions.");
			Real inf=std::numeric_limits<Real>::infinity();
			min=Vector3r(inf,inf,inf); max=Vector3r(-inf,-inf,-inf);
			FOREACH(const shared_ptr<Body>& b, *rb->bodies){
				max=componentMaxVector(max,b->state->pos);
				min=componentMinVector(min,b->state->pos);
			}
		} else {LOG_DEBUG("Using rootBody's Aabb");}
	} else {
		LOG_DEBUG("No rootBody's Aabb; setting scene in cube (-1,-1,-1)x(1,1,1)");
		min=Vector3r(-1,-1,-1); max=Vector3r(1,1,1);
	}
	LOG_DEBUG("Got scene box min="<<min<<" and max="<<max);
	Vector3r center = (max+min)*0.5;
	Vector3r halfSize = (max-min)*0.5;
	float radius=std::max(halfSize[0],std::max(halfSize[1],halfSize[2])); if(radius<=0) radius=1;
	LOG_DEBUG("Scene center="<<center<<", halfSize="<<halfSize<<", radius="<<radius);
	setSceneCenter(qglviewer::Vec(center[0],center[1],center[2]));
	setSceneRadius(radius*1.5);
	showEntireScene();
}

void GLViewer::draw()
{
#ifdef YADE_GL2PS
	if(!nextFrameSnapshotFilename.empty() && boost::algorithm::ends_with(nextFrameSnapshotFilename,".pdf")){
		gl2psStream=fopen(nextFrameSnapshotFilename.c_str(),"wb");
		if(!gl2psStream){ int err=errno; throw runtime_error(string("Error opening file ")+nextFrameSnapshotFilename+": "+strerror(err)); }
		LOG_DEBUG("Start saving snapshot to "<<nextFrameSnapshotFilename);
		size_t nBodies=Omega::instance().getScene()->bodies->size();
		int sortAlgo=(nBodies<100 ? GL2PS_BSP_SORT : GL2PS_SIMPLE_SORT);
		gl2psBeginPage(/*const char *title*/"Some title", /*const char *producer*/ "Yade",
			/*GLint viewport[4]*/ NULL,
			/*GLint format*/ GL2PS_PDF, /*GLint sort*/ sortAlgo, /*GLint options*/GL2PS_SIMPLE_LINE_OFFSET|GL2PS_USE_CURRENT_VIEWPORT|GL2PS_TIGHT_BOUNDING_BOX|GL2PS_COMPRESS|GL2PS_OCCLUSION_CULL|GL2PS_NO_BLENDING, 
			/*GLint colormode*/ GL_RGBA, /*GLint colorsize*/0, 
			/*GL2PSrgba *colortable*/NULL, 
			/*GLint nr*/0, /*GLint ng*/0, /*GLint nb*/0, 
			/*GLint buffersize*/4096*4096 /* 16MB */, /*FILE *stream*/ gl2psStream,
			/*const char *filename*/NULL);
	}
#endif

	qglviewer::Vec vd=camera()->viewDirection(); renderer->viewDirection=Vector3r(vd[0],vd[1],vd[2]);
	if(Omega::instance().getScene()){
		int selection = selectedName();
		if(selection!=-1 && (*(Omega::instance().getScene()->bodies)).exists(selection) && isMoving){
			static int last(-1);
			if(last == selection) // delay by one redraw, so the body will not jump into 0,0,0 coords
			{
				Quaternionr& q = (*(Omega::instance().getScene()->bodies))[selection]->state->ori;
				Vector3r&    v = (*(Omega::instance().getScene()->bodies))[selection]->state->pos;
				float v0,v1,v2; manipulatedFrame()->getPosition(v0,v1,v2);v[0]=v0;v[1]=v1;v[2]=v2;
				double q0,q1,q2,q3; manipulatedFrame()->getOrientation(q0,q1,q2,q3);	q.x()=q0;q.y()=q1;q.z()=q2;q.w()=q3;
			}
			(*(Omega::instance().getScene()->bodies))[selection]->userForcedDisplacementRedrawHook();	
			last=selection;
		}
		if(manipulatedClipPlane>=0){
			assert(manipulatedClipPlane<renderer->numClipPlanes);
			float v0,v1,v2; manipulatedFrame()->getPosition(v0,v1,v2);
			double q0,q1,q2,q3; manipulatedFrame()->getOrientation(q0,q1,q2,q3);
			Se3r newSe3(Vector3r(v0,v1,v2),Quaternionr(q0,q1,q2,q3)); newSe3.orientation.normalize();
			const Se3r& oldSe3=renderer->clipPlaneSe3[manipulatedClipPlane];
			FOREACH(int planeId, boundClipPlanes){
				if(planeId>=renderer->numClipPlanes || !renderer->clipPlaneActive[planeId] || planeId==manipulatedClipPlane) continue;
				Se3r& boundSe3=renderer->clipPlaneSe3[planeId];
				Quaternionr relOrient=oldSe3.orientation.conjugate()*boundSe3.orientation; relOrient.normalize();
				Vector3r relPos=oldSe3.orientation.conjugate()*(boundSe3.position-oldSe3.position);
				boundSe3.position=newSe3.position+newSe3.orientation*relPos;
				boundSe3.orientation=newSe3.orientation*relOrient;
				boundSe3.orientation.normalize();
			}
			renderer->clipPlaneSe3[manipulatedClipPlane]=newSe3;
		}
		renderer->render(Omega::instance().getScene(), selectedName());
	}
}

void GLViewer::drawWithNames(){
	qglviewer::Vec vd=camera()->viewDirection(); renderer->viewDirection=Vector3r(vd[0],vd[1],vd[2]);
	if(Omega::instance().getScene() && Omega::instance().getScene()->bodies->size()<renderer->selectBodyLimit) renderer->renderWithNames(Omega::instance().getScene());
}

// new object selected.
// set frame coordinates, and isDynamic=false;
void GLViewer::postSelection(const QPoint& point) 
{
	LOG_DEBUG("Selection is "<<selectedName());
	int selection = selectedName();
	if(selection<0){
		if(isMoving){
			displayMessage("Moving finished"); mouseMovesCamera(); isMoving=false;
			Omega::instance().getScene()->selectedBody = -1;
		}
		return;
	}
	if(selection>=0 && (*(Omega::instance().getScene()->bodies)).exists(selection)){
		resetManipulation();
		if(Body::byId(body_id_t(selection))->isClumpMember()){ // select clump (invisible) instead of its member
			LOG_DEBUG("Clump member #"<<selection<<" selected, selecting clump instead.");
			selection=Body::byId(body_id_t(selection))->clumpId;
		}
		setSelectedName(selection);
		LOG_DEBUG("New selection "<<selection);
		displayMessage("Selected body #"+lexical_cast<string>(selection)+(Body::byId(selection)->isClump()?" (clump)":""));
		Quaternionr& q = Body::byId(selection)->state->ori;
		Vector3r&    v = Body::byId(selection)->state->pos;
		manipulatedFrame()->setPositionAndOrientation(qglviewer::Vec(v[0],v[1],v[2]),qglviewer::Quaternion(q.x(),q.y(),q.z(),q.w()));
		Omega::instance().getScene()->selectedBody = selection;
			PyGILState_STATE gstate;
			gstate = PyGILState_Ensure();
				python::object main=python::import("__main__");
				python::object global=main.attr("__dict__");
				// the try/catch block must be properly nested inside PyGILState_Ensure and PyGILState_Release
				try{
					python::eval(string("onBodySelect("+lexical_cast<string>(selection)+")").c_str(),global,global);
				} catch (python::error_already_set const &) {
					LOG_DEBUG("unable to call onBodySelect. Not defined?");
				}
			PyGILState_Release(gstate);
			// see https://svn.boost.org/trac/boost/ticket/2781 for exception handling
	}
}

// maybe new object will be selected.
// if so, then set isDynamic of previous selection, to old value
void GLViewer::endSelection(const QPoint &point){
	manipulatedClipPlane=-1;
	//int old = selectedName();
	QGLViewer::endSelection(point);
	// if(old != -1 && old!=selectedName() && (*(Omega::instance().getScene()->bodies)).exists(old)) Body::byId(old)->isDynamic = wasDynamic;
}

qglviewer::Vec GLViewer::displayedSceneCenter(){
	return camera()->unprojectedCoordinatesOf(qglviewer::Vec(width()/2 /* pixels */ ,height()/2 /* pixels */, /*middle between near plane and far plane*/ .5));
}

float GLViewer::displayedSceneRadius(){
	return (camera()->unprojectedCoordinatesOf(qglviewer::Vec(width()/2,height()/2,.5))-camera()->unprojectedCoordinatesOf(qglviewer::Vec(0,0,.5))).norm();
}

void GLViewer::postDraw(){
	Real wholeDiameter=QGLViewer::camera()->sceneRadius()*2;

	renderer->viewInfo.sceneRadius=QGLViewer::camera()->sceneRadius();
	qglviewer::Vec c=QGLViewer::camera()->sceneCenter();
	renderer->viewInfo.sceneCenter=Vector3r(c[0],c[1],c[2]);

	Real dispDiameter=min(wholeDiameter,max((Real)displayedSceneRadius()*2,wholeDiameter/1e3)); // limit to avoid drawing 1e5 lines with big zoom level
	//qglviewer::Vec center=QGLViewer::camera()->sceneCenter();
	Real gridStep=pow(10,(floor(log10(dispDiameter)-.7)));
	Real scaleStep=pow(10,(floor(log10(displayedSceneRadius()*2)-.7))); // unconstrained
	int nSegments=((int)(wholeDiameter/gridStep))+1;
	Real realSize=nSegments*gridStep;
	//LOG_TRACE("nSegments="<<nSegments<<",gridStep="<<gridStep<<",realSize="<<realSize);
	glPushMatrix();

	nSegments *= 2; // there's an error in QGLViewer::drawGrid(), so we need to mitigate it by '* 2'
	// XYZ grids
	glLineWidth(.5);
	if(drawGridXYZ[0]) {glColor3f(0.6,0.3,0.3); glPushMatrix(); glRotated(90.,0.,1.,0.); QGLViewer::drawGrid(realSize,nSegments); glPopMatrix();}
	if(drawGridXYZ[1]) {glColor3f(0.3,0.6,0.3); glPushMatrix(); glRotated(90.,1.,0.,0.); QGLViewer::drawGrid(realSize,nSegments); glPopMatrix();}
	if(drawGridXYZ[2]) {glColor3f(0.3,0.3,0.6); glPushMatrix(); /*glRotated(90.,0.,1.,0.);*/ QGLViewer::drawGrid(realSize,nSegments); glPopMatrix();}
	if(grid_subdivision){
	if(drawGridXYZ[0]) {glColor3f(0.4,0.1,0.1); glPushMatrix(); glRotated(90.,0.,1.,0.); QGLViewer::drawGrid(realSize,nSegments*10); glPopMatrix();}
	if(drawGridXYZ[1]) {glColor3f(0.1,0.4,0.1); glPushMatrix(); glRotated(90.,1.,0.,0.); QGLViewer::drawGrid(realSize,nSegments*10); glPopMatrix();}
	if(drawGridXYZ[2]) {glColor3f(0.1,0.1,0.4); glPushMatrix(); /*glRotated(90.,0.,1.,0.);*/ QGLViewer::drawGrid(realSize,nSegments*10); glPopMatrix();}
	}
	
	// scale
	if(drawScale){
		Real segmentSize=scaleStep;
		qglviewer::Vec screenDxDy[3]; // dx,dy for x,y,z scale segments
		int extremalDxDy[2]={0,0};
		for(int axis=0; axis<3; axis++){
			qglviewer::Vec delta(0,0,0); delta[axis]=segmentSize;
			qglviewer::Vec center=displayedSceneCenter();
			screenDxDy[axis]=camera()->projectedCoordinatesOf(center+delta)-camera()->projectedCoordinatesOf(center);
			for(int xy=0;xy<2;xy++)extremalDxDy[xy]=(axis>0 ? min(extremalDxDy[xy],(int)screenDxDy[axis][xy]) : screenDxDy[axis][xy]);
		}
		//LOG_DEBUG("Screen offsets for axes: "<<" x("<<screenDxDy[0][0]<<","<<screenDxDy[0][1]<<") y("<<screenDxDy[1][0]<<","<<screenDxDy[1][1]<<") z("<<screenDxDy[2][0]<<","<<screenDxDy[2][1]<<")");
		int margin=10; // screen pixels
		int scaleCenter[2]; scaleCenter[0]=abs(extremalDxDy[0])+margin; scaleCenter[1]=abs(extremalDxDy[1])+margin;
		//LOG_DEBUG("Center of scale "<<scaleCenter[0]<<","<<scaleCenter[1]);
		//displayMessage(QString().sprintf("displayed scene radius %g",displayedSceneRadius()));
		startScreenCoordinatesSystem();
			glDisable(GL_LIGHTING);
			glDisable(GL_DEPTH_TEST);
			glLineWidth(3.0);
			for(int axis=0; axis<3; axis++){
				Vector3r color(.4,.4,.4); color[axis]=.9;
				glColor3v(color);
				glBegin(GL_LINES);
				glVertex2f(scaleCenter[0],scaleCenter[1]);
				glVertex2f(scaleCenter[0]+screenDxDy[axis][0],scaleCenter[1]+screenDxDy[axis][1]);
				glEnd();
			}
			glLineWidth(1.);
			glEnable(GL_DEPTH_TEST);
			QGLViewer::drawText(scaleCenter[0],scaleCenter[1],QString().sprintf("%.3g",(double)scaleStep));
		stopScreenCoordinatesSystem();
	}

	// cutting planes (should be moved to OpenGLRenderingEngine perhaps?)
	// only painted if one of those is being manipulated
	if(manipulatedClipPlane>=0){
		for(int planeId=0; planeId<renderer->numClipPlanes; planeId++){
			if(!renderer->clipPlaneActive[planeId] && planeId!=manipulatedClipPlane) continue;
			glPushMatrix();
				const Se3r& se3=renderer->clipPlaneSe3[planeId];
				AngleAxisr aa(angleAxisFromQuat(se3.orientation));	
				glTranslatef(se3.position[0],se3.position[1],se3.position[2]);
				glRotated(aa.angle()*Mathr::RAD_TO_DEG,aa.axis()[0],aa.axis()[1],aa.axis()[2]);
				Real cff=1;
				if(!renderer->clipPlaneActive[planeId]) cff=.4;
				glColor3f(max((Real)0.,cff*cos(planeId)),max((Real)0.,cff*sin(planeId)),planeId==manipulatedClipPlane); // variable colors
				QGLViewer::drawGrid(realSize,2*nSegments);
				drawArrow(wholeDiameter/6);
			glPopMatrix();
		}
	}
	
	Scene* rb=Omega::instance().getScene().get();
	#define _W3 setw(3)<<setfill('0')
	#define _W2 setw(2)<<setfill('0')
	if(timeDispMask!=0){
		const int lineHt=13;
		unsigned x=10,y=height()-3-lineHt*2;
		glColor3v(Vector3r(1,1,1));
		if(timeDispMask & GLViewer::TIME_VIRT){
			ostringstream oss;
			Real t=Omega::instance().getSimulationTime();
			unsigned min=((unsigned)t/60),sec=(((unsigned)t)%60),msec=((unsigned)(1e3*t))%1000,usec=((unsigned long)(1e6*t))%1000,nsec=((unsigned long)(1e9*t))%1000;
			if(min>0) oss<<_W2<<min<<":"<<_W2<<sec<<"."<<_W3<<msec<<"m"<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else if (sec>0) oss<<_W2<<sec<<"."<<_W3<<msec<<"m"<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else if (msec>0) oss<<_W3<<msec<<"m"<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else if (usec>0) oss<<_W3<<usec<<"u"<<_W3<<nsec<<"n";
			else oss<<_W3<<nsec<<"ns";
			QGLViewer::drawText(x,y,oss.str());
			y-=lineHt;
		}
		glColor3v(Vector3r(0,.5,.5));
		if(timeDispMask & GLViewer::TIME_REAL){
			QGLViewer::drawText(x,y,getRealTimeString() /* virtual, since player gets that from db */);
			y-=lineHt;
		}
		if(timeDispMask & GLViewer::TIME_ITER){
			ostringstream oss;
			oss<<"#"<<Omega::instance().getCurrentIteration();
			if(rb->stopAtIteration>rb->currentIteration) oss<<" ("<<setiosflags(ios::fixed)<<setw(3)<<setprecision(1)<<setfill('0')<<(100.*rb->currentIteration)/rb->stopAtIteration<<"%)";
			QGLViewer::drawText(x,y,oss.str());
			y-=lineHt;
		}
		if(drawGridXYZ[0] || drawGridXYZ[1] || drawGridXYZ[2]){
			glColor3v(Vector3r(1,1,0));
			ostringstream oss;
			oss<<"grid: "<<setprecision(4)<<gridStep;
			if(grid_subdivision) oss<<" (minor "<<setprecision(4)<<gridStep*.1<<")";
			QGLViewer::drawText(x,y,oss.str());
			y-=lineHt;
		}
	}
	QGLViewer::postDraw();
	if(!nextFrameSnapshotFilename.empty()){
		#ifdef YADE_GL2PS
			if(boost::algorithm::ends_with(nextFrameSnapshotFilename,".pdf")){
				gl2psEndPage();
				LOG_DEBUG("Finished saving snapshot to "<<nextFrameSnapshotFilename);
				fclose(gl2psStream);
			} else
	#endif
		{
			// save the snapshot
			saveSnapshot(QString(nextFrameSnapshotFilename),/*overwrite*/ true);
		}
		// notify the caller that it is done already (probably not an atomic op :-|, though)
		nextFrameSnapshotFilename.clear();
	}
}

string GLViewer::getRealTimeString(){
	ostringstream oss;
	time_duration t=Omega::instance().getComputationDuration();
	unsigned d=t.hours()/24,h=t.hours()%24,m=t.minutes(),s=t.seconds();
	oss<<"clock ";
	if(d>0) oss<<d<<"days "<<_W2<<h<<":"<<_W2<<m<<":"<<_W2<<s;
	else if(h>0) oss<<_W2<<h<<":"<<_W2<<m<<":"<<_W2<<s;
	else oss<<_W2<<m<<":"<<_W2<<s;
	return oss.str();
}
#undef _W2
#undef _W3

void GLViewer::closeEvent(QCloseEvent *e){
	//emit closeSignal(viewId);
	YadeQtMainWindow::self->closeView(this);
}

void GLViewer::mouseMoveEvent(QMouseEvent *e){
	last_user_event = boost::posix_time::second_clock::local_time();
	QGLViewer::mouseMoveEvent(e);
}

void GLViewer::mousePressEvent(QMouseEvent *e){
	last_user_event = boost::posix_time::second_clock::local_time();
	QGLViewer::mousePressEvent(e);
}

/* Handle double-click event; if clipping plane is manipulated, align it with the global coordinate system.
 * Otherwise pass the event to QGLViewer to handle it normally.
 *
 * mostly copied over from ManipulatedFrame::mouseDoubleClickEvent
 */
void GLViewer::mouseDoubleClickEvent(QMouseEvent *event){
	last_user_event = boost::posix_time::second_clock::local_time();

	if(manipulatedClipPlane<0) { /* LOG_DEBUG("Double click not on clipping plane"); */ QGLViewer::mouseDoubleClickEvent(event); return; }
#if QT_VERSION >= 0x040000
	if (event->modifiers() == Qt::NoModifier)
#else
	if (event->state() == Qt::NoButton)
#endif
	switch (event->button()){
		case Qt::LeftButton:  manipulatedFrame()->alignWithFrame(NULL,true); LOG_DEBUG("Aligning cutting plane"); break;
		// case Qt::RightButton: projectOnLine(camera->position(), camera->viewDirection()); break;
		default: break; // avoid warning
	}
}

void GLViewer::wheelEvent(QWheelEvent* event){
	last_user_event = boost::posix_time::second_clock::local_time();

	if(manipulatedClipPlane<0){ QGLViewer::wheelEvent(event); return; }
	assert(manipulatedClipPlane<renderer->numClipPlanes);
	float distStep=1e-3*sceneRadius();
	//const float wheelSensitivityCoef = 8E-4f;
	//Vec trans(0.0, 0.0, -event->delta()*wheelSensitivity()*wheelSensitivityCoef*(camera->position()-position()).norm());
	float dist=event->delta()*manipulatedFrame()->wheelSensitivity()*distStep;
	Vector3r normal=renderer->clipPlaneSe3[manipulatedClipPlane].orientation*Vector3r(0,0,1);
	qglviewer::Vec newPos=manipulatedFrame()->position()+qglviewer::Vec(normal[0],normal[1],normal[2])*dist;
	manipulatedFrame()->setPosition(newPos);
	renderer->clipPlaneSe3[manipulatedClipPlane].position=Vector3r(newPos[0],newPos[1],newPos[2]);
	updateGL();
	/* in draw, bound cutting planes will be moved as well */
}

// cut&paste from QGLViewer::domElement documentation
QDomElement GLViewer::domElement(const QString& name, QDomDocument& document) const{
	QDomElement de=document.createElement("gridXYZ");
	string val; if(drawGridXYZ[0])val+="x"; if(drawGridXYZ[1])val+="y"; if(drawGridXYZ[2])val+="z";
	de.setAttribute("normals",val);
	QDomElement de2=document.createElement("timeDisplay"); de2.setAttribute("mask",timeDispMask);
	QDomElement res=QGLViewer::domElement(name,document);
	res.appendChild(de);
	res.appendChild(de2);
	return res;
}

// cut&paste from QGLViewer::initFromDomElement documentation
void GLViewer::initFromDOMElement(const QDomElement& element){
	QGLViewer::initFromDOMElement(element);
	QDomElement child=element.firstChild().toElement();
	while (!child.isNull()){
		if (child.tagName()=="gridXYZ" && child.hasAttribute("normals")){
			string val=child.attribute("normals").lower();
			drawGridXYZ[0]=false; drawGridXYZ[1]=false; drawGridXYZ[2]=false;
			if(val.find("x")!=string::npos)drawGridXYZ[0]=true; if(val.find("y")!=string::npos)drawGridXYZ[1]=true; if(val.find("z")!=string::npos)drawGridXYZ[2]=true;
		}
		if(child.tagName()=="timeDisplay" && child.hasAttribute("mask")) timeDispMask=atoi(child.attribute("mask").ascii());
		child = child.nextSibling().toElement();
	}
}

boost::posix_time::ptime GLViewer::getLastUserEvent(){return last_user_event;};


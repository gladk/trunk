// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
#pragma once

#include<yade/lib-multimethods/DynLibDispatcher.hpp>
#include<yade/core/Dispatcher.hpp>
#include<yade/core/Body.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>

#include<yade/pkg-common/GLDrawFunctors.hpp>

struct GlExtraDrawer: public Serializable{
	Scene* scene;
	virtual void render();
	YADE_CLASS_BASE_DOC_ATTRS(GlExtraDrawer,Serializable,"Performing arbitrary OpenGL drawing commands; called from :yref:`OpenGLRenderer` (see :yref:`OpenGLRenderer.extraDrawers`) once regular rendering routines will have finished.\n\nThis class itself does not render anything, derived classes should override the *render* method.",
		((bool,dead,false,,"Deactivate the object (on error/exception)."))
	);
};
REGISTER_SERIALIZABLE(GlExtraDrawer);

class OpenGLRenderer : public Serializable
{
	public:
		int _nothing; // remove later, only target for deprecated selectBodyLimit

		static const int numClipPlanes=3;

		bool pointClipped(const Vector3r& p);
		vector<Vector3r> clipPlaneNormals;
		void setBodiesDispInfo();
		static bool initDone;
		Vector3r viewDirection; // updated from GLViewer regularly
		GLViewInfo viewInfo; // update from GLView regularly
		Vector3r highlightEmission0;
		Vector3r highlightEmission1;

		// normalized saw signal with given periodicity, with values ∈ 〈0,1〉 */
		Real normSaw(Real t, Real period){ Real xi=(t-period*((int)(t/period)))/period; /* normalized value, (0-1〉 */ return (xi<.5?2*xi:2-2*xi); }
		Real normSquare(Real t, Real period){ Real xi=(t-period*((int)(t/period)))/period; /* normalized value, (0-1〉 */ return (xi<.5?0:1); }

		void drawPeriodicCell();

		void setBodiesRefSe3();

		struct BodyDisp{
			Vector3r pos;
			Quaternionr ori;
			bool isDisplayed;
		};
		//! display data for individual bodies
		vector<BodyDisp> bodyDisp;

		virtual ~OpenGLRenderer();

	private:
		void resetSpecularEmission();

		GlBoundDispatcher boundDispatcher;
		GlIGeomDispatcher geomDispatcher;
		GlIPhysDispatcher physDispatcher;
		GlShapeDispatcher shapeDispatcher;
		// GlStateDispatcher stateDispatcher;


		vector<string>
			// stateFunctorNames,
			boundFunctorNames,
			shapeFunctorNames, 
			geomFunctorNames,
			physFunctorNames;

		DECLARE_LOGGER;

	public :
		// updated after every call to render
		shared_ptr<Scene> scene; 

		void init();
		void initgl();
		void render(const shared_ptr<Scene>& scene, Body::id_t selection=Body::id_t(-1));
		void pyRender(){render(Omega::instance().getScene());}
	
		void renderDOF_ID();
		void renderIPhys();
		void renderIGeom();
		void renderBound();
		// called also to render selectable entitites;
		void renderShape();
		void renderAllInteractionsWire();

	YADE_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(OpenGLRenderer,Serializable,"Class responsible for rendering scene on OpenGL devices.",
		((Vector3r,dispScale,((void)"disable scaling",Vector3r::Ones()),,"Artificially enlarge (scale) dispalcements from bodies' :yref:`reference positions<State.refPos>` by this relative amount, so that they become better visible (independently in 3 dimensions). Disbled if (1,1,1)."))
		((Real,rotScale,((void)"disable scaling",1.),,"Artificially enlarge (scale) rotations of bodies relative to their :yref:`reference orientation<State.refOri>`, so the they are better visible."))
		((Vector3r,lightPos,Vector3r(75,130,0),,"Position of OpenGL light source in the scene."))
		((Vector3r,bgColor,Vector3r(.2,.2,.2),,"Color of the backgroud canvas (RGB)"))
		((bool,wire,false,,"Render all bodies with wire only (faster)"))
		((bool,dof,false,,"Show which degrees of freedom are blocked for each body"))
		((bool,id,false,,"Show body id's"))
		((bool,bound,false,,"Render body :yref:`Bound`"))
		((bool,shape,true,,"Render body :yref:`Shape`"))
		((bool,intrWire,false,,"If rendering interactions, use only wires to represent them."))
		((bool,intrGeom,false,,"Render :yref:`Interaction::geom` objects."))
		((bool,intrPhys,false,,"Render :yref:`Interaction::phys` objects"))
		((int,mask,((void)"draw everything",~0),,"Bitmask for showing only bodies where ((mask & :yref:`Body::mask`)!=0)"))
		((Body::id_t,selId,Body::ID_NONE,,"Id of particle that was selected by the user."))
		((vector<Se3r>,clipPlaneSe3,vector<Se3r>(numClipPlanes,Se3r(Vector3r::Zero(),Quaternionr::Identity())),,"Position and orientation of clipping planes"))
		((vector<bool>,clipPlaneActive,vector<bool>(numClipPlanes,false),,"Activate/deactivate respective clipping planes"))
		((vector<shared_ptr<GlExtraDrawer> >,extraDrawers,,,"Additional rendering components (:yref:`GlExtraDrawer`)."))
		((bool,intrAllWire,false,,"Draw wire for all interactions, blue for potential and green for real ones (mostly for debugging)")),
		/*deprec*/
			((Light_position,lightPos,))
			((Background_color,bgColor,))
			((Body_wire,wire,))
			((Show_DOF,dof,))
			((Show_ID,id,))
			((Body_bounding_volume,bound,))
			((Body_interacting_geom,shape,))
			((Interaction_wire,intrWire,))
			((Interaction_geometry,intrGeom,))
			((Interaction_physics,intrPhys,))
			((Draw_mask,mask,))
			((displacementScale,dispScale,))
			((rotationScale,rotScale,))
			((scaleDisplacements,id,"! This option is no longer necessary, set dispScale to value different from (1,1,1) to have scaling applied."))
			((scaleRotations,id,"! This option is no longer necessary, set rotScale to value different from 1.0 to have scaling applied."))
			((selectBodyLimit,_nothing,"There is no longer a limit on selecting particles; the routine is called only with shift-click and has no performance impact otherwise."))
		,
		/*init*/,
		/*ctor*/,
		/*py*/
		.def("setRefSe3",&OpenGLRenderer::setBodiesRefSe3,"Make current positions and orientation reference for scaleDisplacements and scaleRotations.")
		.def("render",&OpenGLRenderer::pyRender,"Render the scene in the current OpenGL context.")
	);
};
REGISTER_SERIALIZABLE(OpenGLRenderer);




/*************************************************************************
*  Copyright (C) 2012 by François Kneib   francois.kneib@gmail.com       *
*  Copyright (C) 2012 by Bruno Chareyre   bruno.chareyre@hmg.inpg.fr     *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/* TABLE OF CONTENT, and the minimum you need to understand.
- 2 new shapes for grids : GridNode (vertices) and GridConnection (edges)
- 2 new contact geometries :
    * GridNodeGeom6D to handle GridNode-GridNode contacts (the internal behaviour of the grid)
    * ScGridCoGeom to handle Sphere-GridConnection contacts (the interaction between the grid and an external sphere)
        Note : the Sphere-Grid contacts are always handled by the GridConnections, but the forces are applied on the related GridNodes.
        Note : there is no contact between two GridConnections, they must be linked with GridNodes.
- The 2 related Ig2 :
    * Ig2_GridNode_GridNode_GridNodeGeom6D (doing almost the same than Ig2_Sphere_Sphere_ScGeom6D)
    * Ig2_Sphere_GridConnection_ScGridCoGeom (the biggest part of the code, it handles contact detection and history when a sphere is sliding on the grid over consecutive GridConnections)
- The Law2_ScGridCoGeom_FrictPhys_CundallStrack who handles the elastic frictional Sphere-GridConnection contact. The GridNode-GridNode law is Law2_ScGeom6D_CohFrictPhys_CohesionMoment by inheritance.
*/



#pragma once
#include "Sphere.hpp"
#include<yade/pkg/dem/FrictPhys.hpp>
#include <yade/pkg/dem/ScGeom.hpp>
#include <yade/core/Body.hpp>
#include<yade/pkg/common/Dispatching.hpp>
#include<yade/pkg/dem/Ig2_Sphere_Sphere_ScGeom.hpp>
#ifdef YADE_OPENGL
	#include<yade/pkg/common/GLDrawFunctors.hpp>
#endif



//!##################	SHAPES   #####################

class GridConnection: public Sphere{
	public:
		virtual ~GridConnection();
		Real getLength();
		Vector3r getSegment();
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(GridConnection,Sphere,"GridConnection shape. Component of a grid designed to link two :yref:`GridNodes<GridNode>`. It's highly recommanded to use utils.gridConnection(...) to generate correct :yref:`GridConnections<GridConnection>`.",
		((shared_ptr<Body> , node1 , ,,"First :yref:`Body` the GridConnection is connected to."))
		((shared_ptr<Body> , node2 , ,,"Second :yref:`Body` the GridConnection is connected to."))
		((bool, periodic, false,,"true if two nodes from different periods are connected."))
		((Vector3i , cellDist , Vector3i(0,0,0),,"missing doc :(")),
		createIndex(); /*ctor*/
	);
	REGISTER_CLASS_INDEX(GridConnection,Sphere);
};
REGISTER_SERIALIZABLE(GridConnection);


class GridNode: public Sphere{
	public:
		virtual ~GridNode();
		void addConnection(shared_ptr<Body> GC);
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(GridNode,Sphere,"GridNode shape, component of a grid.\nTo create a Grid, place the nodes first, they will define the spacial discretisation of it. It's highly recommanded to use utils.gridNode(...) to generate correct :yref:`GridNodes<GridNode>`. Note that the GridNodes should only be in an Interaction with other GridNodes. The Sphere-Grid contact is only handled by the :yref:`GridConnections<GridConnection>`.",
		((vector<shared_ptr<Body> >,ConnList,,,"List of :yref:`GridConnections<GridConnection>` the GridNode is connected to.")),
		/*ctor*/
		createIndex();,
		/*py*/
		.def("addConnection",&GridNode::addConnection,(python::arg("Body")),"Add a GridConnection to the GridNode.")
	);
	REGISTER_CLASS_INDEX(GridNode,Sphere);
};
REGISTER_SERIALIZABLE(GridNode);


//!##################	Contact Geometry   #####################

//!			O-O
class GridNodeGeom6D: public ScGeom6D {
	public:
		virtual ~GridNodeGeom6D();
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(GridNodeGeom6D,ScGeom6D,"Geometry of a :yref:`GridNode`-:yref:`GridNode` contact. Inherits almost everything from :yref:`ScGeom6D`.",
		((shared_ptr<Body>, connectionBody,,,"Reference to the :yref:`GridNode` :yref:`Body` who is linking the two :yref:`GridNodes<GridNode>`.")),
		/* extra initializers */,
		/* ctor */ createIndex();,
		/* py */
	);
	REGISTER_CLASS_INDEX(GridNodeGeom6D,ScGeom6D);
};
REGISTER_SERIALIZABLE(GridNodeGeom6D);

//!			O/
class ScGridCoGeom: public ScGeom {
	public:
		/// Emulate a sphere whose position is the projection of sphere's center on cylinder sphere, and with motion linearly interpolated between nodes
		State fictiousState;
		virtual ~ScGridCoGeom ();
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(ScGridCoGeom,ScGeom,"Geometry of a :yref:`GridConnection`-:yref:`Sphere` contact.",
		((int,isDuplicate,0,,"this flag is turned true (1) automatically if the contact is shared between two Connections. A duplicated interaction will be skipped once by the constitutive law, so that only one contact at a time is effective. If isDuplicate=2, it means one of the two duplicates has no longer geometric interaction, and should be erased by the constitutive laws."))
		((int,trueInt,-1,,"Defines the body id of the :yref:`GridConnection` where the contact is real, when :yref:`ScGridCoGeom::isDuplicate`>0."))
		((int,id3,0,,"id of the first :yref:`GridNode`. |yupdate|"))
		((int,id4,0,,"id of the second :yref:`GridNode`. |yupdate|"))
		((Real,relPos,0,,"position of the contact on the connection (0: node-, 1:node+) |yupdate|")),
		createIndex(); /*ctor*/
	);
	REGISTER_CLASS_INDEX(ScGridCoGeom,ScGeom);
};
REGISTER_SERIALIZABLE(ScGridCoGeom);

//!##################	IGeom Functors   #####################

//!			O-O
class Ig2_GridNode_GridNode_GridNodeGeom6D: public Ig2_Sphere_Sphere_ScGeom{
	public:
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_GridNode_GridNode_GridNodeGeom6D,Ig2_Sphere_Sphere_ScGeom,"Create/update a :yref:`GridNodeGeom6D` instance representing the geometry of a contact point between two :yref:`GridNode<GridNode>`, including relative rotations.",
		((bool,updateRotations,true,,"Precompute relative rotations. Turning this false can speed up simulations when rotations are not needed in constitutive laws (e.g. when spheres are compressed without cohesion and moment in early stage of a triaxial test), but is not foolproof. Change this value only if you know what you are doing."))
		((bool,creep,false,,"Substract rotational creep from relative rotation. The rotational creep :yref:`ScGeom6D::twistCreep` is a quaternion and has to be updated inside a constitutive law, see for instance :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`."))
	);
	FUNCTOR2D(GridNode,GridNode);
	// needed for the dispatcher, even if it is symmetric
	DEFINE_FUNCTOR_ORDER_2D(GridNode,GridNode);
};
REGISTER_SERIALIZABLE(Ig2_GridNode_GridNode_GridNodeGeom6D);

//!			O/
class Ig2_Sphere_GridConnection_ScGridCoGeom: public IGeomFunctor{
	public:
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_Sphere_GridConnection_ScGridCoGeom,IGeomFunctor,"Create/update a :yref:`ScGridCoGeom6D` instance representing the geometry of a contact point between a :yref:`GricConnection` and a :yref:`Sphere` including relative rotations.",
		((Real,interactionDetectionFactor,1,,"Enlarge both radii by this factor (if >1), to permit creation of distant interactions."))
	);
	FUNCTOR2D(Sphere,GridConnection);
	DEFINE_FUNCTOR_ORDER_2D(Sphere,GridConnection);
};
REGISTER_SERIALIZABLE(Ig2_Sphere_GridConnection_ScGridCoGeom);


//!##################	Laws   #####################

//!			O/
class Law2_ScGridCoGeom_FrictPhys_CundallStrack: public LawFunctor{
	public:
		virtual void go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I);
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGridCoGeom_FrictPhys_CundallStrack,LawFunctor,"Law between a frictional :yref:`GridConnection` and a frictional :yref:`Sphere`. Almost the same than :yref:`Law2_ScGeom_FrictPhys_CundallStrack`, but the force is divided and applied on the two :yref:`GridNodes<GridNode>` only.",
		((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
		((bool,traceEnergy,false,Attr::hidden,"Define the total energy dissipated in plastic slips at all contacts."))
		((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
		((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
		,,
	);
	FUNCTOR2D(ScGridCoGeom,FrictPhys);
};
REGISTER_SERIALIZABLE(Law2_ScGridCoGeom_FrictPhys_CundallStrack);



//!##################	Bounds   #####################

class Bo1_GridConnection_Aabb : public BoundFunctor
{
	public :
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r&, const Body*);
	FUNCTOR1D(GridConnection);
	YADE_CLASS_BASE_DOC_ATTRS(Bo1_GridConnection_Aabb,BoundFunctor,"Functor creating :yref:`Aabb` from a :yref:`GridConnection`.",
		((Real,aabbEnlargeFactor,((void)"deactivated",-1),,"Relative enlargement of the bounding box; deactivated if negative."))
	);
};
REGISTER_SERIALIZABLE(Bo1_GridConnection_Aabb);

//!##################	Rendering   #####################
#ifdef YADE_OPENGL
class Gl1_GridConnection : public GlShapeFunctor{
	private:
		//static int glCylinderList;
		//void subdivideTriangle(Vector3r& v1,Vector3r& v2,Vector3r& v3, int depth);
		void drawCylinder(bool wire, Real radius, Real length, const Quaternionr& shift=Quaternionr::Identity());
		//void initGlLists(void);
	public:
		virtual void go(const shared_ptr<Shape>&, const shared_ptr<State>&,bool,const GLViewInfo&);
		void out( Quaternionr q );
	YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_GridConnection,GlShapeFunctor,"Renders :yref:`Cylinder` object",
		((bool,wire,false,,"Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
		((bool,glutNormalize,true,,"Fix normals for non-wire rendering"))
		((int,glutSlices,8,,"Number of sphere slices."))
		((int,glutStacks,4,,"Number of sphere stacks."))
	);
	RENDERS(GridConnection);
};
REGISTER_SERIALIZABLE(Gl1_GridConnection);
#endif










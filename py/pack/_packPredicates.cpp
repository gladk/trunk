// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#include<boost/python.hpp>
#include<yade/extra/boost_python_len.hpp>
#include<yade/lib-base/Logging.hpp>
#include<yade/lib-base/Math.hpp>
#include<yade/lib-pyutil/doc_opts.hpp>

using namespace boost;
using namespace std;
#ifdef YADE_LOG4CXX
	log4cxx::LoggerPtr logger=log4cxx::Logger::getLogger("yade.pack.predicates");
#endif

/*
This file contains various predicates that say whether a given point is within the solid,
or, not closer than "pad" to its boundary, if pad is nonzero
Besides the (point,pad) operator, each predicate defines aabb() method that returns
(min,max) tuple defining minimum and maximum point of axis-aligned bounding box 
for the predicate.

These classes are primarily used for yade.pack.* functions creating packings.
See examples/regular-sphere-pack/regular-sphere-pack.py for an example.

*/

// aux functions
python::tuple vec2tuple(const Vector3r& v){return boost::python::make_tuple(v[0],v[1],v[2]);}
python::tuple vec2tuple(const Vector2r& v){return boost::python::make_tuple(v[0],v[1]);}
Vector3r tuple2vec(const python::tuple& t){return Vector3r(python::extract<double>(t[0])(),python::extract<double>(t[1])(),python::extract<double>(t[2])());}
Vector2r tuple2vec2d(const python::tuple& t){return Vector2r(python::extract<double>(t[0])(),python::extract<double>(t[1])());}
//void ttuple2vvec(const python::tuple& t, Vector3r& v1, Vector3r& v2){ v1=tuple2vec(python::extract<python::tuple>(t[0])()); v2=tuple2vec(python::extract<python::tuple>(t[1])()); }
void ttuple2vvec(const python::tuple& t, Vector3r& v1, Vector3r& v2){ v1=python::extract<Vector3r>(t[0])(); v2=python::extract<Vector3r>(t[1])(); }
python::tuple vvec2ttuple(const Vector3r&v1, const Vector3r&v2){ return python::make_tuple(v1,v2); }

struct Predicate{
	public:
		virtual bool operator() (const Vector3r& pt,Real pad=0.) const = 0;
		virtual python::tuple aabb() const = 0;
		Vector3r dim() const { Vector3r mn,mx; ttuple2vvec(aabb(),mn,mx); return mx-mn; }
		Vector3r center() const { Vector3r mn,mx; ttuple2vvec(aabb(),mn,mx); return .5*(mn+mx); }
};
// make the pad parameter optional
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PredicateCall_overloads,operator(),1,2);

/* Since we want to make Predicate::operator() and Predicate::aabb() callable from c++ on python::object
with the right virtual method resolution, we have to wrap the class in the following way. See 
http://www.boost.org/doc/libs/1_38_0/libs/python/doc/tutorial/doc/html/python/exposing.html for documentation
on exposing virtual methods.

This makes it possible to derive a python class from Predicate, override its aabb() method, for instance,
and use it in PredicateUnion, which will call the python implementation of aabb() as it should. This
approach is used in the inGtsSurface class defined in pack.py.

See scripts/test/gts-operators.py for an example.

NOTE: you still have to call base class ctor in your class' ctor derived in python, e.g.
super(inGtsSurface,self).__init__() so that virtual methods work as expected.
*/
struct PredicateWrap: Predicate, python::wrapper<Predicate>{
	bool operator()(const Vector3r& pt, Real pad=0.) const { return this->get_override("__call__")(pt,pad);}
	python::tuple aabb() const { return this->get_override("aabb")(); }
};
// make the pad parameter optional
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PredicateWrapCall_overloads,operator(),1,2);

/*********************************************************************************
****************** Boolean operations on predicates ******************************
*********************************************************************************/

const Predicate& obj2pred(python::object obj){ return python::extract<const Predicate&>(obj)();}

class PredicateBoolean: public Predicate{
	protected:
		const python::object A,B;
	public:
		PredicateBoolean(const python::object _A, const python::object _B): A(_A), B(_B){}
		const python::object getA(){ return A;}
		const python::object getB(){ return B;}
};

// http://www.linuxtopia.org/online_books/programming_books/python_programming/python_ch16s03.html
class PredicateUnion: public PredicateBoolean{
	public:
		PredicateUnion(const python::object _A, const python::object _B): PredicateBoolean(_A,_B){}
		virtual bool operator()(const Vector3r& pt,Real pad) const {return obj2pred(A)(pt,pad)||obj2pred(B)(pt,pad);}
		virtual python::tuple aabb() const { Vector3r minA,maxA,minB,maxB; ttuple2vvec(obj2pred(A).aabb(),minA,maxA); ttuple2vvec(obj2pred(B).aabb(),minB,maxB); return vvec2ttuple(componentMinVector(minA,minB),componentMaxVector(maxA,maxB));}
};
PredicateUnion makeUnion(const python::object& A, const python::object& B){ return PredicateUnion(A,B);}

class PredicateIntersection: public PredicateBoolean{
	public:
		PredicateIntersection(const python::object _A, const python::object _B): PredicateBoolean(_A,_B){}
		virtual bool operator()(const Vector3r& pt,Real pad) const {return obj2pred(A)(pt,pad) && obj2pred(B)(pt,pad);}
		virtual python::tuple aabb() const { Vector3r minA,maxA,minB,maxB; ttuple2vvec(obj2pred(A).aabb(),minA,maxA); ttuple2vvec(obj2pred(B).aabb(),minB,maxB); return vvec2ttuple(componentMaxVector(minA,minB),componentMinVector(maxA,maxB));}
};
PredicateIntersection makeIntersection(const python::object& A, const python::object& B){ return PredicateIntersection(A,B);}

class PredicateDifference: public PredicateBoolean{
	public:
		PredicateDifference(const python::object _A, const python::object _B): PredicateBoolean(_A,_B){}
		virtual bool operator()(const Vector3r& pt,Real pad) const {return obj2pred(A)(pt,pad) && !obj2pred(B)(pt,-pad);}
		virtual python::tuple aabb() const { return obj2pred(A).aabb(); }
};
PredicateDifference makeDifference(const python::object& A, const python::object& B){ return PredicateDifference(A,B);}

class PredicateSymmetricDifference: public PredicateBoolean{
	public:
		PredicateSymmetricDifference(const python::object _A, const python::object _B): PredicateBoolean(_A,_B){}
		virtual bool operator()(const Vector3r& pt,Real pad) const {bool inA=obj2pred(A)(pt,pad), inB=obj2pred(B)(pt,pad); return (inA && !inB) || (!inA && inB);}
		virtual python::tuple aabb() const { Vector3r minA,maxA,minB,maxB; ttuple2vvec(obj2pred(A).aabb(),minA,maxA); ttuple2vvec(obj2pred(B).aabb(),minB,maxB); return vvec2ttuple(componentMinVector(minA,minB),componentMaxVector(maxA,maxB));}
};
PredicateSymmetricDifference makeSymmetricDifference(const python::object& A, const python::object& B){ return PredicateSymmetricDifference(A,B);}

/*********************************************************************************
****************************** Primitive predicates ******************************
*********************************************************************************/


/*! Sphere predicate */
class inSphere: public Predicate {
	Vector3r center; Real radius;
public:
	inSphere(const Vector3r& _center, Real _radius){center=_center; radius=_radius;}
	virtual bool operator()(const Vector3r& pt, Real pad=0.) const { return ((pt-center).norm()-pad<=radius-pad); }
	virtual python::tuple aabb() const {return vvec2ttuple(Vector3r(center[0]-radius,center[1]-radius,center[2]-radius),Vector3r(center[0]+radius,center[1]+radius,center[2]+radius));}
};

/*! Axis-aligned box predicate */
class inAlignedBox: public Predicate{
	Vector3r mn, mx;
public:
	inAlignedBox(const Vector3r& _mn, const Vector3r& _mx): mn(_mn), mx(_mx) {}
	virtual bool operator()(const Vector3r& pt, Real pad=0.) const {
		return
			mn[0]+pad<=pt[0] && mx[0]-pad>=pt[0] &&
			mn[1]+pad<=pt[1] && mx[1]-pad>=pt[1] &&
			mn[2]+pad<=pt[2] && mx[2]-pad>=pt[2];
	}
	virtual python::tuple aabb() const { return vvec2ttuple(mn,mx); }
};

/*! Arbitrarily oriented cylinder predicate */
class inCylinder: public Predicate{
	Vector3r c1,c2,c12; Real radius,ht;
public:
	inCylinder(const Vector3r& _c1, const Vector3r& _c2, Real _radius){c1=_c1; c2=_c2; c12=c2-c1; radius=_radius; ht=c12.norm(); }
	bool operator()(const Vector3r& pt, Real pad=0.) const {
		Real u=(pt.dot(c12)-c1.dot(c12))/(ht*ht); // normalized coordinate along the c1--c2 axis
		if((u*ht<0+pad) || (u*ht>ht-pad)) return false; // out of cylinder along the axis
		Real axisDist=((pt-c1).cross(pt-c2)).norm()/ht;
		if(axisDist>radius-pad) return false;
		return true;
	}
	python::tuple aabb() const {
		// see http://www.gamedev.net/community/forums/topic.asp?topic_id=338522&forum_id=20&gforum_id=0 for the algorithm
		const Vector3r& A(c1); const Vector3r& B(c2); 
		Vector3r k(
			sqrt((pow(A[1]-B[1],2)+pow(A[2]-B[2],2)))/ht,
			sqrt((pow(A[0]-B[0],2)+pow(A[2]-B[2],2)))/ht,
			sqrt((pow(A[0]-B[0],2)+pow(A[1]-B[1],2)))/ht);
		Vector3r mn=componentMinVector(A,B), mx=componentMaxVector(A,B);
		return vvec2ttuple(mn-radius*k,mx+radius*k);
	}
};

/*! Oriented hyperboloid predicate (cylinder as special case).

See http://mathworld.wolfram.com/Hyperboloid.html for the parametrization and meaning of symbols
*/
class inHyperboloid: public Predicate{
	Vector3r c1,c2,c12; Real R,a,ht,c;
public:
	inHyperboloid(const Vector3r& _c1, const Vector3r& _c2, Real _R, Real _r){
		c1=_c1; c2=_c2; R=_R; a=_r;
		c12=c2-c1; ht=c12.norm();
		Real uMax=sqrt(pow(R/a,2)-1); c=ht/(2*uMax);
	}
	// WARN: this is not accurate, since padding is taken as perpendicular to the axis, not the the surface
	bool operator()(const Vector3r& pt, Real pad=0.) const {
		Real v=(pt.dot(c12)-c1.dot(c12))/(ht*ht); // normalized coordinate along the c1--c2 axis
		if((v*ht<0+pad) || (v*ht>ht-pad)) return false; // out of cylinder along the axis
		Real u=(v-.5)*ht/c; // u from the wolfram parametrization; u is 0 in the center
		Real rHere=a*sqrt(1+u*u); // pad is taken perpendicular to the axis, not to the surface (inaccurate)
		Real axisDist=((pt-c1).cross(pt-c2)).norm()/ht;
		if(axisDist>rHere-pad) return false;
		return true;
	}
	python::tuple aabb() const {
		// the lazy way
		return inCylinder(c1,c2,R).aabb();
	}
};

/*! Axis-aligned ellipsoid predicate */
class inEllipsoid: public Predicate{
	Vector3r c, abc;
public:
	inEllipsoid(const Vector3r& _c, const Vector3r& _abc) {c=_c; abc=_abc;}
	bool operator()(const Vector3r& pt, Real pad=0.) const {
		//Define the ellipsoid X-coordinate of given Y and Z
		Real x = sqrt((1-pow((pt[1]-c[1]),2)/((abc[1]-pad)*(abc[1]-pad))-pow((pt[2]-c[2]),2)/((abc[2]-pad)*(abc[2]-pad)))*((abc[0]-pad)*(abc[0]-pad)))+c[0]; 
		Vector3r edgeEllipsoid(x,pt[1],pt[2]); // create a vector of these 3 coordinates
		//check whether given coordinates lie inside ellipsoid or not
		if ((pt-c).norm()<=(edgeEllipsoid-c).norm()) return true;
		else return false;
	}
	python::tuple aabb() const {
		const Vector3r& center(c); const Vector3r& ABC(abc);
		return vvec2ttuple(Vector3r(center[0]-ABC[0],center[1]-ABC[1],center[2]-ABC[2]),Vector3r(center[0]+ABC[0],center[1]+ABC[1],center[2]+ABC[2]));
	}
};

/*! Negative notch predicate.

Use intersection (& operator) of another predicate with notInNotch to create notched solid.


		
		geometry explanation:
		
			c: the center
			normalHalfHt (in constructor): A-C
			inside: perpendicular to notch edge, points inside the notch (unit vector)
			normal: perpendicular to inside, perpendicular to both notch planes
			edge: unit vector in the direction of the edge

		          ↑ distUp        A
		-------------------------
		                        | C
		         inside(unit) ← * → distInPlane
		                        |
		-------------------------
		          ↓ distDown      B

*/
class notInNotch: public Predicate{
	Vector3r c, edge, normal, inside; Real aperture;
public:
	notInNotch(const Vector3r& _c, const Vector3r& _edge, const Vector3r& _normal, Real _aperture){
		c=_c;
		edge=_edge; edge.normalize();
		normal=_normal; normal-=edge*edge.dot(normal); normal.normalize();
		inside=edge.cross(normal);
		aperture=_aperture;
		// LOG_DEBUG("edge="<<edge<<", normal="<<normal<<", inside="<<inside<<", aperture="<<aperture);
	}
	bool operator()(const Vector3r& pt, Real pad=0.) const {
		Real distUp=normal.dot(pt-c)-aperture/2, distDown=-normal.dot(pt-c)-aperture/2, distInPlane=-inside.dot(pt-c);
		// LOG_DEBUG("pt="<<pt<<", distUp="<<distUp<<", distDown="<<distDown<<", distInPlane="<<distInPlane);
		if(distInPlane>=pad) return true;
		if(distUp     >=pad) return true;
		if(distDown   >=pad) return true;
		if(distInPlane<0) return false;
		if(distUp  >0) return sqrt(pow(distInPlane,2)+pow(distUp,2))>=pad;
		if(distDown>0) return sqrt(pow(distInPlane,2)+pow(distUp,2))>=pad;
		// between both notch planes, closer to the edge than pad (distInPlane<pad)
		return false;
	}
	// This predicate is not bounded, return infinities
	python::tuple aabb() const {
		Real inf=std::numeric_limits<Real>::infinity();
		return vvec2ttuple(Vector3r(-inf,-inf,-inf),Vector3r(inf,inf,inf));
	}
};

#ifdef YADE_GTS
extern "C" {
// HACK
#include"../3rd-party/pygts-0.3.1/pygts.h"

}
/* Helper function for inGtsSurface::aabb() */
static void vertex_aabb(GtsVertex *vertex, pair<Vector3r,Vector3r> *bb)
{
	GtsPoint *_p=GTS_POINT(vertex);
	Vector3r p(_p->x,_p->y,_p->z);
	bb->first=componentMinVector(bb->first,p);
	bb->second=componentMaxVector(bb->second,p);
}

/*
This class plays tricks getting around pyGTS to get GTS objects and cache bb tree to speed
up point inclusion tests. For this reason, we have to link with _gts.so (see corresponding
SConscript file), which is at the same time the python module.
*/
class inGtsSurface: public Predicate{
	python::object pySurf; // to hold the reference so that surf is valid
	GtsSurface *surf;
	bool is_open, noPad, noPadWarned;
	GNode* tree;
public:
	inGtsSurface(python::object _surf, bool _noPad=false): pySurf(_surf), noPad(_noPad), noPadWarned(false) {
		if(!pygts_surface_check(_surf.ptr())) throw invalid_argument("Ctor must receive a gts.Surface() instance."); 
		surf=PYGTS_SURFACE_AS_GTS_SURFACE(PYGTS_SURFACE(_surf.ptr()));
	 	if(!gts_surface_is_closed(surf)) throw invalid_argument("Surface is not closed.");
		is_open=gts_surface_volume(surf)<0.;
		if((tree=gts_bb_tree_surface(surf))==NULL) throw runtime_error("Could not create GTree.");
	}
	~inGtsSurface(){g_node_destroy(tree);}
	python::tuple aabb() const {
		Real inf=std::numeric_limits<Real>::infinity();
		pair<Vector3r,Vector3r> bb; bb.first=Vector3r(inf,inf,inf); bb.second=Vector3r(-inf,-inf,-inf);
		gts_surface_foreach_vertex(surf,(GtsFunc)vertex_aabb,&bb);
		return vvec2ttuple(bb.first,bb.second);
	}
	bool ptCheck(const Vector3r& pt) const{
		GtsPoint gp; gp.x=pt[0]; gp.y=pt[1]; gp.z=pt[2];
		return (bool)gts_point_is_inside_surface(&gp,tree,is_open);
	}
	bool operator()(const Vector3r& pt, Real pad=0.) const {
		if(noPad){
			if(pad!=0. && noPadWarned) LOG_WARN("inGtsSurface constructed with noPad; requested non-zero pad set to zero.");
			return ptCheck(pt);
		}
		return ptCheck(pt) && ptCheck(pt-Vector3r(pad,0,0)) && ptCheck(pt+Vector3r(pad,0,0)) && ptCheck(pt-Vector3r(0,pad,0))&& ptCheck(pt+Vector3r(0,pad,0)) && ptCheck(pt-Vector3r(0,0,pad))&& ptCheck(pt+Vector3r(0,0,pad));
	}
	python::object surface() const {return pySurf;}
};

#endif

BOOST_PYTHON_MODULE(_packPredicates){
	python::scope().attr("__doc__")="Spatial predicates for volumes (defined analytically or by triangulation).";
	YADE_SET_DOCSTRING_OPTS;
	// base predicate class
	python::class_<PredicateWrap,/* necessary, as methods are pure virtual*/ boost::noncopyable>("Predicate")
		.def("__call__",python::pure_virtual(&Predicate::operator()))
		.def("aabb",python::pure_virtual(&Predicate::aabb))
		.def("dim",&Predicate::dim)
		.def("center",&Predicate::center)
		.def("__or__",makeUnion).def("__and__",makeIntersection).def("__sub__",makeDifference).def("__xor__",makeSymmetricDifference);
	// boolean operations
	python::class_<PredicateBoolean,python::bases<Predicate>,boost::noncopyable>("PredicateBoolean","Boolean operation on 2 predicates (abstract class)",python::no_init)
		.add_property("A",&PredicateBoolean::getA).add_property("B",&PredicateBoolean::getB);
	python::class_<PredicateUnion,python::bases<PredicateBoolean> >("PredicateUnion","Union (non-exclusive disjunction) of 2 predicates. A point has to be inside any of the two predicates to be inside. Can be constructed using the ``|`` operator on predicates: ``pred1 | pred2``.",python::init<python::object,python::object>());
	python::class_<PredicateIntersection,python::bases<PredicateBoolean> >("PredicateIntersection","Intersection (conjunction) of 2 predicates. A point has to be inside both predicates. Can be constructed using the ``&`` operator on predicates: ``pred1 & pred2``.",python::init<python::object,python::object >());
	python::class_<PredicateDifference,python::bases<PredicateBoolean> >("PredicateDifference","Difference (conjunction with negative predicate) of 2 predicates. A point has to be inside the first and outside the second predicate. Can be constructed using the ``-`` operator on predicates: ``pred1 - pred2``.",python::init<python::object,python::object >());
	python::class_<PredicateSymmetricDifference,python::bases<PredicateBoolean> >("PredicateSymmetricDifference","SymmetricDifference (exclusive disjunction) of 2 predicates. A point has to be in exactly one predicate of the two. Can be constructed using the ``^`` operator on predicates: ``pred1 ^ pred2``.",python::init<python::object,python::object >());
	// primitive predicates
	python::class_<inSphere,python::bases<Predicate> >("inSphere","Sphere predicate.",python::init<const Vector3r&,Real>(python::args("center","radius"),"Ctor taking center (as a 3-tuple) and radius"));
	python::class_<inAlignedBox,python::bases<Predicate> >("inAlignedBox","Axis-aligned box predicate",python::init<const Vector3r&,const Vector3r&>(python::args("minAABB","maxAABB"),"Ctor taking minumum and maximum points of the box (as 3-tuples)."));
	python::class_<inCylinder,python::bases<Predicate> >("inCylinder","Cylinder predicate",python::init<const Vector3r&,const Vector3r&,Real>(python::args("centerBottom","centerTop","radius"),"Ctor taking centers of the lateral walls (as 3-tuples) and radius."));
	python::class_<inHyperboloid,python::bases<Predicate> >("inHyperboloid","Hyperboloid predicate",python::init<const Vector3r&,const Vector3r&,Real,Real>(python::args("centerBottom","centerTop","radius","skirt"),"Ctor taking centers of the lateral walls (as 3-tuples), radius at bases and skirt (middle radius)."));
	python::class_<inEllipsoid,python::bases<Predicate> >("inEllipsoid","Ellipsoid predicate",python::init<const Vector3r&,const Vector3r&>(python::args("centerPoint","abc"),"Ctor taking center of the ellipsoid (3-tuple) and its 3 radii (3-tuple)."));
	python::class_<notInNotch,python::bases<Predicate> >("notInNotch","Outside of infinite, rectangle-shaped notch predicate",python::init<const Vector3r&,const Vector3r&,const Vector3r&,Real>(python::args("centerPoint","edge","normal","aperture"),"Ctor taking point in the symmetry plane, vector pointing along the edge, plane normal and aperture size.\nThe side inside the notch is edge×normal.\nNormal is made perpendicular to the edge.\nAll vectors are normalized at construction time.")); 
	#ifdef YADE_GTS
		python::class_<inGtsSurface,python::bases<Predicate> >("inGtsSurface","GTS surface predicate",python::init<python::object,python::optional<bool> >(python::args("surface","noPad"),"Ctor taking a gts.Surface() instance, which must not be modified during instance lifetime.\nThe optional noPad can disable padding (if set to True), which speeds up calls several times.\nNote: padding checks inclusion of 6 points along +- cardinal directions in the pad distance from given point, which is not exact."))
			.add_property("surf",&inGtsSurface::surface,"The associated gts.Surface object.");
	#endif
}

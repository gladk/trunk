#include<yade/pkg-dem/Shop.hpp>
#include<boost/python.hpp>
#include<yade/extra/boost_python_len.hpp>
#include<yade/core/Scene.hpp>
#include<yade/core/Omega.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/DemXDofGeom.hpp>
#include<yade/pkg-common/Facet.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/NormShearPhys.hpp>
#include<yade/lib-computational-geometry/Hull2d.hpp>
#include<yade/lib-pyutil/doc_opts.hpp>
#include<cmath>
#include<yade/pkg-dem/ViscoelasticPM.hpp>

#include<numpy/ndarrayobject.h>

using namespace std;
namespace py = boost::python;

#ifdef YADE_LOG4CXX
	log4cxx::LoggerPtr logger=log4cxx::Logger::getLogger("yade.utils");
#endif

py::tuple vec2tuple(const Vector3r& v){return py::make_tuple(v[0],v[1],v[2]);}
Vector3r tuple2vec(const py::tuple& t){return Vector3r(py::extract<double>(t[0])(),py::extract<double>(t[1])(),py::extract<double>(t[2])());}
bool isInBB(Vector3r p, Vector3r bbMin, Vector3r bbMax){return p[0]>bbMin[0] && p[0]<bbMax[0] && p[1]>bbMin[1] && p[1]<bbMax[1] && p[2]>bbMin[2] && p[2]<bbMax[2];}

bool ptInAABB(py::tuple p, py::tuple bbMin, py::tuple bbMax){return isInBB(tuple2vec(p),tuple2vec(bbMin),tuple2vec(bbMax));}

/* \todo implement groupMask */
py::tuple aabbExtrema(Real cutoff=0.0, bool centers=false){
	if(cutoff<0. || cutoff>1.) throw invalid_argument("Cutoff must be >=0 and <=1.");
	Real inf=std::numeric_limits<Real>::infinity();
	Vector3r minimum(inf,inf,inf),maximum(-inf,-inf,-inf);
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		shared_ptr<Sphere> s=dynamic_pointer_cast<Sphere>(b->shape); if(!s) continue;
		Vector3r rrr(s->radius,s->radius,s->radius);
		minimum=componentMinVector(minimum,Vector3r(b->state->pos-(centers?Vector3r::Zero():rrr)));
		maximum=componentMaxVector(maximum,Vector3r(b->state->pos+(centers?Vector3r::Zero():rrr)));
	}
	Vector3r dim=maximum-minimum;
	return py::make_tuple(Vector3r(minimum+.5*cutoff*dim),Vector3r(maximum-.5*cutoff*dim));
}

py::tuple negPosExtremeIds(int axis, Real distFactor=1.1){
	py::tuple extrema=aabbExtrema();
	Real minCoord=py::extract<double>(extrema[0][axis])(),maxCoord=py::extract<double>(extrema[1][axis])();
	py::list minIds,maxIds;
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		shared_ptr<Sphere> s=dynamic_pointer_cast<Sphere>(b->shape); if(!s) continue;
		if(b->state->pos[axis]-s->radius*distFactor<=minCoord) minIds.append(b->getId());
		if(b->state->pos[axis]+s->radius*distFactor>=maxCoord) maxIds.append(b->getId());
	}
	return py::make_tuple(minIds,maxIds);
}
BOOST_PYTHON_FUNCTION_OVERLOADS(negPosExtremeIds_overloads,negPosExtremeIds,1,2);

py::tuple coordsAndDisplacements(int axis,py::tuple Aabb=py::tuple()){
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero()); bool useBB=py::len(Aabb)>0;
	if(useBB){bbMin=py::extract<Vector3r>(Aabb[0])();bbMax=py::extract<Vector3r>(Aabb[1])();}
	py::list retCoord,retDispl;
	FOREACH(const shared_ptr<Body>&b, *Omega::instance().getScene()->bodies){
		if(useBB && !isInBB(b->state->pos,bbMin,bbMax)) continue;
		retCoord.append(b->state->pos[axis]);
		retDispl.append(b->state->pos[axis]-b->state->refPos[axis]);
	}
	return py::make_tuple(retCoord,retDispl);
}

void setRefSe3(){
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		b->state->refPos=b->state->pos;
		b->state->refOri=b->state->ori;
	}
}

Real PWaveTimeStep(){return Shop::PWaveTimeStep();};
Real getSpheresVolume(){ return Shop::getSpheresVolume(); }

Real elasticEnergyInAABB(py::tuple Aabb){
	Vector3r bbMin=py::extract<Vector3r>(Aabb[0])(), bbMax=py::extract<Vector3r>(Aabb[1])();
	shared_ptr<Scene> rb=Omega::instance().getScene();
	Real E=0;
	FOREACH(const shared_ptr<Interaction>&i, *rb->interactions){
		if(!i->interactionPhysics) continue;
		shared_ptr<NormShearPhys> bc=dynamic_pointer_cast<NormShearPhys>(i->interactionPhysics); if(!bc) continue;
		shared_ptr<Dem3DofGeom> geom=dynamic_pointer_cast<Dem3DofGeom>(i->interactionGeometry); if(!geom){LOG_ERROR("NormShearPhys contact doesn't have Dem3DofGeom associated?!"); continue;}
		const shared_ptr<Body>& b1=Body::byId(i->getId1(),rb), b2=Body::byId(i->getId2(),rb);
		bool isIn1=isInBB(b1->state->pos,bbMin,bbMax), isIn2=isInBB(b2->state->pos,bbMin,bbMax);
		if(!isIn1 && !isIn2) continue;
		LOG_DEBUG("Interaction #"<<i->getId1()<<"--#"<<i->getId2());
		Real weight=1.;
		if((!isIn1 && isIn2) || (isIn1 && !isIn2)){
			//shared_ptr<Body> bIn=isIn1?b1:b2, bOut=isIn2?b2:b1;
			Vector3r vIn=(isIn1?b1:b2)->state->pos, vOut=(isIn2?b1:b2)->state->pos;
			#define _WEIGHT_COMPONENT(axis) if(vOut[axis]<bbMin[axis]) weight=min(weight,abs((vOut[axis]-bbMin[axis])/(vOut[axis]-vIn[axis]))); else if(vOut[axis]>bbMax[axis]) weight=min(weight,abs((vOut[axis]-bbMax[axis])/(vOut[axis]-vIn[axis])));
			_WEIGHT_COMPONENT(0); _WEIGHT_COMPONENT(1); _WEIGHT_COMPONENT(2);
			assert(weight>=0 && weight<=1);
			//cerr<<"Interaction crosses Aabb boundary, weight is "<<weight<<endl;
			//LOG_DEBUG("Interaction crosses Aabb boundary, weight is "<<weight);
		} else {assert(isIn1 && isIn2); /* cerr<<"Interaction inside, weight is "<<weight<<endl;*/ /*LOG_DEBUG("Interaction inside, weight is "<<weight);*/}
		E+=geom->refLength*weight*(.5*bc->kn*pow(geom->strainN(),2)+.5*bc->ks*pow(geom->strainT().norm(),2));
	}
	return E;
}

/* return histogram ([bin1Min,bin2Min,…],[value1,value2,…]) from projections of interactions
 * to the plane perpendicular to axis∈[0,1,2]; the number of bins can be specified and they cover
 * the range (0,π), since interactions are bidirecional, hence periodically the same on (π,2π).
 *
 * Only contacts using ScGeom are considered.
 * Both bodies must be in the mask (except if mask is 0, when all bodies are considered)
 * If the projection is shorter than minProjLen, it is skipped.
 *
 * If both bodies are _outside_ the aabb (if specified), the interaction is skipped.
 *
 */
py::tuple interactionAnglesHistogram(int axis, int mask=0, size_t bins=20, py::tuple aabb=py::tuple(), Real minProjLen=1e-6){
	if(axis<0||axis>2) throw invalid_argument("Axis must be from {0,1,2}=x,y,z.");
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero()); bool useBB=py::len(aabb)>0; if(useBB){bbMin=py::extract<Vector3r>(aabb[0])();bbMax=py::extract<Vector3r>(aabb[1])();}
	Real binStep=Mathr::PI/bins; int axis2=(axis+1)%3, axis3=(axis+2)%3;
	vector<Real> cummProj(bins,0.);
	shared_ptr<Scene> rb=Omega::instance().getScene();
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions){
		if(!i->isReal()) continue;
		const shared_ptr<Body>& b1=Body::byId(i->getId1(),rb), b2=Body::byId(i->getId2(),rb);
		if(!b1->maskOk(mask) || !b2->maskOk(mask)) continue;
		if(useBB && !isInBB(b1->state->pos,bbMin,bbMax) && !isInBB(b2->state->pos,bbMin,bbMax)) continue;
		GenericSpheresContact* geom=dynamic_cast<GenericSpheresContact*>(i->interactionGeometry.get());
		if(!geom) continue;
		Vector3r n(geom->normal); n[axis]=0.; Real nLen=n.norm();
		if(nLen<minProjLen) continue; // this interaction is (almost) exactly parallel to our axis; skip that one
		Real theta=acos(n[axis2]/nLen)*(n[axis3]>0?1:-1); if(theta<0) theta+=Mathr::PI;
		int binNo=theta/binStep;
		cummProj[binNo]+=nLen;
	}
	py::list val,binMid;
	for(size_t i=0; i<(size_t)bins; i++){ val.append(cummProj[i]); binMid.append(i*binStep);}
	return py::make_tuple(binMid,val);
}
BOOST_PYTHON_FUNCTION_OVERLOADS(interactionAnglesHistogram_overloads,interactionAnglesHistogram,1,4);

py::tuple bodyNumInteractionsHistogram(py::tuple aabb=py::tuple()){
	Vector3r bbMin(Vector3r::Zero()), bbMax(Vector3r::Zero()); bool useBB=py::len(aabb)>0; if(useBB){bbMin=py::extract<Vector3r>(aabb[0])();bbMax=py::extract<Vector3r>(aabb[1])();}
	const shared_ptr<Scene>& rb=Omega::instance().getScene();
	vector<int> bodyNumIntr; bodyNumIntr.resize(rb->bodies->size(),0);
	int maxIntr=0;
	FOREACH(const shared_ptr<Interaction>& i, *rb->interactions){
		if(!i->isReal()) continue;
		const body_id_t id1=i->getId1(), id2=i->getId2(); const shared_ptr<Body>& b1=Body::byId(id1,rb), b2=Body::byId(id2,rb);
		if((useBB && isInBB(b1->state->pos,bbMin,bbMax)) || !useBB) bodyNumIntr[id1]+=1;
		if((useBB && isInBB(b2->state->pos,bbMin,bbMax)) || !useBB) bodyNumIntr[id2]+=1;
		maxIntr=max(max(maxIntr,bodyNumIntr[b1->getId()]),bodyNumIntr[b2->getId()]);
	}
	vector<int> bins; bins.resize(maxIntr+1,0);
	for(size_t id=0; id<bodyNumIntr.size(); id++){
		if(bodyNumIntr[id]>0) bins[bodyNumIntr[id]]+=1;
		// 0 is handled specially: add body to the 0 bin only if it is inside the bb requested (if applicable)
		// otherwise don't do anything, since it is outside the volume of interest
		else if((useBB && isInBB(Body::byId(id,rb)->state->pos,bbMin,bbMax)) || !useBB) bins[0]+=1;
	}
	py::list count,num;
	for(size_t n=0; n<bins.size(); n++){
		if(bins[n]==0) continue;
		num.append(n); count.append(bins[n]);
	}
	return py::make_tuple(num,count);
}
BOOST_PYTHON_FUNCTION_OVERLOADS(bodyNumInteractionsHistogram_overloads,bodyNumInteractionsHistogram,0,1);

Vector3r inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2)
{
	return Shop::inscribedCircleCenter(v0,v1,v2);
}
py::dict getViscoelasticFromSpheresInteraction(Real m, Real tc, Real en, Real es)
{
	shared_ptr<ViscElMat> b = shared_ptr<ViscElMat>(new ViscElMat());
	Shop::getViscoelasticFromSpheresInteraction(m,tc,en,es,b);
	py::dict d;
	d["kn"]=b->kn;
	d["cn"]=b->cn;
	d["ks"]=b->ks;
	d["cs"]=b->cs;
	return d;
}
/* reset highlight of all bodies */
void highlightNone(){
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		if(!b->shape) continue;
		b->shape->highlight=false;
	}
}

/*!Sum moments acting on given bodies
 *
 * @param ids is the calculated bodies ids
 * @param axis is the direction of axis with respect to which the moment is calculated.
 * @param axisPt is a point on the axis.
 *
 * The computation is trivial: moment from force is is by definition r×F, where r
 * is position relative to axisPt; moment from moment is m; such moment per body is
 * projected onto axis.
 */
Real sumTorques(py::tuple ids, const Vector3r& axis, const Vector3r& axisPt){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	rb->forces.sync();
	Real ret=0;
	size_t len=py::len(ids);
	for(size_t i=0; i<len; i++){
		const Body* b=(*rb->bodies)[py::extract<int>(ids[i])].get();
		const Vector3r& m=rb->forces.getTorque(b->getId());
		const Vector3r& f=rb->forces.getForce(b->getId());
		Vector3r r=b->state->pos-axisPt;
		ret+=axis.dot(m+r.cross(f));
	}
	return ret;
}
/* Sum forces acting on bodies within mask.
 *
 * @param ids list of ids
 * @param direction direction in which forces are summed
 *
 */
Real sumForces(py::tuple ids, const Vector3r& direction){
	shared_ptr<Scene> rb=Omega::instance().getScene();
	rb->forces.sync();
	Real ret=0;
	size_t len=py::len(ids);
	for(size_t i=0; i<len; i++){
		body_id_t id=py::extract<int>(ids[i]);
		const Vector3r& f=rb->forces.getForce(id);
		ret+=direction.dot(f);
	}
	return ret;
}

/* Sum force acting on facets given by their ids in the sense of their respective normals.
   If axis is given, it will sum forces perpendicular to given axis only (not the the facet normals).
*/
Real sumFacetNormalForces(vector<body_id_t> ids, int axis=-1){
	shared_ptr<Scene> rb=Omega::instance().getScene(); rb->forces.sync();
	Real ret=0;
	FOREACH(const body_id_t id, ids){
		Facet* f=YADE_CAST<Facet*>(Body::byId(id,rb)->shape.get());
		if(axis<0) ret+=rb->forces.getForce(id).dot(f->nf);
		else {
			Vector3r ff=rb->forces.getForce(id); ff[axis]=0;
			ret+=ff.dot(f->nf);
		}
	}
	return ret;
}

/* Set wire display of all/some/none bodies depending on the filter. */
void wireSome(string filter){
	enum{none,all,noSpheres,unknown};
	int mode=(filter=="none"?none:(filter=="all"?all:(filter=="noSpheres"?noSpheres:unknown)));
	if(mode==unknown) { LOG_WARN("Unknown wire filter `"<<filter<<"', using noSpheres instead."); mode=noSpheres; }
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		if(!b->shape) return;
		bool wire;
		switch(mode){
			case none: wire=false; break;
			case all: wire=true; break;
			case noSpheres: wire=!(bool)(dynamic_pointer_cast<Sphere>(b->shape)); break;
			default: throw logic_error("No such case possible");
		}
		b->shape->wire=wire;
	}
}
void wireAll(){wireSome("all");}
void wireNone(){wireSome("none");}
void wireNoSpheres(){wireSome("noSpheres");}


/* Tell us whether a point lies in polygon given by array of points.
 *  @param xy is the point that is being tested
 *  @param vertices is Numeric.array (or list or tuple) of vertices of the polygon.
 *         Every row of the array is x and y coordinate, numer of rows is >= 3 (triangle).
 *
 * Copying the algorithm from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
 * is gratefully acknowledged:
 *
 * License to Use:
 * Copyright (c) 1970-2003, Wm. Randolph Franklin
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
 *   2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
 *   3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * http://numpy.scipy.org/numpydoc/numpy-13.html told me how to use Numeric.array from c
 */
bool pointInsidePolygon(py::tuple xy, py::object vertices){
	Real testx=py::extract<double>(xy[0])(),testy=py::extract<double>(xy[1])();
	char** vertData; int rows, cols; PyArrayObject* vert=(PyArrayObject*)vertices.ptr();
	int result=PyArray_As2D((PyObject**)&vert /* is replaced */ ,&vertData,&rows,&cols,PyArray_DOUBLE);
	if(result!=0) throw invalid_argument("Unable to cast vertices to 2d array");
	if(cols!=2 || rows<3) throw invalid_argument("Vertices must have 2 columns (x and y) and at least 3 rows.");
	int i /*current node*/, j/*previous node*/; bool inside=false;
	for(i=0,j=rows-1; i<rows; j=i++){
		double vx_i=*(double*)(vert->data+i*vert->strides[0]), vy_i=*(double*)(vert->data+i*vert->strides[0]+vert->strides[1]), vx_j=*(double*)(vert->data+j*vert->strides[0]), vy_j=*(double*)(vert->data+j*vert->strides[0]+vert->strides[1]);
		if (((vy_i>testy)!=(vy_j>testy)) && (testx < (vx_j-vx_i) * (testy-vy_i) / (vy_j-vy_i) + vx_i) ) inside=!inside;
	}
	Py_DECREF(vert);
	return inside;
}

#if 0
/* Compute convex hull of given points, given as python list of Vector2 */
py::list convexHull2d(const py::list& pts){
	size_t l=py::len(pts);
	py::list ret;
	std::list<Vector2r> pts2;
	for(size_t i=0; i<l; i++){
		pts2.push_back(py::extract<Vector2r>(pts[i]));
		cerr<<*pts2.rbegin()<<endl;
	}
	ConvexHull2d ch2d(pts2);
	vector<Vector2r> hull=ch2d();
	FOREACH(const Vector2r& v, hull) ret.append(v);
	cout<<pts2.size()<<" "<<hull.size()<<"@@@"<<endl;
	return ret;
}
#endif 

/* Compute area of convex hull when when taking (swept) spheres crossing the plane at coord, perpendicular to axis.

	All spheres that touch the plane are projected as hexagons on their circumference to the plane.
	Convex hull from this cloud is computed.
	The area of the hull is returned.

*/
Real approxSectionArea(Real coord, int axis){
	std::list<Vector2r> cloud;
	if(axis<0 || axis>2) throw invalid_argument("Axis must be ∈ {0,1,2}");
	const int ax1=(axis+1)%3, ax2=(axis+2)%3;
	const Real sqrt3=sqrt(3);
	Vector2r mm,mx; int i=0;
	FOREACH(const shared_ptr<Body>& b, *Omega::instance().getScene()->bodies){
		Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
		if(!s) continue;
		const Vector3r& pos(b->state->pos); const Real r(s->radius);
		if((pos[axis]>coord && (pos[axis]-r)>coord) || (pos[axis]<coord && (pos[axis]+r)<coord)) continue;
		Vector2r c(pos[ax1],pos[ax2]);
		cloud.push_back(c+Vector2r(r,0.)); cloud.push_back(c+Vector2r(-r,0.));
		cloud.push_back(c+Vector2r( r/2., sqrt3*r)); cloud.push_back(c+Vector2r( r/2.,-sqrt3*r));
		cloud.push_back(c+Vector2r(-r/2., sqrt3*r)); cloud.push_back(c+Vector2r(-r/2.,-sqrt3*r));
		if(i++==0){ mm=c, mx=c;}
		mm=Vector2r(min(c[0]-r,mm[0]),min(c[1]-r,mm[1]));
		mx=Vector2r(max(c[0]+r,mx[0]),max(c[1]+r,mx[1]));
	}
	if(cloud.size()==0) return 0;
	ConvexHull2d ch2d(cloud);
	vector<Vector2r> hull=ch2d();
	return simplePolygonArea2d(hull);
}
/* Find all interactions deriving from NormShearPhys that cross plane given by a point and normal
	(the normal may not be normalized in this case, though) and sum forces (both normal and shear) on them.
	
	Return a 3-tuple with the components along global x,y,z axes.

	(This could be easily extended to return sum of only normal forces or only of shear forces.)
*/
Vector3r forcesOnPlane(const Vector3r& planePt, const Vector3r&  normal){
	Vector3r ret(Vector3r::Zero());
	Scene* rootBody=Omega::instance().getScene().get();
	FOREACH(const shared_ptr<Interaction>&I, *rootBody->interactions){
		if(!I->isReal()) continue;
		NormShearPhys* nsi=dynamic_cast<NormShearPhys*>(I->interactionPhysics.get());
		if(!nsi) continue;
		Vector3r pos1,pos2;
		Dem3DofGeom* d3dg=dynamic_cast<Dem3DofGeom*>(I->interactionGeometry.get()); // Dem3DofGeom has copy of se3 in itself, otherwise we have to look up the bodies
		if(d3dg){ pos1=d3dg->se31.position; pos2=d3dg->se32.position; }
		else{ pos1=Body::byId(I->getId1(),rootBody)->state->pos; pos2=Body::byId(I->getId2(),rootBody)->state->pos; }
		Real dot1=(pos1-planePt).dot(normal), dot2=(pos2-planePt).dot(normal);
		if(dot1*dot2>0) continue; // both interaction points on the same side of the plane
		// if pt1 is on the negative plane side, d3dg->normal.Dot(normal)>0, the force is well oriented;
		// otherwise, reverse its contribution
		ret+=(dot1<0.?1.:-1.)*(nsi->normalForce+nsi->shearForce);
	}
	return ret;
}

/* Less general than forcesOnPlane, computes force on plane perpendicular to axis, passing through coordinate coord. */
Vector3r forcesOnCoordPlane(Real coord, int axis){
	Vector3r planePt(Vector3r::Zero()); planePt[axis]=coord;
	Vector3r normal(Vector3r::Zero()); normal[axis]=1;
	return forcesOnPlane(planePt,normal);
}


py::tuple spiralProject(const Vector3r& pt, Real dH_dTheta, int axis=2, Real periodStart=std::numeric_limits<Real>::quiet_NaN(), Real theta0=0){
	Real r,h,theta;
	boost::tie(r,h,theta)=Shop::spiralProject(pt,dH_dTheta,axis,periodStart,theta0);
	return py::make_tuple(py::make_tuple(r,h),theta);
}
//BOOST_PYTHON_FUNCTION_OVERLOADS(spiralProject_overloads,spiralProject,2,5);

shared_ptr<Interaction> Shop__createExplicitInteraction(body_id_t id1, body_id_t id2){ return Shop::createExplicitInteraction(id1,id2,/*force*/true);}

//Vector3r Shop__scalarOnColorScale(Real scalar){ return Shop::scalarOnColorScale(scalar);}

Real Shop__kineticEnergy(){return Shop::kineticEnergy();}
Real Shop__unbalancedForce(bool useMaxForce /*false by default*/){return Shop::unbalancedForce(useMaxForce);}
py::tuple Shop__totalForceInVolume(){Real stiff; Vector3r ret=Shop::totalForceInVolume(stiff); return py::make_tuple(ret,stiff); }

BOOST_PYTHON_MODULE(_utils){
	// http://numpy.scipy.org/numpydoc/numpy-13.html mentions this must be done in module init, otherwise we will crash
	import_array();

	YADE_SET_DOCSTRING_OPTS;


	py::def("PWaveTimeStep",PWaveTimeStep,"Get timestep accoring to the velocity of P-Wave propagation; computed from sphere radii, rigidities and masses.");
	py::def("getSpheresVolume",getSpheresVolume,"Compute the total volume of spheres in the simulation (might crash for now if dynamic bodies are not spheres)");
	py::def("aabbExtrema",aabbExtrema,(py::arg("cutoff")=0.0,py::arg("centers")=false),"Return coordinates of box enclosing all bodies\n\n:param bool centers: do not take sphere radii in account, only their centroids\n:param float∈〈0…1〉 cutoff: relative dimension by which the box will be cut away at its boundaries.\n\n\n:return: (lower corner, upper corner) as (Vector3,Vector3)\n\n");
	py::def("ptInAABB",ptInAABB,"Return True/False whether the point (3-tuple) p is within box given by its min (3-tuple) and max (3-tuple) corners");
	py::def("negPosExtremeIds",negPosExtremeIds,negPosExtremeIds_overloads(py::args("axis","distFactor"),"Return list of ids for spheres (only) that are on extremal ends of the specimen along given axis; distFactor multiplies their radius so that sphere that do not touch the boundary coordinate can also be returned."));
	py::def("approxSectionArea",approxSectionArea,"Compute area of convex hull when when taking (swept) spheres crossing the plane at coord, perpendicular to axis.");
	#if 0
		py::def("convexHull2d",convexHull2d,"Return 2d convex hull of list of 2d points, as list of polygon vertices.");
	#endif
	py::def("coordsAndDisplacements",coordsAndDisplacements,(python::arg("axis"),python::arg("Aabb")=python::tuple()),"Return tuple of 2 same-length lists for coordinates and displacements (coordinate minus reference coordinate) along given axis (1st arg); if the Aabb=((x_min,y_min,z_min),(x_max,y_max,z_max)) box is given, only bodies within this box will be considered.");
	py::def("setRefSe3",setRefSe3,"Set reference :yref:`positions<State::refPos>` and :yref:`orientations<State::refOri>` of all :yref:`bodies<Body>` equal to their current :yref:`positions<State::pos>` and :yref:`orientations<State::ori>`.");
	py::def("interactionAnglesHistogram",interactionAnglesHistogram,interactionAnglesHistogram_overloads(py::args("axis","mask","bins","aabb")));
	py::def("bodyNumInteractionsHistogram",bodyNumInteractionsHistogram,bodyNumInteractionsHistogram_overloads(py::args("aabb")));
	py::def("elasticEnergy",elasticEnergyInAABB);
	py::def("inscribedCircleCenter",inscribedCircleCenter,(py::arg("v1"),py::arg("v2"),py::arg("v3")),"Return center of inscribed circle for triangle given by its vertices *v1*, *v2*, *v3*.");
	py::def("unbalancedForce",&Shop__unbalancedForce,(py::args("useMaxForce")=false),"Compute the ratio of mean (or maximum, if *useMaxForce*) summary force on bodies and maximum force magnitude on interactions. For perfectly static equilibrium, summary force on all bodies is zero (since forces from interactions cancel out and induce no acceleration of particles); this ratio will tend to zero as simulation stabilizes, though zero is never reached because of finite precision computation. Sufficiently small value can be e.g. 1e-2 or smaller, depending on how much equilibrium it should be.");
	py::def("kineticEnergy",Shop__kineticEnergy,"Compute overall kinetic energy of the simulation as\n\n.. math:: \\sum\\frac{1}{2}\\left(m_i\\vec{v}_i^2+\\vec{\\omega}(\\mat{I}\\vec{\\omega}^T)\\right).\n\n.. warning::\n\n\tNo transformation of inertia tensor (in local frame) $\\mat{I}$ is done, although it is multiplied by angular velocity $\\vec{\\omega}$ (in global frame); the value will not be accurate for aspherical particles.\n");
	py::def("sumForces",sumForces,(py::arg("ids"),py::arg("direction")),"Return summary force on bodies with given *ids*, projected on the *direction* vector.");
	py::def("sumTorques",sumTorques,(py::arg("ids"),py::arg("axis"),py::arg("axisPt")),"Sum forces and torques on bodies given in *ids* with respect to axis specified by a point *axisPt* and its direction *axis*.");
	py::def("sumFacetNormalForces",sumFacetNormalForces,(py::arg("ids"),py::arg("axis")=-1),"Sum force magnitudes on given bodies (must have :yref:`shape<Body.shape>` of the :yref:`Facet` type), considering only part of forces perpendicular to each :yref:`facet's<Facet>` face; if *axis* has positive value, then the specified axis (0=x, 1=y, 2=z) will be used instead of facet's normals.");
	py::def("forcesOnPlane",forcesOnPlane,(py::arg("planePt"),py::arg("normal")),"Find all interactions deriving from :yref:`NormShearPhys` that cross given plane and sum forces (both normal and shear) on them.\n\n:param Vector3 planePt: a point on the plane\n:param Vector3 normal: plane normal (will be normalized).\n");
	py::def("forcesOnCoordPlane",forcesOnCoordPlane);
	py::def("totalForceInVolume",Shop__totalForceInVolume,"Return summed forces on all interactions and average isotropic stiffness, as tuple (Vector3,float)");
	py::def("createInteraction",Shop__createExplicitInteraction,(py::arg("id1"),py::arg("id2")),"Create interaction between given bodies by hand.\n\nCurrent engines are searched for :yref:`InteractionGeometryDispatcher` and :yref:`InteractionPhysicsDispatcher` (might be both hidden in :yref:`InteractionDispatchers`). Geometry is created using ``force`` parameter of the :yref:`geometry dispatcher<InteractionGeometryDispatcher>`, wherefore the interaction will exist even if bodies do not spatially overlap and the functor would return ``false`` under normal circumstances. \n\n.. warning::\n\tThis function will very likely behave incorrectly for periodic simulations (though it could be extended it to handle it farily easily).");
	py::def("spiralProject",spiralProject,(py::arg("pt"),py::arg("dH_dTheta"),py::arg("axis")=2,py::arg("periodStart")=std::numeric_limits<Real>::quiet_NaN(),py::arg("theta0")=0));
	py::def("pointInsidePolygon",pointInsidePolygon);
	py::def("scalarOnColorScale",Shop::scalarOnColorScale);
	py::def("highlightNone",highlightNone,"Reset :yref:`highlight<Shape::highlight>` on all bodies.");
	py::def("wireAll",wireAll,"Set :yref:`Shape::wire` on all bodies to True, rendering them with wireframe only.");
	py::def("wireNone",wireNone,"Set :yref:`Shape::wire` on all bodies to False, rendering them as solids.");
	py::def("wireNoSpheres",wireNoSpheres,"Set :yref:`Shape::wire` to True on non-spherical bodies (:yref:`Facets<Facet>`, :yref:`Walls<Wall>`).");
	py::def("flipCell",&Shop::flipCell,(py::arg("flip")=Matrix3r(Matrix3r::Zero())),"Flip periodic cell so that angles between $R^3$ axes and transformed axes are as small as possible. This function relies on the fact that periodic cell defines by repetition or its corners regular grid of points in $R^3$; however, all cells generating identical grid are equivalent and can be flipped one over another. This necessiatates adjustment of :yref:`Interaction.cellDist` for interactions that cross boundary and didn't before (or vice versa), and re-initialization of collider. The *flip* argument can be used to specify desired flip: integers, each column for one axis; if zero matrix, best fit (minimizing the angles) is computed automatically.\n\nIn c++, this function is accessible as ``Shop::flipCell``.\n\n.. warning::\n\t This function is currently broken and should not be used.");
	py::def("getViscoelasticFromSpheresInteraction",getViscoelasticFromSpheresInteraction,(py::arg("m"),py::arg("tc"),py::arg("en"),py::arg("es")),"Get viscoelastic interaction parameters from analytical solution of a pair spheres collision problem. \n\n:Parameters:\n\t`m` : float\n\t\tsphere mass\n\t`tc` : float\n\t\tcollision time\n\t`en` : float\n\t\tnormal restitution coefficient\n\t`es` : float\n\t\ttangential restitution coefficient.\n\n\n:return: \n\tdict with keys:\n\n\tkn : float\n\t\tnormal elastic coefficient computed as:\n\n.. math:: k_n=\\frac{m}{t_c^2}\\left(\\pi^2+(\\ln e_n)^2\\right)\n\n\tcn : float\n\t\tnormal viscous coefficient computed as:\n\n.. math:: c_n=-\\frac{2m}{t_c}\\ln e_n\n\n\n\tkt : float\n\t\ttangential elastic coefficient computed as:\n\n.. math:: k_t=\\frac27\\frac{m}{t_c^2}\\left(\\pi^2+(\\ln e_t)^2\\right)\n\n\tct : float\n\t\ttangential viscous coefficient computed as:\n\n.. math:: c_t=-\\frac27\\frac{m}{t_c}\\ln e_t.\n\n\nFor details see [Pournin2001]_. ");
}

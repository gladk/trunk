/*************************************************************************
*  Copyright (C) 2008 by Bruno Chareyre                                  *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<yade/lib-pyutil/numpy.hpp>
#include<boost/python.hpp>
#include<yade/extra/boost_python_len.hpp>
#include<yade/pkg-dem/Shop.hpp>
#include"TesselationWrapper.hpp"

YADE_PLUGIN((TesselationWrapper));
YADE_REQUIRE_FEATURE(CGAL)
CREATE_LOGGER(TesselationWrapper);

//spatial sort traits to use with a pair of CGAL::sphere pointers and integer.
//template<class _Triangulation>
struct RTraits_for_spatial_sort : public CGT::RTriangulation::Geom_traits {
	//typedef typename _Triangulation::Geom_traits Gt;
	typedef CGT::RTriangulation::Geom_traits Gt;
	typedef std::pair<const CGT::Sphere*,body_id_t> Point_3;

	struct Less_x_3 {
		bool operator()(const Point_3& p,const Point_3& q) const {
			return Gt::Less_x_3()(* (p.first),* (q.first));
		}
	};
	struct Less_y_3 {
		bool operator()(const Point_3& p,const Point_3& q) const {
			return Gt::Less_y_3()(* (p.first),* (q.first));
		}
	};
	struct Less_z_3 {
		bool operator()(const Point_3& p,const Point_3& q) const {
			return Gt::Less_z_3()(* (p.first),* (q.first));
		}
	};
	Less_x_3  less_x_3_object() const {return Less_x_3();}
	Less_y_3  less_y_3_object() const {return Less_y_3();}
	Less_z_3  less_z_3_object() const {return Less_z_3();}
};


//function inserting points into a triangulation (where YADE::Sphere is converted to CGT::Sphere)
//and setting the info field to the bodies id.
//Possible improvements : use bodies pointers to avoid one copy, use aabb's lists to replace the shuffle/sort part
// template <class Triangulation>
void build_triangulation_with_ids(const shared_ptr<BodyContainer>& bodies, TesselationWrapper &TW)
{
	TW.clear();
	CGT::Tesselation& Tes = *(TW.Tes);
	CGT::RTriangulation& T = Tes.Triangulation();
	std::vector<CGT::Sphere> spheres;
	std::vector<std::pair<const CGT::Sphere*,body_id_t> > pointsPtrs;
	spheres.reserve(bodies->size());
	pointsPtrs.reserve(bodies->size());

	BodyContainer::iterator biBegin    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();
	BodyContainer::iterator bi = biBegin;

	body_id_t Ng = 0;
	body_id_t MaxId=0;
	TW.mean_radius = 0;
	for (; bi!=biEnd ; ++bi) {
		if ((*bi)->isDynamic) { //then it is a sphere (not a wall) FIXME : need test if isSphere
			const Sphere* s = YADE_CAST<Sphere*> ((*bi)->shape.get());
			const Vector3r& pos = (*bi)->state->pos;
			const Real rad = s->radius;
			CGT::Sphere sp(CGT::Point(pos[0],pos[1],pos[2]),rad*rad);
			spheres.push_back(sp);
			pointsPtrs.push_back(std::make_pair(&(spheres[Ng]/*.point()*/),(*bi)->getId()));
			TW.Pmin = CGT::Point(min(TW.Pmin.x(),pos.x()-rad),min(TW.Pmin.y(), pos.y()-rad),min(TW.Pmin.z(), pos.z()-rad));
			TW.Pmax = CGT::Point(max(TW.Pmax.x(),pos.x()+rad),max(TW.Pmax.y(),pos.y()+rad),max(TW.Pmax.z(),pos.z()+rad));
			Ng++; TW.mean_radius += rad;
			MaxId = max(MaxId,(*bi)->getId());
		}
	}
	TW.mean_radius /= Ng; TW.rad_divided = true;
	spheres.resize(Ng);
	pointsPtrs.resize(Ng);
	Tes.vertexHandles.resize(MaxId+1);
	Tes.redirected = 1;
	std::random_shuffle(pointsPtrs.begin(), pointsPtrs.end());
	spatial_sort(pointsPtrs.begin(),pointsPtrs.end(), RTraits_for_spatial_sort()/*, CGT::RTriangulation::Weighted_point*/);

	CGT::RTriangulation::Cell_handle hint;

	TW.n_spheres = 0;
	for (std::vector<std::pair<const CGT::Sphere*,body_id_t> >::const_iterator
			p = pointsPtrs.begin();p != pointsPtrs.end(); ++p) {
		CGT::RTriangulation::Locate_type lt;
		CGT::RTriangulation::Cell_handle c;
		int li, lj;
		c = T.locate(* (p->first), lt, li, lj, hint);
		CGT::RTriangulation::Vertex_handle v = T.insert(*(p->first),lt,c,li,lj);
		if (v==CGT::RTriangulation::Vertex_handle())
			hint=c;
		else {
			v->info() = (const unsigned int) p->second;
			//Vh->info().isFictious = false;//false is the default
			Tes.max_id = std::max(Tes.max_id,(int) p->second);
			Tes.vertexHandles[p->second]=v;
			hint=v->cell();
			++TW.n_spheres;
		}
	}
	//cerr << " loaded : " << Ng<<", triangulated : "<<TW.n_spheres<<", mean radius = " << TW.mean_radius<<endl;
}


// is this a joke?? #include<limits> std::numeric_limits<double>::infinity();
//__attribute__((unused)) static double inf = 1e10; 
double pminx=0;
double pminy=0;
double pminz=0;
double pmaxx=0;
double pmaxy=0;
double pmaxz=0;
double thickness = 0;

//  Finite_edges_iterator facet_it; //an edge in a triangulation is a facet in corresponding tesselation, remember...
//That explain the name.



TesselationWrapper::~TesselationWrapper() { delete Tes;}

void TesselationWrapper::clear(void)
{
	Tes->Clear();
	Pmin = CGT::Point(inf, inf, inf);
	Pmax = CGT::Point(-inf, -inf, -inf);
	mean_radius = 0;
	n_spheres = 0;
	rad_divided = false;
	bounded = false;
	facet_it = Tes->Triangulation().finite_edges_end();
}

void TesselationWrapper::clear2(void) //for testing purpose
{
	Tes->Clear();

//   Pmin = Point(inf, inf, inf);
//  Pmax = Point(-inf, -inf, -inf);
//  mean_radius = 0;
//  n_spheres = 0;
//  rad_divided = false;
// bounded = false;
//  facet_it = Tes->Triangulation().finite_edges_end ();
}

bool TesselationWrapper::insertSceneSpheres(const Scene* scene)
{
	const shared_ptr<BodyContainer>& bodies = scene->bodies;
	build_triangulation_with_ids(bodies, *this);
	return 1;
}

double TesselationWrapper::Volume(unsigned int id) {return Tes->Volume(id);}

bool TesselationWrapper::insert(double x, double y, double z, double rad, unsigned int id)
{
	using namespace std;
	Pmin = CGT::Point(min(Pmin.x(), x-rad),
					  min(Pmin.y(), y-rad),
					  min(Pmin.z(), z-rad));
	Pmax = CGT::Point(max(Pmax.x(), x+rad),
					  max(Pmax.y(), y+rad),
					  max(Pmax.z(), z+rad));
	mean_radius += rad;
	++n_spheres;
	return (Tes->insert(x,y,z,rad,id)!=NULL);
}

void TesselationWrapper::checkMinMax(double x, double y, double z, double rad)
{
	using namespace std;
	Pmin = CGT::Point(min(Pmin.x(), x-rad),
					  min(Pmin.y(), y-rad),
					  min(Pmin.z(), z-rad));
	Pmax = CGT::Point(max(Pmax.x(), x+rad),
					  max(Pmax.y(), y+rad),
					  max(Pmax.z(), z+rad));
	mean_radius += rad;
	++n_spheres;
}


bool TesselationWrapper::move(double x, double y, double z, double rad, unsigned int id)
{
	using namespace std;

	Pmin = CGT::Point(min(Pmin.x(), x-rad),
					  min(Pmin.y(), y-rad),
					  min(Pmin.z(), z-rad));
	Pmax = CGT::Point(max(Pmax.x(), x+rad),
					  max(Pmax.y(), y+rad),
					  max(Pmax.z(), z+rad));
	mean_radius += rad;

	if (Tes->move(x,y,z,rad,id)!=NULL)
		return true;
	else {
		cerr << "Tes->move(x,y,z,rad,id)==NULL" << endl; return false;
	}
}

void TesselationWrapper::ComputeTesselation(void)
{
	if (!rad_divided) {
		mean_radius /= n_spheres;
		rad_divided = true;
	}
	Tes->Compute();
}

void TesselationWrapper::ComputeTesselation(double pminx, double pmaxx, double pminy, double pmaxy, double pminz, double pmaxz, double dt)
{
	AddBoundingPlanes(pminx, pmaxx,  pminy,  pmaxy, pminz, pmaxz, dt);
	ComputeTesselation();
}

void TesselationWrapper::ComputeVolumes(void)
{
	ComputeTesselation();
	Tes->ComputeVolumes();
}
unsigned int TesselationWrapper::NumberOfFacets(bool initIters)
{
	if (initIters) InitIter();
	return Tes->Triangulation().number_of_finite_edges();
}

void TesselationWrapper::InitIter(void)
{
	facet_begin = Tes->Triangulation().finite_edges_begin();
	facet_end = Tes->Triangulation().finite_edges_end();
	facet_it = facet_begin;
}

bool TesselationWrapper::nextFacet(std::pair<unsigned int,unsigned int>& facet)
{
	if (facet_end==facet_it) return false;
	facet.first = facet_it->first->vertex(facet_it->second)->info().id();
	facet.second = facet_it->first->vertex((facet_it)->third)->info().id();
	++facet_it;
	return true;
}

void TesselationWrapper::AddBoundingPlanes(double pminx, double pmaxx, double pminy, double pmaxy, double pminz, double pmaxz,double dt)
{
	//Not sure this hack form JFJ works in all cases (?)
	if (dt == 0) {
		thickness = -1*pminx;
	}
	if (!bounded) {
		if (!rad_divided) {
			mean_radius /= n_spheres;
			rad_divided = true;
		}
		double FAR = 10000;
		Tes->vertexHandles[0]=Tes->insert(0.5*(pminx+pmaxx),pminy-FAR*(pmaxx-pminx),0.5*(pmaxz+pminz),FAR*(pmaxx-pminx)+thickness,0,true);
		Tes->vertexHandles[1]=Tes->insert(0.5*(pminx+pmaxx), pmaxy+FAR*(pmaxx-pminx),0.5*(pmaxz+pminz),FAR*(pmaxx-pminx)+thickness, 1, true);
		Tes->vertexHandles[2]=Tes->insert(pminx-FAR*(pmaxy-pminy), 0.5*(pmaxy+pminy), 0.5*(pmaxz+pminz), FAR*(pmaxy-pminy)+thickness, 2, true);
		Tes->vertexHandles[3]=Tes->insert(pmaxx+FAR*(pmaxx-pminy), 0.5*(pmaxy+pminy), 0.5*(pmaxz+pminz), FAR*(pmaxy-pminy)+thickness, 3, true);
		Tes->vertexHandles[4]=Tes->insert(0.5*(pminx+pmaxx), 0.5*(pmaxy+pminy), pminz-FAR*(pmaxy-pminy), FAR*(pmaxy-pminy)+thickness, 4, true);
		Tes->vertexHandles[5]=Tes->insert(0.5*(pminx+pmaxx), 0.5*(pmaxy+pminy), pmaxz+FAR*(pmaxy-pminy), FAR*(pmaxy-pminy)+thickness, 5, true);
		bounded = true;
	}
}

void  TesselationWrapper::AddBoundingPlanes(void)
{
	if (!bounded) {
		if (!rad_divided) {
			mean_radius /= n_spheres;
			rad_divided = true;
		}
		double FAR = 10000;
		//Add big bounding spheres with isFictious=true
		Tes->vertexHandles[0]=Tes->insert(0.5*(Pmin.x()+Pmax.x()),Pmin.y()-FAR*(Pmax.x()-Pmin.x()),0.5*(Pmax.z()+Pmin.z()),FAR*(Pmax.x()-Pmin.x()), 0, true);
		Tes->vertexHandles[1]=Tes->insert(0.5*(Pmin.x()+Pmax.x()),Pmax.y()+FAR*(Pmax.x()-Pmin.x()),0.5*(Pmax.z()+Pmin.z()),FAR*(Pmax.x()-Pmin.x()), 1, true);
		Tes->vertexHandles[2]=Tes->insert(Pmin.x()-FAR*(Pmax.y()-Pmin.y()),0.5*(Pmax.y()+Pmin.y()),0.5*(Pmax.z()+Pmin.z()),FAR*(Pmax.y()-Pmin.y()),2, true);
		Tes->vertexHandles[3]=Tes->insert(Pmax.x()+FAR*(Pmax.y()-Pmin.y()),0.5*(Pmax.y()+Pmin.y()),0.5*(Pmax.z()+Pmin.z()),FAR*(Pmax.y()-Pmin.y()),3,true);
		Tes->vertexHandles[4]=Tes->insert(0.5*(Pmin.x()+Pmax.x()),0.5*(Pmax.y()+Pmin.y()),Pmin.z()-FAR*(Pmax.y()-Pmin.y()),FAR*(Pmax.y()-Pmin.y()),4,true);
		Tes->vertexHandles[5]=Tes->insert(0.5*(Pmin.x()+Pmax.x()),0.5*(Pmax.y()+Pmin.y()),Pmax.z()+FAR*(Pmax.y()-Pmin.y()),FAR*(Pmax.y()-Pmin.y()),5, true);
		bounded = true;
	}
}

void  TesselationWrapper::RemoveBoundingPlanes(void)
{
	Tes->remove(0);
	Tes->remove(1);
	Tes->remove(2);
	Tes->remove(3);
	Tes->remove(4);
	Tes->remove(5);
	Pmin = CGT::Point(inf, inf, inf);
	Pmax = CGT::Point(-inf, -inf, -inf);
	mean_radius = 0;
	rad_divided = false;
	bounded = false;
}

void TesselationWrapper::setState (bool state){ mma.setState(state ? 2 : 1);}

python::dict TesselationWrapper::getVolPoroDef(bool deformation)
{
		Scene* scene=Omega::instance().getScene().get();
		delete Tes;
		CGT::TriaxialState* ts;
		if (deformation){//use the final state to compute volumes
			mma.analyser->ComputeParticlesDeformation();
			Tes = &mma.analyser->TS1->tesselation();
			ts = mma.analyser->TS1;
			}
		else {	Tes = &mma.analyser->TS0->tesselation();//no reason to use the final state if we don't want to compute deformations, keep using the initial
			ts = mma.analyser->TS0;}
		CGT::RTriangulation& Tri = Tes->Triangulation();
		Pmin=ts->box.base; Pmax=ts->box.sommet;
		//if (!scene->isPeriodic) AddBoundingPlanes();
		ComputeVolumes();
		int bodiesDim = scene->bodies->size();
		int dim1[]={bodiesDim};
		int dim2[]={bodiesDim,9};
		/// This is the code that needs numpy include
		//numpy_boost<body_id_t,1> id(dim1);
 		numpy_boost<double,1> vol(dim1);
 		numpy_boost<double,1> poro(dim1);
 		numpy_boost<double,2> def(dim2);
 		//FOREACH(const shared_ptr<Body>& b, *scene->bodies){
 		for (CGT::RTriangulation::Finite_vertices_iterator  V_it = Tri.finite_vertices_begin(); V_it !=  Tri.finite_vertices_end(); V_it++) {
 			//id[]=V_it->info().id()
 			//if(!b) continue;
 			const body_id_t id = V_it->info().id();
 			Real sphereVol = 4.188790 * std::pow ( ( V_it->point().weight() ),1.5 );// 4/3*PI*R³ = 4.188...*R³
 			vol[id]=V_it->info().v();			
 			poro[id]=(V_it->info().v() - sphereVol)/V_it->info().v();
			if (deformation) MATRIX3R_TO_NUMPY(mma.analyser->ParticleDeformation[id],def[id]);
 			//cerr << V_it->info().v()<<" "<<ParticleDeformation[id]<<endl;
 		}
 		python::dict ret;
 		ret["vol"]=vol;
 		ret["poro"]=poro;
 		if (deformation) ret["def"]=def;		
 		return ret;
}

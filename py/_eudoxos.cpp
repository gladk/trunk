#if 0
	#include<yade/lib-pyutil/numpy.hpp>
#endif

#include<yade/pkg-dem/ConcretePM.hpp>
#include<boost/python.hpp>
#include<yade/extra/boost_python_len.hpp>
#include<yade/pkg-dem/Shop.hpp>

#ifdef YADE_VTK
	#pragma GCC diagnostic ignored "-Wdeprecated"
		#include<vtkPointLocator.h>
		#include<vtkIdList.h>
		#include<vtkUnstructuredGrid.h>
		#include<vtkPoints.h>
	#pragma GCC diagnostic warning "-Wdeprecated"
#endif

namespace py = boost::python;
using namespace std;
#ifdef YADE_LOG4CXX
	log4cxx::LoggerPtr logger=log4cxx::Logger::getLogger("yade.eudoxos");
#endif

#ifdef YADE_CPM_FULL_MODEL_AVAILABLE
	#include"../../brefcom-mm.hh"
#endif

# if 0
Real elasticEnergyDensityInAABB(python::tuple Aabb){
	Vector3r bbMin=tuple2vec(python::extract<python::tuple>(Aabb[0])()), bbMax=tuple2vec(python::extract<python::tuple>(Aabb[1])()); Vector3r box=bbMax-bbMin;
	shared_ptr<Scene> rb=Omega::instance().getScene();
	Real E=0;
	FOREACH(const shared_ptr<Interaction>&i, *rb->interactions){
		if(!i->phys) continue;
		shared_ptr<CpmPhys> bc=dynamic_pointer_cast<CpmPhys>(i->phys); if(!bc) continue;
		const shared_ptr<Body>& b1=Body::byId(i->getId1(),rb), b2=Body::byId(i->getId2(),rb);
		bool isIn1=isInBB(b1->state->pos,bbMin,bbMax), isIn2=isInBB(b2->state->pos,bbMin,bbMax);
		if(!isIn1 && !isIn2) continue;
		Real weight=1.;
		if((!isIn1 && isIn2) || (isIn1 && !isIn2)){
			//shared_ptr<Body> bIn=isIn1?b1:b2, bOut=isIn2?b2:b1;
			Vector3r vIn=(isIn1?b1:b2)->state->pos, vOut=(isIn2?b1:b2)->state->pos;
			#define _WEIGHT_COMPONENT(axis) if(vOut[axis]<bbMin[axis]) weight=min(weight,abs((vOut[axis]-bbMin[axis])/(vOut[axis]-vIn[axis]))); else if(vOut[axis]>bbMax[axis]) weight=min(weight,abs((vOut[axis]-bbMax[axis])/(vOut[axis]-vIn[axis])));
			_WEIGHT_COMPONENT(0); _WEIGHT_COMPONENT(1); _WEIGHT_COMPONENT(2);
			assert(weight>=0 && weight<=1);
		}
		E+=.5*bc->E*bc->crossSection*pow(bc->epsN,2)+.5*bc->G*bc->crossSection*pow(bc->epsT.norm(),2);
	}
	return E/(box[0]*box[1]*box[2]);
}
#endif

/*! Set velocity of all dynamic particles so that if their motion were unconstrained,
 * axis given by axisPoint and axisDirection would be reached in timeToAxis
 * (or, point at distance subtractDist from axis would be reached).
 *
 * The code is analogous to AxialGravityEngine and is intended to give initial motion
 * to particles subject to axial compaction to speed up the process. */
void velocityTowardsAxis(const Vector3r& axisPoint, const Vector3r& axisDirection, Real timeToAxis, Real subtractDist=0., Real perturbation=0.1){
	FOREACH(const shared_ptr<Body>&b, *(Omega::instance().getScene()->bodies)){
		if(!b->isDynamic()) continue;
		const Vector3r& x0=b->state->pos;
		const Vector3r& x1=axisPoint;
		const Vector3r x2=axisPoint+axisDirection;
		Vector3r closestAxisPoint=(x2-x1) * /* t */ (-(x1-x0).dot(x2-x1))/((x2-x1).squaredNorm());
		Vector3r toAxis=closestAxisPoint-x0;
		if(subtractDist>0) toAxis*=(toAxis.norm()-subtractDist)/toAxis.norm();
		b->state->vel=toAxis/timeToAxis;
		Vector3r ppDiff=perturbation*(1./sqrt(3.))*Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom())*b->state->vel.norm();
		b->state->vel+=ppDiff;
	}
}
BOOST_PYTHON_FUNCTION_OVERLOADS(velocityTowardsAxis_overloads,velocityTowardsAxis,3,5);

void particleConfinement(){
	CpmStateUpdater csu; csu.update();
}

// makes linker error out with monolithic build..
#if 0
python::dict testNumpy(){
	Scene* scene=Omega::instance().getScene().get();
	int dim1[]={scene->bodies->size()};
	int dim2[]={scene->bodies->size(),3};
	numpy_boost<Real,1> mass(dim1);
	numpy_boost<Real,2> vel(dim2);
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if(!b) continue;
		mass[b->getId()]=b->state->mass;
		VECTOR3R_TO_NUMPY(vel[b->getId()],b->state->vel);
	}
	py::dict ret;
	ret["mass"]=mass;
	ret["vel"]=vel;
	return ret;
}
#endif

#if 0
/* Compute stress tensor on each particle */
void particleMacroStress(void){
	Scene* scene=Omega::instance().getScene().get();
	// find interactions of each particle
	std::vector<std::list<Body::id_t> > bIntr(scene->bodies->size());
	FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
		if(!i->isReal) continue;
		// only contacts between 2 spheres
		Sphere* s1=dynamic_cast<Sphere*>(Body::byId(i->getId1(),scene)->shape.get())
		Sphere* s2=dynamic_cast<Sphere*>(Body::byId(i->getId2(),scene)->shape.get())
		if(!s1 || !s2) continue;
		bIntr[i.getId1()].push_back(i.getId2());
		bIntr[i.getId2()].push_back(i.getId1());
	}
	for(Body::id_t id1=0; id1<(Body::id_t)bIntr.size(); id1++){
		if(bIntr[id1].size()==0) continue;
		Matrix3r ss(Matrix3r::Zero()); // stress tensor on current particle
		FOREACH(Body::id_t id2, bIntr[id1]){
			const shared_ptr<Interaction> i=scene->interactions->find(id1,id2);
			assert(i);
			Dem3DofGeom* geom=YADE_CAST<Dem3DofGeom*>(i->geom);
			CpmPhys* phys=YADE_CAST<CpmPhys*>(i->phys);
			Real d=(geom->se31->pos-geom->se32->pos).norm(); // current contact length
			const Vector3r& n=geom->normal;
			const Real& A=phys->crossSection;
			const Vector3r& sigmaT=phys->sigmaT;
			const Real& sigmaN=phys->sigmaN;
			for(int i=0; i<3; i++) for(int j=0;j<3; j++){
				ss[i][j]=d*A*(sigmaN*n[i]*n[j]+.5*(sigmaT[i]*n[j]+sigmaT[j]*n[i]));
			}
		}
		// divide ss by V of the particle
		// FIXME: for now, just use 2*(4/3)*π*r³ (.5 porosity)
		ss*=1/(2*(4/3)*Mathr::PI*);
	}
}
#endif
#include<yade/lib-smoothing/WeightedAverage2d.hpp>
/* Fastly locate interactions within given distance from a point in 2d (projected to plane) */
struct HelixInteractionLocator2d{
	struct FlatInteraction{ Real r,h,theta; shared_ptr<Interaction> i; FlatInteraction(Real _r, Real _h, Real _theta, const shared_ptr<Interaction>& _i): r(_r), h(_h), theta(_theta), i(_i){}; };
	shared_ptr<GridContainer<FlatInteraction> > grid;
	Real thetaSpan;
	int axis;
	HelixInteractionLocator2d(Real dH_dTheta, int _axis, Real periodStart, Real theta0, Real thetaMin, Real thetaMax): axis(_axis){
		Scene* scene=Omega::instance().getScene().get();
		Real inf=std::numeric_limits<Real>::infinity();
		Vector2r lo=Vector2r(inf,inf), hi(-inf,-inf);
		Real minD0(inf),maxD0(-inf);
		Real minTheta(inf), maxTheta(-inf);
		std::list<FlatInteraction> intrs;
		// first pass: find extrema for positions and interaction lengths, build interaction list
		FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
			Dem3DofGeom* ge=dynamic_cast<Dem3DofGeom*>(i->geom.get());
			CpmPhys* ph=dynamic_cast<CpmPhys*>(i->phys.get());
			if(!ge || !ph) continue;
			Real r,h,theta;
			boost::tie(r,h,theta)=Shop::spiralProject(ge->contactPoint,dH_dTheta,axis,periodStart,theta0);
			if(!isnan(thetaMin) && theta<thetaMin) continue;
			if(!isnan(thetaMax) && theta>thetaMax) continue;
			lo=lo.cwise().min(Vector2r(r,h)); hi=hi.cwise().max(Vector2r(r,h));
			minD0=min(minD0,ge->refLength); maxD0=max(maxD0,ge->refLength);
			minTheta=min(minTheta,theta); maxTheta=max(maxTheta,theta);
			intrs.push_back(FlatInteraction(r,h,theta,i));
		}
		// create grid, put interactions inside
		Vector2i nCells=Vector2i(max(10,(int)((hi[0]-lo[0])/(2*minD0))),max(10,(int)((hi[1]-lo[1])/(2*minD0))));
		Vector2r hair=1e-2*Vector2r((hi[0]-lo[0])/nCells[0],(hi[1]-lo[1])/nCells[1]); // avoid rounding issue on boundary, enlarge by 1/100 cell size on each side
		grid=shared_ptr<GridContainer<FlatInteraction> >(new GridContainer<FlatInteraction>(lo-hair,hi+hair,nCells));
		FOREACH(const FlatInteraction& fi, intrs){
			grid->add(fi,Vector2r(fi.r,fi.h));
		}
		thetaSpan=maxTheta-minTheta;
	}
	py::list intrsAroundPt(const Vector2r& pt, Real radius){
		py::list ret;
		FOREACH(const Vector2i& v, grid->circleFilter(pt,radius)){
			FOREACH(const FlatInteraction& fi, grid->grid[v[0]][v[1]]){
				if((pow(fi.r-pt[0],2)+pow(fi.h-pt[1],2))>radius*radius) continue; // too far
				ret.append(fi.i);
			}
		}
		return ret;
	}
	// return macroscopic stress around interactions that project around given point and their average omega
	// stresses are rotated around axis back by theta angle
	python::tuple macroAroundPt(const Vector2r& pt, Real radius){
		Matrix3r ss(Matrix3r::Zero());
		Real omegaCumm=0, kappaCumm=0; int nIntr=0;
		FOREACH(const Vector2i& v, grid->circleFilter(pt,radius)){
			FOREACH(const FlatInteraction& fi, grid->grid[v[0]][v[1]]){
				if((pow(fi.r-pt[0],2)+pow(fi.h-pt[1],2))>radius*radius) continue; // too far
				Dem3DofGeom* geom=YADE_CAST<Dem3DofGeom*>(fi.i->geom.get());
				CpmPhys* phys=YADE_CAST<CpmPhys*>(fi.i->phys.get());
				// transformation matrix, to rotate to the plane
				Vector3r ax(Vector3r::Zero()); ax[axis]=1.;
				Quaternionr q(AngleAxisr(fi.theta,ax)); q=q.conjugate();
				Matrix3r TT=q.toRotationMatrix();
				//
				Real d=(geom->se31.position-geom->se32.position).norm(); // current contact length
				const Vector3r& n=TT*geom->normal;
				const Real& A=phys->crossSection;
				const Vector3r& sigmaT=TT*phys->sigmaT;
				const Real& sigmaN=phys->sigmaN;
				for(int i=0; i<3; i++) for(int j=0;j<3; j++){
					ss(i,j)+=d*A*(sigmaN*n[i]*n[j]+.5*(sigmaT[i]*n[j]+sigmaT[j]*n[i]));
				}
				omegaCumm+=phys->omega; kappaCumm+=phys->kappaD;
				nIntr++;
			}
		}
		// divide by approx spatial volume over which we averaged:
		// spiral cylinder with two half-spherical caps at ends
		ss*=1/((4/3.)*Mathr::PI*pow(radius,3)+Mathr::PI*pow(radius,2)*(thetaSpan*pt[0]-2*radius)); 
		return python::make_tuple(nIntr,ss,omegaCumm/nIntr,kappaCumm/nIntr);
	}
	Vector2r getLo(){ return grid->getLo(); }
	Vector2r getHi(){ return grid->getHi(); }

};

#ifdef YADE_VTK

/* Fastly locate interactions within given distance from a given point. See python docs for details. */
class InteractionLocator{
	// object doing the work, see http://www.vtk.org/doc/release/5.2/html/a01048.html
	vtkPointLocator *locator;
	// used by locator
	vtkUnstructuredGrid* grid;
	vtkPoints* points;
	// maps vtk's id to Interaction objects
	vector<shared_ptr<Interaction> > intrs;
	// axis-aligned bounds of all interactions
	Vector3r mn,mx;
	// count of interactions we hold
	int cnt;
	public:
	InteractionLocator(){
		// traverse all real interactions in the current simulation
		// add them to points, create grid with those points
		// store shared_ptr's to interactions in intrs separately
		Scene* scene=Omega::instance().getScene().get();
		locator=vtkPointLocator::New();
		points=vtkPoints::New();
		int count=0;
		FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
			if(!i->isReal()) continue;
			Dem3DofGeom* ge=dynamic_cast<Dem3DofGeom*>(i->geom.get());
			if(!ge) continue;
			const Vector3r& cp(ge->contactPoint);
			int id=points->InsertNextPoint(cp[0],cp[1],cp[2]);
			if(intrs.size()<=(size_t)id){intrs.resize(id+1000);}
			intrs[id]=i;
			count++;
		}
		if(count==0) throw std::runtime_error("Zero interactions when constructing InteractionLocator!");
		double bb[6];
		points->ComputeBounds(); points->GetBounds(bb);
		mn=Vector3r(bb[0],bb[2],bb[4]); mx=Vector3r(bb[1],bb[3],bb[5]);
		cnt=points->GetNumberOfPoints();

		grid=vtkUnstructuredGrid::New();
		grid->SetPoints(points);
		locator->SetDataSet(grid);
		locator->BuildLocator();
	}
	// cleanup
	~InteractionLocator(){ locator->Delete(); points->Delete(); grid->Delete(); }

	py::list intrsAroundPt(const Vector3r& pt, Real radius){
		vtkIdList *ids=vtkIdList::New();
		locator->FindPointsWithinRadius(radius,(const double*)(&pt),ids);
		int numIds=ids->GetNumberOfIds();
		py::list ret;
		for(int i=0; i<numIds; i++){
			// LOG_TRACE("Add "<<i<<"th id "<<ids->GetId(i));
			ret.append(intrs[ids->GetId(i)]);
		}
		return ret;
	}
	python::tuple macroAroundPt(const Vector3r& pt, Real radius, Real forceVolume=-1){
		Matrix3r ss(Matrix3r::Zero());
		vtkIdList *ids=vtkIdList::New();
		locator->FindPointsWithinRadius(radius,(const double*)(&pt),ids);
		int numIds=ids->GetNumberOfIds();
		Real omegaCumm=0, kappaCumm=0;
		for(int k=0; k<numIds; k++){
			const shared_ptr<Interaction>& I(intrs[ids->GetId(k)]);
			Dem3DofGeom* geom=YADE_CAST<Dem3DofGeom*>(I->geom.get());
			CpmPhys* phys=YADE_CAST<CpmPhys*>(I->phys.get());
			Real d=(geom->se31.position-geom->se32.position).norm(); // current contact length
			const Vector3r& n=geom->normal;
			const Real& A=phys->crossSection;
			const Vector3r& sigmaT=phys->sigmaT;
			const Real& sigmaN=phys->sigmaN;
			for(int i=0; i<3; i++) for(int j=0;j<3; j++){
				ss(i,j)+=d*A*(sigmaN*n[i]*n[j]+.5*(sigmaT[i]*n[j]+sigmaT[j]*n[i]));
			}
			omegaCumm+=phys->omega; kappaCumm+=phys->kappaD;
		}
		Real volume=(forceVolume>0?forceVolume:(4/3.)*Mathr::PI*pow(radius,3));
		ss*=1/volume;
		return py::make_tuple(ss,omegaCumm/numIds,kappaCumm/numIds);
	}
	py::tuple getBounds(){ return py::make_tuple(mn,mx);}
	int getCnt(){ return cnt; }
};
#endif

BOOST_PYTHON_MODULE(_eudoxos){
	YADE_SET_DOCSTRING_OPTS;
	py::def("velocityTowardsAxis",velocityTowardsAxis,velocityTowardsAxis_overloads(py::args("axisPoint","axisDirection","timeToAxis","subtractDist","perturbation")));
	// def("spiralSphereStresses2d",spiralSphereStresses2d,(python::arg("dH_dTheta"),python::arg("axis")=2));
	py::def("particleConfinement",particleConfinement);
	#if 0
		import_array();
		py::def("testNumpy",testNumpy);
	#endif
#ifdef YADE_VTK
	py::class_<InteractionLocator>("InteractionLocator","Locate all (real) interactions in space by their :yref:`contact point<Dem3DofGeom::contactPoint>`. When constructed, all real interactions are spatially indexed (uses `vtkPointLocator <http://www.vtk.org/doc/release/5.4/html/a01247.html>`_ internally). Use instance methods to use those data. \n\n.. note::\n\tData might become inconsistent with real simulation state if simulation is being run between creation of this object and spatial queries.")
		.def("intrsAroundPt",&InteractionLocator::intrsAroundPt,((python::arg("point"),python::arg("maxDist"))),"Return list of real interactions that are not further than *maxDist* from *point*.")
		.def("macroAroundPt",&InteractionLocator::macroAroundPt,((python::arg("point"),python::arg("maxDist"),python::arg("forceVolume")=-1)),"Return tuple of averaged stress tensor (as Matrix3), average omega and average kappa values. *forceVolume* can be used (if positive) rather than the sphere (with *maxDist* radius) volume for the computation. (This is useful if *point* and *maxDist* encompass empty space that you want to avoid.)")
		.add_property("bounds",&InteractionLocator::getBounds,"Return coordinates of lower and uppoer corner of axis-aligned abounding box of all interactions")
		.add_property("count",&InteractionLocator::getCnt,"Number of interactions held")
	;
#endif
	py::class_<HelixInteractionLocator2d>("HelixInteractionLocator2d",
		"Locate all real interactions in 2d plane (reduced by spiral projection from 3d, using ``Shop::spiralProject``, which is the same as :yref:`yade.utils.spiralProject`) using their :yref:`contact points<Dem3DofGeom::contactPoint>`. \n\n.. note::\n\tDo not run simulation while using this object.",
		python::init<Real,int,Real,Real,Real,Real>((python::arg("dH_dTheta"),python::arg("axis")=0,python::arg("periodStart")=NaN,python::arg("theta0")=0,python::arg("thetaMin")=NaN,python::arg("thetaMax")=NaN),":param float dH_dTheta: Spiral inclination, i.e. height increase per 1 radian turn;\n:param int axis: axis of rotation (0=x,1=y,2=z)\n:param float theta: spiral angle at zero height (theta intercept)\n:param float thetaMin: only interactions with $\\theta$≥\\ *thetaMin* will be considered (NaN to deactivate)\n:param float thetaMax: only interactions with $\\theta$≤\\ *thetaMax* will be considered (NaN to deactivate)\n\nSee :yref:`yade.utils.spiralProject`.")
	)
		.def("intrsAroundPt",&HelixInteractionLocator2d::intrsAroundPt,(python::arg("pt2d"),python::arg("radius")),"Return list of interaction objects that are not further from *pt2d* than *radius* in the projection plane")
		.def("macroAroundPt",&HelixInteractionLocator2d::macroAroundPt,(python::arg("pt2d"),python::arg("radius")),"Compute macroscopic stress around given point; the interaction ($n$ and $\\sigma^T$ are rotated to the projection plane by $\\theta$ (as given by :yref:`yade.utils.spiralProject`) first, but no skew is applied). The formula used is\n\n.. math::\n\n    \\sigma_{ij}=\\frac{1}{V}\\sum_{IJ}d^{IJ}A^{IJ}\\left[\\sigma^{N,IJ}n_i^{IJ}n_j^{IJ}+\\frac{1}{2}\\left(\\sigma_i^{T,IJ}n_j^{IJ}+\\sigma_j^{T,IJ}n_i^{IJ}\\right)\\right]\n\nwhere the sum is taken over volume $V$ containing interactions $IJ$ between spheres $I$ and $J$;\n\n* $i$, $j$ indices denote Cartesian components of vectors and tensors,\n* $d^{IJ}$ is current distance between spheres $I$ and $J$,\n* $A^{IJ}$ is area of contact $IJ$,\n* $n$ is ($\\theta$-rotated) interaction normal (unit vector pointing from center of $I$ to the center of $J$)\n* $\\sigma^{N,IJ}$  is normal stress (as scalar) in contact $IJ$,\n* $\\sigma^{T,IJ}$ is shear stress in contact $IJ$ in global coordinates and $\\theta$-rotated. \n\nAdditionally, computes average of :yref:`CpmPhys.omega` ($\\bar\\omega$) and :yref:`CpmPhys.kappaD` ($\\bar\\kappa_D$). *N* is the number of interactions in the volume given.\n\n:return: tuple of (*N*, $\\tens{\\sigma}$, $\\bar\\omega$, $\\bar\\kappa_D$).\n")
		.add_property("lo",&HelixInteractionLocator2d::getLo,"Return lower corner of the rectangle containing all interactions.")
		.add_property("hi",&HelixInteractionLocator2d::getHi,"Return upper corner of the rectangle containing all interactions.");

}

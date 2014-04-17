/*************************************************************************
*  Copyright (C) 2009 by Emanuele Catalano <catalano@grenoble-inp.fr>    *
*  Copyright (C) 2009 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*  Copyright (C) 2012 by Donia Marzougui <donia.marzougui@grenoble-inp.fr>*
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef FLOW_ENGINE

// #define XVIEW
#include "FlowBoundingSphere.hpp"//include after #define XVIEW
#include <iostream>
#include <fstream>
#include <new>
#include <utility>
#include "vector"
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>


#ifdef XVIEW
// #include "Vue3D.h" //FIXME implicit dependencies will look for this class (out of tree) even ifndef XVIEW
#endif

#ifdef YADE_OPENMP
//   #define GS_OPEN_MP //It should never be defined if Yade is not using openmp
#endif

// #define USE_FAST_MATH 1


using namespace std;
namespace CGT
{

typedef vector<double> VectorR;

//! Use this factor, or minLength, to reduce max permeability values (see usage below))
const double minLength = 0.01;//percentage of mean rad

//! Factors including the effect of 1/2 symmetry in hydraulic radii
const Real multSym1 = 1/pow(2,0.25);
const Real multSym2 = 1/pow(4,0.25);

#ifdef XVIEW
Vue3D Vue1;
#endif
template<class Tesselation>
FlowBoundingSphere<Tesselation>::~FlowBoundingSphere()
{
}
template <class Tesselation> 
FlowBoundingSphere<Tesselation>::FlowBoundingSphere()
{
	xMin = 1000.0, xMax = -10000.0, yMin = 1000.0, yMax = -10000.0, zMin = 1000.0, zMax = -10000.0;
	currentTes = 0;
	nOfSpheres = 0;
	sectionArea = 0, Height=0, vTotal=0;
	vtkInfiniteVertices=0, vtkInfiniteCells=0;
	viscosity = 1;
	fluidBulkModulus = 0;
	tessBasedForce = true;
	for (int i=0;i<6;i++) boundsIds[i] = 0;
	minPermLength=-1;
	slipBoundary = false;//no-slip/symmetry conditions on lateral boundaries
	tolerance = 1e-07;
	relax = 1.9;
	ks=0;
	distanceCorrection = true;
	clampKValues = true;
	meanKStat = true; KOptFactor=0;
	noCache=true;
	pressureChanged=false;
	computeAllCells=true;//might be turned false IF the code is reorganized (we can make a separate function to compute unitForceVectors outside compute_Permeability) AND it really matters for CPU time
	debugOut = true;
	rAverage = false; /** use the average between the effective radius (inscribed sphere in facet) and the equivalent (circle surface = facet fluid surface) **/
	OUTPUT_BOUDARIES_RADII = false;
	rAverage = false; /** if true use the average between the effective radius (inscribed sphere in facet) and the equivalent (circle surface = facet fluid surface) **/
// 	areaR2Permeability=true;
	permeabilityMap = false;
	computedOnce=false;
	minKdivKmean=0.0001;
	maxKdivKmean=100.;
	ompThreads=1;
	errorCode=0;
}

template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::resetNetwork() {T[currentTes].Clear();noCache=true;}

template <class Tesselation>
void FlowBoundingSphere<Tesselation>::averageRelativeCellVelocity()
{
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
        Point posAvFacet;
        int numCells = 0;
        double facetFlowRate = 0;
	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
        for ( FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++ ) {
		if (cell->info().isGhost) continue;
		cell->info().averageVelocity() =CGAL::NULL_VECTOR;
                numCells++;
		Real totFlowRate = 0;//used to acount for influxes in elements where pressure is imposed
                for ( int i=0; i<4; i++ ) if (!Tri.is_infinite(cell->neighbor(i))){
				CVector Surfk = cell->info()-cell->neighbor(i)->info();
				Real area = sqrt ( Surfk.squared_length() );
				Surfk = Surfk/area;
                        	CVector branch = cell->vertex ( facetVertices[i][0] )->point() - cell->info();
                        	posAvFacet = (Point) cell->info() + ( branch*Surfk ) *Surfk;
				facetFlowRate = (cell->info().kNorm())[i] * (cell->info().shiftedP() - cell->neighbor (i)->info().shiftedP());
				totFlowRate += facetFlowRate;
				cell->info().averageVelocity() = cell->info().averageVelocity() + (facetFlowRate) * ( posAvFacet-CGAL::ORIGIN );
		}
		//This is the influx term
		if (cell->info().Pcondition) cell->info().averageVelocity() = cell->info().averageVelocity() - (totFlowRate)*((Point) cell->info()-CGAL::ORIGIN );
		//now divide by volume
		cell->info().averageVelocity() = cell->info().averageVelocity() /abs(cell->info().volume());
	}
}



template <class Tesselation> 
bool FlowBoundingSphere<Tesselation>::isOnSolid  (double X, double Y, double Z)
{
  RTriangulation& Tri = T[currentTes].Triangulation();
  FiniteCellsIterator cellEnd = Tri.finiteCellsEnd();
  for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
    for (int i=0; i<4; i++){
      double radius = sqrt(cell->vertex(i)->point().weight());
      if (X < (cell->vertex(i)->point().x()+radius) && X > (cell->vertex(i)->point().x()-radius)){
	if (Y < (cell->vertex(i)->point().y()+radius) && Y > (cell->vertex(i)->point().y()-radius)){
	  if (Z < (cell->vertex(i)->point().z()+radius) && Z > (cell->vertex(i)->point().z()-radius)){
	    return true;}}}}}
      return false;
}
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::averageFluidVelocity()
{
	averageRelativeCellVelocity();
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
	int numVertex = 0;
	FiniteVerticesIterator verticesEnd = Tri.finite_vertices_end();
	for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt !=  verticesEnd; vIt++) {
	  numVertex++;}
	
	vector<Real> volumes;
	vector<CGT::CVector> velocityVolumes;
	velocityVolumes.resize(numVertex);
	volumes.resize(numVertex);
	
	for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt !=  verticesEnd; vIt++) {
	  velocityVolumes[vIt->info().id()]=CGAL::NULL_VECTOR;
	  volumes[vIt->info().id()]=0.f;}
	
	FiniteCellsIterator cellEnd = Tri.finiteCellsEnd();
	for ( FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++ )
	{
	  if (cell->info().fictious()==0){
	    for (int i=0;i<4;i++){
	      velocityVolumes[cell->vertex(i)->info().id()] =  velocityVolumes[cell->vertex(i)->info().id()] + cell->info().averageVelocity()*cell->info().volume();
	      volumes[cell->vertex(i)->info().id()] = volumes[cell->vertex(i)->info().id()] + cell->info().volume();}
	  }}	    
	
	std::ofstream fluid_vel ("Velocity", std::ios::out);
	double Rx = (xMax-xMin) /10;
        double Ry = (yMax-yMin) /12;
	double Rz = (zMax-zMin) /20;
	CellHandle cellula;
	
	CVector velocity = CGAL::NULL_VECTOR;
	int i=0;
	for(double X=xMin+Rx;X<xMax;X+=Rx){
	  for (double Y=yMin+Ry;Y<yMax;Y+=Ry){
	    velocity = CGAL::NULL_VECTOR; i=0;
	    for (double Z=zMin+Rz;Z<zMax;Z+=Rz){
	      cellula = Tri.locate(Point(X,Y,Z));
	      for (int y=0;y<4;y++) {if (!cellula->vertex(y)->info().isFictious) {velocity = velocity + (velocityVolumes[cellula->vertex(y)->info().id()]/volumes[cellula->vertex(y)->info().id()]);i++;}}
	    }velocity = velocity/i;
	    fluid_vel << X << " " << Y << " " << velocity << endl;
	  }}
}
template <class Tesselation> 
vector<Real> FlowBoundingSphere<Tesselation>::averageFluidVelocityOnSphere(unsigned int Id_sph)
{
	averageRelativeCellVelocity();
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();	
	Real volumes; CGT::CVector velocityVolumes;
	vector<Real> result;
	result.resize(3);	
	velocityVolumes=CGAL::NULL_VECTOR;
	volumes=0.f;
	
	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
	for ( FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++ )
	{
	  if (cell->info().fictious()==0){
	    for (unsigned int i=0;i<4;i++){
	      if (cell->vertex(i)->info().id()==Id_sph){
		velocityVolumes = velocityVolumes + cell->info().averageVelocity()*cell->info().volume();
		volumes = volumes + cell->info().volume();}}}}
		
	for (int i=0;i<3;i++) result[i] += velocityVolumes[i]/volumes;
	return result;
}
template <class Tesselation> 
double FlowBoundingSphere<Tesselation>::getPorePressure (double X, double Y, double Z)
{
	if (noCache && T[!currentTes].Max_id()<=0) return 0;//the engine never solved anything
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
	CellHandle cell = Tri.locate(Point(X,Y,Z));
	return cell->info().p();
}

template <class Tesselation>
double FlowBoundingSphere<Tesselation>::getCell (double X, double Y, double Z)
{
	if (noCache) {cerr<<"Triangulation does not exist. Waht did you do?!"<<endl; return -1;}
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
	CellHandle cell = Tri.locate(Point(X,Y,Z));
	return cell->info().id;
}

template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::measurePressureProfile(double WallUpy, double WallDowny)
{  
	if (noCache && T[!currentTes].Max_id()<=0) return;//the engine never solved anything
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
        CellHandle permeameter;
	std::ofstream capture ("Pressure_profile", std::ios::app);
        int intervals = 5;
	int captures = 6;
        double Rz = (zMax-zMin)/intervals;
	double Ry = (WallUpy-WallDowny)/captures;
	double X=(xMax+xMin)/2;
	double Y = WallDowny;
	double pressure = 0.f;
	int cell=0;
	for (int i=0; i<captures; i++){
        for (double Z=min(zMin,zMax); Z<=max(zMin,zMax); Z+=abs(Rz)) {
		permeameter = Tri.locate(Point(X, Y, Z));
		pressure+=permeameter->info().p();
		cell++;
        }
        Y += Ry;
        capture  << pressure/cell << endl;}
	
}
template <class Tesselation> 
double FlowBoundingSphere<Tesselation>::averageSlicePressure(double Y)
{
  RTriangulation& Tri = T[currentTes].Triangulation();
  double P_ave = 0.f;
  int n = 0;
  double Ry = (yMax-yMin)/30;
  double Rx = (xMax-xMin)/30;
  double Rz = (zMax-zMin)/30;
  for (double X=xMin; X<=xMax+Ry/10; X=X+Rx) {
	for (double Z=zMin; Z<=zMax+Ry/10; Z=Z+Rz) {
	  P_ave+=Tri.locate(Point(X, Y, Z))->info().p();
	  n++;
	}
  }
  P_ave/=n;
  return P_ave;
}
template <class Tesselation> 
double FlowBoundingSphere<Tesselation>::averagePressure()
{
  RTriangulation& Tri = T[currentTes].Triangulation();
  double P = 0.f, Ppond=0.f, Vpond=0.f;
  int n = 0;
  for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); cell++) {
	P+=cell->info().p();
	n++;
	Ppond+=cell->info().p()*cell->info().volume();
	Vpond+=cell->info().volume();
  }
  P/=n;
  Ppond/=Vpond;
  return Ppond;
}


template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::computeFacetForcesWithCache(bool onlyCache)
{
	RTriangulation& Tri = T[currentTes].Triangulation();
	CVector nullVect(0,0,0);
	//reset forces
	if (!onlyCache) for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) v->info().forces=nullVect;

	#ifdef parallel_forces
	if (noCache) {
		perVertexUnitForce.clear(); perVertexPressure.clear();
		perVertexUnitForce.resize(T[currentTes].maxId+1);
		perVertexPressure.resize(T[currentTes].maxId+1);}
	#endif
	CellHandle neighbourCell;
	VertexHandle mirrorVertex;
	CVector tempVect;
	//FIXME : Ema, be carefull with this (noCache), it needs to be turned true after retriangulation
	if (noCache) {for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
			CellHandle& cell = *cellIt;
			//reset cache
			for (int k=0;k<4;k++) cell->info().unitForceVectors[k]=nullVect;

			for (int j=0; j<4; j++) if (!Tri.is_infinite(cell->neighbor(j))) {
					neighbourCell = cell->neighbor(j);
					const CVector& Surfk = cell->info().facetSurfaces[j];
					//FIXME : later compute that fluidSurf only once in hydraulicRadius, for now keep full surface not modified in cell->info for comparison with other forces schemes
					//The ratio void surface / facet surface
					Real area = sqrt(Surfk.squared_length()); if (area<=0) cerr <<"AREA <= 0!!"<<endl;
					CVector facetNormal = Surfk/area;
					const std::vector<CVector>& crossSections = cell->info().facetSphereCrossSections;
					CVector fluidSurfk = cell->info().facetSurfaces[j]*cell->info().facetFluidSurfacesRatio[j];
					/// handle fictious vertex since we can get the projected surface easily here
					if (cell->vertex(j)->info().isFictious) {
						Real projSurf=abs(Surfk[boundary(cell->vertex(j)->info().id()).coordinate]);
						tempVect=-projSurf*boundary(cell->vertex(j)->info().id()).normal;
						cell->vertex(j)->info().forces = cell->vertex(j)->info().forces+tempVect*cell->info().p();
						//define the cached value for later use with cache*p
						cell->info().unitForceVectors[j]=cell->info().unitForceVectors[j]+ tempVect;
					}
					/// Apply weighted forces f_k=sqRad_k/sumSqRad*f
					CVector facetUnitForce = -fluidSurfk*cell->info().solidSurfaces[j][3];
					CVector facetForce = cell->info().p()*facetUnitForce;
										
					for (int y=0; y<3;y++) {
						cell->vertex(facetVertices[j][y])->info().forces = cell->vertex(facetVertices[j][y])->info().forces + facetForce*cell->info().solidSurfaces[j][y];
						//add to cached value
						cell->info().unitForceVectors[facetVertices[j][y]]=cell->info().unitForceVectors[facetVertices[j][y]]+facetUnitForce*cell->info().solidSurfaces[j][y];
						//uncomment to get total force / comment to get only viscous forces (Bruno)
						if (!cell->vertex(facetVertices[j][y])->info().isFictious) {
							cell->vertex(facetVertices[j][y])->info().forces = cell->vertex(facetVertices[j][y])->info().forces -facetNormal*cell->info().p()*crossSections[j][y];
							//add to cached value
							cell->info().unitForceVectors[facetVertices[j][y]]=cell->info().unitForceVectors[facetVertices[j][y]]-facetNormal*crossSections[j][y];
						}
					}
					#ifdef parallel_forces
					perVertexUnitForce[cell->vertex(j)->info().id()].push_back(&(cell->info().unitForceVectors[j]));
					perVertexPressure[cell->vertex(j)->info().id()].push_back(&(cell->info().p()));
					#endif
			}
		}
		noCache=false;//cache should always be defined after execution of this function
		if (onlyCache) return;
	} else {//use cached values
		#ifndef parallel_forces
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			for (int yy=0;yy<4;yy++) cell->vertex(yy)->info().forces = cell->vertex(yy)->info().forces + cell->info().unitForceVectors[yy]*cell->info().p();}
			
		#else
		#pragma omp parallel for num_threads(ompThreads)
		for (int vn=0; vn<= T[currentTes].maxId; vn++) {
			if (T[currentTes].vertexHandles[vn]==NULL) continue;
			VertexHandle& v = T[currentTes].vertexHandles[vn];
			const int& id =  v->info().id();
			CVector tf (0,0,0);
			int k=0;
			for (vector<const Real*>::iterator c = perVertexPressure[id].begin(); c != perVertexPressure[id].end(); c++)
				tf = tf + (*(perVertexUnitForce[id][k++]))*(**c);
			v->info().forces = tf;
		}
		#endif
	}
	if (debugOut) {
		CVector totalForce = nullVect;
		for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v)	{
			if (!v->info().isFictious) totalForce = totalForce + v->info().forces;
			else if (boundary(v->info().id()).flowCondition==1) totalForce = totalForce + v->info().forces;	}
		cout << "totalForce = "<< totalForce << endl;}
}

template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::applySinusoidalPressure(RTriangulation& Tri, double amplitude, double averagePressure, double loadIntervals)
{
	double step = 1/loadIntervals;
	VectorCell tmpCells;
	tmpCells.resize(10000);
	VCellIterator cellsIt = tmpCells.begin();
	for (double alpha=0; alpha<1.001; alpha+=step)
	{
	  VCellIterator cellsEnd = Tri.incident_cells(T[currentTes].vertexHandles[yMaxId],cellsIt);
	  for (VCellIterator it = tmpCells.begin(); it != cellsEnd; it++)
	  {
	    if(!Tri.is_infinite(*it)){
	      Point& p1 = (*it)->info();
	      CellHandle& cell = *it;
	      if (p1.x()<xMin) cell->info().p() = averagePressure+amplitude;
	      else if (p1.x()>xMax) cell->info().p() = averagePressure-amplitude;
	      else if (p1.x()>(xMin+alpha*(xMax-xMin)) && p1.x()<(xMin+(alpha+step)*(xMax-xMin))) cell->info().p() = averagePressure + (amplitude)*(cos(alpha*M_PI));
	  }
	  }
	}
}
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::interpolate(Tesselation& Tes, Tesselation& NewTes)
{
        CellHandle oldCell;
        RTriangulation& Tri = Tes.Triangulation();
	for (typename VectorCell::iterator cellIt=NewTes.cellHandles.begin(); cellIt!=NewTes.cellHandles.end(); cellIt++){
		CellHandle& newCell = *cellIt;
		if (newCell->info().Pcondition || newCell->info().isGhost) continue;
		CVector center ( 0,0,0 );
		if (newCell->info().fictious()==0) for ( int k=0;k<4;k++ ) center= center + 0.25* (Tes.vertex(newCell->vertex(k)->info().id())->point()-CGAL::ORIGIN);
		else {
			Real boundPos=0; int coord=0;
			for ( int k=0;k<4;k++ ) {
				if (!newCell->vertex (k)->info().isFictious) center= center+0.3333333333*(Tes.vertex(newCell->vertex(k)->info().id())->point()-CGAL::ORIGIN);
				else {
					coord=boundary (newCell->vertex(k)->info().id()).coordinate;
					boundPos=boundary (newCell->vertex(k)->info().id()).p[coord];
				}
			}
			center=CVector(coord==0?boundPos:center[0],coord==1?boundPos:center[1],coord==2?boundPos:center[2]);
		}
                oldCell = Tri.locate(Point(center[0],center[1],center[2]));
		newCell->info().getInfo(oldCell->info());
//                 newCell->info().p() = oldCell->info().shiftedP();
        }
//  	Tes.Clear();//Don't reset to avoid segfault when getting pressure in scripts just after interpolation
}

template <class Tesselation> 
Real FlowBoundingSphere<Tesselation>::checkSphereFacetOverlap(const Sphere& v0, const Sphere& v1, const Sphere& v2)
{
	//First, check that v0 projection fall between v1 and v2...
	Real dist=(v0-v1)*(v2-v1);
	if (dist<0) return 0;
	Real v1v2=(v2-v1).squared_length();
	if (dist>v1v2) return 0;
	//... then, check distance
	Real m=(cross_product(v0-v1,v2-v1)).squared_length()/v1v2;
	if (m<v0.weight()) {
		Real d=2*sqrt((v0.weight()-m));
		Real teta=2*acos(sqrt(m/v0.weight()));
		return 0.5*(teta*v0.weight()-d*sqrt(m));//this is S0, we use crossSection to avoid computing an "asin"
// 		return crossSection-m*d;
	} else return 0;
}

template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::computePermeability()
{
	if (debugOut)  cout << "----Computing_Permeability------" << endl;
	RTriangulation& Tri = T[currentTes].Triangulation();
	VSolidTot = 0, Vtotalissimo = 0, vPoral = 0, sSolidTot = 0, vTotalPorosity=0, vPoralPorosity=0;
	FiniteCellsIterator cellEnd = Tri.finite_cells_end();

	CellHandle neighbourCell;

	double k=0, distance = 0, radius = 0;
	int surfneg=0;
	int NEG=0, POS=0, pass=0;

	bool ref = Tri.finite_cells_begin()->info().isvisited;
	Real meanK=0, STDEV=0, meanRadius=0, meanDistance=0;
	Real infiniteK=1e10;

	for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
		CellHandle& cell = *cellIt;
		Point& p1 = cell->info();
		for (int j=0; j<4; j++) {
			neighbourCell = cell->neighbor(j);
			Point& p2 = neighbourCell->info();
			if (!Tri.is_infinite(neighbourCell) && (neighbourCell->info().isvisited==ref || computeAllCells)) {
				//compute and store the area of sphere-facet intersections for later use
				VertexHandle W [3];
				for (int kk=0; kk<3; kk++) {
					W[kk] = cell->vertex(facetVertices[j][kk]);
				}
				Sphere& v0 = W[0]->point();
				Sphere& v1 = W[1]->point();
				Sphere& v2 = W[2]->point();

				cell->info().facetSphereCrossSections[j]=CVector(
				   W[0]->info().isFictious ? 0 : 0.5*v0.weight()*acos((v1-v0)*(v2-v0)/sqrt((v1-v0).squared_length()*(v2-v0).squared_length())),
				   W[1]->info().isFictious ? 0 : 0.5*v1.weight()*acos((v0-v1)*(v2-v1)/sqrt((v1-v0).squared_length()*(v2-v1).squared_length())),
				   W[2]->info().isFictious ? 0 : 0.5*v2.weight()*acos((v0-v2)*(v1-v2)/sqrt((v1-v2).squared_length()*(v2-v0).squared_length())));

				pass+=1;
				CVector l = p1 - p2;
				distance = sqrt(l.squared_length());
				if (!rAverage) radius = 2* computeHydraulicRadius(cell, j);
				else radius = (computeEffectiveRadius(cell, j)+computeEquivalentRadius(cell,j))*0.5;
				if (radius<0) NEG++;
				else POS++;
				if (radius==0) {
					cout << "INS-INS PROBLEM!!!!!!!" << endl;
				}
				Real fluidArea=0;
				if (distance!=0) {
					if (minPermLength>0 && distanceCorrection) distance=max(minPermLength,distance);
					const CVector& Surfk = cell->info().facetSurfaces[j];
					Real area = sqrt(Surfk.squared_length());
					const CVector& crossSections = cell->info().facetSphereCrossSections[j];
					Real S0=0;
					S0=checkSphereFacetOverlap(v0,v1,v2);
					if (S0==0) S0=checkSphereFacetOverlap(v1,v2,v0);
					if (S0==0) S0=checkSphereFacetOverlap(v2,v0,v1);
					//take absolute value, since in rare cases the surface can be negative (overlaping spheres)
					fluidArea=abs(area-crossSections[0]-crossSections[1]-crossSections[2]+S0);
					cell->info().facetFluidSurfacesRatio[j]=fluidArea/area;
					k=(fluidArea * pow(radius,2)) / (8*viscosity*distance);
					 meanDistance += distance;
					 meanRadius += radius;
					 meanK +=  k*kFactor;

				if (k<0 && debugOut) {surfneg+=1;
				cout<<"__ k<0 __"<<k<<" "<<" fluidArea "<<fluidArea<<" area "<<area<<" "<<crossSections[0]<<" "<<crossSections[1]<<" "<<crossSections[2] <<" "<<W[0]->info().id()<<" "<<W[1]->info().id()<<" "<<W[2]->info().id()<<" "<<p1<<" "<<p2<<" test "<<endl;}				     
				} else  {cout <<"infinite K1!"<<endl; k = infiniteK;}//Will be corrected in the next loop

				(cell->info().kNorm())[j]= k*kFactor;
				if (!neighbourCell->info().isGhost) (neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]= (cell->info().kNorm())[j];
			}
		}
		cell->info().isvisited = !ref;
	}
	if (debugOut) cout<<"surfneg est "<<surfneg<<endl;
	meanK /= pass;
	meanRadius /= pass;
	meanDistance /= pass;
	Real globalK=kFactor*meanDistance*vPoral/(sSolidTot*8.*viscosity);//An approximate value of macroscopic permeability, for clamping local values below
	if (debugOut) {
		cout << "PassCompK = " << pass << endl;
		cout << "meanK = " << meanK << endl;
		cout << "globalK = " << globalK << endl;
		cout << "maxKdivKmean*globalK = " << maxKdivKmean*globalK << endl;
		cout << "minKdivKmean*globalK = " << minKdivKmean*globalK << endl;
		cout << "meanTubesRadius = " << meanRadius << endl;
		cout << "meanDistance = " << meanDistance << endl;
	}
	ref = Tri.finite_cells_begin()->info().isvisited;
	pass=0;

	if (clampKValues) for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
		CellHandle& cell = *cellIt;
		for (int j=0; j<4; j++) {
			neighbourCell = cell->neighbor(j);
			if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited==ref) {
				pass++;
				(cell->info().kNorm())[j] = max(minKdivKmean*globalK ,min((cell->info().kNorm())[j], maxKdivKmean*globalK));
				(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]=(cell->info().kNorm())[j];
			}
		}
	}
	if (debugOut) cout << "PassKcorrect = " << pass << endl;
	if (debugOut) cout << "POS = " << POS << " NEG = " << NEG << " pass = " << pass << endl;

	// A loop to compute the standard deviation of the local K distribution, and use it to include/exclude K values higher then (meanK +/- K_opt_factor*STDEV)
	if (meanKStat)
	{
		std::ofstream k_opt_file("k_stdev.txt" ,std::ios::out);
		ref = Tri.finite_cells_begin()->info().isvisited;
		pass=0;
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			for (int j=0; j<4; j++) {
				neighbourCell = cell->neighbor(j);
				if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited==ref) {
					pass++;
					STDEV += pow(((cell->info().kNorm())[j]-meanK),2);
				}
			}cell->info().isvisited = !ref;
		}
		STDEV = sqrt(STDEV/pass);
		if (debugOut) cout << "PassSTDEV = " << pass << endl << "STATISTIC K" << endl;
		double k_min = 0, k_max = meanK + KOptFactor*STDEV;
		cout << "Kmoy = " << meanK << " Standard Deviation = " << STDEV << endl<< "kmin = " << k_min << " kmax = " << k_max << endl;
		ref = Tri.finite_cells_begin()->info().isvisited;
		pass=0;
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			for (int j=0; j<4; j++) {
				neighbourCell = cell->neighbor(j);
				if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited==ref) {
					pass+=1;
					if ((cell->info().kNorm())[j]>k_max) {
						(cell->info().kNorm())[j]=k_max;
						(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]= (cell->info().kNorm())[j];
					}
					k_opt_file << KOptFactor << " " << (cell->info().kNorm())[j] << endl;
				}
			}cell->info().isvisited=!ref;
		}
		if (debugOut) cout << "PassKopt = " << pass << endl;
	}
	if (debugOut) {
		FiniteVerticesIterator verticesEnd = Tri.finite_vertices_end();
		Real Vgrains = 0;
		int grains=0;
		for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt !=  verticesEnd; vIt++) {
			if (!vIt->info().isFictious && !vIt->info().isGhost) {
				grains +=1;
				Vgrains += 1.33333333 * M_PI * pow(vIt->point().weight(),1.5);}}
		cout<<grains<<"grains - " <<"vTotal = " << vTotal << " Vgrains = " << Vgrains << " vPoral1 = " << (vTotal-Vgrains) << endl;
		cout << "Vtotalissimo = " << Vtotalissimo/2 << " VSolidTot = " << VSolidTot/2 << " vPoral2 = " << vPoral/2  << " sSolidTot = " << sSolidTot << endl<< endl;
		if (!rAverage) cout << "------Hydraulic Radius is used for permeability computation------" << endl << endl;
		else cout << "------Average Radius is used for permeability computation------" << endl << endl;
		cout << "-----computed_Permeability-----" << endl;}
}

template <class Tesselation> 
vector<double> FlowBoundingSphere<Tesselation>::getConstrictions()
{
	RTriangulation& Tri = T[currentTes].Triangulation();
	vector<double> constrictions;
	for (FiniteFacetsIterator f_it=Tri.finite_facets_begin(); f_it != Tri.finite_facets_end();f_it++){
		//in the periodic case, we skip facets with lowest id out of the base period
		if ( ((f_it->first->info().index <= f_it->first->neighbor(f_it->second)->info().index) && f_it->first->info().isGhost)
		||  ((f_it->first->info().index >= f_it->first->neighbor(f_it->second)->info().index) && f_it->first->neighbor(f_it->second)->info().isGhost)
		|| f_it->first->info().index == 0 || f_it->first->neighbor(f_it->second)->info().index == 0) continue;
		constrictions.push_back(computeEffectiveRadius(f_it->first, f_it->second));
	}
	return constrictions;
}

template <class Tesselation>
vector<Constriction> FlowBoundingSphere<Tesselation>::getConstrictionsFull()
{
	RTriangulation& Tri = T[currentTes].Triangulation();
	vector<Constriction> constrictions;
	for (FiniteFacetsIterator f_it=Tri.finite_facets_begin(); f_it != Tri.finite_facets_end();f_it++){
		//in the periodic case, we skip facets with lowest id out of the base period
 		 if ( ((f_it->first->info().index <= f_it->first->neighbor(f_it->second)->info().index) && f_it->first->info().isGhost)
		||  ((f_it->first->info().index >= f_it->first->neighbor(f_it->second)->info().index) && f_it->first->neighbor(f_it->second)->info().isGhost)
		|| f_it->first->info().index == 0 || f_it->first->neighbor(f_it->second)->info().index == 0) continue;
		vector<double> rn;
		const CVector& normal = f_it->first->info().facetSurfaces[f_it->second];
		if (!normal[0] && !normal[1] && !normal[2]) continue;
		rn.push_back(computeEffectiveRadius(f_it->first, f_it->second));
		rn.push_back(normal[0]);
		rn.push_back(normal[1]);
		rn.push_back(normal[2]);
		Constriction cons (pair<int,int>(f_it->first->info().index,f_it->first->neighbor(f_it->second)->info().index),rn);
		constrictions.push_back(cons);
	}
	return constrictions;
}

template <class Tesselation> 
double FlowBoundingSphere<Tesselation>::computeEffectiveRadius(CellHandle cell, int j)
{
	RTriangulation& Tri = T[currentTes].Triangulation();
        if (Tri.is_infinite(cell->neighbor(j))) return 0;

	CVector B = cell->vertex(facetVertices[j][1])->point().point()-cell->vertex(facetVertices[j][0])->point().point();
	CVector x = B/sqrt(B.squared_length());
	CVector C = cell->vertex(facetVertices[j][2])->point().point()-cell->vertex(facetVertices[j][0])->point().point();
	CVector z = CGAL::cross_product(x,C);
	CVector y = CGAL::cross_product(x,z);
	y = y/sqrt(y.squared_length());

	double b1[2]; b1[0] = B*x; b1[1] = B*y;
	double c1[2]; c1[0] = C*x; c1[1] = C*y;

	double rA = sqrt(cell->vertex(facetVertices[j][0])->point().weight());
	double rB = sqrt(cell->vertex(facetVertices[j][1])->point().weight());
	double rC = sqrt(cell->vertex(facetVertices[j][2])->point().weight());

	double A = ((pow(rA,2))*(1-c1[0]/b1[0])+((pow(rB,2)*c1[0])/b1[0])-pow(rC,2)+pow(c1[0],2)+pow(c1[1],2)-((pow(b1[0],2)+pow(b1[1],2))*c1[0]/b1[0]))/(2*c1[1]-2*b1[1]*c1[0]/b1[0]);
	double BB = (rA-rC-((rA-rB)*c1[0]/b1[0]))/(c1[1]-b1[1]*c1[0]/b1[0]);
	double CC = (pow(rA,2)-pow(rB,2)+pow(b1[0],2)+pow(b1[1],2))/(2*b1[0]);
	double D = (rA-rB)/b1[0];
	double E = b1[1]/b1[0];
	double F = pow(CC,2)+pow(E,2)*pow(A,2)-2*CC*E*A;

	double c = -F-pow(A,2)+pow(rA,2);
	double b = 2*rA-2*(D-BB*E)*(CC-E*A)-2*A*BB;
	double a = 1-pow((D-BB*E),2)-pow(BB,2);

	if ((pow(b,2)-4*a*c)<0){cout << "NEGATIVE DETERMINANT" << endl; }
	double reff = (-b+sqrt(pow(b,2)-4*a*c))/(2*a);
	if (reff<0) return 0;//happens very rarely, with bounding spheres most probably
	//if the facet involves one ore more bounding sphere, we return R with a minus sign
	if (cell->vertex(facetVertices[j][2])->info().isFictious || cell->vertex(facetVertices[j][1])->info().isFictious || cell->vertex(facetVertices[j][2])->info().isFictious) return -reff;
	else return reff;
}

template <class Tesselation> 
double FlowBoundingSphere<Tesselation>::computeEquivalentRadius(CellHandle cell, int j)
{
	Real fluidSurf = sqrt(cell->info().facetSurfaces[j].squared_length())*cell->info().facetFluidSurfacesRatio[j];
	return sqrt(fluidSurf/M_PI);
}
template <class Tesselation> 
double FlowBoundingSphere<Tesselation>::computeHydraulicRadius(CellHandle cell, int j)
{
	RTriangulation& Tri = T[currentTes].Triangulation();
        if (Tri.is_infinite(cell->neighbor(j))) return 0;
	double Vpore = this->volumePoreVoronoiFraction(cell, j);
	double Ssolid = this->surfaceSolidPore(cell, j, slipBoundary, /*reuse the same facet data*/ true);

	//handle symmetry (tested ok)
	if (slipBoundary && facetNFictious>0) {
		//! Include a multiplier so that permeability will be K/2 or K/4 in symmetry conditions
		Real mult= facetNFictious==1 ? multSym1 : multSym2;
		return Vpore/Ssolid*mult;}
	return Vpore/Ssolid;
}
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::initializePressure( double pZero )
{
        RTriangulation& Tri = T[currentTes].Triangulation();
        FiniteCellsIterator cellEnd = Tri.finite_cells_end();

        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++){
		cell->info().p() = pZero; cell->info().dv()=0;}

        for (int bound=0; bound<6;bound++) {
                int& id = *boundsIds[bound];
		boundingCells[bound].clear();
		if (id<0) continue;
                Boundary& bi = boundary(id);
                if (!bi.flowCondition) {
                        VectorCell tmpCells;
                        tmpCells.resize(10000);
                        VCellIterator cells_it = tmpCells.begin();
                        VCellIterator cells_end = Tri.incident_cells(T[currentTes].vertexHandles[id],cells_it);
                        for (VCellIterator it = tmpCells.begin(); it != cells_end; it++){
				(*it)->info().p() = bi.value;(*it)->info().Pcondition=true;
				boundingCells[bound].push_back(*it);
			}
                }
        }
        IPCells.clear();
        for (unsigned int n=0; n<imposedP.size();n++) {
		CellHandle cell=Tri.locate(imposedP[n].first);
		//check redundancy
		for (unsigned int kk=0;kk<IPCells.size();kk++){
			if (cell==IPCells[kk]) cerr<<"Two imposed pressures fall in the same cell."<<endl;
			else if  (cell->info().Pcondition) cerr<<"Imposed pressure fall in a boundary condition."<<endl;}
		IPCells.push_back(cell);
		cell->info().p()=imposedP[n].second;
		cell->info().Pcondition=true;}
	pressureChanged=false;

	IFCells.clear();
	for (unsigned int n=0; n<imposedF.size();n++) {
		CellHandle cell=Tri.locate(imposedF[n].first);
		//check redundancy
		for (unsigned int kk=0;kk<IPCells.size();kk++){
			if (cell==IPCells[kk]) cerr<<"Both flux and pressure are imposed in the same cell."<<endl;
			else if  (cell->info().Pcondition) cerr<<"Imposed flux fall in a pressure boundary condition."<<endl;}
		IFCells.push_back(cell);
		cell->info().Pcondition=false;}

}

template <class Tesselation> 
bool FlowBoundingSphere<Tesselation>::reApplyBoundaryConditions()
{
	if (!pressureChanged) return false;
        for (int bound=0; bound<6;bound++) {
                int& id = *boundsIds[bound];
		if (id<0) continue;
                Boundary& bi = boundary(id);
                if (!bi.flowCondition) {
                        for (VCellIterator it = boundingCells[bound].begin(); it != boundingCells[bound].end(); it++){
			(*it)->info().p() = bi.value; (*it)->info().Pcondition=true;
			}
                }
        }
        for (unsigned int n=0; n<imposedP.size();n++) {
		IPCells[n]->info().p()=imposedP[n].second;
		IPCells[n]->info().Pcondition=true;}
	pressureChanged=false;
	return true;
}

template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::gaussSeidel(Real dt)
{
	reApplyBoundaryConditions();
	RTriangulation& Tri = T[currentTes].Triangulation();
	int j = 0;
	double m, n, dp_max, p_max, sum_p, p_moy, dp, sum_dp;
	double compFlowFactor=0;
	vector<Real> previousP;
	previousP.resize(Tri.number_of_finite_cells());
	const int num_threads=1;
	bool compressible= (fluidBulkModulus>0);
#ifdef GS_OPEN_MP
	omp_set_num_threads(num_threads);
#endif

       if(debugOut){ cout << "tolerance = " << tolerance << endl;
        cout << "relax = " << relax << endl;}
			vector<Real> t_sum_p, t_dp_max, t_sum_dp, t_p_max;
			t_sum_dp.resize(num_threads);
			t_dp_max.resize(num_threads);
			t_p_max.resize(num_threads);
			t_sum_p.resize(num_threads);
        FiniteCellsIterator cellEnd = Tri.finite_cells_end();
	#ifdef GS_OPEN_MP
		vector<FiniteCellsIterator> cells_its;
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) if ( !cell->info().Pcondition ) cells_its.push_back(cell);
		int numCells=cells_its.size();
		cout<<"cells_its.size() "<<cells_its.size();
	#endif
// 	#pragma omp parallel shared(t_sum_dp, t_dp_max, sum_p, sum_dp,cells_its, j, Tri, relax)
	{
        do {
                int cell2=0; dp_max = 0;p_max = 0;p_moy=0;sum_p=0;sum_dp=0;
		#ifdef GS_OPEN_MP
		cell2=numCells;
		for (int ii=0;ii<num_threads;ii++) t_p_max[ii] =0;
		for (int ii=0;ii<num_threads;ii++) t_dp_max[ii] =0;
		for (int ii=0;ii<num_threads;ii++) t_sum_p[ii]=0;
                for (int ii=0;ii<num_threads;ii++) t_sum_dp[ii]=0;
		int kk=0;
		const int numCells2 = numCells;
		#pragma omp parallel for private(dp, m, n, kk) shared(tolerance, t_sum_dp, t_dp_max, sum_p, sum_dp,cells_its, j, Tri, relax) schedule(dynamic, 1000)
		for (kk=0; kk<numCells2; kk++) {
			const FiniteCellsIterator& cell = cells_its[kk];
			{
		#else
		int bb=-1;
                for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			bb++;
			if ( !cell->info().Pcondition ) {
		                cell2++;
		#endif
				if (compressible && j==0) { previousP[bb]=cell->info().p(); }
				m=0, n=0;
				for (int j2=0; j2<4; j2++) {
					if (!Tri.is_infinite(cell->neighbor(j2))) { 
						/// COMPRESSIBLE: 
						if ( compressible ) {
							compFlowFactor = fluidBulkModulus*dt*cell->info().invVoidVolume();
							m += compFlowFactor*(cell->info().kNorm())[j2] * cell->neighbor(j2)->info().p();
							if (j==0) n +=compFlowFactor*(cell->info().kNorm())[j2];
						} else {							
						/// INCOMPRESSIBLE 
							m += (cell->info().kNorm())[j2] * cell->neighbor(j2)->info().p();
							if ( isinf(m) && j<10 ) cout << "(cell->info().kNorm())[j2] = " << (cell->info().kNorm())[j2] << " cell->neighbor(j2)->info().p() = " << cell->neighbor(j2)->info().p() << endl;
							if (j==0) n += (cell->info().kNorm())[j2];
						}  
					}
				}
				dp = cell->info().p();
				if (n!=0 || j!=0) {
					if (j==0) { if (compressible) cell->info().invSumK=1/(1+n); else cell->info().invSumK=1/n; }
					if ( compressible ) {
					/// COMPRESSIBLE cell->info().p() = ( (previousP - compFlowFactor*cell->info().dv()) + m ) / n ;
						cell->info().p() = ( ((previousP[bb] - ((fluidBulkModulus*dt*cell->info().invVoidVolume())*(cell->info().dv()))) + m) * cell->info().invSumK - cell->info().p()) * relax + cell->info().p();
					} else {
					/// INCOMPRESSIBLE cell->info().p() =   - ( cell->info().dv() - m ) / ( n ) = ( -cell.info().dv() + m ) / n ;
						cell->info().p() = (- (cell->info().dv() - m) * cell->info().invSumK - cell->info().p()) * relax + cell->info().p();
					}
					#ifdef GS_OPEN_MP
					#endif
				}
                                dp -= cell->info().p();
                                #ifdef GS_OPEN_MP
                                const int tn=omp_get_thread_num();
				t_sum_dp[tn] += std::abs(dp);
				t_dp_max[tn]=max(t_dp_max[tn], std::abs(dp));
				t_p_max[tn]= max(t_p_max[tn], std::abs(cell->info().p()));
				t_sum_p[tn]+= std::abs(cell->info().p());
				#else
                                dp_max = max(dp_max, std::abs(dp));
                                p_max = max(p_max, std::abs(cell->info().p()));
                                sum_p += std::abs(cell->info().p());
                                sum_dp += std::abs(dp);
				#endif
                        }
                }
                #ifdef GS_OPEN_MP

                for (int ii=0;ii<num_threads;ii++) p_max =max(p_max, t_p_max[ii]);
		for (int ii=0;ii<num_threads;ii++) dp_max =max(dp_max, t_dp_max[ii]);
		for (int ii=0;ii<num_threads;ii++) sum_p+=t_sum_p[ii];
                for (int ii=0;ii<num_threads;ii++) sum_dp+=t_sum_dp[ii];
                #endif
		p_moy = sum_p/cell2;
 
		#ifdef GS_OPEN_MP
		#pragma omp master
		#endif
		j++;
	#ifdef GS_OPEN_MP
	} while (j<1500);
	#else
	} while ((dp_max/p_max) > tolerance /*&& j<4000*/ /*&& ( dp_max > tolerance )*//* &&*/ /*( j<50 )*/);
	#endif
	}

        if (debugOut) {cout << "pmax " << p_max << "; pmoy : " << p_moy << endl;
        cout << "iteration " << j <<"; erreur : " << dp_max/p_max << endl;}
	computedOnce=true;
}

template <class Tesselation>
double FlowBoundingSphere<Tesselation>::boundaryFlux(unsigned int boundaryId)
{
	RTriangulation& Tri = T[currentTes].Triangulation();
	double Q1=0;

	VectorCell tmpCells;
	tmpCells.resize(10000);
	VCellIterator cells_it = tmpCells.begin();

	VCellIterator cell_up_end = Tri.incident_cells(T[currentTes].vertexHandles[boundaryId],cells_it);
	for (VCellIterator it = tmpCells.begin(); it != cell_up_end; it++)
	{
		const CellHandle& cell = *it;
		if (cell->info().isGhost) continue;
		Q1 -= cell->info().dv();
		for (int j2=0; j2<4; j2++)
			Q1 += (cell->info().kNorm())[j2]* (cell->neighbor(j2)->info().shiftedP()-cell->info().shiftedP());
	}
	return Q1;
}

template <class Tesselation> 
double FlowBoundingSphere<Tesselation>::permeameter(double PInf, double PSup, double Section, double DeltaY, const char *file)
{
  RTriangulation& Tri = T[currentTes].Triangulation();
  std::ofstream kFile(file, std::ios::out);
  double Q2=0, Q1=0;
  int cellQ1=0, cellQ2=0;
  double p_out_max=-10000000, p_out_min=10000000, p_in_max=-100000000, p_in_min=10000000,p_out_moy=0, p_in_moy=0;

  VectorCell tmpCells;
  tmpCells.resize(10000);
  VCellIterator cells_it = tmpCells.begin();

  VCellIterator cell_up_end = Tri.incident_cells(T[currentTes].vertexHandles[yMaxId],cells_it);
  for (VCellIterator it = tmpCells.begin(); it != cell_up_end; it++)
  {
    CellHandle& cell = *it;
    for (int j2=0; j2<4; j2++) {
      if (!cell->neighbor(j2)->info().Pcondition){
	Q1 = Q1 + (cell->neighbor(j2)->info().kNorm())[Tri.mirror_index(cell, j2)]* (cell->neighbor(j2)->info().p()-cell->info().p());
	cellQ1+=1;
	p_out_max = max(cell->neighbor(j2)->info().p(), p_out_max);
	p_out_min = min(cell->neighbor(j2)->info().p(), p_out_min);
	p_out_moy += cell->neighbor(j2)->info().p();}
  }}

  VectorCell tmpCells2;
  tmpCells2.resize(10000);
  VCellIterator cells_it2 = tmpCells2.begin();

  VCellIterator cell_down_end = Tri.incident_cells(T[currentTes].vertexHandles[yMinId],cells_it2);
  for (VCellIterator it = tmpCells2.begin(); it != cell_down_end; it++)
  {
    CellHandle& cell = *it;
    for (int j2=0; j2<4; j2++){
      if (!cell->neighbor(j2)->info().Pcondition){
	Q2 = Q2 + (cell->neighbor(j2)->info().kNorm())[Tri.mirror_index(cell, j2)]* (cell->info().p()-cell->neighbor(j2)->info().p());
	cellQ2+=1;
	p_in_max = max(cell->neighbor(j2)->info().p(), p_in_max);
	p_in_min = min(cell->neighbor(j2)->info().p(), p_in_min);
	p_in_moy += cell->neighbor(j2)->info().p();}
  }}

	double density = 1;
        double viscosity = viscosity;
        double gravity = 1;
        double Vdarcy = Q1/Section;
	double DeltaP = abs(PInf-PSup);
	double DeltaH = DeltaP/ (density*gravity);
	double k = viscosity*Vdarcy*DeltaY / DeltaP; /**m²**/
	double Ks = k*(density*gravity)/viscosity; /**m/s**/
	
	if (debugOut){
	cout << "the maximum superior pressure is = " << p_out_max << " the min is = " << p_out_min << endl;
	cout << "the maximum inferior pressure is = " << p_in_max << " the min is = " << p_in_min << endl;
	cout << "superior average pressure is " << p_out_moy/cellQ1 << endl;
        cout << "inferior average pressure is " << p_in_moy/cellQ2 << endl;
        cout << "celle comunicanti in basso = " << cellQ2 << endl;
        cout << "celle comunicanti in alto = " << cellQ1 << endl;
        cout << "The incoming average flow rate is = " << Q2 << " m^3/s " << endl;
        cout << "The outgoing average flow rate is = " << Q1 << " m^3/s " << endl;
        cout << "The gradient of charge is = " << DeltaH/DeltaY << " [-]" << endl;
        cout << "Darcy's velocity is = " << Vdarcy << " m/s" <<endl;
        cout << "The permeability of the sample is = " << k << " m^2" <<endl;
	cout << endl << "The hydraulic conductivity of the sample is = " << Ks << " m/s" << endl << endl;
	}
	kFile << "yMax id = "<<yMaxId<< "yMin id = "<<yMinId<<endl;
	kFile << "the maximum superior pressure is = " << p_out_max << " the min is = " << p_out_min << endl;
	kFile << "the maximum inferior pressure is = " << p_in_max << " the min is = " << p_in_min << endl;
        kFile << "superior average pressure is " << p_out_moy/cellQ2 << endl;
        kFile << "inferior average pressure is " << p_in_moy/cellQ1 << endl;
        kFile << "celle comunicanti in basso = " << cellQ2 << endl;
        kFile << "celle comunicanti in alto = " << cellQ1 << endl;
	kFile << "The incoming average flow rate is = " << Q2 << " m^3/s " << endl;
        kFile << "The outgoing average flow rate is = " << Q1 << " m^3/s " << endl;
        kFile << "The gradient of charge is = " << DeltaH/DeltaY << " [-]" << endl;
        kFile << "Darcy's velocity is = " << Vdarcy << " m/s" <<endl;
        kFile << "The hydraulic conductivity of the sample is = " << Ks << " m/s" <<endl;
        kFile << "The permeability of the sample is = " << k << " m^2" <<endl;
	return k;
}
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::displayStatistics()
{
	RTriangulation& Tri = T[currentTes].Triangulation();
        int Zero =0, Inside=0, Fictious=0;
        FiniteCellsIterator cellEnd = Tri.finite_cells_end();
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
                int zeros =0;
                for (int j=0; j!=4; j++) {
                        if ((cell->info().kNorm())[j]==0) {
                                zeros+=1;
                        }
                }
                if (zeros==4) {
                        Zero+=1;
                }
                if (!cell->info().fictious()) {
                        Inside+=1;
                } else {
                        Fictious+=1;
                }
        }
        int fict=0, real=0;
        for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
                if (v->info().isFictious) fict+=1;
                else real+=1;
        }
        long Vertices = Tri.number_of_vertices();
	long Cells = Tri.number_of_finite_cells();
	long Facets = Tri.number_of_finite_facets();
        if(debugOut) {cout << "zeros = " << Zero << endl;
	cout << "There are " << Vertices << " vertices, dont " << fict << " fictious et " << real << " reeeeeel" << std::endl;
        cout << "There are " << Cells << " cells " << std::endl;
        cout << "There are " << Facets << " facets " << std::endl;
        cout << "There are " << Inside << " cells INSIDE." << endl;
        cout << "There are " << Fictious << " cells FICTIOUS." << endl;}

	num_particles = real;
}
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::saveVtk(const char* folder)
{
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
        static unsigned int number = 0;
        char filename[80];
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(filename,"%s/out_%d.vtk",folder,number++);
	int firstReal=-1;

	//count fictious vertices and cells
	vtkInfiniteVertices=vtkInfiniteCells=0;
 	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (!isDrawable) vtkInfiniteCells+=1;
	}
	for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
                if (!v->info().isReal()) vtkInfiniteVertices+=1;
                else if (firstReal==-1) firstReal=vtkInfiniteVertices;}

        basicVTKwritter vtkfile((unsigned int) Tri.number_of_vertices()-vtkInfiniteVertices, (unsigned int) Tri.number_of_finite_cells()-vtkInfiniteCells);

        vtkfile.open(filename,"test");

        vtkfile.begin_vertices();
        double x,y,z;
        for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
		if (v->info().isReal()){
		x = (double)(v->point().point()[0]);
                y = (double)(v->point().point()[1]);
                z = (double)(v->point().point()[2]);
                vtkfile.write_point(x,y,z);}
        }
        vtkfile.end_vertices();

        vtkfile.begin_cells();
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
        	if (isDrawable){vtkfile.write_cell(cell->vertex(0)->info().id()-firstReal, cell->vertex(1)->info().id()-firstReal, cell->vertex(2)->info().id()-firstReal, cell->vertex(3)->info().id()-firstReal);}
        }
        vtkfile.end_cells();

	if (permeabilityMap){
	vtkfile.begin_data("Permeability",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().s);}
	}
	vtkfile.end_data();}
	else{
	vtkfile.begin_data("Pressure",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().p());}
	}
	vtkfile.end_data();}

	if (1){
	averageRelativeCellVelocity();
	vtkfile.begin_data("Velocity",CELL_DATA,VECTORS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().averageVelocity()[0],cell->info().averageVelocity()[1],cell->info().averageVelocity()[2]);}
	}
	vtkfile.end_data();}
}

#ifdef XVIEW
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::dessineTriangulation(Vue3D &Vue, RTriangulation &T)
{
        double* Segments = NULL;
        long N_seg = newListeEdges(T, &Segments);
        Vue.Dessine_Segment(Segments, N_seg);
        deleteListeEdges(&Segments, N_seg);
}
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::dessineShortTesselation(Vue3D &Vue, Tesselation &Tes)
{
        if (!Tes.computed()) Tes.compute();
        double* Segments = NULL;
        long N_seg = Tes.newListeShortEdges(&Segments);
        Vue.Dessine_Segment(Segments, N_seg);
        deleteListeEdges(&Segments, N_seg);
}
#endif
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::generateVoxelFile( )
{
	RTriangulation& Tri = T[currentTes].Triangulation();
        double l = 1;
        int dx = 200;
        double eps = l/dx;

        std::ofstream voxelfile("MATRIX",std::ios::out);
        bool solid=false;

        for (double y=0; y<=l; y+=eps) {
                for (double z=0; z<=l; z+=eps) {
                        for (double x=0; x<=l; x+=eps) {
                                solid=false;

                                for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != Tri.finite_vertices_end(); vIt++) {
                                        double radius = sqrt(vIt->point().weight());
                                        if ((sqrt(pow((x- (vIt->point()[0])),2) +pow((y- (vIt->point()[1])),2) +pow((z- (vIt->point()[2])),2))) <= radius) solid=true;
                                }
                                if (solid) voxelfile << 1;
                                else voxelfile << 0;
                        }
                        voxelfile << endl;
                }
        }
}

template <class Tesselation>
double FlowBoundingSphere<Tesselation>::samplePermeability(double& xMin,double& xMax ,double& yMin,double& yMax,double& zMin,double& zMax/*, string key*/)
{
        double Section = (xMax-xMin) * (zMax-zMin);
        double DeltaY = yMax-yMin;
        boundary(yMinId).flowCondition=0;
        boundary(yMaxId).flowCondition=0;
        boundary(yMinId).value=0;
        boundary(yMaxId).value=1;
	double pZero = abs((boundary(yMinId).value-boundary(yMaxId).value)/2);
	initializePressure( pZero );
	gaussSeidel();
	const char *kk = "Permeability";
        return permeameter(boundary(yMinId).value, boundary(yMaxId).value, Section, DeltaY, kk);
}
template <class Tesselation> 
bool FlowBoundingSphere<Tesselation>::isInsideSphere(double& x, double& y, double& z)
{
	RTriangulation& Tri = T[currentTes].Triangulation();
        for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != Tri.finite_vertices_end(); vIt++) {
                double radius = vIt->point().weight();
                if (pow((x- (vIt->point()[0])),2) +pow((y- (vIt->point()[1])),2) +pow((z- (vIt->point()[2])),2)   <= radius) return true;
        }
        return false;
}
template <class Tesselation> 
void FlowBoundingSphere<Tesselation>::sliceField(const char *filename)
{
        /** Pressure field along one cutting plane **/
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
        CellHandle permeameter;

        std::ofstream consFile(filename,std::ios::out);

        int intervals = 400;
        double Ry = (yMax-yMin) /intervals;
        double Rz = (zMax-zMin) /intervals;
	double X=0.5;
                for (double Y=min(yMax,yMin); Y<=max(yMax,yMin); Y=Y+abs(Ry)) {
                        for (double Z=min(zMin,zMax); Z<=max(zMin,zMax); Z=Z+abs(Rz)) {
			  permeameter = Tri.locate(Point(X, Y, Z));
			  consFile << permeameter->info().p() <<" ";
                        }
                        consFile << endl;}
                consFile << endl;
        consFile.close();
}

template <class Tesselation>
void  FlowBoundingSphere<Tesselation>::computeEdgesSurfaces()
{
	RTriangulation& Tri = T[currentTes].Triangulation();
	//first, copy interacting pairs and normal lub forces form prev. triangulation in a sorted structure for initializing the new lub. Forces
	vector<vector<pair<unsigned int,Real> > > lubPairs;
	lubPairs.resize(Tri.number_of_vertices()+1);
	for (unsigned int k=0; k<edgeNormalLubF.size(); k++)
		lubPairs[min(edgeIds[k].first->id(),edgeIds[k].second->id())].push_back(
			pair<int,Real> (max(edgeIds[k].first->id(),edgeIds[k].second->id()),edgeNormalLubF[k]));

	//Now we reset the containers and initialize them
	edgeSurfaces.clear(); edgeIds.clear(); edgeNormalLubF.clear();
	FiniteEdgesIterator ed_it;
	for ( FiniteEdgesIterator ed_it = Tri.finite_edges_begin(); ed_it!=Tri.finite_edges_end();ed_it++ )
	{
		const VertexInfo& vi1=(ed_it->first)->vertex(ed_it->second)->info();
		const VertexInfo& vi2=(ed_it->first)->vertex(ed_it->third)->info();
	  
		//We eliminate edges that would be periodic replications or involving two bounding objects, i.e. the min id must be non-ghost and non-fictious
		if (vi1.id()<vi2.id()) {if (vi1.isFictious || vi2.isGhost) continue;}
		else if (vi2.isFictious || vi2.isGhost) continue;
		double area = T[currentTes].computeVFacetArea(ed_it);
		edgeSurfaces.push_back(area);
		unsigned int id1 = vi1.id();
		unsigned int id2 = vi2.id();
		edgeIds.push_back(pair<const VertexInfo*,const VertexInfo*>(&vi1,&vi2));
		//For persistant edges, we must transfer the lub. force value from the older triangulation structure
		if (id1>id2) swap(id1,id2);
		unsigned int i=0;
		//Look for the pair (id1,id2) in lubPairs
		while (i<lubPairs[id1].size()) {
			if (lubPairs[id1][i].first == id2) {
				//it's found, we copy the lub force
				edgeNormalLubF.push_back(lubPairs[id1][i].second);
				break;}
			++i;
		}
		// not found, we initialize with zero lub force
    		if (i==lubPairs[id1].size()) edgeNormalLubF.push_back(0);
	}
}

template <class Tesselation> 
Vector3r FlowBoundingSphere<Tesselation>::computeViscousShearForce(const Vector3r& deltaV, const int& edge_id, const Real& Rh)
{
    Vector3r tau = deltaV*viscosity/Rh;
    return tau * edgeSurfaces[edge_id];
}

template <class Tesselation> 
Vector3r FlowBoundingSphere<Tesselation>::computeShearLubricationForce(const Vector3r& deltaShearV, const Real& dist, const int& edge_id, const Real& eps, const Real& centerDist, const Real& meanRad )
{
    Real d = max(dist,0.) + 2.*eps*meanRad;
    Vector3r viscLubF = 0.5 * Mathr::PI * viscosity * (-2*meanRad + centerDist*log(centerDist/d)) * deltaShearV;
    return viscLubF;
}

template <class Tesselation> 
Vector3r FlowBoundingSphere<Tesselation>::computePumpTorque(const Vector3r& deltaShearAngV, const Real& dist, const int& edge_id, const Real& eps, const Real& meanRad )
{
    Real d = max(dist,0.) + 2.*eps*meanRad;
    Vector3r viscPumpC = Mathr::PI * viscosity * pow(meanRad,3) *(3./20. * log(meanRad/d) + 63./500. * (d/meanRad) * log(meanRad/d)) * deltaShearAngV;
    return viscPumpC;
}

template <class Tesselation> 
Vector3r FlowBoundingSphere<Tesselation>::computeTwistTorque(const Vector3r& deltaNormAngV, const Real& dist, const int& edge_id, const Real& eps, const Real& meanRad )
{
    Real d = max(dist,0.) + 2.*eps*meanRad;
    Vector3r twistC = Mathr::PI * viscosity * pow(meanRad,2) * d * log(meanRad/d) * deltaNormAngV;
    return twistC;
}


template <class Tesselation> 
Real FlowBoundingSphere<Tesselation>::computeNormalLubricationForce(const Real& deltaNormV, const Real& dist, const int& edge_id, const Real& eps, const Real& stiffness, const Real& dt, const Real& meanRad)
{
	//FIXME: here introduce elasticity
	Real d = max(dist,0.) + 2.*eps*meanRad;//account for grains roughness
	if (stiffness>0) {
		const Real k = stiffness*meanRad;
		Real prevForce = edgeNormalLubF[edge_id];
		Real instantVisc = 1.5*Mathr::PI*pow(meanRad,2)*viscosity/(d-prevForce/k);
		Real normLubF = instantVisc*(deltaNormV + prevForce/(k*dt))/(1+instantVisc/(k*dt));
		edgeNormalLubF[edge_id]=normLubF;
		return normLubF;
	} else {
		Real normLubF = (1.5*Mathr::PI*pow(meanRad,2)* viscosity* deltaNormV)/d;
		return normLubF;
	}
}

} //namespace CGT

#endif //FLOW_ENGINE


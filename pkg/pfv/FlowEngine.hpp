/*************************************************************************
*  Copyright (C) 2009 by Emanuele Catalano <catalano@grenoble-inp.fr>    *
*  Copyright (C) 2009 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/*!

FlowEngine is an interface between Yade and the fluid solvers defined in lib/triangulation and using the PFV scheme.
There are also a few non-trivial functions defined here, such has building triangulation and computating elements volumes.
The code strongly relies on CGAL library for building triangulations on the top of sphere packings.
CGAL's trangulations introduce the data associated to each element of the triangulation through template parameters (Cell_info and Vertex_info).  
FlowEngine is thus defined via the template class TemplateFlowEngine(Cellinfo,VertexInfo), for easier addition of variables in various couplings.
PeriodicFlowEngine is a variant for periodic boundary conditions.

Which solver will be actually used internally to obtain pore pressure will depend partly on compile time flags (libcholmod available and #define LINSOLV),
and on runtime settings (useSolver=0: Gauss-Seidel (iterative), useSolver=1: CHOLMOD if available (direct sparse cholesky)).

The files defining lower level classes are in yade/lib. The code uses CGAL::Triangulation3 for managing the mesh and storing data. Eigen3::Sparse, suitesparse::cholmod, and metis are used for solving the linear systems with a direct method (Cholesky). An iterative method (Gauss-Seidel) is implemented directly in Yade and can be used as a fallback (see FlowEngine::useSolver).

Most classes in lib/triangulation are templates, and are therefore completely defined in header files.
A pseudo hpp/cpp split is reproduced for clarity, with hpp/ipp extensions (but again, in the end they are all include files).
The files are
- RegularTriangulation.h : declaration of info types embeded in the mesh (template parameters of CGAL::Triangulation3), and instanciation of RTriangulation classes
- Tesselation.h/ipp : encapsulate RegularTriangulation and adds functions to manipulate the dual (Voronoi) graph of the triangulation
- Network.hpp/ipp : specialized for PFV model (the two former are used independently by TesselationWrapper), a set of functions to determine volumes and surfaces of intersections between spheres and various subdomains. Contains two triangulations for smooth transitions while remeshing - e.g. interpolating values in the new mesh using the previous one.
- FlowBoundingSphere.hpp/ipp and PeriodicFlow.hpp/ipp + LinSolv variants: implement the solver in itself (mesh, boundary conditions, solving, defining fluid-particles interactions)
- FlowEngine.hpp/ipp/cpp (this file)

Variants for periodic boundary conditions are also present.

*/

#pragma once
#include<yade/core/PartialEngine.hpp>
#include<yade/pkg/dem/TriaxialCompressionEngine.hpp>
#include<yade/pkg/dem/TesselationWrapper.hpp>
#include<yade/lib/triangulation/FlowBoundingSphere.hpp>
#include<yade/lib/triangulation/PeriodicFlow.hpp>

/// Frequently used:
typedef CGT::CVector CVector;
typedef CGT::Point Point;

/// Converters for Eigen and CGAL vectors
inline CVector makeCgVect ( const Vector3r& yv ) {return CVector ( yv[0],yv[1],yv[2] );}
inline Point makeCgPoint ( const Vector3r& yv ) {return Point ( yv[0],yv[1],yv[2] );}
inline Vector3r makeVector3r ( const Point& yv ) {return Vector3r ( yv[0],yv[1],yv[2] );}
inline Vector3r makeVector3r ( const CVector& yv ) {return Vector3r ( yv[0],yv[1],yv[2] );}

/// C++ templates and YADE_CLASS_... macro family are not very compatible, this #define is a hack to make it work
/// TEMPLATE_FLOW_NAME will be the name of a serializable TemplateFlowEngine<...> instance, which can in turn be
/// inherited from. The instance itself will be useless for actual simulations.
#ifndef TEMPLATE_FLOW_NAME
#error You must define TEMPLATE_FLOW_NAME in your source file before including FlowEngine.hpp
#endif

#ifdef LINSOLV
#define DEFAULTSOLVER CGT::FlowBoundingSphereLinSolv<_Tesselation>
#else
#define DEFAULTSOLVER CGT::FlowBoundingSphere<_Tesselation>
#endif
	
template<class _CellInfo, class _VertexInfo, class _Tesselation=CGT::_Tesselation<CGT::TriangulationTypes<_VertexInfo,_CellInfo> >, class solverT=DEFAULTSOLVER >
class TemplateFlowEngine : public PartialEngine
{	
	public :
		typedef solverT									FlowSolver;
		typedef FlowSolver								Solver;//FIXME: useless alias, search/replace then remove
		typedef typename FlowSolver::Tesselation					Tesselation;
		typedef typename FlowSolver::RTriangulation					RTriangulation;
		typedef typename FlowSolver::FiniteVerticesIterator                  	  	FiniteVerticesIterator;
		typedef typename FlowSolver::FiniteCellsIterator				FiniteCellsIterator;
		typedef typename FlowSolver::CellHandle						CellHandle;
		typedef typename FlowSolver::FiniteEdgesIterator				FiniteEdgesIterator;
		typedef typename FlowSolver::VertexHandle                    			VertexHandle;
		typedef typename RTriangulation::Triangulation_data_structure::Cell::Info       CellInfo;
		typedef typename RTriangulation::Triangulation_data_structure::Vertex::Info     VertexInfo;
		
	protected:
		shared_ptr<FlowSolver> solver;
		shared_ptr<FlowSolver> backgroundSolver;
		volatile bool backgroundCompleted;
		Cell cachedCell;
		struct posData {Body::id_t id; Vector3r pos; Real radius; bool isSphere; bool exists; posData(){exists=0;}};
		vector<posData> positionBufferCurrent;//reflect last known positions before we start computations
		vector<posData> positionBufferParallel;//keep the positions from a given step for multithread factorization
		//copy positions in a buffer for faster and/or parallel access
		void setPositionsBuffer(bool current);
		virtual void trickPermeability() {};

	public :
		int retriangulationLastIter;
		enum {wall_xmin, wall_xmax, wall_ymin, wall_ymax, wall_zmin, wall_zmax};
		Vector3r normal [6];
		bool currentTes;
		bool metisForced;
		int idOffset;
		double epsVolCumulative;
		int ReTrg;
		int ellapsedIter;
		void initSolver (FlowSolver& flow);
		#ifdef LINSOLV
		void setForceMetis (bool force);
		bool getForceMetis ();
		#endif
		void triangulate (Solver& flow);
		void addBoundary (Solver& flow);
		void buildTriangulation (double pZero, Solver& flow);
		void buildTriangulation (Solver& flow);
		void updateVolumes (Solver& flow);
		void initializeVolumes (Solver& flow);
		void boundaryConditions(Solver& flow);
		void updateBCs ( Solver& flow ) {
			if (flow.T[flow.currentTes].maxId>0) boundaryConditions(flow);//avoids crash at iteration 0, when the packing is not bounded yet
			else LOG_ERROR("updateBCs not applied");
			flow.pressureChanged=true;}

		void imposeFlux(Vector3r pos, Real flux);
		unsigned int imposePressure(Vector3r pos, Real p);
		void setImposedPressure(unsigned int cond, Real p);
		void clearImposedPressure();
		void clearImposedFlux();
		void computeViscousForces(Solver& flow);
		Real getCellFlux(unsigned int cond);
		Real getBoundaryFlux(unsigned int boundary) {return solver->boundaryFlux(boundary);}
		Vector3r fluidForce(unsigned int idSph) {
			const CGT::CVector& f=solver->T[solver->currentTes].vertex(idSph)->info().forces; return Vector3r(f[0],f[1],f[2]);}
			
		Vector3r shearLubForce(unsigned int id_sph) {
			return (solver->shearLubricationForces.size()>id_sph)?solver->shearLubricationForces[id_sph]:Vector3r::Zero();}
		Vector3r shearLubTorque(unsigned int id_sph) {
			return (solver->shearLubricationTorques.size()>id_sph)?solver->shearLubricationTorques[id_sph]:Vector3r::Zero();}
		Vector3r pumpLubTorque(unsigned int id_sph) {
			return (solver->pumpLubricationTorques.size()>id_sph)?solver->pumpLubricationTorques[id_sph]:Vector3r::Zero();}
		Vector3r twistLubTorque(unsigned int id_sph) {
			return (solver->twistLubricationTorques.size()>id_sph)?solver->twistLubricationTorques[id_sph]:Vector3r::Zero();}
		Vector3r normalLubForce(unsigned int id_sph) {
			return (solver->normalLubricationForce.size()>id_sph)?solver->normalLubricationForce[id_sph]:Vector3r::Zero();}
		Matrix3r bodyShearLubStress(unsigned int id_sph) {
			return (solver->shearLubricationBodyStress.size()>id_sph)?solver->shearLubricationBodyStress[id_sph]:Matrix3r::Zero();}
		Matrix3r bodyNormalLubStress(unsigned int id_sph) {
			return (solver->normalLubricationBodyStress.size()>id_sph)?solver->normalLubricationBodyStress[id_sph]:Matrix3r::Zero();}	
		Vector3r shearVelocity(unsigned int interaction) {
			return (solver->deltaShearVel[interaction]);}
		Vector3r normalVelocity(unsigned int interaction) {
			return (solver->deltaNormVel[interaction]);}
		Matrix3r normalStressInteraction(unsigned int interaction) {
			return (solver->normalStressInteraction[interaction]);}
		Matrix3r shearStressInteraction(unsigned int interaction) {
			return (solver->shearStressInteraction[interaction]);}
		Vector3r normalVect(unsigned int interaction) {
			return (solver->normalV[interaction]);}
		Real surfaceDistanceParticle(unsigned int interaction) {
			return (solver->surfaceDistance[interaction]);}
		Real edgeSize() {
			return (solver->edgeIds.size());}
		Real OSI() {
			return (solver->onlySpheresInteractions.size());}
		int onlySpheresInteractions(unsigned int interaction) {
			return (solver->onlySpheresInteractions[interaction]);}
		python::list getConstrictions(bool all) {
			vector<Real> csd=solver->getConstrictions(); python::list pycsd;
			for (unsigned int k=0;k<csd.size();k++) if ((all && csd[k]!=0) || csd[k]>0) pycsd.append(csd[k]); return pycsd;}
		python::list getConstrictionsFull(bool all) {
			vector<Constriction> csd=solver->getConstrictionsFull(); python::list pycsd;
			for (unsigned int k=0;k<csd.size();k++) if ((all && csd[k].second[0]!=0) || csd[k].second[0]>0) {
				python::list cons;
				cons.append(csd[k].first.first);
				cons.append(csd[k].first.second);
				cons.append(csd[k].second[0]);
				cons.append(csd[k].second[1]);
				cons.append(csd[k].second[2]);
				cons.append(csd[k].second[3]);
				pycsd.append(cons);}
			return pycsd;}
		
		template<class Cellhandle>
		Real volumeCellSingleFictious (Cellhandle cell);
		template<class Cellhandle>
		Real volumeCellDoubleFictious (Cellhandle cell);
		template<class Cellhandle>
		Real volumeCellTripleFictious (Cellhandle cell);
		template<class Cellhandle>
		Real volumeCell (Cellhandle cell);
		void Oedometer_Boundary_Conditions();
		void averageRealCellVelocity();
		void saveVtk(const char* folder) {solver->saveVtk(folder);}
		vector<Real> avFlVelOnSph(unsigned int idSph) {return solver->averageFluidVelocityOnSphere(idSph);}

// 		void setBoundaryVel(Vector3r vel) {topBoundaryVelocity=vel; updateTriangulation=true;}
		void pressureProfile(double wallUpY, double wallDownY) {return solver->measurePressureProfile(wallUpY,wallDownY);}
		double getPorePressure(Vector3r pos){return solver->getPorePressure(pos[0], pos[1], pos[2]);}
		int getCell(double posX, double posY, double posZ){return solver->getCell(posX, posY, posZ);}
		unsigned int nCells(){return solver->T[solver->currentTes].cellHandles.size();}
		python::list getVertices(unsigned int id){
			python::list ids;
			if (id>=solver->T[solver->currentTes].cellHandles.size()) {LOG_ERROR("id out of range, max value is "<<solver->T[solver->currentTes].cellHandles.size()); return ids;}			
			for (unsigned int i=0;i<4;i++) ids.append(solver->T[solver->currentTes].cellHandles[id]->vertex(i)->info().id());
			return ids;
		}
		double averageSlicePressure(double posY){return solver->averageSlicePressure(posY);}
		double averagePressure(){return solver->averagePressure();}

		void emulateAction(){
			scene = Omega::instance().getScene().get();
			action();}
		#ifdef LINSOLV
		void	exportMatrix(string filename) {if (useSolver==3) solver->exportMatrix(filename.c_str());
				else cerr<<"available for Cholmod solver (useSolver==3)"<<endl;}
		void	exportTriplets(string filename) {if (useSolver==3) solver->exportTriplets(filename.c_str());
				else cerr<<"available for Cholmod solver (useSolver==3)"<<endl;}
		void	cholmodStats() {
					cerr << cholmod_print_common((char*)string("PFV Cholmod factorization").c_str(),&(solver->eSolver.cholmod()))<<endl;
					cerr << "cholmod method:" << solver->eSolver.cholmod().selected<<endl;
					cerr << "METIS called:"<<solver->eSolver.cholmod().called_nd<<endl;}
		bool	metisUsed() {return bool(solver->eSolver.cholmod().called_nd);}
		#endif

		virtual ~TemplateFlowEngine();
		virtual void action();
		virtual void backgroundAction();
		
		//commodities
		void compTessVolumes() {
			solver->T[solver->currentTes].compute();
			solver->T[solver->currentTes].computeVolumes();
		}
		Real getVolume (Body::id_t id) {
			if (solver->T[solver->currentTes].Max_id() <= 0) {emulateAction(); /*LOG_WARN("Not triangulated yet, emulating action");*/}
			if (solver->T[solver->currentTes].Volume(id) == -1) {compTessVolumes();/* LOG_WARN("Computing all volumes now, as you did not request it explicitely.");*/}
			return (solver->T[solver->currentTes].Max_id() >= id) ? solver->T[solver->currentTes].Volume(id) : -1;}

		YADE_CLASS_PYCLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(TemplateFlowEngine,TEMPLATE_FLOW_NAME,PartialEngine,"An engine to solve flow problem in saturated granular media. Model description can be found in [Chareyre2012a]_ and [Catalano2014a]_. See the example script FluidCouplingPFV/oedometer.py. More documentation to come.\n\n.. note::Multi-threading seems to work fine for Cholesky decomposition, but it fails for the solve phase in which -j1 is the fastest, here we specify thread numbers independently using :yref:`FlowEngine::numFactorizeThreads` and :yref:`FlowEngine::numSolveThreads`. These multhreading settings are only impacting the behaviour of openblas library and are relatively independant of :yref:`FlowEngine::multithread`. However, the settings have to be globally consistent. For instance, :yref:`multithread<FlowEngine::multithread>` =True with  yref:`numFactorizeThreads<FlowEngine::numFactorizeThreads>` = yref:`numSolveThreads<FlowEngine::numSolveThreads>` = 4 implies that openblas will mobilize 8 processors at some point. If the system does not have so many procs. it will hurt performance.",
		((bool,isActivated,true,,"Activates Flow Engine"))
		((bool,first,true,,"Controls the initialization/update phases"))
		((bool,doInterpolate,false,,"Force the interpolation of cell's info while remeshing. By default, interpolation would be done only for compressible fluids. It can be forced with this flag."))
		((double, fluidBulkModulus, 0.,,"Bulk modulus of fluid (inverse of compressibility) K=-dP*V/dV [Pa]. Flow is compressible if fluidBulkModulus > 0, else incompressible."))
		((Real, dt, 0,,"timestep [s]"))
		((bool,permeabilityMap,false,,"Enable/disable stocking of average permeability scalar in cell infos."))
		((bool, slipBoundary, true,, "Controls friction condition on lateral walls"))
		((bool,waveAction, false,, "Allow sinusoidal pressure condition to simulate ocean waves"))
		((double, sineMagnitude, 0,, "Pressure value (amplitude) when sinusoidal pressure is applied (p )"))
		((double, sineAverage, 0,,"Pressure value (average) when sinusoidal pressure is applied"))
		((bool, debug, false,,"Activate debug messages"))
		((double, wallThickness,0.001,,"Walls thickness"))
		((double,pZero,0,,"The value used for initializing pore pressure. It is useless for incompressible fluid, but important for compressible model."))
		((double,tolerance,1e-06,,"Gauss-Seidel tolerance"))
		((double,relax,1.9,,"Gauss-Seidel relaxation"))
		((bool, updateTriangulation, 0,,"If true the medium is retriangulated. Can be switched on to force retriangulation after some events (else it will be true periodicaly based on :yref:`FlowEngine::defTolerance` and :yref:`FlowEngine::meshUpdateInterval`. Of course, it costs CPU time."))
		((int,meshUpdateInterval,1000,,"Maximum number of timesteps between re-triangulation events. See also :yref:`FlowEngine::defTolerance`."))
		((double, epsVolMax, 0,(Attr::readonly),"Maximal absolute volumetric strain computed at each iteration. |yupdate|"))
		((double, defTolerance,0.05,,"Cumulated deformation threshold for which retriangulation of pore space is performed. If negative, the triangulation update will occure with a fixed frequency on the basis of :yref:`FlowEngine::meshUpdateInterval`"))
		((double, porosity, 0,(Attr::readonly),"Porosity computed at each retriangulation |yupdate|"))
		((bool,meanKStat,false,,"report the local permeabilities' correction"))
		((bool,clampKValues,true,,"If true, clamp local permeabilities in [minKdivKmean,maxKdivKmean]*globalK. This clamping can avoid singular values in the permeability matrix and may reduce numerical errors in the solve phase. It will also hide junk values if they exist, or bias all values in very heterogeneous problems. So, use this with care."))
		((Real,minKdivKmean,0.0001,,"define the min K value (see :yref:`FlowEngine::clampKValues`)"))
		((Real,maxKdivKmean,100,,"define the max K value (see :yref:`FlowEngine::clampKValues`)"))
		((double,permeabilityFactor,1.0,,"permability multiplier"))
		((double,viscosity,1.0,,"viscosity of the fluid"))
		((double,stiffness, 10000,,"equivalent contact stiffness used in the lubrication model"))
		((int, useSolver, 0,, "Solver to use 0=G-Seidel, 1=Taucs, 2-Pardiso, 3-CHOLMOD"))
		((int, xmin,0,(Attr::readonly),"Index of the boundary $x_{min}$. This index is not equal the the id of the corresponding body in general, it may be used to access the corresponding attributes (e.g. flow.bndCondValue[flow.xmin], flow.wallId[flow.xmin],...)."))
		((int, xmax,1,(Attr::readonly),"See :yref:`FlowEngine::xmin`."))
		((int, ymin,2,(Attr::readonly),"See :yref:`FlowEngine::xmin`."))
		((int, ymax,3,(Attr::readonly),"See :yref:`FlowEngine::xmin`."))
		((int, zmin,4,(Attr::readonly),"See :yref:`FlowEngine::xmin`."))
		((int, zmax,5,(Attr::readonly),"See :yref:`FlowEngine::xmin`."))

		((vector<bool>, bndCondIsPressure, vector<bool>(6,false),,"defines the type of boundary condition for each side. True if pressure is imposed, False for no-flux. Indexes can be retrieved with :yref:`FlowEngine::xmin` and friends."))
		((vector<double>, bndCondValue, vector<double>(6,0),,"Imposed value of a boundary condition. Only applies if the boundary condition is imposed pressure, else the imposed flux is always zero presently (may be generalized to non-zero imposed fluxes in the future)."))
		//FIXME: to be implemented:
		((vector<Vector3r>, boundaryVelocity, vector<Vector3r>(6,Vector3r::Zero()),, "velocity on top boundary, only change it using :yref:`FlowEngine::setBoundaryVel`"))
		((int, ignoredBody,-1,,"Id of a sphere to exclude from the triangulation.)"))
		((vector<int>, wallIds,vector<int>(6),,"body ids of the boundaries (default values are ok only if aabbWalls are appended before spheres, i.e. numbered 0,...,5)"))
		((vector<bool>, boundaryUseMaxMin, vector<bool>(6,true),,"If true (default value) bounding sphere is added as function of max/min sphere coord, if false as function of yade wall position"))
// 					((bool, display_force, false,,"display the lubrication force applied on particles"))
// 					((bool, create_file, false,,"create file of velocities"))
		((bool, viscousShear, false,,"compute viscous shear terms as developped by Donia Marzougui (FIXME: ref.)"))
		((bool, shearLubrication, false,,"compute shear lubrication force as developped by Brule (FIXME: ref.) "))
		((bool, pumpTorque, false,,"Compute pump torque applied on particles "))
		((bool, twistTorque, false,,"Compute twist torque applied on particles "))
		((double, eps, 0.00001,,"roughness defined as a fraction of particles size, giving the minimum distance between particles in the lubrication model."))
		((bool, pressureForce, true,,"compute the pressure field and associated fluid forces. WARNING: turning off means fluid flow is not computed at all."))
		((bool, normalLubrication, false,,"compute normal lubrication force as developped by Brule"))
		((bool, viscousNormalBodyStress, false,,"compute normal viscous stress applied on each body"))
		((bool, viscousShearBodyStress, false,,"compute shear viscous stress applied on each body"))
		((bool, multithread, false,,"Build triangulation and factorize in the background (multi-thread mode)"))
		#ifdef EIGENSPARSE_LIB
		((int, numSolveThreads, 1,,"number of openblas threads in the solve phase."))
		((int, numFactorizeThreads, 1,,"number of openblas threads in the factorization phase"))
		#endif
		((vector<Real>, boundaryPressure,vector<Real>(),,"values defining pressure along x-axis for the top surface. See also :yref:`TEMPLATE_FLOW_NAME::boundaryXPos`"))
		((vector<Real>, boundaryXPos,vector<Real>(),,"values of the x-coordinate for which pressure is defined. See also :yref:`TEMPLATE_FLOW_NAME::boundaryPressure`"))
		,
		/*deprec*/
		((meanK_opt,clampKValues,"the name changed"))
		,,
		timingDeltas=shared_ptr<TimingDeltas>(new TimingDeltas);
		for (int i=0; i<6; ++i){normal[i]=Vector3r::Zero(); wallIds[i]=i;}
		normal[wall_ymin].y()=normal[wall_xmin].x()=normal[wall_zmin].z()=1;
		normal[wall_ymax].y()=normal[wall_xmax].x()=normal[wall_zmax].z()=-1;
		solver = shared_ptr<FlowSolver> (new FlowSolver);
		first=true;
		epsVolMax=epsVolCumulative=retriangulationLastIter=0;
		ReTrg=1;
		backgroundCompleted=true;
		ellapsedIter=0;
		metisForced=false;
		,
		.def("imposeFlux",&TemplateFlowEngine::imposeFlux,(python::arg("pos"),python::arg("p")),"Impose incoming flux in boundary cell of location 'pos'.")
		.def("imposePressure",&TemplateFlowEngine::imposePressure,(python::arg("pos"),python::arg("p")),"Impose pressure in cell of location 'pos'. The index of the condition is returned (for multiple imposed pressures at different points).")
		.def("setImposedPressure",&TemplateFlowEngine::setImposedPressure,(python::arg("cond"),python::arg("p")),"Set pressure value at the point indexed 'cond'.")
		.def("clearImposedPressure",&TemplateFlowEngine::clearImposedPressure,"Clear the list of points with pressure imposed.")
		.def("clearImposedFlux",&TemplateFlowEngine::clearImposedFlux,"Clear the list of points with flux imposed.")
		.def("getCellFlux",&TemplateFlowEngine::getCellFlux,(python::arg("cond")),"Get influx in cell associated to an imposed P (indexed using 'cond').")
		.def("getBoundaryFlux",&TemplateFlowEngine::getBoundaryFlux,(python::arg("boundary")),"Get total flux through boundary defined by its body id.\n\n.. note:: The flux may be not zero even for no-flow condition. This artifact comes from cells which are incident to two or more boundaries (along the edges of the sample, typically). Such flux evaluation on impermeable boundary is just irrelevant, it does not imply that the boundary condition is not applied properly.")
		.def("getConstrictions",&TemplateFlowEngine::getConstrictions,(python::arg("all")=true),"Get the list of constriction radii (inscribed circle) for all finite facets (if all==True) or all facets not incident to a virtual bounding sphere (if all==False).  When all facets are returned, negative radii denote facet incident to one or more fictious spheres.")
		.def("getConstrictionsFull",&TemplateFlowEngine::getConstrictionsFull,(python::arg("all")=true),"Get the list of constrictions (inscribed circle) for all finite facets (if all==True), or all facets not incident to a fictious bounding sphere (if all==False). When all facets are returned, negative radii denote facet incident to one or more fictious spheres. The constrictions are returned in the format {{cell1,cell2}{rad,nx,ny,nz}}")
		.def("edgeSize",&TemplateFlowEngine::edgeSize,"Return the number of interactions.")
		.def("OSI",&TemplateFlowEngine::OSI,"Return the number of interactions only between spheres.")
		
		.def("saveVtk",&TemplateFlowEngine::saveVtk,(python::arg("folder")="./VTK"),"Save pressure field in vtk format. Specify a folder name for output.")
		.def("avFlVelOnSph",&TemplateFlowEngine::avFlVelOnSph,(python::arg("idSph")),"compute a sphere-centered average fluid velocity")
		.def("fluidForce",&TemplateFlowEngine::fluidForce,(python::arg("idSph")),"Return the fluid force on sphere idSph.")
		.def("shearLubForce",&TemplateFlowEngine::shearLubForce,(python::arg("idSph")),"Return the shear lubrication force on sphere idSph.")
		.def("shearLubTorque",&TemplateFlowEngine::shearLubTorque,(python::arg("idSph")),"Return the shear lubrication torque on sphere idSph.")
		.def("normalLubForce",&TemplateFlowEngine::normalLubForce,(python::arg("idSph")),"Return the normal lubrication force on sphere idSph.")
		.def("bodyShearLubStress",&TemplateFlowEngine::bodyShearLubStress,(python::arg("idSph")),"Return the shear lubrication stress on sphere idSph.")
		.def("bodyNormalLubStress",&TemplateFlowEngine::bodyNormalLubStress,(python::arg("idSph")),"Return the normal lubrication stress on sphere idSph.")
		.def("shearVelocity",&TemplateFlowEngine::shearVelocity,(python::arg("idSph")),"Return the shear velocity of the interaction.")
		.def("normalVelocity",&TemplateFlowEngine::normalVelocity,(python::arg("idSph")),"Return the normal velocity of the interaction.")
		.def("normalVect",&TemplateFlowEngine::normalVect,(python::arg("idSph")),"Return the normal vector between particles.")
		.def("surfaceDistanceParticle",&TemplateFlowEngine::surfaceDistanceParticle,(python::arg("interaction")),"Return the distance between particles.")
		.def("onlySpheresInteractions",&TemplateFlowEngine::onlySpheresInteractions,(python::arg("interaction")),"Return the id of the interaction only between spheres.")
		.def("pressureProfile",&TemplateFlowEngine::pressureProfile,(python::arg("wallUpY"),python::arg("wallDownY")),"Measure pore pressure in 6 equally-spaced points along the height of the sample")
		.def("getPorePressure",&TemplateFlowEngine::getPorePressure,(python::arg("pos")),"Measure pore pressure in position pos[0],pos[1],pos[2]")
		.def("averageSlicePressure",&TemplateFlowEngine::averageSlicePressure,(python::arg("posY")),"Measure slice-averaged pore pressure at height posY")
		.def("averagePressure",&TemplateFlowEngine::averagePressure,"Measure averaged pore pressure in the entire volume")
		.def("updateBCs",&TemplateFlowEngine::updateBCs,"tells the engine to update it's boundary conditions before running (especially useful when changing boundary pressure - should not be needed for point-wise imposed pressure)")
		.def("emulateAction",&TemplateFlowEngine::emulateAction,"get scene and run action (may be used to manipulate an engine outside the timestepping loop).")
		.def("getCell",&TemplateFlowEngine::getCell,(python::arg("pos")),"get id of the cell containing (X,Y,Z).")
		.def("nCells",&TemplateFlowEngine::nCells,"get the total number of finite cells in the triangulation.")
		.def("getVertices",&TemplateFlowEngine::getVertices,(python::arg("id")),"get the vertices of a cell")
		#ifdef LINSOLV
		.def("exportMatrix",&TemplateFlowEngine::exportMatrix,(python::arg("filename")="matrix"),"Export system matrix to a file with all entries (even zeros will displayed).")
		.def("exportTriplets",&TemplateFlowEngine::exportTriplets,(python::arg("filename")="triplets"),"Export system matrix to a file with only non-zero entries.")
		.def("cholmodStats",&TemplateFlowEngine::cholmodStats,"get statistics of cholmod solver activity")
		.def("metisUsed",&TemplateFlowEngine::metisUsed,"check wether metis lib is effectively used")
		.add_property("forceMetis",&TemplateFlowEngine::getForceMetis,&TemplateFlowEngine::setForceMetis,"If true, METIS is used for matrix preconditioning, else Cholmod is free to choose the best method (which may be METIS to, depending on the matrix). See ``nmethods`` in Cholmod documentation")
		#endif
		.def("compTessVolumes",&TemplateFlowEngine::compTessVolumes,"Like TesselationWrapper::computeVolumes()")
		.def("volume",&TemplateFlowEngine::getVolume,(python::arg("id")=0),"Returns the volume of Voronoi's cell of a sphere.")
		)
};
// Definition of functions in a separate file for clarity 
#include<yade/pkg/pfv/FlowEngine.ipp>

class FlowCellInfo : public CGT::SimpleCellInfo {
	public:
	//For vector storage of all cells
	unsigned int index;
	int volumeSign;
	bool Pcondition;
	Real invVoidV;
	Real t;
	int fict;
 	Real volumeVariation;
	double pression;
	 //average relative (fluid - facet translation) velocity defined for a single cell as 1/Volume * SUM_ON_FACETS(x_average_facet*average_facet_flow_rate)
	CVector averageCellVelocity;
	// Surface vectors of facets, pointing from outside toward inside the cell
	std::vector<CVector> facetSurfaces;
	//Ratio between fluid surface and facet surface 
	std::vector<Real> facetFluidSurfacesRatio;
	// Reflects the geometrical property of the cell, so that the force by cell fluid on grain "i" is pressure*unitForceVectors[i]
	std::vector<CVector> unitForceVectors;
	// Store the area of triangle-sphere intersections for each facet (used in forces definition)
	std::vector<CVector> facetSphereCrossSections;
	std::vector<CVector> cellForce;
	std::vector<double> rayHydr;
	std::vector<double> modulePermeability;
	// Partial surfaces of spheres in the double-tetrahedron linking two voronoi centers. [i][j] is for sphere facet "i" and sphere facetVertices[i][j]. Last component for 1/sum_surfaces in the facet.doInterpolate
	double solidSurfaces [4][4];

	FlowCellInfo (void)
	{
		modulePermeability.resize(4, 0);
		cellForce.resize(4,CGAL::NULL_VECTOR);
		facetSurfaces.resize(4,CGAL::NULL_VECTOR);
		facetFluidSurfacesRatio.resize(4,0);
		facetSphereCrossSections.resize(4,CGAL::NULL_VECTOR);
		unitForceVectors.resize(4,CGAL::NULL_VECTOR);
		for (int k=0; k<4;k++) for (int l=0; l<3;l++) solidSurfaces[k][l]=0;
		rayHydr.resize(4, 0);
		invSumK=index=volumeSign=s=volumeVariation=pression=invVoidV=fict=0;
		isFictious=false; Pcondition = false; isGhost = false;
		isvisited = false; 		
		isGhost=false;
	}	
	bool isGhost;
	double invSumK;
	bool isvisited;
	
	inline Real& volume (void) {return t;}
	inline const Real& invVoidVolume (void) const {return invVoidV;}
	inline Real& invVoidVolume (void) {return invVoidV;}
	inline Real& dv (void) {return volumeVariation;}
	inline int& fictious (void) {return fict;}
	inline double& p (void) {return pression;}
	inline const double shiftedP (void) const {return pression;} //For compatibility with the periodic case
	inline const std::vector<double>& kNorm (void) const {return modulePermeability;}
	inline std::vector<double>& kNorm (void) {return modulePermeability;}
	inline std::vector< CVector >& facetSurf (void) {return facetSurfaces;}
	inline std::vector<CVector>& force (void) {return cellForce;}
	inline std::vector<double>& Rh (void) {return rayHydr;}
	inline CVector& averageVelocity (void) {return averageCellVelocity;}
	//used for transfering values between two triangulations, overload with more variables in derived classes (see e.g. SoluteFlow)
	inline void getInfo(const FlowCellInfo& otherCellInfo) {p()=otherCellInfo.shiftedP();} 
};

class FlowVertexInfo : public CGT::SimpleVertexInfo {
	CVector grainVelocity;
	Real volumeIncidentCells;
public:
	CVector forces;
	bool isGhost;
	FlowVertexInfo (void) {isGhost=false;}
	inline CVector& force (void) {return forces;}
	inline CVector& vel (void) {return grainVelocity;}
	inline Real& volCells (void) {return volumeIncidentCells;}
	inline const CVector ghostShift (void) const {return CGAL::NULL_VECTOR;}
};




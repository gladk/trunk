/*************************************************************************
*  Copyright (C) 2009 by Emanuele Catalano                               *
*  emanuele.catalano@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include<yade/core/PartialEngine.hpp>
#include<yade/pkg/dem/TriaxialCompressionEngine.hpp>
#include<yade/lib/triangulation/FlowBoundingSphere.hpp>
#include<yade/pkg/dem/TesselationWrapper.hpp>
#include<yade/lib/triangulation/PeriodicFlow.hpp>

class Flow;
class TesselationWrapper;
#ifdef LINSOLV
#define _FlowSolver CGT::FlowBoundingSphereLinSolv<CGT::FlowBoundingSphere<FlowTesselation> >
#else
#define _FlowSolver CGT::FlowBoundingSphere<FlowTesselation>
#endif


class FlowEngine : public PartialEngine
{
	public :
	typedef _FlowSolver							FlowSolver;
	typedef FlowTesselation							Tesselation;
	typedef FlowSolver::RTriangulation					RTriangulation;
	typedef FlowSolver::Finite_vertices_iterator                    	Finite_vertices_iterator;
	typedef FlowSolver::Finite_cells_iterator				Finite_cells_iterator;
	typedef FlowSolver::Cell_handle						Cell_handle;
	typedef RTriangulation::Finite_edges_iterator				Finite_edges_iterator;
	
	protected:
		shared_ptr<FlowSolver> solver;

	public :
		int retriangulationLastIter;
		enum {wall_left=0, wall_right, wall_bottom, wall_top, wall_back, wall_front};
		Vector3r normal [6];
		bool currentTes;
		int id_offset;
		double Eps_Vol_Cumulative;
		int ReTrg;
		template<class Solver>
		void Triangulate (Solver& flow);
		template<class Solver>
		void AddBoundary (Solver& flow);
		template<class Solver>
		void Build_Triangulation (double P_zero, Solver& flow);
		template<class Solver>
		void Build_Triangulation (Solver& flow);
		template<class Solver>
		void UpdateVolumes (Solver& flow);
		template<class Solver>
		void Initialize_volumes (Solver& flow);
		template<class Cellhandle>
		Real Volume_cell_single_fictious (Cellhandle cell);
		template<class Cellhandle>
		Real Volume_cell_double_fictious (Cellhandle cell);
		template<class Cellhandle>
		Real Volume_cell_triple_fictious (Cellhandle cell);
		template<class Cellhandle>
		Real Volume_cell (Cellhandle cell);
		void Oedometer_Boundary_Conditions();
		template<class Solver>
		void BoundaryConditions(Solver& flow);
		template<class Solver>
		void imposeFlux(Vector3r pos, Real flux,Solver& flow);
		template<class Solver>
		unsigned int imposePressure(Vector3r pos, Real p,Solver& flow);
		template<class Solver>
		void setImposedPressure(unsigned int cond, Real p,Solver& flow);
		template<class Solver>
		void clearImposedPressure(Solver& flow);
		void Average_real_cell_velocity();
		template<class Solver>
		void ApplyViscousForces(Solver& flow);
		template<class Solver>
		Real getFlux(unsigned int cond,Solver& flow);
		void saveVtk() {solver->saveVtk();}
		vector<Real> AvFlVelOnSph(unsigned int id_sph) {return solver->Average_Fluid_Velocity_On_Sphere(id_sph);}
		python::list getConstrictions() {
			vector<Real> csd=solver->getConstrictions(); python::list pycsd;
			for (unsigned int k=0;k<csd.size();k++) pycsd.append(csd[k]); return pycsd;}
		Vector3r fluidForce(unsigned int id_sph) {const CGT::Vecteur& f=solver->T[solver->currentTes].vertex(id_sph)->info().forces; return Vector3r(f[0],f[1],f[2]);}
		template<class Solver>
		Vector3r fluidShearForce(unsigned int id_sph, Solver& flow) {return (flow->viscousShearForces.size()>id_sph)?flow->viscousShearForces[id_sph]:Vector3r::Zero();}
		void setBoundaryVel(Vector3r vel) {topBoundaryVelocity=vel; Update_Triangulation=true;}
		void PressureProfile(double wallUpY, double wallDownY) {return solver->MeasurePressureProfile(wallUpY,wallDownY);}
		double MeasurePorePressure(double posX, double posY, double posZ){return solver->MeasurePorePressure(posX, posY, posZ);}
		double MeasureAveragedPressure(double posY){return solver->MeasureAveragedPressure(posY);}
		double MeasureTotalAveragedPressure(){return solver->MeasureTotalAveragedPressure();}
		
		//Instanciation of templates for python binding
		Vector3r 	_fluidShearForce(unsigned int id_sph) {return fluidShearForce(id_sph,solver);}
		void 		_imposeFlux(Vector3r pos, Real flux) {return imposeFlux(pos,flux,*solver);}
		unsigned int 	_imposePressure(Vector3r pos, Real p) {return imposePressure(pos,p,solver);}	
		void 		_setImposedPressure(unsigned int cond, Real p) {setImposedPressure(cond,p,solver);}
		void 		_clearImposedPressure() {clearImposedPressure(solver);}
		Real 		_getFlux(unsigned int cond) {return getFlux(cond,solver);}

		virtual ~FlowEngine();

		virtual void action();

		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(FlowEngine,PartialEngine,"An engine to solve flow problem in saturated granular media",
					((bool,isActivated,true,,"Activates Flow Engine"))
					((bool,first,true,,"Controls the initialization/update phases"))
					((double, fluidBulkModulus, 0.,,"Bulk modulus of fluid (inverse of compressibility) K=-dP*V/dV [Pa]. Flow is compressible if fluidBulkModulus > 0, else incompressible."))
					((Real, dt, 0,,"timestep [s]"))
// 					((bool,save_vtk,false,,"Enable/disable vtk files creation for visualization"))
					((bool,permeability_map,false,,"Enable/disable stocking of average permeability scalar in cell infos."))
					((bool, save_mgpost, false,,"Enable/disable mgpost file creation"))
					((bool, slice_pressures, false, ,"Enable/Disable slice pressure measurement"))
					((bool, velocity_profile, false, ,"Enable/Disable slice velocity measurement"))
					((bool, consolidation,false,,"Enable/Disable storing consolidation files"))
					((bool, slip_boundary, true,, "Controls friction condition on lateral walls"))
					((bool,WaveAction, false,, "Allow sinusoidal pressure condition to simulate ocean waves"))
					((double, Sinus_Amplitude, 0,, "Pressure value (amplitude) when sinusoidal pressure is applied"))
					((double, Sinus_Average, 0,,"Pressure value (average) when sinusoidal pressure is applied"))
					((bool, CachedForces, true,,"Des/Activate the cached forces calculation"))
					((bool, Debug, false,,"Activate debug messages"))
					((double, wall_thickness,0.001,,"Walls thickness"))
					((double,P_zero,0,,"Initial internal pressure for oedometer test"))
					((double,Tolerance,1e-06,,"Gauss-Seidel Tolerance"))
					((double,Relax,1.9,,"Gauss-Seidel relaxation"))
					((bool, Update_Triangulation, 0,,"If true the medium is retriangulated"))
					((int,PermuteInterval,100000,,"Pore space re-triangulation period"))
					((double, eps_vol_max, 0,,"Maximal absolute volumetric strain computed at each iteration"))
					((double, EpsVolPercent_RTRG,0.01,,"Percentuage of cumulate eps_vol at which retriangulation of pore space is performed"))
					((double, porosity, 0,,"Porosity computed at each retriangulation"))
					((bool,compute_K,false,,"Activates permeability measure within a granular sample"))
					((bool,meanK_correction,true,,"Local permeabilities' correction through meanK threshold"))
					((bool,meanK_opt,false,,"Local permeabilities' correction through an optimized threshold"))
					((double,permeability_factor,0.0,,"permability multiplier"))
					((double,viscosity,1.0,,"viscosity of fluid"))
					((Real,loadFactor,1.1,,"Load multiplicator for oedometer test"))
					((double, K, 0,, "Permeability of the sample"))
					((int, useSolver, 0,, "Solver to use 0=G-Seidel, 1=Taucs, 2-Pardiso"))
// 					((std::string,key,"",,"A string appended at the output files, use it to name simulations."))
					((double, V_d, 0,,"darcy velocity of fluid in sample"))
					((bool, Flow_imposed_TOP_Boundary, true,, "if false involve pressure imposed condition"))
					((bool, Flow_imposed_BOTTOM_Boundary, true,, "if false involve pressure imposed condition"))
					((bool, Flow_imposed_FRONT_Boundary, true,, "if false involve pressure imposed condition"))
					((bool, Flow_imposed_BACK_Boundary, true,, "if false involve pressure imposed condition"))
					((bool, Flow_imposed_LEFT_Boundary, true,, "if false involve pressure imposed condition"))
					((bool, Flow_imposed_RIGHT_Boundary, true,,"if false involve pressure imposed condition"))
					((double, Pressure_TOP_Boundary, 0,, "Pressure imposed on top boundary"))
					((double, Pressure_BOTTOM_Boundary,  0,, "Pressure imposed on bottom boundary"))
					((double, Pressure_FRONT_Boundary,  0,, "Pressure imposed on front boundary"))
					((double, Pressure_BACK_Boundary,  0,,"Pressure imposed on back boundary"))
					((double, Pressure_LEFT_Boundary,  0,, "Pressure imposed on left boundary"))
					((double, Pressure_RIGHT_Boundary,  0,, "Pressure imposed on right boundary"))
					((Vector3r, topBoundaryVelocity, Vector3r::Zero(),, "velocity on top boundary, only change it using :yref:`FlowEngine::setBoundaryVel`"))
					((int, wallTopId,3,,"Id of top boundary (default value is ok if aabbWalls are appended BEFORE spheres.)"))
					((int, wallBottomId,2,,"Id of bottom boundary (default value is ok if aabbWalls are appended BEFORE spheres.)"))
					((int, wallFrontId,5,,"Id of front boundary (default value is ok if aabbWalls are appended BEFORE spheres.)"))
					((int, wallBackId,4,,"Id of back boundary (default value is ok if aabbWalls are appended BEFORE spheres.)"))
					((int, wallLeftId,0,,"Id of left boundary (default value is ok if aabbWalls are appended BEFORE spheres.)"))
					((int, wallRightId,1,,"Id of right boundary (default value is ok if aabbWalls are appended BEFORE spheres.)"))
					((Vector3r, id_force, 0,, "Fluid force acting on sphere with id=flow.id_sphere"))
					((bool, BOTTOM_Boundary_MaxMin, 1,,"If true bounding sphere is added as function fo max/min sphere coord, if false as function of yade wall position"))
					((bool, TOP_Boundary_MaxMin, 1,,"If true bounding sphere is added as function fo max/min sphere coord, if false as function of yade wall position"))
					((bool, RIGHT_Boundary_MaxMin, 1,,"If true bounding sphere is added as function fo max/min sphere coord, if false as function of yade wall position"))
					((bool, LEFT_Boundary_MaxMin, 1,,"If true bounding sphere is added as function fo max/min sphere coord, if false as function of yade wall position"))
					((bool, FRONT_Boundary_MaxMin, 1,,"If true bounding sphere is added as function fo max/min sphere coord, if false as function of yade wall position"))
					((bool, BACK_Boundary_MaxMin, 1,,"If true bounding sphere is added as function fo max/min sphere coord, if false as function of yade wall position"))
					((bool, areaR2Permeability, 1,,"Use corrected formula for permeabilities calculation in flowboundingsphere (areaR2permeability variable)"))
					((bool, viscousShear, false,,"Compute viscous shear terms as developped by Donia Marzougui"))
					,,
					timingDeltas=shared_ptr<TimingDeltas>(new TimingDeltas);
					for (int i=0; i<6; ++i){normal[i]=Vector3r::Zero();}
					normal[wall_bottom].y()=normal[wall_left].x()=normal[wall_back].z()=1;
					normal[wall_top].y()=normal[wall_right].x()=normal[wall_front].z()=-1;					
					solver = shared_ptr<FlowSolver> (new FlowSolver);
					first=true;
					Update_Triangulation=false;
					eps_vol_max=Eps_Vol_Cumulative=retriangulationLastIter=0;
					ReTrg=1;
					,
					.def("imposeFlux",&FlowEngine::_imposeFlux,(python::arg("pos"),python::arg("p")),"Impose incoming flux in boundary cell of location 'pos'.")
					.def("imposePressure",&FlowEngine::_imposePressure,(python::arg("pos"),python::arg("p")),"Impose pressure in cell of location 'pos'. The index of the condition is returned (for multiple imposed pressures at different points).")
					.def("setImposedPressure",&FlowEngine::_setImposedPressure,(python::arg("cond"),python::arg("p")),"Set pressure value at the point indexed 'cond'.")
					.def("clearImposedPressure",&FlowEngine::_clearImposedPressure,"Clear the list of points with pressure imposed.")
					.def("getFlux",&FlowEngine::_getFlux,(python::arg("cond")),"Get influx in cell associated to an imposed P (indexed using 'cond').")
					.def("getConstrictions",&FlowEngine::getConstrictions,"Get the list of constrictions (inscribed circle) for all finite facets.")
					.def("saveVtk",&FlowEngine::saveVtk,"Save pressure field in vtk format.")
					.def("AvFlVelOnSph",&FlowEngine::AvFlVelOnSph,(python::arg("Id_sph")),"Compute a sphere-centered average fluid velocity")
					.def("fluidForce",&FlowEngine::fluidForce,(python::arg("Id_sph")),"Return the fluid force on sphere Id_sph.")
					.def("fluidShearForce",&FlowEngine::_fluidShearForce,(python::arg("Id_sph")),"Return the viscous shear force on sphere Id_sph.")
					.def("setBoundaryVel",&FlowEngine::setBoundaryVel,(python::arg("vel")),"Change velocity on top boundary.")
					.def("PressureProfile",&FlowEngine::PressureProfile,(python::arg("wallUpY"),python::arg("wallDownY")),"Measure pore pressure in 6 equally-spaced points along the height of the sample")
					.def("MeasurePorePressure",&FlowEngine::MeasurePorePressure,(python::arg("posX"),python::arg("posY"),python::arg("posZ")),"Measure pore pressure in position pos[0],pos[1],pos[2]")
					.def("MeasureAveragedPressure",&FlowEngine::MeasureAveragedPressure,(python::arg("posY")),"Measure slice-averaged pore pressure at height posY")
					.def("MeasureTotalAveragedPressure",&FlowEngine::MeasureTotalAveragedPressure,"Measure averaged pore pressure in the entire volume")
					)
		DECLARE_LOGGER;
};


template<class Solver>
unsigned int FlowEngine::imposePressure(Vector3r pos, Real p,Solver& flow)
{
	if (!flow) LOG_ERROR("no flow defined yet, run at least one iter");
	flow->imposedP.push_back( pair<CGT::Point,Real>(CGT::Point(pos[0],pos[1],pos[2]),p) );
	//force immediate update of boundary conditions
	Update_Triangulation=true;
	return flow->imposedP.size()-1;
}



REGISTER_SERIALIZABLE(FlowEngine);

#ifdef LINSOLV
#define _PeriFlowSolver CGT::PeriodicFlowLinSolv
#else
#define _PeriFlowSolver CGT::PeriodicFlow
#endif

class PeriodicFlowEngine : public FlowEngine
{
	public :
		public :
		typedef _PeriFlowSolver							FlowSolver;
		typedef PeriFlowTesselation						Tesselation;
		typedef FlowSolver::RTriangulation					RTriangulation;
		typedef FlowSolver::Finite_vertices_iterator                    	Finite_vertices_iterator;
		typedef FlowSolver::Finite_cells_iterator				Finite_cells_iterator;
		typedef FlowSolver::Cell_handle						Cell_handle;
		typedef RTriangulation::Finite_edges_iterator				Finite_edges_iterator;
		typedef RTriangulation::Vertex_handle					Vertex_handle;
		
		shared_ptr<FlowSolver> solver;
		
		void Triangulate ();
// 		void AddBoundary ();
		void Build_Triangulation (Real pzero);
		void Initialize_volumes ();
		void UpdateVolumes ();
		Real Volume_cell (Cell_handle cell);

		Real Volume_cell_single_fictious (Cell_handle cell);
		inline void locateCell(Cell_handle baseCell, unsigned int& index, unsigned int count=0);
		Vector3r meanVelocity();
		
		virtual ~PeriodicFlowEngine();

		virtual void action();
		
		//Instanciation of templates for python binding
		Vector3r 	_fluidShearForce(unsigned int id_sph) {return fluidShearForce(id_sph,solver);}
		void 		_imposeFlux(Vector3r pos, Real flux) {return imposeFlux(pos,flux,*solver);}
		unsigned int 	_imposePressure(Vector3r pos, Real p) {return imposePressure(pos,p,this->solver);}	
// 		void 		_setImposedPressure(unsigned int cond, Real p) {setImposedPressure(cond,p,solver);}
// 		void 		_clearImposedPressure() {clearImposedPressure(solver);}
// 		Real 		_getFlux(unsigned int cond) {getFlux(cond,solver);}

		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(PeriodicFlowEngine,FlowEngine,"An engine to solve flow problem in saturated granular media",
			((Real,duplicateThreshold, 0.06,,"distance from cell borders that will triger periodic duplication in the triangulation |yupdate|"))
			((Vector3r, gradP, Vector3r::Zero(),,"Macroscopic pressure gradient"))
			,,
			wallTopId=wallBottomId=wallFrontId=wallBackId=wallLeftId=wallRightId=-1;
// 			FlowEngine::solver=shared_ptr<FlowSolver>;
			solver = shared_ptr<FlowSolver> (new FlowSolver);
			Update_Triangulation=false;
			eps_vol_max=Eps_Vol_Cumulative=retriangulationLastIter=0;
			ReTrg=1;
			,
			.def("meanVelocity",&PeriodicFlowEngine::meanVelocity,"measure the mean velocity in the period")
// 			.def("imposeFlux",&FlowEngine::_imposeFlux,(python::arg("pos"),python::arg("p")),"Impose incoming flux in boundary cell of location 'pos'.")
			.def("imposePressure",&PeriodicFlowEngine::_imposePressure,(python::arg("pos"),python::arg("p")),"Impose pressure in cell of location 'pos'. The index of the condition is returned (for multiple imposed pressures at different points).")
// 			.def("setImposedPressure",&FlowEngine::_setImposedPressure,(python::arg("cond"),python::arg("p")),"Set pressure value at the point indexed 'cond'.")
// 			.def("clearImposedPressure",&FlowEngine::_clearImposedPressure,"Clear the list of points with pressure imposed.")
// 			.def("getFlux",&FlowEngine::_getFlux,(python::arg("cond")),"Get influx in cell associated to an imposed P (indexed using 'cond').")
		)
		DECLARE_LOGGER;


};

REGISTER_SERIALIZABLE(PeriodicFlowEngine);


// #endif

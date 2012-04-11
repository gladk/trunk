/*************************************************************************
*  Copyright (C) 2010 by Emanuele Catalano <catalano@grenoble-inp.fr>    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef FLOW_ENGINE

#ifndef _FLOWBOUNDINGSPHERE_H
#define _FLOWBOUNDINGSPHERE_H

#include "Network.hpp"
#include "Timer.h"
#include "basicVTKwritter.hpp"
#include "Timer.h"

#ifdef XVIEW
#include "Vue3D.h" //FIXME implicit dependencies will look for this class (out of tree) even ifndef XVIEW
#endif

using namespace std;
namespace CGT {

template<class _Tesselation>
class FlowBoundingSphere : public Network<_Tesselation>
{
	public:
		typedef _Tesselation			Tesselation;
		typedef Network<Tesselation>		_N;
		DECLARE_TESSELATION_TYPES(Network<Tesselation>)
		
		//painfull, but we need that for templates inheritance...
		using _N::T; using _N::x_min; using _N::x_max; using _N::y_min; using _N::y_max; using _N::z_min; using _N::z_max; using _N::Rmoy; using _N::SectionArea; using _N::Height; using _N::Vtotale; using _N::currentTes; using _N::DEBUG_OUT; using _N::nOfSpheres; using _N::x_min_id; using _N::x_max_id; using _N::y_min_id; using _N::y_max_id; using _N::z_min_id; using _N::z_max_id; using _N::boundsIds; using _N::Corner_min; using _N::Corner_max;  using _N::Vsolid_tot; using _N::Vtotalissimo; using _N::Vporale; using _N::Ssolid_tot; using _N::V_porale_porosity; using _N::V_totale_porosity; using _N::boundaries; using _N::id_offset; using _N::vtk_infinite_vertices; using _N::vtk_infinite_cells; using _N::num_particles; using _N::fictious_vertex; using _N::boundingCells;
		//same for functions
		using _N::Define_fictious_cells; using _N::AddBoundingPlanes; using _N::boundary;
		
		virtual ~FlowBoundingSphere();
 		FlowBoundingSphere();

		bool SLIP_ON_LATERALS;
		bool areaR2Permeability;
		double TOLERANCE;
		double RELAX;
		double ks; //Hydraulic Conductivity
		bool meanK_LIMIT, meanK_STAT, distance_correction;
		bool OUTPUT_BOUDARIES_RADII;
		bool noCache;//flag for checking if cached values cell->unitForceVectors have been defined
		//Handling imposed pressures on elements in the form of {point,value} pairs, IPCells contains the cell handles corresponding to point
		vector<pair<Point,Real> > imposedP;
		vector<Cell_handle> IPCells;
		
		void initNewTri () {noCache=true; /*isLinearSystemSet=false; areCellsOrdered=false;*/}//set flags after retriangulation
		bool permeability_map;

		bool computeAllCells;//exececute computeHydraulicRadius for all facets and all spheres (double cpu time but needed for now in order to define crossSections correctly)
		double K_opt_factor;
		int Iterations;

		bool RAVERAGE;
		int walls_id[6];
		vector <double> Edge_Surfaces;
		vector <pair<int,int> > Edge_ids;
		vector <Real> Edge_HydRad;
		vector <Vector3r> Edge_normal;
		vector <Vector3r> viscousShearForces;
	
		void mplot ( char *filename);
		void Localize();
		void Compute_Permeability();
		virtual void GaussSeidel (Real dt=0);
		virtual void ResetNetwork();

		double MeasurePorePressure (double X, double Y, double Z);
		void Fictious_cells ( );

		double k_factor; //permeability moltiplicator
		std::string key; //to give to consolidation files a name with iteration number
		std::vector<double> Pressures; //for automatic write maximum pressures during consolidation
		bool tess_based_force; //allow the force computation method to be chosen from FlowEngine
		Real minPermLength; //min branch length for Poiseuille

		double P_SUP, P_INF, P_INS, VISCOSITY;
		bool compressible;
		double fluidBulkModulus;
		
		Tesselation& Compute_Action ( );
		Tesselation& Compute_Action ( int argc, char *argv[ ], char *envp[ ] );
		Tesselation& LoadPositions(int argc, char *argv[ ], char *envp[ ]);
		void SpheresFileCreator ();
		void DisplayStatistics();
		void Initialize_pressures ( double P_zero );
		void reApplyBoundaryConditions ();
		/// Define forces using the same averaging volumes as for permeability
		void ComputeTetrahedralForces();
		/// Define forces spliting drag and buoyancy terms
		void ComputeFacetForces();
		void ComputeFacetForcesWithCache();
		void saveVtk ( );
		void MGPost ( );
#ifdef XVIEW
		void Dessine_Triangulation ( Vue3D &Vue, RTriangulation &T );
		void Dessine_Short_Tesselation ( Vue3D &Vue, Tesselation &Tes );
#endif
		double Permeameter ( double P_Inf, double P_Sup, double Section, double DeltaY, const char *file );
		double Sample_Permeability( double& x_Min,double& x_Max ,double& y_Min,double& y_Max,double& z_Min,double& z_Max);
		double Compute_HydraulicRadius (Cell_handle cell, int j );

		double dotProduct ( Vecteur x, Vecteur y );
		double Compute_EffectiveRadius(Cell_handle cell, int j);
		double Compute_EquivalentRadius(Cell_handle cell, int j);
		//return the list of constriction values
		vector<double> getConstrictions();

		void GenerateVoxelFile ( );
		
		void ComputeEdgesSurfaces();
		Vector3r ComputeViscousForce(Vector3r deltaV, int edge_id);

		RTriangulation& Build_Triangulation ( Real x, Real y, Real z, Real radius, unsigned const id );

		bool isInsideSphere ( double& x, double& y, double& z );

		void SliceField (const char *filename);
		void ComsolField();

		void Interpolate ( Tesselation& Tes, Tesselation& NewTes );
		void Average_Relative_Cell_Velocity();
		void Average_Fluid_Velocity();
		void ApplySinusoidalPressure(RTriangulation& Tri, double Amplitude, double Average_Pressure, double load_intervals);
		bool isOnSolid  (double X, double Y, double Z);
		void MeasurePressureProfile(double Wall_up_y, double Wall_down_y);
		double MeasureAveragedPressure(double Y);
		double MeasureTotalAveragedPressure();
		
		vector<Real> Average_Fluid_Velocity_On_Sphere(unsigned int Id_sph);
		//Solver?
		int useSolver;//(0 : GaussSeidel, 1 : TAUCS, 2 : PARDISO)
};

} //namespace CGT

#ifdef LINSOLV
#include "FlowBoundingSphereLinSolv.hpp"
#endif

/// _____ Template Implementation ____
#include "FlowBoundingSphere.ipp"

#endif //FLOW_ENGINE

#endif

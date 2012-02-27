#ifndef _FLOWBOUNDINGSPHERE_H
#define _FLOWBOUNDINGSPHERE_H

#include "Operations.h"
#include "Timer.h"
#include "Tesselation.h"
#include "basicVTKwritter.hpp"
#include "Timer.h"
#include "stdafx.h"
#include "Empilement.h"

#ifdef XVIEW
//#include "Vue3D.h" //FIXME implicit dependencies will look for this class (out of tree) even ifndef XVIEW
#endif

#ifdef FLOW_ENGINE


namespace CGT{

struct Boundary
{
	Point p;
	Vecteur normal;
	int coordinate;
	bool flowCondition;//flowCondition=0, pressure is imposed // flowCondition=1, flow is imposed
	Real value;
};

class FlowBoundingSphere
{
	public:
		virtual ~FlowBoundingSphere();
 		FlowBoundingSphere();
		
		int x_min_id, x_max_id, y_min_id, y_max_id, z_min_id, z_max_id;
		int boundsIds [6];
		bool currentTes;
		bool SLIP_ON_LATERALS;
		double TOLERANCE;
		double RELAX;
		double ks; //Hydraulic Conductivity
		bool meanK_LIMIT, meanK_STAT, distance_correction;
		double K_opt_factor;
		int Iterations;
		
		Boundary boundaries [6];
		int walls_id[6];
		short id_offset;
 		Boundary& boundary (int b) {return boundaries[b-id_offset];}
		
		void mplot (RTriangulation& Tri, char *filename);
		void Localize ();
		void Compute_Permeability();
		
		void DisplayStatistics();
		void GaussSeidel ( );
		void Compute_Forces ();
		double surface_external_triple_fictious ( Real position[3], Cell_handle cell, Boundary b );
		void Fictious_cells ( );

		Tesselation T [2];
		
		double x_min, x_max, y_min, y_max, z_min, z_max, Rmoy;
		Real Vsolid_tot, Vtotalissimo, Vporale, Ssolid_tot;
		double k_factor; //permeability moltiplicator
		std::string key; //to give to consolidation files a name with iteration number
		std::vector<double> Pressures; //for automatic write maximum pressures during consolidation
		bool tess_based_force; //allow the force computation method to be chosen from FlowEngine
		Real minPermLength; //min branch length for Poiseuille
		
		double P_SUP, P_INF, P_INS;
		
		void AddBoundingPlanes ( Tesselation& Tes, double x_Min,double x_Max ,double y_Min,double y_Max,double z_Min,double z_Max );
		void AddBoundingPlanes(int y_min_id, int y_max_id, int x_min_id, int x_max_id, int z_min_id, int z_max_id);
		void AddBoundingPlanes();
		void AddBoundingPlanes(Real center[3], Real Extents[3], int id);
		
		void Compute_Action ( );
		void Compute_Action ( int argc, char *argv[ ], char *envp[ ] );
		void DisplayStatistics ( RTriangulation& Tri );
		Vecteur external_force_single_fictious ( Cell_handle cell );
		void SpheresFileCreator ();
// 		void Analytical_Consolidation ( );
		
		void Boundary_Conditions ( RTriangulation& Tri );
		void Initialize_pressures ( double P_zero );
		/// Define forces using the same averaging volumes as for permeability
		void ComputeTetrahedralForces();
		void save_vtk_file ( RTriangulation &T );
		void MGPost ( RTriangulation& Tri );
#ifdef XVIEW
		void Dessine_Triangulation ( Vue3D &Vue, RTriangulation &T );
		void Dessine_Short_Tesselation ( Vue3D &Vue, Tesselation &Tes );
#endif
		double Permeameter ( RTriangulation& Tri, double P_Inf, double P_Sup, double Section, double DeltaY, char *file );
		double Sample_Permeability ( RTriangulation& Tri, double x_Min,double x_Max ,double y_Min,double y_Max,double z_Min,double z_Max, std::string key );
		double Compute_HydraulicRadius ( RTriangulation& Tri, Cell_handle cell, int j );
		double PermeameterCurve ( RTriangulation& Tri, char *filename, Real time, int intervals );

		double dotProduct ( Vecteur x, Vecteur y );

// 		double crossProduct ( double x[3], double y[3] );

		double surface_solid_facet ( Sphere ST1, Sphere ST2, Sphere ST3 );
		Vecteur surface_double_fictious_facet ( Vertex_handle fSV1, Vertex_handle fSV2, Vertex_handle SV3 );
		Vecteur surface_single_fictious_facet ( Vertex_handle fSV1, Vertex_handle SV2, Vertex_handle SV3 );
		
		double surface_solid_fictious_facet ( Vertex_handle ST1, Vertex_handle ST2, Vertex_handle ST3 );
		
		double surface_external_triple_fictious (Cell_handle cell, Boundary b );
		
		double surface_external_double_fictious ( Cell_handle cell, Boundary b );
		
		double surface_external_single_fictious ( Cell_handle cell, Boundary b );

		void GenerateVoxelFile ( RTriangulation& Tri );

		RTriangulation& Build_Triangulation ( Real x, Real y, Real z, Real radius, unsigned const id );

		void Build_Tessalation ( RTriangulation& Tri );

		double spherical_triangle_area ( Sphere STA1, Sphere STA2, Sphere STA3, Point PTA1 );
		
		double fast_spherical_triangle_area ( const Sphere& STA1, const Point& STA2, const Point& STA3, const Point& PTA1 );
// 		Real solid_angle ( const Point& STA1, const Point& STA2, const Point& STA3, const Point& PTA1 );
		double spherical_triangle_volume ( const Sphere& ST1, const Point& PT1, const Point& PT2, const Point& PT3 );
		Real fast_solid_angle ( const Point& STA1, const Point& PTA1, const Point& PTA2, const Point& PTA3 );
		
		bool isInsideSphere ( RTriangulation& Tri, double x, double y, double z );
		
		void SliceField ( RTriangulation& Tri );

		void Interpolate ( Tesselation& Tes, Tesselation& NewTes );
		
		double volume_single_fictious_pore ( Vertex_handle SV1, Vertex_handle SV2, Vertex_handle SV3, Point PV1 );
		//Fast version, assign surface of facet for future forces calculations (pointing from PV2 to PV1)
		double volume_single_fictious_pore ( const Vertex_handle& SV1, const Vertex_handle& SV2, const Vertex_handle& SV3, const Point& PV1,  const Point& PV2, Vecteur& facetSurface);
		double volume_double_fictious_pore ( Vertex_handle SV1, Vertex_handle SV2, Vertex_handle SV3, Point PV1 );
		//Fast version, assign surface of facet for future forces calculations (pointing from PV2 to PV1)
		double volume_double_fictious_pore (Vertex_handle SV1, Vertex_handle SV2, Vertex_handle SV3, Point& PV1, Point& PV2, Vecteur& facetSurface);
};

} //namespace CGT
#endif //FLOW_ENGINE

#endif

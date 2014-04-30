/*************************************************************************
*  Copyright (C) 2009 by Emanuele Catalano <catalano@grenoble-inp.fr>    *
*  Copyright (C) 2009 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*  Copyright (C) 2012 by Donia Marzougui <donia.marzougui@grenoble-inp.fr>*
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef FLOW_ENGINE

#pragma once

#include "Network.hpp"
#include "Timer.h"
#include "basicVTKwritter.hpp"
#include "Timer.h"

typedef pair<pair<int,int>, vector<double> > Constriction;

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
		using _N::T; using _N::xMin; using _N::xMax; using _N::yMin; using _N::yMax; using _N::zMin; using _N::zMax; using _N::Rmoy; using _N::sectionArea; using _N::Height; using _N::vTotal; using _N::currentTes; using _N::debugOut; using _N::nOfSpheres; using _N::xMinId; using _N::xMaxId; using _N::yMinId; using _N::yMaxId; using _N::zMinId; using _N::zMaxId; using _N::boundsIds; using _N::cornerMin; using _N::cornerMax;  using _N::VSolidTot; using _N::Vtotalissimo; using _N::vPoral; using _N::sSolidTot; using _N::vPoralPorosity; using _N::vTotalPorosity; using _N::boundaries; using _N::idOffset; using _N::vtkInfiniteVertices; using _N::vtkInfiniteCells; using _N::num_particles; using _N::boundingCells; using _N::facetVertices; using _N::facetNFictious;
		//same for functions
		using _N::defineFictiousCells; using _N::addBoundingPlanes; using _N::boundary;

		virtual ~FlowBoundingSphere();
 		FlowBoundingSphere();

		bool slipBoundary;
		double tolerance;
		double relax;
		double ks; //Hydraulic Conductivity
		bool clampKValues, meanKStat, distanceCorrection;
		bool OUTPUT_BOUDARIES_RADII;
		bool noCache;//flag for checking if cached values cell->unitForceVectors have been defined
		bool computedOnce;//flag for checking if current triangulation has been computed at least once
		bool pressureChanged;//are imposed pressures modified (on python side)? When it happens, we have to reApplyBoundaryConditions
		int errorCode;
		
		//Handling imposed pressures/fluxes on elements in the form of {point,value} pairs, IPCells contains the cell handles corresponding to point
		vector<pair<Point,Real> > imposedP;
		vector<CellHandle> IPCells;
		vector<pair<Point,Real> > imposedF;
		vector<CellHandle> IFCells;
		//Pointers to vectors used for user defined boundary pressure
		vector<Real> *pxpos, *ppval;
		
		void initNewTri () {noCache=true; /*isLinearSystemSet=false; areCellsOrdered=false;*/}//set flags after retriangulation
		bool permeabilityMap;

		bool computeAllCells;//exececute computeHydraulicRadius for all facets and all spheres (double cpu time but needed for now in order to define crossSections correctly)
		double KOptFactor;
		double minKdivKmean;
		double maxKdivKmean;
		int Iterations;

		bool rAverage;
		int walls_id[6];
		#define parallel_forces
		#ifdef parallel_forces
		int ompThreads;
		vector< vector<const CVector*> > perVertexUnitForce;
		vector< vector<const Real*> > perVertexPressure;
		#endif
		vector <double> edgeSurfaces;
		vector <pair<const VertexInfo*,const VertexInfo*> > edgeIds;
		vector <Real> edgeNormalLubF;
		vector <Vector3r> shearLubricationForces;
		vector <Vector3r> shearLubricationTorques;
		vector <Vector3r> pumpLubricationTorques;
		vector <Vector3r> twistLubricationTorques;
		vector <Vector3r> normalLubricationForce;
		vector <Matrix3r> shearLubricationBodyStress;
		vector <Matrix3r> normalLubricationBodyStress;
		vector <Vector3r> deltaNormVel;
		vector <Vector3r> deltaShearVel;
		vector <Vector3r> normalV;
		vector <Real> surfaceDistance;
		vector <int> onlySpheresInteractions;
		vector <Matrix3r> shearStressInteraction;
		vector <Matrix3r> normalStressInteraction;
		
		void Localize();
		void computePermeability();
		virtual void gaussSeidel (Real dt=0);
		virtual void resetNetwork();


		double kFactor; //permeability moltiplicator
		std::string key; //to give to consolidation files a name with iteration number
// 		std::vector<double> pressures; //for automatic write maximum pressures during consolidation
		bool tessBasedForce; //allow the force computation method to be chosen from FlowEngine
		Real minPermLength; //min branch length for Poiseuille

		double viscosity;
		double fluidBulkModulus;
		
		void displayStatistics();
		void initializePressure ( double pZero );
		bool reApplyBoundaryConditions ();
		void computeFacetForcesWithCache(bool onlyCache=false);
		void saveVtk (const char* folder);
#ifdef XVIEW
		void dessineTriangulation ( Vue3D &Vue, RTriangulation &T );
		void dessineShortTesselation ( Vue3D &Vue, Tesselation &Tes );
#endif
		double permeameter ( double PInf, double PSup, double Section, double DeltaY, const char *file );
		double samplePermeability( double& xMin,double& xMax ,double& yMin,double& yMax,double& zMin,double& zMax);
		double computeHydraulicRadius (CellHandle cell, int j );
		Real checkSphereFacetOverlap(const Sphere& v0, const Sphere& v1, const Sphere& v2);

		double dotProduct ( CVector x, CVector y );
		double computeEffectiveRadius(CellHandle cell, int j);
		double computeEquivalentRadius(CellHandle cell, int j);
		//return the list of constriction values
		vector<double> getConstrictions();
		vector<Constriction> getConstrictionsFull();

		void generateVoxelFile ( );
		
		void computeEdgesSurfaces();
		Vector3r computeViscousShearForce(const Vector3r& deltaV, const int& edge_id, const Real& Rh);
		Real computeNormalLubricationForce(const Real& deltaNormV, const Real& dist, const int& edge_id, const Real& eps, const Real& stiffness, const Real& dt, const Real& meanRad);
		Vector3r computeShearLubricationForce(const Vector3r& deltaShearV, const Real& dist, const int& edge_id, const Real& eps, const Real& centerDist, const Real& meanRad);
		Vector3r computePumpTorque(const Vector3r& deltaShearAngV, const Real& dist, const int& edge_id, const Real& eps, const Real& meanRad );
		Vector3r computeTwistTorque(const Vector3r& deltaNormAngV, const Real& dist, const int& edge_id, const Real& eps, const Real& meanRad );



		RTriangulation& buildTriangulation ( Real x, Real y, Real z, Real radius, unsigned const id );

		bool isInsideSphere ( double& x, double& y, double& z );

		void sliceField (const char *filename);
		void comsolField();

		void interpolate ( Tesselation& Tes, Tesselation& NewTes );
		virtual void averageRelativeCellVelocity();
		void averageFluidVelocity();
		void applySinusoidalPressure(RTriangulation& Tri, double amplitude, double averagePressure, double loadIntervals);
		void applyUserDefinedPressure(RTriangulation& Tri, vector<Real>& xpos, vector<Real>& pval);
		bool isOnSolid  (double X, double Y, double Z);
		double getPorePressure (double X, double Y, double Z);
		void measurePressureProfile(double WallUpy, double WallDowny);
		double averageSlicePressure(double Y);
		double averagePressure();
		double getCell (double X,double Y,double Z);
		double boundaryFlux(unsigned int boundaryId);
		
		vector<Real> averageFluidVelocityOnSphere(unsigned int Id_sph);
		//Solver?
		int useSolver;//(0 : GaussSeidel, 1:CHOLMOD)
};

} //namespace CGT
#include <yade/lib/triangulation/FlowBoundingSphere.ipp>
#ifdef LINSOLV
#include "yade/lib/triangulation/FlowBoundingSphereLinSolv.hpp"
#endif

/// _____ Template Implementation ____
// #include "yade/lib/triangulation/FlowBoundingSphereLinSolv.ipp"

#endif //FLOW_ENGINE

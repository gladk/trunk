/*************************************************************************
*  Copyright (C) 2006 by Bruno Chareyre,  bruno.chareyre@hmg.inpg.fr     *
*  Copyright (C) 2010 by Emanuele Catalano, emanuele.catalano@hmg.inpg.fr*
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

//Define basic types from CGAL templates

#ifndef _Def_types
#define _Def_types

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/circulator.h>
#include <CGAL/number_utils.h>
#include <boost/static_assert.hpp>

namespace CGT {
//Robust kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//A bit faster, but gives crash eventualy
// typedef CGAL::Cartesian<double> K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   				Traits;
typedef K::Point_3									Point;
typedef Traits::Vector_3 								Vecteur;
typedef Traits::Segment_3								Segment;
#ifndef NO_REAL_CHECK
/** compilation inside yade: check that Real in yade is the same as Real we will define; otherwise it might make things go wrong badly (perhaps) **/
BOOST_STATIC_ASSERT(sizeof(Traits::RT)==sizeof(Real));
#endif
typedef Traits::RT									Real; //Dans cartesian, RT = FT
typedef Traits::Weighted_point								Sphere;
typedef Traits::Line_3									Droite;
typedef Traits::Plane_3									Plan;
typedef Traits::Triangle_3								Triangle;
typedef Traits::Tetrahedron_3								Tetraedre;

class SimpleCellInfo : public Point {
	public:
	Real s;
	bool isFictious;
	SimpleCellInfo (void) {isFictious=false; s=0;}
	SimpleCellInfo& operator= (const Point &p) { Point::operator= (p); return *this; }
	SimpleCellInfo& operator= (const float &scalar) { s=scalar; return *this; }
	inline Real x (void) {return Point::x();}
	inline Real y (void) {return Point::y();}
	inline Real z (void) {return Point::z();}
	inline Real& f (void) {return s;}
};

class SimpleVertexInfo : public Vecteur {
protected:
	Real s;
	unsigned int i;
	Real vol;
public:
	bool isFictious;
	SimpleVertexInfo& operator= (const Vecteur &u) { Vecteur::operator= (u); return *this; }
	SimpleVertexInfo& operator= (const float &scalar) { s=scalar; return *this; }
	SimpleVertexInfo& operator= (const unsigned int &id) { i= id; return *this; }
	inline Real ux (void) {return Vecteur::x();}
	inline Real uy (void) {return Vecteur::y();}
	inline Real uz (void) {return Vecteur::z();}
	inline Real& f (void) {return s;}
	inline Real& v (void) {return vol;}
	inline const unsigned int& id (void) const {return i;}
	SimpleVertexInfo (void) {isFictious=false; s=0; i=0;}
};


class FlowCellInfo : public SimpleCellInfo {

	public:
	//For vector storage of all cells
	unsigned int index;
	int volumeSign;
	bool Pcondition;
	Real invVoidV;
	Real t;
	int fict;
 	Real VolumeVariation;
	double pression;
	 //average relative (fluid - facet translation) velocity defined for a single cell as 1/Volume * SUM_ON_FACETS(x_average_facet*average_facet_flow_rate)
	Vecteur Average_Cell_Velocity;
	// Surface vectors of facets, pointing from outside toward inside the cell
	std::vector<Vecteur> facetSurfaces;
	//Ratio between fluid surface and facet surface 
	std::vector<Real> facetFluidSurfacesRatio;
	// Reflects the geometrical property of the cell, so that the force by cell fluid on grain "i" is pressure*unitForceVectors[i]
	std::vector<Vecteur> unitForceVectors;
	// Store the area of triangle-sphere intersections for each facet (used in forces definition)
	std::vector<Vecteur> facetSphereCrossSections;
	std::vector<Vecteur> cell_force;
	std::vector<double> RayHydr;
// 	std::vector<double> flow_rate;
	std::vector<double> module_permeability;
	// Partial surfaces of spheres in the double-tetrahedron linking two voronoi centers. [i][j] is for sphere facet "i" and sphere facetVertices[i][j]. Last component for 1/sum_surfaces in the facet.
	double solidSurfaces [4][4];

	FlowCellInfo (void)
	{
		module_permeability.resize(4, 0);
		cell_force.resize(4);
		facetSurfaces.resize(4);
		facetFluidSurfacesRatio.resize(4);
		facetSphereCrossSections.resize(4);
		unitForceVectors.resize(4);
		for (int k=0; k<4;k++) for (int l=0; l<3;l++) solidSurfaces[k][l]=0;
		RayHydr.resize(4, 0);
		isInside = false;
		inv_sum_k=0;
		isFictious=false; Pcondition = false; isInferior = false; isSuperior = false; isLateral = false; isvisited = false; isExternal=false;
		index=0;
		volumeSign=0;
		s=0;
		VolumeVariation=0;
		pression=0;
		invVoidV=0;
	}	

	double inv_sum_k;
	bool isInside;
	bool isInferior;
	bool isSuperior;
	bool isLateral;
	bool isvisited;
	bool isExternal;
	
	FlowCellInfo& operator= (const std::vector<double> &v) { for (int i=0; i<4;i++) module_permeability[i]= v[i]; return *this; }
	FlowCellInfo& operator= (const Point &p) { Point::operator= (p); return *this; }
	FlowCellInfo& operator= (const float &scalar) { s=scalar; return *this; }
	
	inline Real& volume (void) {return t;}
	inline Real& invVoidVolume (void) {return invVoidV;}
	inline Real& dv (void) {return VolumeVariation;}
	inline int& fictious (void) {return fict;}
	inline double& p (void) {return pression;}
	inline std::vector<double>& k_norm (void) {return module_permeability;}
	inline std::vector< Vecteur >& facetSurf (void) {return facetSurfaces;}
	inline std::vector<Vecteur>& force (void) {return cell_force;}
	inline std::vector<double>& Rh (void) {return RayHydr;}

	inline Vecteur& av_vel (void) {return Average_Cell_Velocity;}
};

class FlowVertexInfo : public SimpleVertexInfo {
	Vecteur Grain_Velocity;
	Real volume_incident_cells;
public:
	FlowVertexInfo& operator= (const Vecteur &u) { Vecteur::operator= (u); return *this; }
	FlowVertexInfo& operator= (const float &scalar) { s=scalar; return *this; }
	FlowVertexInfo& operator= (const unsigned int &id) { i= id; return *this; }
	Vecteur forces;
	inline Vecteur& force (void) {return forces;}
	inline Vecteur& vel (void) {return Grain_Velocity;}
	inline Real& vol_cells (void) {return volume_incident_cells;}
};

class PeriodicCellInfo : public FlowCellInfo
{	
	public:
	static Vecteur gradP;
	int period[3];
	static Vecteur hSize[3];
	static Vecteur deltaP;
	int ghost;
	Real* _pression;
	bool isGhost;
	PeriodicCellInfo (void)
	{
		_pression=&pression;
		period[0]=period[1]=period[2]=0;
		isGhost=false;
	}
	~PeriodicCellInfo (void) {}
	PeriodicCellInfo& operator= (const Point &p) { Point::operator= (p); return *this; }
	PeriodicCellInfo& operator= (const float &scalar) { s=scalar; return *this; }
	
	inline const double shiftedP (void) const {return isGhost? (*_pression)+pShift() :(*_pression) ;}
	inline const double pShift (void) const {return deltaP[0]*period[0] + deltaP[1]*period[1] +deltaP[2]*period[2];}
// 	inline const double p (void) {return shiftedP();}
	inline void setP (const Real& p) {pression=p;}
};

class PeriodicVertexInfo : public FlowVertexInfo {
	public:
	PeriodicVertexInfo& operator= (const Vecteur &u) { Vecteur::operator= (u); return *this; }
	PeriodicVertexInfo& operator= (const float &scalar) { s=scalar; return *this; }
	PeriodicVertexInfo& operator= (const unsigned int &id) { i= id; return *this; }
	bool isGhost;
	int period[3];
	inline const Vecteur ghostShift (void) {return isGhost? (period[0]*PeriodicCellInfo::hSize[0]+period[1]*PeriodicCellInfo::hSize[1]+period[2]*PeriodicCellInfo::hSize[2]) : CGAL::NULL_VECTOR;}
};


} // namespace CGT

#endif
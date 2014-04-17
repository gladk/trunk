#ifdef FLOW_ENGINE

#pragma once

#include<yade/lib/triangulation/FlowBoundingSphere.hpp>//include after #define XVIEW
// #include "Timer.h"
// #include "PeriodicTesselation.h"
// #include "basicVTKwritter.hpp"
// #include "Timer.h"
// #include "Network.hpp"
// 
// #ifdef XVIEW
// #include "Vue3D.h" //FIXME implicit dependencies will look for this class (out of tree) even ifndef XVIEW
// #endif


namespace CGT{

// 	typedef CGT::FlowBoundingSphere<PeriFlowTesselation> PeriodicFlowBoundingSphere;
	template<class _Tesselation>
	class PeriodicFlow : public CGT::FlowBoundingSphere<_Tesselation>
	{
		public:
		typedef _Tesselation			Tesselation;
		typedef Network<Tesselation>		_N;
		DECLARE_TESSELATION_TYPES(Network<Tesselation>)
		typedef CGT::FlowBoundingSphere<_Tesselation> BaseFlowSolver;
		
		//painfull, but we need that for templates inheritance...
		using _N::T; using _N::xMin; using _N::xMax; using _N::yMin; using _N::yMax; using _N::zMin; using _N::zMax; using _N::Rmoy; using _N::sectionArea; using _N::Height; using _N::vTotal; using _N::currentTes; using _N::debugOut; using _N::nOfSpheres; using _N::xMinId; using _N::xMaxId; using _N::yMinId; using _N::yMaxId; using _N::zMinId; using _N::zMaxId; using _N::boundsIds; using _N::cornerMin; using _N::cornerMax;  using _N::VSolidTot; using _N::Vtotalissimo; using _N::vPoral; using _N::sSolidTot; using _N::vPoralPorosity; using _N::vTotalPorosity; using _N::boundaries; using _N::idOffset; using _N::vtkInfiniteVertices; using _N::vtkInfiniteCells; using _N::num_particles; using _N::boundingCells; using _N::facetVertices; using _N::facetNFictious;
		using BaseFlowSolver::noCache; using BaseFlowSolver::rAverage; using BaseFlowSolver::distanceCorrection; using BaseFlowSolver::minPermLength; using BaseFlowSolver::checkSphereFacetOverlap; using BaseFlowSolver::viscosity; using BaseFlowSolver::kFactor; using BaseFlowSolver::permeabilityMap; using BaseFlowSolver::maxKdivKmean; using BaseFlowSolver::clampKValues; using BaseFlowSolver::KOptFactor; using BaseFlowSolver::meanKStat; using BaseFlowSolver::fluidBulkModulus; using BaseFlowSolver::relax; using BaseFlowSolver::tolerance; using BaseFlowSolver::minKdivKmean;
		
		//same for functions
		using _N::defineFictiousCells; using _N::addBoundingPlanes; using _N::boundary;
		
		void interpolate(Tesselation& Tes, Tesselation& NewTes);
		void computeFacetForcesWithCache(bool onlyCache=false);
		void computePermeability();
		void gaussSeidel(Real dt=0);
		void displayStatistics();
		#ifdef EIGENSPARSE_LIB
		//Eigen's sparse matrix for forces computation
// 		Eigen::SparseMatrix<double> FIntegral;
// 		Eigen::SparseMatrix<double> PshiftsInts;
// 		Eigen::SparseMatrix<double> FRHS;
// 		Eigen::VectorXd forces;
		#endif
	};
	
template<class _Tesselation>	
void PeriodicFlow<_Tesselation>::interpolate(Tesselation& Tes, Tesselation& NewTes)
{
        CellHandle oldCell;
        RTriangulation& Tri = Tes.Triangulation();
	for (VCellIterator cellIt=NewTes.cellHandles.begin(); cellIt!=NewTes.cellHandles.end(); cellIt++){
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
		//FIXME: should use getInfo
                newCell->info().p() = oldCell->info().shiftedP();
        }
//  	Tes.Clear();//Don't reset to avoid segfault when getting pressure in scripts just after interpolation
}
	

template<class _Tesselation>
void PeriodicFlow<_Tesselation>::computeFacetForcesWithCache(bool onlyCache)
{
	RTriangulation& Tri = T[currentTes].Triangulation();
	CVector nullVect(0,0,0);
	static vector<CVector> oldForces;
	if (oldForces.size()<=Tri.number_of_vertices()) oldForces.resize(Tri.number_of_vertices()+1);
	//reset forces
	for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
		if (noCache) {oldForces[v->info().id()]=nullVect; v->info().forces=nullVect;}
		else {oldForces[v->info().id()]=v->info().forces; v->info().forces=nullVect;}
	}
	CellHandle neighbourCell;
	VertexHandle mirrorVertex;
	CVector tempVect;
	
	//FIXME : Ema, be carefull with this (noCache), it needs to be turned true after retriangulation
	if (noCache) {
		for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
		CellHandle& cell = *cellIt;
			for (int k=0;k<4;k++) cell->info().unitForceVectors[k]=nullVect;
			for (int j=0; j<4; j++) if (!Tri.is_infinite(cell->neighbor(j))) {
// 				#ifdef EIGENSPARSE_LIB
// 				if (!cell->info().Pcondition) ++nPCells;
// 				#endif
				neighbourCell = cell->neighbor(j);
				const CVector& Surfk = cell->info().facetSurfaces[j];
				//FIXME : later compute that fluidSurf only once in hydraulicRadius, for now keep full surface not modified in cell->info for comparison with other forces schemes
				//The ratio void surface / facet surface
				Real area = sqrt(Surfk.squared_length());
				CVector facetNormal = Surfk/area;
				const std::vector<CVector>& crossSections = cell->info().facetSphereCrossSections;
				CVector fluidSurfk = cell->info().facetSurfaces[j]*cell->info().facetFluidSurfacesRatio[j];
				/// handle fictious vertex since we can get the projected surface easily here
				if (cell->vertex(j)->info().isFictious) {
					Real projSurf=abs(Surfk[boundary(cell->vertex(j)->info().id()).coordinate]);
					tempVect=-projSurf*boundary(cell->vertex(j)->info().id()).normal;
					//define the cached value for later use with cache*p
					cell->info().unitForceVectors[j]=cell->info().unitForceVectors[j]+ tempVect;
				}
				/// Apply weighted forces f_k=sqRad_k/sumSqRad*f
				CVector facetUnitForce = -fluidSurfk*cell->info().solidSurfaces[j][3];
				for (int y=0; y<3;y++) {
					//add to cached value
					cell->info().unitForceVectors[facetVertices[j][y]]=cell->info().unitForceVectors[facetVertices[j][y]]+facetUnitForce*cell->info().solidSurfaces[j][y];
					//uncomment to get total force / comment to get only viscous forces (Bruno)
					if (!cell->vertex(facetVertices[j][y])->info().isFictious) {
						//add to cached value
						cell->info().unitForceVectors[facetVertices[j][y]]=cell->info().unitForceVectors[facetVertices[j][y]]-facetNormal*crossSections[j][y];
					}
				}
			}
	}

	noCache=false;//cache should always be defined after execution of this function
	if (onlyCache) return;
	}// end if(noCache)
	
	//use cached values
	//First define products that will be used below for all cells:
	Real pDeltas [3];
	for (unsigned int k=0; k<3;k++) pDeltas[k]=CellInfo::hSize[k]*CellInfo::gradP;
	//Then compute the forces
	for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
		const CellHandle& cell = *cellIt;
		for (int yy=0;yy<4;yy++) {
			VertexInfo& vhi = cell->vertex(yy)->info();
			Real unshiftedP = cell->info().p();
			//the pressure translated to a ghost cell adjacent to the non-ghost vertex
			unshiftedP -= pDeltas[0]*vhi.period[0] + pDeltas[1]*vhi.period[1] +pDeltas[2]*vhi.period[2];
			T[currentTes].vertexHandles[vhi.id()]->info().forces=T[currentTes].vertexHandles[vhi.id()]->info().forces + cell->info().unitForceVectors[yy]*unshiftedP;
		}
	}
	if (debugOut) {
		CVector totalForce = nullVect;
		for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); v++)
		{
			if (!v->info().isFictious /*&& !v->info().isGhost*/ ){
				totalForce = totalForce + v->info().forces;}
			else /*if (!v->info().isGhost)*/{
				if (boundary(v->info().id()).flowCondition==1) totalForce = totalForce + v->info().forces;
			}
		}
		cout << "totalForce = "<< totalForce << endl;}
}

template<class _Tesselation>
void PeriodicFlow<_Tesselation>::computePermeability()
{
	if (debugOut)  cout << "----Computing_Permeability (Periodic)------" << endl;
	RTriangulation& Tri = T[currentTes].Triangulation();
	VSolidTot = 0, Vtotalissimo = 0, vPoral = 0, sSolidTot = 0, vTotalPorosity=0, vPoralPorosity=0;
	CellHandle neighbourCell;

	double k=0, distance = 0, radius = 0;
	int surfneg=0; int NEG=0, POS=0, pass=0;
	Real meanK=0, STDEV=0, meanRadius=0, meanDistance=0;
	Real infiniteK=1e3;
	double volume_sub_pore = 0.f;
	VectorCell& cellHandles= T[currentTes].cellHandles;
	for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
			CellHandle& cell = *cellIt;
			Point& p1 = cell->info();
			for (int j=0; j<4; j++){
				neighbourCell = cell->neighbor(j);
				Point& p2 = neighbourCell->info();
				if (!Tri.is_infinite(neighbourCell) /*&& (neighbour_cell->info().isvisited==ref || computeAllCells)*/) {
					//compute and store the area of sphere-facet intersections for later use
					VertexHandle W [3];
					for (int kk=0; kk<3; kk++) {
						W[kk] = cell->vertex(facetVertices[j][kk]);
					}
					Sphere& v0 = W[0]->point();
					Sphere& v1 = W[1]->point();
					Sphere& v2 = W[2]->point();
#ifdef USE_FAST_MATH
					//FIXME : code not compiling,, do the same as in "else"
					assert((W[permut3[jj][1]]->point()-W[permut3[jj][0]]->point())*(W[permut3[jj][2]]->point()-W[permut3[jj][0]]->point())>=0 && (W[permut3[jj][1]]->point()-W[permut3[jj][0]]->point())*(W[permut3[jj][2]]->point()-W[permut3[jj][0]]->point())<=1);
					for (int jj=0;jj<3;jj++)
						cell->info().facetSphereCrossSections[j][jj]=0.5*W[jj]->point().weight()*Wm3::FastInvCos1((W[permut3[jj][1]]->point()-W[permut3[jj][0]]->point())*(W[permut3[jj][2]]->point()-W[permut3[jj][0]]->point()));
#else
					cell->info().facetSphereCrossSections[j]=CVector(
					W[0]->info().isFictious ? 0 : 0.5*v0.weight()*acos((v1-v0)*(v2-v0)/sqrt((v1-v0).squared_length()*(v2-v0).squared_length())),
					W[1]->info().isFictious ? 0 : 0.5*v1.weight()*acos((v0-v1)*(v2-v1)/sqrt((v1-v0).squared_length()*(v2-v1).squared_length())),
					W[2]->info().isFictious ? 0 : 0.5*v2.weight()*acos((v0-v2)*(v1-v2)/sqrt((v1-v2).squared_length()*(v2-v0).squared_length())));
#endif
					pass+=1;
					CVector l = p1 - p2;
					distance = sqrt(l.squared_length());
					if (!rAverage) radius = 2* this->computeHydraulicRadius(cell, j);
					else radius = (this->computeEffectiveRadius(cell, j)+this->computeEquivalentRadius(cell,j))*0.5;
					if (radius<0) NEG++;
					else POS++;
					if (radius==0) {
						cout << "INS-INS PROBLEM!!!!!!!" << endl;
					}
					Real fluidArea=0;
					int test=0;
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
						meanK += k*kFactor;

					if (k<0 && debugOut) {surfneg+=1;
					cout<<"__ k<0 __"<<k<<" "<<" fluidArea "<<fluidArea<<" area "<<area<<" "<<crossSections[0]<<" "<<crossSections[1]<<" "<<crossSections[2] <<" "<<W[0]->info().id()<<" "<<W[1]->info().id()<<" "<<W[2]->info().id()<<" "<<p1<<" "<<p2<<" test "<<test<<endl;}
					     
					} else  cout <<"infinite K2! surfaces will be missing (FIXME)"<<endl; k = infiniteK;
					//Will be corrected in the next loop
					(cell->info().kNorm())[j]= k*kFactor;
					if (!neighbourCell->info().isGhost) (neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]= (cell->info().kNorm())[j];
					//The following block is correct but not very usefull, since all values are clamped below with MIN and MAX, skip for now
// 					else {//find the real neighbor connected to our cell through periodicity
// 						CellHandle true_neighbour_cell = cellHandles[neighbour_cell->info().baseIndex];
// 						for (int ii=0;ii<4;ii++)
// 							if (true_neighbour_cell->neighbor(ii)->info().index == cell->info().index){
// 								(true_neighbour_cell->info().kNorm())[ii]=(cell->info().kNorm())[j]; break;
// 							}
// 					}
					if(permeabilityMap){
						CellHandle c = cell;
						cell->info().s = cell->info().s + k*distance/fluidArea*this->volumePoreVoronoiFraction (c,j);
						volume_sub_pore += this->volumePoreVoronoiFraction (c,j);
					}
				}
			}
// 			cell->info().isvisited = !ref;
			if(permeabilityMap){
				cell->info().s = cell->info().s/volume_sub_pore;
				volume_sub_pore = 0.f;
			}
// 		}
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
	pass=0;
	if (clampKValues) for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
		CellHandle& cell = *cellIt;
		for (int j=0; j<4; j++) {
			neighbourCell = cell->neighbor(j);
			if (!Tri.is_infinite(neighbourCell) /*&& neighbour_cell->info().isvisited==ref*/) {
				pass++;
				(cell->info().kNorm())[j] = max(minKdivKmean*globalK, min((cell->info().kNorm())[j], maxKdivKmean*globalK));
				(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]=(cell->info().kNorm())[j];
				if (!neighbourCell->info().isGhost) (neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]= (cell->info().kNorm())[j];
					else {//find the real neighbor connected to our cell through periodicity, as we want exactly the same permeability without rounding errors
						CellHandle& true_neighbourCell = cellHandles[neighbourCell->info().baseIndex];
						for (int ii=0;ii<4;ii++)
							if (true_neighbourCell->neighbor(ii)->info().index == cell->info().index){
								(true_neighbourCell->info().kNorm())[ii]=(cell->info().kNorm())[j]; break;
							}
					}

			}
		}
	}
	if (debugOut) cout << "PassKcorrect = " << pass << endl;

	if (debugOut) cout << "POS = " << POS << " NEG = " << NEG << " pass = " << pass << endl;

// A loop to compute the standard deviation of the local K distribution, and use it to include/exclude K values higher then (meanK +/- K_opt_factor*STDEV)
	if (meanKStat)
	{
		std::ofstream k_opt_file("k_stdev.txt" ,std::ios::out);
		pass=0;
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			for (int j=0; j<4; j++) {
				neighbourCell = cell->neighbor(j);
				if (!Tri.is_infinite(neighbourCell) /*&& neighbour_cell->info().isvisited==ref*/) {
					pass++;
					STDEV += pow(((cell->info().kNorm())[j]-meanK),2);
				}
			}
		}
		STDEV = sqrt(STDEV/pass);
		if (debugOut) cout << "PassSTDEV = " << pass << endl;
		cout << "STATISTIC K" << endl;
		double k_min = 0, k_max = meanK + KOptFactor*STDEV;
		cout << "Kmoy = " << meanK << " Standard Deviation = " << STDEV << endl;
		cout << "kmin = " << k_min << " kmax = " << k_max << endl;
		pass=0;
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			for (int j=0; j<4; j++) {
				neighbourCell = cell->neighbor(j);
				if (!Tri.is_infinite(neighbourCell) /*&& neighbour_cell->info().isvisited==ref*/) {
					pass+=1;
					if ((cell->info().kNorm())[j]>k_max) {
						(cell->info().kNorm())[j]=k_max;
						(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]= (cell->info().kNorm())[j];
					}
					k_opt_file << KOptFactor << " " << (cell->info().kNorm())[j] << endl;
				}
			}
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
				Vgrains += 1.33333333 * M_PI * pow(vIt->point().weight(),1.5);
			}
		}
		cout<<grains<<"grains - " <<"vTotal = " << vTotal << " Vgrains = " << Vgrains << " vPoral1 = " << (vTotal-Vgrains) << endl;
		cout << "Vtotalissimo = " << Vtotalissimo/2 << " VSolidTot = " << VSolidTot/2 << " vPoral2 = " << vPoral/2  << " sSolidTot = " << sSolidTot << endl<< endl;

		if (!rAverage) cout << "------Hydraulic Radius is used for permeability computation------" << endl << endl;
		else cout << "------Average Radius is used for permeability computation------" << endl << endl;
		cout << "-----computed_Permeability Periodic-----" << endl;
	}
}




template<class _Tesselation>
void PeriodicFlow<_Tesselation>::gaussSeidel(Real dt)
{
    RTriangulation& Tri = T[currentTes].Triangulation();
    int j = 0;
    double m, n, dp_max, p_max, sum_p, dp, sum_dp;
    double compFlowFactor=0;
    vector<Real> previousP;
    previousP.resize(Tri.number_of_finite_cells());
    const int num_threads=1;
    bool compressible= fluidBulkModulus>0;

    vector<Real> t_sum_p, t_dp_max, t_sum_dp, t_p_max;
    t_sum_dp.resize(num_threads);
    t_dp_max.resize(num_threads);
    t_p_max.resize(num_threads);
    t_sum_p.resize(num_threads);
    FiniteCellsIterator cellEnd = Tri.finite_cells_end();
    do {
        int cell2=0;
        dp_max = 0;
        p_max = 0;
        sum_p=0;
        sum_dp=0;
        int bb=-1;
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
            bb++;
            if ( !cell->info().Pcondition && !cell->info().isGhost) {
		cell2++;
		if (compressible && j==0) previousP[bb]=cell->info().shiftedP();
		m=0, n=0;
		for (int j2=0; j2<4; j2++) {
		  if (!Tri.is_infinite(cell->neighbor(j2))) {
			if ( compressible ) {
				compFlowFactor = fluidBulkModulus*dt*cell->info().invVoidVolume();
				m += compFlowFactor*(cell->info().kNorm())[j2]*cell->neighbor(j2)->info().shiftedP();
				if (j==0) n += compFlowFactor*(cell->info().kNorm())[j2];
			} else {
				m += (cell->info().kNorm())[j2]*cell->neighbor(j2)->info().shiftedP();
				if ( isinf(m) && j<10 ) cout << "(cell->info().kNorm())[j2] = " << (cell->info().kNorm())[j2] << " cell->neighbor(j2)->info().shiftedP() = " << cell->neighbor(j2)->info().shiftedP() << endl;
				if (j==0) n += (cell->info().kNorm())[j2];
			} 
		  }
		}
		dp = cell->info().p();
		if (n!=0 || j!=0) {
			if (j==0) { if (compressible) cell->info().invSumK=1/(1+n); else cell->info().invSumK=1/n; }
			if ( compressible ) { 
				cell->info().setP( ( ((previousP[bb] - ((fluidBulkModulus*dt*cell->info().invVoidVolume())*(cell->info().dv()))) + m) * cell->info().invSumK - cell->info().shiftedP()) * relax + cell->info().shiftedP());
			} else {
				cell->info().setP((-(cell->info().dv()-m)*cell->info().invSumK-cell->info().p())*relax+cell->info().shiftedP());
			}
		}
		dp -= cell->info().p();
		dp_max = max(dp_max, std::abs(dp));
		p_max = max(p_max, std::abs(cell->info().shiftedP()));
		sum_p += cell->info().shiftedP();
		sum_dp += std::abs(dp);
            }
        }
        j++;

	if (j>=40000) cerr<<"\r GS not converged after 40k iterations, break";

    } while ((dp_max/p_max) > tolerance && j<40000 /*&& ( dp_max > tolerance )*//* &&*/ /*( j<50 )*/);

    int cel=0;
    double Pav=0;
    for (FiniteCellsIterator cell = Tri.finite_cells_begin();
            cell != cellEnd;
            cell++) {
        cel++;
        Pav+=cell->info().shiftedP();
    }
    Pav/=cel;
}

template<class _Tesselation>
void PeriodicFlow<_Tesselation>::displayStatistics()
{
	RTriangulation& Tri = T[currentTes].Triangulation();
	int Zero =0, Inside=0, Fictious=0, ghostC=0,realC=0, ghostV=0, realV=0;
	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
		int zeros =0;
		for (int j=0; j!=4; j++) 
			  if ((cell->info().kNorm())[j]==0) zeros+=1;
		if (zeros==4) Zero+=1;
		if (!cell->info().fictious()) Inside+=1; else Fictious+=1;
		if (cell->info().isGhost)  ghostC+=1; else realC+=1;
	}
	int fict=0, real=0;
	for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
		if (v->info().isFictious) fict+=1;
		else real+=1;
	}
	long Vertices = Tri.number_of_vertices();
	long Cells = Tri.number_of_finite_cells();
	long Facets = Tri.number_of_finite_facets();
	if(debugOut) {
		cout << "zeros = " << Zero << endl;
		cout << "There are " << Vertices << " vertices, dont " << fict << " fictious et " << real << " reeeeeel" << std::endl;
		cout << "There are " << ghostV+realV << " vertices, dont " << ghostV << " ghost et " << realV << " reeeeeel" << std::endl;
		cout << "There are " << ghostC+realC << " cells, dont " << ghostC << " ghost et " << realC << " reeeeeel" << std::endl;
		cout << "There are " << Cells << " cells " << std::endl;
		cout << "There are " << Facets << " facets " << std::endl;
		cout << "There are " << Inside << " cells INSIDE." << endl;
		cout << "There are " << Fictious << " cells FICTIOUS." << endl;
	}
	vtkInfiniteVertices = fict;
	vtkInfiniteCells = Fictious;
	num_particles = real;
}
	
} //END NAMESPACE

#ifdef LINSOLV
#include "PeriodicFlowLinSolv.hpp"
#endif


#endif //FLOW_ENGINE

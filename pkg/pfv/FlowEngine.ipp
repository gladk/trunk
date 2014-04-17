/*************************************************************************
*  Copyright (C) 2009 by Emanuele Catalano <catalano@grenoble-inp.fr>    *
*  Copyright (C) 2009 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef YADE_CGAL

#ifdef FLOW_ENGINE
#include<yade/core/Scene.hpp>
#include<yade/lib/base/Math.hpp>
#include<yade/pkg/dem/TesselationWrapper.hpp>
#include<yade/pkg/common/Sphere.hpp>
#include<yade/pkg/common/Wall.hpp>
#include<yade/pkg/common/Box.hpp>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#ifdef LINSOLV
#include <cholmod.h>
#endif

#include "FlowEngine.hpp"

template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::~TemplateFlowEngine() {}

// YADE_PLUGIN((TFlowEng));

template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
unsigned int TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::imposePressure(Vector3r pos, Real p)
{
// 	if (!flow) LOG_ERROR("no flow defined yet, run at least one iter");
	solver->imposedP.push_back( pair<CGT::Point,Real>(CGT::Point(pos[0],pos[1],pos[2]),p) );
	//force immediate update of boundary conditions
	updateTriangulation=true;
	return solver->imposedP.size()-1;
}

template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::action()
{
       if ( !isActivated ) return;
        timingDeltas->start();
	setPositionsBuffer(true);
	timingDeltas->checkpoint ( "Position buffer" );
        if (first) {
	  if (multithread) setPositionsBuffer(false);
	  buildTriangulation(pZero,*solver);
	  initializeVolumes(*solver);
	  backgroundSolver=solver;
	  backgroundCompleted=true;
	}
	solver->ompThreads = ompThreads>0? ompThreads : omp_get_max_threads();

        timingDeltas->checkpoint ( "Triangulating" );
	updateVolumes ( *solver );
        timingDeltas->checkpoint ( "Update_Volumes" );
	
        epsVolCumulative += epsVolMax;
	retriangulationLastIter++;
	if (!updateTriangulation) updateTriangulation = // If not already set true by another function of by the user, check conditions
		(defTolerance>0 && epsVolCumulative > defTolerance) || retriangulationLastIter>meshUpdateInterval;

        ///compute flow and and forces here
	if (pressureForce){
		solver->gaussSeidel(scene->dt);
		timingDeltas->checkpoint ( "Gauss-Seidel (includes matrix construct and factorization in single-thread mode)" );
		solver->computeFacetForcesWithCache();}
        timingDeltas->checkpoint ( "compute_Forces" );
        ///Application of vicscous forces
        scene->forces.sync();
	timingDeltas->checkpoint ( "forces.sync()" );
	computeViscousForces ( *solver );
	timingDeltas->checkpoint ( "viscous forces" );
	Vector3r force;
	Vector3r torque;
        FiniteVerticesIterator verticesEnd = solver->T[solver->currentTes].Triangulation().finite_vertices_end();
        for ( FiniteVerticesIterator vIt = solver->T[solver->currentTes].Triangulation().finite_vertices_begin(); vIt !=  verticesEnd; vIt++ ) {
		force = pressureForce ? Vector3r ( vIt->info().forces[0],vIt->info().forces[1],vIt->info().forces[2] ): Vector3r(0,0,0);
		torque = Vector3r(0,0,0);
                if (shearLubrication || viscousShear){
			force = force + solver->shearLubricationForces[vIt->info().id()];
			torque = torque + solver->shearLubricationTorques[vIt->info().id()];
			if (pumpTorque)
				torque = torque + solver->pumpLubricationTorques[vIt->info().id()];
		}
		if (twistTorque)
			torque = torque + solver->twistLubricationTorques[vIt->info().id()];
		if (normalLubrication)
			force = force + solver-> normalLubricationForce[vIt->info().id()];
		scene->forces.addForce ( vIt->info().id(), force);
		scene->forces.addTorque ( vIt->info().id(), torque);
        }
        ///End compute flow and forces
        timingDeltas->checkpoint ( "Applying Forces" );
	int sleeping = 0;
	if (multithread && !first) {
		while (updateTriangulation && !backgroundCompleted) { /*cout<<"sleeping..."<<sleeping++<<endl;*/
		  sleeping++;
		boost::this_thread::sleep(boost::posix_time::microseconds(1000));}
		if (debug && sleeping) cerr<<"sleeping..."<<sleeping<<endl;
		if (updateTriangulation || (ellapsedIter>(0.5*meshUpdateInterval) && backgroundCompleted)) {
			if (debug) cerr<<"switch flow solver"<<endl;
			if (useSolver==0) LOG_ERROR("background calculations not available for Gauss-Seidel");
			if (fluidBulkModulus>0 || doInterpolate) solver->interpolate (solver->T[solver->currentTes], backgroundSolver->T[backgroundSolver->currentTes]);
			solver=backgroundSolver;
			backgroundSolver = shared_ptr<FlowSolver> (new FlowSolver);
			if (metisForced) {backgroundSolver->eSolver.cholmod().nmethods=1; backgroundSolver->eSolver.cholmod().method[0].ordering=CHOLMOD_METIS;}
			//Copy imposed pressures/flow from the old solver
			backgroundSolver->imposedP = vector<pair<CGT::Point,Real> >(solver->imposedP);
			backgroundSolver->imposedF = vector<pair<CGT::Point,Real> >(solver->imposedF);
			if (debug) cerr<<"switched"<<endl;
			setPositionsBuffer(false);//set "parallel" buffer for background calculation 
			backgroundCompleted=false;
			retriangulationLastIter=ellapsedIter;
			updateTriangulation=false;
			epsVolCumulative=0;
			ellapsedIter=0;
			boost::thread workerThread(&TemplateFlowEngine::backgroundAction,this);
			workerThread.detach();
			if (debug) cerr<<"backgrounded"<<endl;
			initializeVolumes(*solver);
			computeViscousForces(*solver);
			if (debug) cerr<<"volumes initialized"<<endl;
		}
		else {
			if (debug && !backgroundCompleted) cerr<<"still computing solver in the background, ellapsedIter="<<ellapsedIter<<endl;
			ellapsedIter++;
		}
	} else {
	        if (updateTriangulation && !first) {
			buildTriangulation (pZero, *solver);
			initializeVolumes(*solver);
			computeViscousForces(*solver);
               		updateTriangulation = false;
			epsVolCumulative=0;
			retriangulationLastIter=0;
			ReTrg++;}
        }
        first=false;
        timingDeltas->checkpoint ( "triangulate + init volumes" );
}

template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::backgroundAction()
{
	if (useSolver<1) {LOG_ERROR("background calculations not available for Gauss-Seidel"); return;}
        buildTriangulation ( pZero,*backgroundSolver );
	//FIXME: GS is computing too much, we need only matrix factorization in fact
	backgroundSolver->gaussSeidel(scene->dt);
	//FIXME(2): and here we need only cached variables, not forces
	backgroundSolver->computeFacetForcesWithCache(/*onlyCache?*/ true);
// 	boost::this_thread::sleep(boost::posix_time::seconds(5));
 	backgroundCompleted = true;
}

template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::boundaryConditions ( Solver& flow )
{
	for (int k=0;k<6;k++)	{
		flow.boundary (wallIds[k]).flowCondition=!bndCondIsPressure[k];
                flow.boundary (wallIds[k]).value=bndCondValue[k];
                flow.boundary (wallIds[k]).velocity = boundaryVelocity[k];//FIXME: needs correct implementation, maybe update the cached pos/vel?
	}
}

template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::setImposedPressure ( unsigned int cond, Real p)
{
        if ( cond>=solver->imposedP.size() ) LOG_ERROR ( "Setting p with cond higher than imposedP size." );
        solver->imposedP[cond].second=p;
        //force immediate update of boundary conditions
	solver->pressureChanged=true;
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::imposeFlux ( Vector3r pos, Real flux){
        solver->imposedF.push_back ( pair<CGT::Point,Real> ( CGT::Point ( pos[0],pos[1],pos[2] ), flux ) );
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::clearImposedPressure () { solver->imposedP.clear(); solver->IPCells.clear();}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::clearImposedFlux () { solver->imposedF.clear(); solver->IFCells.clear();}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
Real TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::getCellFlux ( unsigned int cond)
{
	if ( cond>=solver->imposedP.size() ) {LOG_ERROR ( "Getting flux with cond higher than imposedP size." ); return 0;}
        double flux=0;
        typename Solver::CellHandle& cell= solver->IPCells[cond];
        for ( int ngb=0;ngb<4;ngb++ ) {
                /*if (!cell->neighbor(ngb)->info().Pcondition)*/
                flux+= cell->info().kNorm() [ngb]* ( cell->info().p()-cell->neighbor ( ngb )->info().p() );
        }
        return flux+cell->info().dv();
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::initSolver ( FlowSolver& flow )
{
       	flow.Vtotalissimo=0; flow.VSolidTot=0; flow.vPoral=0; flow.sSolidTot=0;
        flow.slipBoundary=slipBoundary;
        flow.kFactor = permeabilityFactor;
        flow.debugOut = debug;
        flow.useSolver = useSolver;
	#ifdef EIGENSPARSE_LIB
	flow.numSolveThreads = numSolveThreads;
	flow.numFactorizeThreads = numFactorizeThreads;
	#endif
	flow.meanKStat = meanKStat;
        flow.viscosity = viscosity;
        flow.tolerance=tolerance;
        flow.relax=relax;
        flow.clampKValues = clampKValues;
	flow.maxKdivKmean = maxKdivKmean;
	flow.minKdivKmean = minKdivKmean;
        flow.meanKStat = meanKStat;
        flow.permeabilityMap = permeabilityMap;
        flow.fluidBulkModulus = fluidBulkModulus;
//         flow.T[flow.currentTes].Clear();
        flow.T[flow.currentTes].maxId=-1;
        flow.xMin = 1000.0, flow.xMax = -10000.0, flow.yMin = 1000.0, flow.yMax = -10000.0, flow.zMin = 1000.0, flow.zMax = -10000.0;
}

#ifdef LINSOLV
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::setForceMetis ( bool force )
{
        if (force) {
		metisForced=true;
		solver->eSolver.cholmod().nmethods=1;
		solver->eSolver.cholmod().method[0].ordering=CHOLMOD_METIS;
	} else {cholmod_defaults(&(solver->eSolver.cholmod())); metisForced=false;}
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
bool TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::getForceMetis () {return (solver->eSolver.cholmod().nmethods==1);}
#endif
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::buildTriangulation ( Solver& flow )
{
        buildTriangulation ( 0.f,flow );
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::buildTriangulation ( double pZero, Solver& flow )
{
 	if (first) flow.currentTes=0;
        else {
                flow.currentTes=!flow.currentTes;
                if (debug) cout << "--------RETRIANGULATION-----------" << endl;
        }
	flow.resetNetwork();
	initSolver(flow);

        addBoundary ( flow );
        triangulate ( flow );
        if ( debug ) cout << endl << "Tesselating------" << endl << endl;
        flow.T[flow.currentTes].compute();

        flow.defineFictiousCells();
	// For faster loops on cells define this vector
	flow.T[flow.currentTes].cellHandles.clear();
	flow.T[flow.currentTes].cellHandles.reserve(flow.T[flow.currentTes].Triangulation().number_of_finite_cells());
	FiniteCellsIterator cell_end = flow.T[flow.currentTes].Triangulation().finite_cells_end();
	int k=0;
	for ( FiniteCellsIterator cell = flow.T[flow.currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ){
		flow.T[flow.currentTes].cellHandles.push_back(cell);
		cell->info().id=k++;}//define unique numbering now, corresponds to position in cellHandles
        flow.displayStatistics ();
        flow.computePermeability();
	//This virtual function does nothing yet, derived class may overload it to make permeability different (see DFN engine)
	trickPermeability();
        porosity = flow.vPoralPorosity/flow.vTotalPorosity;

        boundaryConditions ( flow );
        flow.initializePressure ( pZero );
	
        if ( !first && !multithread && (useSolver==0 || fluidBulkModulus>0 || doInterpolate)) flow.interpolate ( flow.T[!flow.currentTes], flow.T[flow.currentTes] );
        if ( waveAction ) flow.applySinusoidalPressure ( flow.T[flow.currentTes].Triangulation(), sineMagnitude, sineAverage, 30 );
        if (normalLubrication || shearLubrication || viscousShear) flow.computeEdgesSurfaces();
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::setPositionsBuffer(bool current)
{
	vector<posData>& buffer = current? positionBufferCurrent : positionBufferParallel;
	buffer.clear();
	buffer.resize(scene->bodies->size());
	shared_ptr<Sphere> sph ( new Sphere );
        const int Sph_Index = sph->getClassIndexStatic();
	FOREACH ( const shared_ptr<Body>& b, *scene->bodies ) {
                if (!b || ignoredBody==b->getId()) continue;
                posData& dat = buffer[b->getId()];
		dat.id=b->getId();
		dat.pos=b->state->pos;
		dat.isSphere= (b->shape->getClassIndex() ==  Sph_Index);
		if (dat.isSphere) dat.radius = YADE_CAST<Sphere*>(b->shape.get())->radius;
		dat.exists=true;
	}
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::addBoundary ( Solver& flow )
{
	vector<posData>& buffer = multithread ? positionBufferParallel : positionBufferCurrent;
        solver->xMin = Mathr::MAX_REAL, solver->xMax = -Mathr::MAX_REAL, solver->yMin = Mathr::MAX_REAL, solver->yMax = -Mathr::MAX_REAL, solver->zMin = Mathr::MAX_REAL, solver->zMax = -Mathr::MAX_REAL;
        FOREACH ( const posData& b, buffer ) {
                if ( !b.exists ) continue;
                if ( b.isSphere ) {
                        const Real& rad = b.radius;
                        const Real& x = b.pos[0];
                        const Real& y = b.pos[1];
                        const Real& z = b.pos[2];
                        flow.xMin = min ( flow.xMin, x-rad );
                        flow.xMax = max ( flow.xMax, x+rad );
                        flow.yMin = min ( flow.yMin, y-rad );
                        flow.yMax = max ( flow.yMax, y+rad );
                        flow.zMin = min ( flow.zMin, z-rad );
                        flow.zMax = max ( flow.zMax, z+rad );
                }
        }
	//FIXME idOffset must be set correctly, not the case here (always 0), then we need walls first or it will fail
        idOffset = flow.T[flow.currentTes].maxId+1;
        flow.idOffset = idOffset;
        flow.sectionArea = ( flow.xMax - flow.xMin ) * ( flow.zMax-flow.zMin );
        flow.vTotal = ( flow.xMax-flow.xMin ) * ( flow.yMax-flow.yMin ) * ( flow.zMax-flow.zMin );
        flow.yMinId=wallIds[ymin];
        flow.yMaxId=wallIds[ymax];
        flow.xMaxId=wallIds[xmax];
        flow.xMinId=wallIds[xmin];
        flow.zMinId=wallIds[zmin];
        flow.zMaxId=wallIds[zmax];

        //FIXME: Id's order in boundsIds is done according to the enumeration of boundaries from TXStressController.hpp, line 31. DON'T CHANGE IT!
        flow.boundsIds[0]= &flow.xMinId;
        flow.boundsIds[1]= &flow.xMaxId;
        flow.boundsIds[2]= &flow.yMinId;
        flow.boundsIds[3]= &flow.yMaxId;
        flow.boundsIds[4]= &flow.zMinId;
        flow.boundsIds[5]= &flow.zMaxId;

	for (int k=0;k<6;k++) flow.boundary ( *flow.boundsIds[k] ).useMaxMin = boundaryUseMaxMin[k];

        flow.cornerMin = CGT::Point ( flow.xMin, flow.yMin, flow.zMin );
        flow.cornerMax = CGT::Point ( flow.xMax, flow.yMax, flow.zMax );
 
        //assign BCs types and values
        boundaryConditions ( flow );

        double center[3];
        for ( int i=0; i<6; i++ ) {
                if ( *flow.boundsIds[i]<0 ) continue;
                CGT::CVector Normal ( normal[i].x(), normal[i].y(), normal[i].z() );
                if ( flow.boundary ( *flow.boundsIds[i] ).useMaxMin ) flow.addBoundingPlane(Normal, *flow.boundsIds[i] );
                else {
			for ( int h=0;h<3;h++ ) center[h] = buffer[*flow.boundsIds[i]].pos[h];
// 			cerr << "id="<<*flow.boundsIds[i] <<" center="<<center[0]<<","<<center[1]<<","<<center[2]<<endl;
                        flow.addBoundingPlane ( center, wallThickness, Normal,*flow.boundsIds[i] );
                }
        }
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::triangulate ( Solver& flow )
{
///Using Tesselation wrapper (faster)
// 	TesselationWrapper TW;
// 	if (TW.Tes) delete TW.Tes;
// 	TW.Tes = &(flow.T[flow.currentTes]);//point to the current Tes we have in Flowengine
// 	TW.insertSceneSpheres();//TW is now really inserting in TemplateFlowEngine, using the faster insert(begin,end)
// 	TW.Tes = NULL;//otherwise, Tes would be deleted by ~TesselationWrapper() at the end of the function.
///Using one-by-one insertion
	vector<posData>& buffer = multithread ? positionBufferParallel : positionBufferCurrent;
	FOREACH ( const posData& b, buffer ) {
                if ( !b.exists ) continue;
                if ( b.isSphere ) {
			if (b.id==ignoredBody) continue;
                        flow.T[flow.currentTes].insert ( b.pos[0], b.pos[1], b.pos[2], b.radius, b.id );}
        }
	flow.T[flow.currentTes].redirected=true;//By inserting one-by-one, we already redirected
	flow.shearLubricationForces.resize ( flow.T[flow.currentTes].maxId+1 );
	flow.shearLubricationTorques.resize ( flow.T[flow.currentTes].maxId+1 );
	flow.pumpLubricationTorques.resize ( flow.T[flow.currentTes].maxId+1 );
	flow.twistLubricationTorques.resize ( flow.T[flow.currentTes].maxId+1 );
	flow.shearLubricationBodyStress.resize ( flow.T[flow.currentTes].maxId+1 );
	flow.normalLubricationForce.resize ( flow.T[flow.currentTes].maxId+1 );
	flow.normalLubricationBodyStress.resize ( flow.T[flow.currentTes].maxId+1 );
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::initializeVolumes ( Solver& flow )
{
	typedef typename Solver::FiniteVerticesIterator FiniteVerticesIterator;
	
	FiniteVerticesIterator vertices_end = flow.T[flow.currentTes].Triangulation().finite_vertices_end();
	CGT::CVector Zero(0,0,0);
	for (FiniteVerticesIterator V_it = flow.T[flow.currentTes].Triangulation().finite_vertices_begin(); V_it!= vertices_end; V_it++) V_it->info().forces=Zero;

	FOREACH(CellHandle& cell, flow.T[flow.currentTes].cellHandles)
	{
		switch ( cell->info().fictious() )
		{
			case ( 0 ) : cell->info().volume() = volumeCell ( cell ); break;
			case ( 1 ) : cell->info().volume() = volumeCellSingleFictious ( cell ); break;
			case ( 2 ) : cell->info().volume() = volumeCellDoubleFictious ( cell ); break;
			case ( 3 ) : cell->info().volume() = volumeCellTripleFictious ( cell ); break;
			default: break; 
		}
		if (flow.fluidBulkModulus>0) { cell->info().invVoidVolume() = 1 / ( abs(cell->info().volume()) - flow.volumeSolidPore(cell) ); }
	}
	if (debug) cout << "Volumes initialised." << endl;
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::averageRealCellVelocity()
{
        solver->averageRelativeCellVelocity();
        Vector3r Vel ( 0,0,0 );
        //AVERAGE CELL VELOCITY
        FiniteCellsIterator cell_end = solver->T[solver->currentTes].Triangulation().finite_cells_end();
        for ( FiniteCellsIterator cell = solver->T[solver->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {
                for ( int g=0;g<4;g++ ) {
                        if ( !cell->vertex ( g )->info().isFictious ) {
                                const shared_ptr<Body>& sph = Body::byId ( cell->vertex ( g )->info().id(), scene );
                                for ( int i=0;i<3;i++ ) Vel[i] = Vel[i] + sph->state->vel[i]/4;
                        }
                }
                RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
                CGT::Point pos_av_facet;
                double volume_facet_translation = 0;
                CGT::CVector Vel_av ( Vel[0], Vel[1], Vel[2] );
                for ( int i=0; i<4; i++ ) {
                        volume_facet_translation = 0;
                        if ( !Tri.is_infinite ( cell->neighbor ( i ) ) ) {
                                CGT::CVector Surfk = cell->info()-cell->neighbor ( i )->info();
                                Real area = sqrt ( Surfk.squared_length() );
                                Surfk = Surfk/area;
                                CGT::CVector branch = cell->vertex ( facetVertices[i][0] )->point() - cell->info();
                                pos_av_facet = ( CGT::Point ) cell->info() + ( branch*Surfk ) *Surfk;
                                volume_facet_translation += Vel_av*cell->info().facetSurfaces[i];
                                cell->info().averageVelocity() = cell->info().averageVelocity() - volume_facet_translation/cell->info().volume() * ( pos_av_facet-CGAL::ORIGIN );
                        }
                }
        }
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::updateVolumes ( Solver& flow )
{
        if ( debug ) cout << "Updating volumes.............." << endl;
        Real invDeltaT = 1/scene->dt;
        epsVolMax=0;
        Real totVol=0; Real totDVol=0;
	#ifdef YADE_OPENMP
	const long size=flow.T[flow.currentTes].cellHandles.size();
	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for(long i=0; i<size; i++){
		CellHandle& cell = flow.T[flow.currentTes].cellHandles[i];
	#else
	FOREACH(CellHandle& cell, flow.T[flow.currentTes].cellHandles){
	#endif
		double newVol, dVol;
                switch ( cell->info().fictious() ) {
                	case ( 3 ) : newVol = volumeCellTripleFictious ( cell ); break;
               		case ( 2 ) : newVol = volumeCellDoubleFictious ( cell ); break;
                	case ( 1 ) : newVol = volumeCellSingleFictious ( cell ); break;
			case ( 0 ) : newVol = volumeCell (cell ); break;
                	default: newVol = 0; break;}
                dVol=cell->info().volumeSign* ( newVol - cell->info().volume() );
		cell->info().dv() = dVol*invDeltaT;
                cell->info().volume() = newVol;
		if (defTolerance>0) { //if the criterion is not used, then we skip these updates a save a LOT of time when Nthreads > 1
			#pragma omp atomic
			totVol+=newVol;
			#pragma omp atomic
                	totDVol+=abs(dVol);}
        }
	if (defTolerance>0)  epsVolMax = totDVol/totVol;
	for (unsigned int n=0; n<flow.imposedF.size();n++) {
		flow.IFCells[n]->info().dv()+=flow.imposedF[n].second;
		flow.IFCells[n]->info().Pcondition=false;}
        if ( debug ) cout << "Updated volumes, total =" <<totVol<<", dVol="<<totDVol<<endl;
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
template<class Cellhandle>
Real TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::volumeCellSingleFictious ( Cellhandle cell )
{
        Vector3r V[3];
        int b=0;
        int w=0;
        cell->info().volumeSign=1;
        Real Wall_coordinate=0;

        for ( int y=0;y<4;y++ ) {
                if ( ! ( cell->vertex ( y )->info().isFictious ) ) {
                        V[w]=positionBufferCurrent[cell->vertex ( y )->info().id()].pos;
			w++;
                } else {
                        b = cell->vertex ( y )->info().id();
                        const shared_ptr<Body>& wll = Body::byId ( b , scene );
                        if ( !solver->boundary ( b ).useMaxMin ) Wall_coordinate = wll->state->pos[solver->boundary ( b ).coordinate]+ ( solver->boundary ( b ).normal[solver->boundary ( b ).coordinate] ) *wallThickness/2.;
                        else Wall_coordinate = solver->boundary ( b ).p[solver->boundary ( b ).coordinate];
                }
        }
        Real Volume = 0.5* ( ( V[0]-V[1] ).cross ( V[0]-V[2] ) ) [solver->boundary ( b ).coordinate] * ( 0.33333333333* ( V[0][solver->boundary ( b ).coordinate]+ V[1][solver->boundary ( b ).coordinate]+ V[2][solver->boundary ( b ).coordinate] ) - Wall_coordinate );
        return abs ( Volume );
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
template<class Cellhandle>
Real TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::volumeCellDoubleFictious ( Cellhandle cell )
{
        Vector3r A=Vector3r::Zero(), AS=Vector3r::Zero(),B=Vector3r::Zero(), BS=Vector3r::Zero();

        cell->info().volumeSign=1;
        int b[2];
        int coord[2];
        Real Wall_coordinate[2];
        int j=0;
        bool first_sph=true;

        for ( int g=0;g<4;g++ ) {
                if ( cell->vertex ( g )->info().isFictious ) {
                        b[j] = cell->vertex ( g )->info().id();
                        coord[j]=solver->boundary ( b[j] ).coordinate;
                        if ( !solver->boundary ( b[j] ).useMaxMin ) Wall_coordinate[j] = positionBufferCurrent[b[j]].pos[coord[j]] + ( solver->boundary ( b[j] ).normal[coord[j]] ) *wallThickness/2.;
                        else Wall_coordinate[j] = solver->boundary ( b[j] ).p[coord[j]];
                        j++;
                } else if ( first_sph ) {
                        A=AS=/*AT=*/ positionBufferCurrent[cell->vertex(g)->info().id()].pos;
                        first_sph=false;
                } else {
                        B=BS=/*BT=*/ positionBufferCurrent[cell->vertex(g)->info().id()].pos;;
                }
        }
        AS[coord[0]]=BS[coord[0]] = Wall_coordinate[0];

        //first pyramid with triangular base (A,B,BS)
        Real Vol1=0.5* ( ( A-BS ).cross ( B-BS ) ) [coord[1]]* ( 0.333333333* ( 2*B[coord[1]]+A[coord[1]] )-Wall_coordinate[1] );
        //second pyramid with triangular base (A,AS,BS)
        Real Vol2=0.5* ( ( AS-BS ).cross ( A-BS ) ) [coord[1]]* ( 0.333333333* ( B[coord[1]]+2*A[coord[1]] )-Wall_coordinate[1] );
        return abs ( Vol1+Vol2 );
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
template<class Cellhandle>
Real TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::volumeCellTripleFictious ( Cellhandle cell )
{
        Vector3r A;

        int b[3];
        int coord[3];
        Real Wall_coordinate[3];
        int j=0;
        cell->info().volumeSign=1;

        for ( int g=0;g<4;g++ ) {
                if ( cell->vertex ( g )->info().isFictious ) {
                        b[j] = cell->vertex ( g )->info().id();
                        coord[j]=solver->boundary ( b[j] ).coordinate;
                        const shared_ptr<Body>& wll = Body::byId ( b[j] , scene );
                        if ( !solver->boundary ( b[j] ).useMaxMin ) Wall_coordinate[j] = wll->state->pos[coord[j]] + ( solver->boundary ( b[j] ).normal[coord[j]] ) *wallThickness/2.;
                        else Wall_coordinate[j] = solver->boundary ( b[j] ).p[coord[j]];
                        j++;
                } else {
                        const shared_ptr<Body>& sph = Body::byId ( cell->vertex ( g )->info().id(), scene );
                        A= ( sph->state->pos );
                }
        }
        Real Volume = ( A[coord[0]]-Wall_coordinate[0] ) * ( A[coord[1]]-Wall_coordinate[1] ) * ( A[coord[2]]-Wall_coordinate[2] );
        return abs ( Volume );
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
template<class Cellhandle>
Real TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::volumeCell ( Cellhandle cell )
{
	static const Real inv6 = 1/6.;
	const Vector3r& p0 = positionBufferCurrent[cell->vertex ( 0 )->info().id()].pos;
	const Vector3r& p1 = positionBufferCurrent[cell->vertex ( 1 )->info().id()].pos;
	const Vector3r& p2 = positionBufferCurrent[cell->vertex ( 2 )->info().id()].pos;
	const Vector3r& p3 = positionBufferCurrent[cell->vertex ( 3 )->info().id()].pos;
	Real volume = inv6 * ((p0-p1).cross(p0-p2)).dot(p0-p3);
        if ( ! ( cell->info().volumeSign ) ) cell->info().volumeSign= ( volume>0 ) ?1:-1;
        return volume;
}
template< class _CellInfo, class _VertexInfo, class _Tesselation, class solverT >
void TemplateFlowEngine<_CellInfo,_VertexInfo,_Tesselation,solverT>::computeViscousForces ( Solver& flow )
{
	if (normalLubrication || shearLubrication || viscousShear){
		if ( debug ) cout << "Application of viscous forces" << endl;
		if ( debug ) cout << "Number of edges = " << flow.edgeIds.size() << endl;
		for ( unsigned int k=0; k<flow.shearLubricationForces.size(); k++ ) flow.shearLubricationForces[k]=Vector3r::Zero();
		for ( unsigned int k=0; k<flow.shearLubricationTorques.size(); k++ ) flow.shearLubricationTorques[k]=Vector3r::Zero();
		for ( unsigned int k=0; k<flow.pumpLubricationTorques.size(); k++ ) flow.pumpLubricationTorques[k]=Vector3r::Zero();
		for ( unsigned int k=0; k<flow.twistLubricationTorques.size(); k++ ) flow.twistLubricationTorques[k]=Vector3r::Zero();
		for ( unsigned int k=0; k<flow.shearLubricationBodyStress.size(); k++) flow.shearLubricationBodyStress[k]=Matrix3r::Zero();
		for ( unsigned int k=0; k<flow.normalLubricationForce.size(); k++ ) flow.normalLubricationForce[k]=Vector3r::Zero();
		for ( unsigned int k=0; k<flow.normalLubricationBodyStress.size(); k++) flow.normalLubricationBodyStress[k]=Matrix3r::Zero();

		typedef typename Solver::Tesselation Tesselation;
		const Tesselation& Tes = flow.T[flow.currentTes];
		flow.deltaShearVel.clear(); flow.normalV.clear(); flow.deltaNormVel.clear(); flow.surfaceDistance.clear(); flow.onlySpheresInteractions.clear(); flow.normalStressInteraction.clear(); flow.shearStressInteraction.clear();


		for ( int i=0; i< ( int ) flow.edgeIds.size(); i++ ) {
			const VertexInfo& vi1 = *flow.edgeIds[i].first;
			const VertexInfo& vi2 = *flow.edgeIds[i].second;
			const int& id1 = vi1.id();
			const int& id2 = vi2.id();
			
			int hasFictious= Tes.vertex ( id1 )->info().isFictious +  Tes.vertex ( id2 )->info().isFictious;
			if (hasFictious>0 or id1==id2) continue;
			const shared_ptr<Body>& sph1 = Body::byId ( id1, scene );
			const shared_ptr<Body>& sph2 = Body::byId ( id2, scene );
			Sphere* s1=YADE_CAST<Sphere*> ( sph1->shape.get() );
			Sphere* s2=YADE_CAST<Sphere*> ( sph2->shape.get() );
			const Real& r1 = s1->radius;
			const Real& r2 = s2->radius;
			Vector3r deltaV; Real deltaNormV; Vector3r deltaShearV;
			Vector3r O1O2Vector; Real O1O2; Vector3r normal; Real surfaceDist; Vector3r O1CVector; Vector3r O2CVector;Real meanRad ;Real Rh; Vector3r deltaAngVel; Vector3r deltaShearAngVel;
			Vector3r shearLubF; Vector3r normaLubF; Vector3r pumpT; Vector3r deltaAngNormVel; Vector3r twistT; Vector3r angVel1; Vector3r angVel2; 
		//FIXME: if periodic and velGrad!=0, then deltaV should account for velGrad, not the case currently
			if ( !hasFictious ){
				O1O2Vector = sph2->state->pos + makeVector3r(vi2.ghostShift()) - sph1->state->pos - makeVector3r(vi1.ghostShift());
				O1O2 = O1O2Vector.norm(); 
				normal= (O1O2Vector/O1O2);
				surfaceDist = O1O2 - r2 - r1;
				O1CVector = (O1O2/2. + (pow(r1,2) - pow(r2,2)) / (2.*O1O2))*normal;
				O2CVector = -(O1O2Vector - O1CVector);
				meanRad = (r2 + r1)/2.;
				Rh = (r1 < r2)? surfaceDist + 0.45 * r1 : surfaceDist + 0.45 * r2;
				deltaV = (sph2->state->vel + sph2->state->angVel.cross(-r2 * normal)) - (sph1->state->vel+ sph1->state->angVel.cross(r1 * normal));
				angVel1 = sph1->state->angVel;
				angVel2 = sph2->state->angVel;
				deltaAngVel = sph2->state->angVel - sph1->state->angVel;

			} else {
				if ( hasFictious==1 ) {//for the fictious sphere, use velocity of the boundary, not of the body
					bool v1fictious = Tes.vertex ( id1 )->info().isFictious;
					int bnd = v1fictious? id1 : id2;
					int coord = flow.boundary(bnd).coordinate;
					O1O2 = v1fictious ? abs((sph2->state->pos + makeVector3r(Tes.vertex(id2)->info().ghostShift()))[coord] - flow.boundary(bnd).p[coord]) : abs((sph1->state->pos + makeVector3r(Tes.vertex(id1)->info().ghostShift()))[coord] - flow.boundary(bnd).p[coord]);
					if(v1fictious)
						normal = makeVector3r(flow.boundary(id1).normal);
					else
						normal = -makeVector3r(flow.boundary(id2).normal);
					O1O2Vector = O1O2 * normal;
					meanRad = v1fictious ? r2:r1;
					surfaceDist = O1O2 - meanRad;
					if (v1fictious){
						O1CVector = Vector3r::Zero();
						O2CVector = - O1O2Vector;}
					else{
						O1CVector =  O1O2Vector;
						O2CVector = Vector3r::Zero();}
				
					Rh = surfaceDist + 0.45 * meanRad;
					Vector3r v1 = ( Tes.vertex ( id1 )->info().isFictious ) ? flow.boundary ( id1 ).velocity:sph1->state->vel + sph1->state->angVel.cross(r1 * normal);
					Vector3r v2 = ( Tes.vertex ( id2 )->info().isFictious ) ? flow.boundary ( id2 ).velocity:sph2->state->vel + sph2->state->angVel.cross(-r2 * (normal));
					deltaV = v2-v1;
					angVel1 = ( Tes.vertex ( id1 )->info().isFictious ) ? Vector3r::Zero() : sph1->state->angVel;
					angVel2 = ( Tes.vertex ( id2 )->info().isFictious ) ? Vector3r::Zero() : sph2->state->angVel;
					deltaAngVel = angVel2 - angVel1;
				}
			}
			deltaShearV = deltaV - ( normal.dot ( deltaV ) ) *normal;
			deltaShearAngVel = deltaAngVel - ( normal.dot ( deltaAngVel ) ) *normal;
			flow.deltaShearVel.push_back(deltaShearV);
			flow.normalV.push_back(normal);
			flow.surfaceDistance.push_back(max(surfaceDist, 0.) + eps*meanRad);

			/// Compute the  shear Lubrication force and torque on each particle
			
			if (shearLubrication)
				shearLubF = flow.computeShearLubricationForce(deltaShearV,surfaceDist,i,eps,O1O2,meanRad);
			else if (viscousShear) 
				shearLubF = flow.computeViscousShearForce ( deltaShearV, i , Rh);
				
			if (viscousShear || shearLubrication){

				flow.shearLubricationForces[id1]+=shearLubF;
				flow.shearLubricationForces[id2]+=(-shearLubF);
				flow.shearLubricationTorques[id1]+=O1CVector.cross(shearLubF);
				flow.shearLubricationTorques[id2]+=O2CVector.cross(-shearLubF);
				
				/// Compute the  pump Lubrication torque on each particle
				
				if (pumpTorque){
					pumpT = flow.computePumpTorque(deltaShearAngVel, surfaceDist, i, eps, meanRad );
					flow.pumpLubricationTorques[id1]+=(-pumpT);
					flow.pumpLubricationTorques[id2]+=pumpT;}
				
				/// Compute the  twist Lubrication torque on each particle
				
				if (twistTorque){
					deltaAngNormVel = (normal.dot(deltaAngVel))*normal ;
					twistT = flow.computeTwistTorque(deltaAngNormVel, surfaceDist, i, eps, meanRad );
					flow.twistLubricationTorques[id1]+=(-twistT);
					flow.twistLubricationTorques[id2]+=twistT;
				}
			}		
			/// Compute the viscous shear stress on each particle
			
			if (viscousShearBodyStress){
				flow.shearLubricationBodyStress[id1] += shearLubF * O1CVector.transpose()/ (4.0/3.0 *3.14* pow(r1,3));
				flow.shearLubricationBodyStress[id2] += (-shearLubF) * O2CVector.transpose()/ (4.0/3.0 *3.14* pow(r2,3));
				flow.shearStressInteraction.push_back(shearLubF * O1O2Vector.transpose()/(4.0/3.0 *3.14* pow(r1,3)));
				}

			/// Compute the normal lubrication force applied on each particle
			
			if (normalLubrication){
				deltaNormV = normal.dot(deltaV);
				flow.deltaNormVel.push_back(deltaNormV * normal);
				normaLubF = flow.computeNormalLubricationForce (deltaNormV, surfaceDist, i,eps,stiffness,scene->dt,meanRad)*normal;
				flow.normalLubricationForce[id1]+=normaLubF;
				flow.normalLubricationForce[id2]+=(-normaLubF);

				/// Compute the normal lubrication stress on each particle
				
				if (viscousNormalBodyStress){
					flow.normalLubricationBodyStress[id1] += normaLubF * O1CVector.transpose()/ (4.0/3.0 *3.14* pow(r1,3));
					flow.normalLubricationBodyStress[id2] += (-normaLubF) *O2CVector.transpose() / (4.0/3.0 *3.14* pow(r2,3));
					flow.normalStressInteraction.push_back(normaLubF * O1O2Vector.transpose()/(4.0/3.0 *3.14* pow(r1,3)));
				}
			}
			
			if (!hasFictious)
				flow.onlySpheresInteractions.push_back(i);
				
		}
	}
}

#endif //FLOW_ENGINE

#endif /* YADE_CGAL */


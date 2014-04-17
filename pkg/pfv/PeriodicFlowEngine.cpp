/*************************************************************************
*  Copyright (C) 2009 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*  Copyright (C) 2012 by Donia Marzougui <donia.marzougui@grenoble-inp.fr>*
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef YADE_CGAL
#ifdef FLOW_ENGINE

/// The periodic variant of FlowEngine is defined here. It should become a template class for more flexibility.
/// It is a bit more complicated as for FlowEngine, though, because we need template inheriting from template, which breaks YADE_CLASS_XXX logic_error
/// See below the commented exemple, for a possible solution

#define TEMPLATE_FLOW_NAME FlowEngine_PeriodicInfo
#include <yade/pkg/pfv/FlowEngine.hpp>

class PeriodicCellInfo : public FlowCellInfo
{	
	public:
	static CVector gradP;
	//for real cell, baseIndex is the rank of the cell in cellHandles. For ghost cells, it is the baseIndex of the corresponding real cell.
	//Unlike ordinary index, baseIndex is also indexing cells with imposed pressures
	int baseIndex;
	int period[3];
	static CVector hSize[3];
	static CVector deltaP;
	int ghost;
	Real* _pression;
	PeriodicCellInfo (void){
		_pression=&pression;
		period[0]=period[1]=period[2]=0;
		baseIndex=-1;
		volumeSign=0;}
	~PeriodicCellInfo (void) {}

	inline const double shiftedP (void) const {return isGhost? (*_pression)+pShift() :(*_pression) ;}
	inline const double pShift (void) const {return deltaP[0]*period[0] + deltaP[1]*period[1] +deltaP[2]*period[2];}
// 	inline const double p (void) {return shiftedP();}
	inline void setP (const Real& p) {pression=p;}
	bool isReal (void) {return !(isFictious || isGhost);}
};


class PeriodicVertexInfo : public FlowVertexInfo {
	public:
	PeriodicVertexInfo& operator= (const CVector &u) { CVector::operator= (u); return *this; }
	PeriodicVertexInfo& operator= (const float &scalar) { s=scalar; return *this; }
	PeriodicVertexInfo& operator= (const unsigned int &id) { i= id; return *this; }
	int period[3];
	//FIXME: the name is misleading, even non-ghost can be out of the period and therefore they need to be shifted as well
	inline const CVector ghostShift (void) const {
		return period[0]*PeriodicCellInfo::hSize[0]+period[1]*PeriodicCellInfo::hSize[1]+period[2]*PeriodicCellInfo::hSize[2];}
	PeriodicVertexInfo (void) {isFictious=false; s=0; i=0; period[0]=period[1]=period[2]=0; isGhost=false;}
	bool isReal (void) {return !(isFictious || isGhost);}
};

typedef CGT::TriangulationTypes<PeriodicVertexInfo,PeriodicCellInfo>			PeriFlowTriangulationTypes;
typedef CGT::PeriodicTesselation<CGT::_Tesselation<PeriFlowTriangulationTypes> >	PeriFlowTesselation;
#ifdef LINSOLV
#define _PeriFlowSolver CGT::PeriodicFlowLinSolv<PeriFlowTesselation>
#else
#define _PeriFlowSolver CGT::PeriodicFlow<PeriFlowTesselation>
#endif

typedef TemplateFlowEngine<	PeriodicCellInfo,
				PeriodicVertexInfo,
				CGT::PeriodicTesselation<CGT::_Tesselation<CGT::TriangulationTypes<PeriodicVertexInfo,PeriodicCellInfo> > >,
				_PeriFlowSolver>
				FlowEngine_PeriodicInfo;

REGISTER_SERIALIZABLE(FlowEngine_PeriodicInfo);
YADE_PLUGIN((FlowEngine_PeriodicInfo));


class PeriodicFlowEngine : public FlowEngine_PeriodicInfo
{
	public :
		void triangulate (FlowSolver& flow);
		void buildTriangulation (Real pzero, FlowSolver& flow);
		void initializeVolumes (FlowSolver&  flow);
		void updateVolumes (FlowSolver&  flow);
		Real volumeCell (CellHandle cell);

		Real volumeCellSingleFictious (CellHandle cell);
		inline void locateCell(CellHandle baseCell, unsigned int& index, int& baseIndex, FlowSolver& flow, unsigned int count=0);
		Vector3r meanVelocity();

		virtual ~PeriodicFlowEngine();

		virtual void action();
		//Cache precomputed values for pressure shifts, based on current hSize and pGrad
		void preparePShifts();
		
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(PeriodicFlowEngine,FlowEngine_PeriodicInfo,"A variant of :yref:`FlowEngine` implementing periodic boundary conditions. The API is very similar.",
		((Real,duplicateThreshold, 0.06,,"distance from cell borders that will triger periodic duplication in the triangulation |yupdate|"))
		((Vector3r, gradP, Vector3r::Zero(),,"Macroscopic pressure gradient"))
		,,
		wallIds=vector<int>(6,-1);
		solver = shared_ptr<FlowSolver> (new FlowSolver);
		epsVolMax=epsVolCumulative=retriangulationLastIter=0;
		ReTrg=1;
		first=true;
		,
		//nothing special to define, we re-use FlowEngine methods
		//.def("meanVelocity",&PeriodicFlowEngine::meanVelocity,"measure the mean velocity in the period")
		)
		DECLARE_LOGGER;


};
REGISTER_SERIALIZABLE(PeriodicFlowEngine);


CVector PeriodicCellInfo::hSize[]={CVector(),CVector(),CVector()};
CVector PeriodicCellInfo::deltaP=CVector();
CVector PeriodicCellInfo::gradP=CVector();

CREATE_LOGGER ( PeriodicFlowEngine );

PeriodicFlowEngine::~PeriodicFlowEngine(){}

void PeriodicFlowEngine:: action()
{
        if ( !isActivated ) return;
	timingDeltas->start();
	preparePShifts();
	setPositionsBuffer(true);
	if (first) {
		if (multithread) setPositionsBuffer(false);
		cachedCell= Cell(*(scene->cell));
		buildTriangulation(pZero,*solver);
		if (solver->errorCode>0) {LOG_INFO("triangulation error, pausing"); Omega::instance().pause(); return;}
		initializeVolumes(*solver); backgroundSolver=solver; backgroundCompleted=true;}
//         if ( first ) {buildTriangulation ( pZero ); updateTriangulation = false; initializeVolumes();}
	
	timingDeltas->checkpoint("Triangulating");
        updateVolumes (*solver);
        epsVolCumulative += epsVolMax;
	retriangulationLastIter++;
	if (!updateTriangulation) updateTriangulation = // If not already set true by another function of by the user, check conditions
		(defTolerance>0 && epsVolCumulative > defTolerance) || retriangulationLastIter>meshUpdateInterval;

	timingDeltas->checkpoint("Update_Volumes");

	///compute flow and and forces here
	if (pressureForce){
		solver->gaussSeidel(scene->dt);
		timingDeltas->checkpoint("Gauss-Seidel");
		solver->computeFacetForcesWithCache();}
	timingDeltas->checkpoint("compute_Pressure_Forces");

        ///compute vicscous forces
        scene->forces.sync();
        computeViscousForces(*solver);
	timingDeltas->checkpoint("compute_Viscous_Forces");
	Vector3r force;
	Vector3r torque;
	const Tesselation& Tes = solver->T[solver->currentTes];
	for (int id=0; id<=Tes.maxId; id++){
		assert (Tes.vertexHandles[id] != NULL);
		const Tesselation::VertexInfo& v_info = Tes.vertexHandles[id]->info();
		force =(pressureForce) ? Vector3r ( ( v_info.forces ) [0],v_info.forces[1],v_info.forces[2] ) : Vector3r(0,0,0);
		torque = Vector3r(0,0,0);
                if (shearLubrication || viscousShear){
			force = force +solver->shearLubricationForces[v_info.id()];
			torque = torque +solver->shearLubricationTorques[v_info.id()];
			if (pumpTorque)
				torque = torque +solver->pumpLubricationTorques[v_info.id()];
			if (twistTorque)
				torque = torque +solver->twistLubricationTorques[v_info.id()];
		}
		
		if (normalLubrication)
			force = force + solver->normalLubricationForce[v_info.id()];
		scene->forces.addForce ( v_info.id(), force);
		scene->forces.addTorque ( v_info.id(), torque);
	}
        ///End Compute flow and forces
	timingDeltas->checkpoint("Applying Forces");
	if (multithread && !first) {
		while (updateTriangulation && !backgroundCompleted) { /*cout<<"sleeping..."<<sleeping++<<endl;*/ 	boost::this_thread::sleep(boost::posix_time::microseconds(1000));}
		if (updateTriangulation || ellapsedIter>(0.5*meshUpdateInterval)) {
			if (useSolver==0) LOG_ERROR("background calculations not available for Gauss-Seidel");
			if (fluidBulkModulus>0 || doInterpolate) solver->interpolate (solver->T[solver->currentTes], backgroundSolver->T[backgroundSolver->currentTes]);
			solver=backgroundSolver;
			backgroundSolver = shared_ptr<FlowSolver> (new FlowSolver);
			//Copy imposed pressures/flow from the old solver
			backgroundSolver->imposedP = vector<pair<CGT::Point,Real> >(solver->imposedP);
			backgroundSolver->imposedF = vector<pair<CGT::Point,Real> >(solver->imposedF);
			setPositionsBuffer(false);
			cachedCell= Cell(*(scene->cell));
			backgroundCompleted=false;
			retriangulationLastIter=ellapsedIter;
			ellapsedIter=0;
			epsVolCumulative=0;
			boost::thread workerThread(&PeriodicFlowEngine::backgroundAction,this);
			workerThread.detach();
			initializeVolumes(*solver);
			computeViscousForces(*solver);
		}
		else if (debug && !first) {
			if (debug && !backgroundCompleted) cerr<<"still computing solver in the background"<<endl;
			ellapsedIter++;}
	} else {
	        if (updateTriangulation && !first) {
			cachedCell= Cell(*(scene->cell));
			buildTriangulation (pZero, *solver);
			initializeVolumes(*solver);
			computeViscousForces(*solver);
               		updateTriangulation = false;
			epsVolCumulative=0;
                	retriangulationLastIter=0;
                	ReTrg++;}
        }
        first=false;
	timingDeltas->checkpoint("Ending");
}


// void PeriodicFlowEngine::backgroundAction()
// {
// 	if (useSolver<1) {LOG_ERROR("background calculations not available for Gauss-Seidel"); return;}
//         buildTriangulation (pZero,*backgroundSolver);
// 	//FIXME: GS is computing too much, we need only matrix factorization in fact
// 	backgroundSolver->gaussSeidel(scene->dt);
// 	backgroundSolver->computeFacetForcesWithCache(/*onlyCache?*/ true);
// // 	boost::this_thread::sleep(boost::posix_time::seconds(10));
// 	backgroundCompleted = true;
// }

void PeriodicFlowEngine::triangulate( FlowSolver& flow )
{
        Tesselation& Tes = flow.T[flow.currentTes];
	vector<posData>& buffer = multithread ? positionBufferParallel : positionBufferCurrent;
	FOREACH ( const posData& b, buffer ) {
                if ( !b.exists || !b.isSphere || b.id==ignoredBody) continue;
                Vector3i period; Vector3r wpos;
		// FIXME: use "sheared" variant if the cell is sheared
		wpos=cachedCell.wrapPt ( b.pos,period );
		const Body::id_t& id = b.id;
		const Real& rad = b.radius;
                const Real& x = wpos[0];
                const Real& y = wpos[1];
                const Real& z = wpos[2];
                VertexHandle vh0=Tes.insert ( x, y, z, rad, id );
//                 VertexHandle vh0=Tes.insert ( b.pos[0], b.pos[1], b.pos[2], b.radius, b.id );
		if (vh0==NULL) {
			flow.errorCode = 2;
			LOG_ERROR("Vh NULL in PeriodicFlowEngine::triangulate(), check input data"); continue;}
		for ( int k=0;k<3;k++ ) vh0->info().period[k]=-period[k];
                const Vector3r cellSize ( cachedCell.getSize() );
		//FIXME: if hasShear, comment in
//                 wpos=scene->cell->unshearPt ( wpos );
                // traverse all periodic cells around the body, to see if any of them touches
                Vector3r halfSize= ( rad+duplicateThreshold ) *Vector3r ( 1,1,1 );
                Vector3r pmin,pmax;
                Vector3i i;
                for ( i[0]=-1; i[0]<=1; i[0]++ )
                        for ( i[1]=-1;i[1]<=1; i[1]++ )
                                for ( i[2]=-1; i[2]<=1; i[2]++ ) {
                                        if ( i[0]!=0 || i[1]!=0 || i[2]!=0 ) { // middle; already rendered above
                                        Vector3r pos2=wpos+Vector3r ( cellSize[0]*i[0],cellSize[1]*i[1],cellSize[2]*i[2] ); // shift, but without shear!
                                        pmin=pos2-halfSize;
                                        pmax=pos2+halfSize;
                                        if ( (pmin[0]<=cellSize[0]) && (pmax[0]>=0) && (pmin[1]<=cellSize[1]) && (pmax[1]>=0) && (pmin[2]<=cellSize[2]) && (pmax[2]>=0) ) {
                                                //with shear:
                                                //Vector3r pt=scene->cell->shearPt ( pos2 );
                                                //without shear:
                                                const Vector3r& pt= pos2;
                                                VertexHandle vh=Tes.insert ( pt[0],pt[1],pt[2],rad,id,false,id );
                                                for ( int k=0;k<3;k++ ) vh->info().period[k]=i[k]-period[k];}}
				}
		//re-assign the original vertex pointer since duplicates may have overwrite it
		Tes.vertexHandles[id]=vh0;
        }
        Tes.redirected=true;//By inserting one-by-one, we already redirected
        flow.shearLubricationForces.resize ( Tes.maxId+1 );
	flow.shearLubricationTorques.resize ( Tes.maxId+1 );
	flow.pumpLubricationTorques.resize ( Tes.maxId+1 );
	flow.twistLubricationTorques.resize ( Tes.maxId+1 );
	flow.shearLubricationBodyStress.resize ( Tes.maxId+1 );
	flow.normalLubricationForce.resize ( Tes.maxId+1 );
	flow.normalLubricationBodyStress.resize ( Tes.maxId+1 );
}


Real PeriodicFlowEngine::volumeCell ( CellHandle cell )
{
	static const Real inv6 = 1/6.;	
	const Vector3r p0 = positionBufferCurrent[cell->vertex(0)->info().id()].pos + makeVector3r(cell->vertex(0)->info().ghostShift());
	const Vector3r p1 = positionBufferCurrent[cell->vertex(1)->info().id()].pos + makeVector3r(cell->vertex(1)->info().ghostShift());
	const Vector3r p2 = positionBufferCurrent[cell->vertex(2)->info().id()].pos + makeVector3r(cell->vertex(2)->info().ghostShift());
	const Vector3r p3 = positionBufferCurrent[cell->vertex(3)->info().id()].pos + makeVector3r(cell->vertex(3)->info().ghostShift());
	Real volume = inv6*((p0-p1).cross(p0-p2)).dot(p0-p3);
        if ( ! ( cell->info().volumeSign ) ) cell->info().volumeSign= ( volume>0 ) ?1:-1;
        return volume;
}

Real PeriodicFlowEngine::volumeCellSingleFictious ( CellHandle cell )
{
        Vector3r V[3];
        int b=0;
        int w=0;
        cell->info().volumeSign=1;
        Real Wall_coordinate=0;

        for ( int y=0;y<4;y++ ) {
                if ( ! ( cell->vertex ( y )->info().isFictious ) ) {
                        const shared_ptr<Body>& sph = Body::byId ( cell->vertex ( y )->info().id(), scene );
                        V[w]=sph->state->pos+ makeVector3r ( cell->vertex ( y )->info().ghostShift() );
                        w++;
                } else {
                        b = cell->vertex ( y )->info().id();
                        const shared_ptr<Body>& wll = Body::byId ( b,scene );
                        if ( !solver->boundary ( b ).useMaxMin ) Wall_coordinate = wll->state->pos[solver->boundary ( b ).coordinate]+ ( solver->boundary ( b ).normal[solver->boundary ( b ).coordinate] ) *wallThickness/2.;
                        else Wall_coordinate = solver->boundary ( b ).p[solver->boundary ( b ).coordinate];
                }
        }
        Real Volume = 0.5* ( ( V[0]-V[1] ).cross ( V[0]-V[2] ) ) [solver->boundary ( b ).coordinate] * ( 0.33333333333* ( V[0][solver->boundary ( b ).coordinate]+ V[1][solver->boundary ( b ).coordinate]+ V[2][solver->boundary ( b ).coordinate] ) - Wall_coordinate );
        return abs ( Volume );
}


void PeriodicFlowEngine::locateCell ( CellHandle baseCell, unsigned int& index, int& baseIndex, FlowSolver& flow, unsigned int count)
{
        if (count>10) {
		LOG_ERROR("More than 10 attempts to locate a cell, duplicateThreshold may be too small, resulting in periodicity inconsistencies.");
		flow.errorCode=1; return;
	}
	PeriFlowTesselation::CellInfo& baseInfo = baseCell->info();
        //already located, return FIXME: is inline working correctly? else move this test outside the function, just before the calls
	if ( baseInfo.index>0 || baseInfo.isGhost ) return;
	RTriangulation& Tri = flow.T[flow.currentTes].Triangulation();
	Vector3r center ( 0,0,0 );
	Vector3i period;

	if (baseCell->info().fictious()==0)
		for ( int k=0;k<4;k++ ) center+= 0.25*makeVector3r (baseCell->vertex(k)->point());
	else {
		
		Real boundPos=0; int coord=0;
		for ( int k=0;k<4;k++ ) {
			if ( !baseCell->vertex ( k )->info().isFictious ) center+= 0.3333333333*makeVector3r ( baseCell->vertex ( k )->point() );
			else {
				coord=flow.boundary ( baseCell->vertex ( k )->info().id() ).coordinate;
				boundPos=flow.boundary ( baseCell->vertex ( k )->info().id() ).p[coord];}
		}
		center[coord]=boundPos;
	}
	Vector3r wdCenter= cachedCell.wrapPt ( center,period );
	if ( period[0]!=0 || period[1]!=0 || period[2]!=0 ) {
		if ( baseCell->info().index>0 ) {
			cout<<"indexed cell is found ghost!"<<baseInfo.index <<endl;
			baseInfo.isGhost=false;
			return;
		}
		CellHandle ch= Tri.locate ( CGT::Point ( wdCenter[0],wdCenter[1],wdCenter[2] )
// 					     ,/*hint*/ v0
					     );
		baseInfo.period[0]=period[0];
		baseInfo.period[1]=period[1];
		baseInfo.period[2]=period[2];
		//call recursively, since the returned cell could be also a ghost (especially if baseCell is a non-periodic type from the external contour
		locateCell ( ch,index,baseIndex,flow,++count );
		if ( ch==baseCell ) cerr<<"WTF!!"<<endl;
		//check consistency
		bool checkC=false;
		for (int kk=0; kk<4;kk++) if ((!baseCell->vertex(kk)->info().isGhost) && ((!baseCell->vertex(kk)->info().isFictious))) checkC = true;
		if (checkC) {
			bool checkV=true;
			for (int kk=0; kk<4;kk++) {
				checkV=false;
				for (int jj=0; jj<4;jj++)
					if (baseCell->vertex(kk)->info().id() == ch->vertex(jj)->info().id()) checkV = true;
				if (!checkV) {cerr <<"periodicity is broken"<<endl;
				for (int jj=0; jj<4;jj++) cerr<<baseCell->vertex(jj)->info().id()<<" ";
				cerr<<" vs. ";
				for (int jj=0; jj<4;jj++) cerr<<ch->vertex(jj)->info().id()<<" ";
				cerr<<endl;}
			}
		} else {
// 			bool checkV=true;
// 			for (int kk=0; kk<4;kk++) {
// 				checkV=false;
// 				for (int jj=0; jj<4;jj++)
// 					if (baseCell->vertex(kk)->info().id() == ch->vertex(jj)->info().id()) checkV = true;
// 				if (!checkV) {cerr <<"periodicity is broken (that's ok probably)"<<endl;
// 				for (int jj=0; jj<4;jj++) cerr<<baseCell->vertex(jj)->info().id()<<" ";
// 				cerr<<" vs. ";
// 				for (int jj=0; jj<4;jj++) cerr<<ch->vertex(jj)->info().id()<<" ";
// 				cerr<<endl;}
// 			}
		}

		baseInfo.isGhost=true;
		baseInfo._pression=& ( ch->info().p() );
		baseInfo.index=ch->info().index;
		baseInfo.baseIndex=ch->info().baseIndex;
		baseInfo.Pcondition=ch->info().Pcondition;
	} else {
		baseInfo.isGhost=false;
		//index is 1-based, if it is zero it is not initialized, we define it here
		if (  baseInfo.baseIndex<0 ){
			baseInfo.baseIndex=++baseIndex;
			if (!baseInfo.Pcondition) baseInfo.index=++index;}
	}
}

Vector3r PeriodicFlowEngine::meanVelocity()
{
        solver->averageRelativeCellVelocity();
        Vector3r meanVel ( 0,0,0 );
        Real volume=0;
        FiniteCellsIterator cell_end = solver->T[solver->currentTes].Triangulation().finite_cells_end();
        for ( FiniteCellsIterator cell = solver->T[solver->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {
		//We could also define velocity using cell's center
//                 if ( !cell->info().isReal() ) continue;
                if ( cell->info().isGhost ) continue;
                for ( int i=0;i<3;i++ )
                        meanVel[i]=meanVel[i]+ ( ( cell->info().averageVelocity() ) [i] * abs ( cell->info().volume() ) );
                volume+=abs ( cell->info().volume() );
        }
        return ( meanVel/volume );
}

void PeriodicFlowEngine::updateVolumes (FlowSolver& flow)
{
        if ( debug ) cout << "Updating volumes.............." << endl;
        Real invDeltaT = 1/scene->dt;
        double newVol, dVol;
        epsVolMax=0;
        Real totVol=0;
        Real totDVol=0;
        Real totVol0=0;
        Real totVol1=0;

	FOREACH(CellHandle& cell, flow.T[flow.currentTes].cellHandles){
                switch ( cell->info().fictious() ) {
                case ( 1 ) :
                        newVol = volumeCellSingleFictious ( cell );
                        totVol1+=newVol;
                        break;
                case ( 0 ) :
                        newVol = volumeCell ( cell );
                        totVol0+=newVol;
                        break;
                default:
                        newVol = 0;
                        break;
                }
                totVol+=newVol;
                dVol=cell->info().volumeSign * ( newVol - cell->info().volume() );
                totDVol+=dVol;
                epsVolMax = max ( epsVolMax, abs ( dVol/newVol ) );
                cell->info().dv() = dVol * invDeltaT;
                cell->info().volume() = newVol;
        }
        if ( debug ) cout << "Updated volumes, total =" <<totVol<<", dVol="<<totDVol<<" "<< totVol0<<" "<< totVol1<<endl;
}


void PeriodicFlowEngine::initializeVolumes (FlowSolver& flow)
{
        FiniteVerticesIterator vertices_end = flow.T[flow.currentTes].Triangulation().finite_vertices_end();
        CGT::CVector Zero ( 0,0,0 );
        for ( FiniteVerticesIterator V_it = flow.T[flow.currentTes].Triangulation().finite_vertices_begin(); V_it!= vertices_end; V_it++ ) V_it->info().forces=Zero;

	FOREACH(CellHandle& cell, flow.T[flow.currentTes].cellHandles){
		switch ( cell->info().fictious() )
		{
			case ( 0 ) : cell->info().volume() = volumeCell ( cell ); break;
			case ( 1 ) : cell->info().volume() = volumeCellSingleFictious ( cell ); break;
			default:  cell->info().volume() = 0; break;
		}
		//FIXME: the void volume is negative sometimes, hence crashing...
		if (flow.fluidBulkModulus>0) { cell->info().invVoidVolume() = 1. / (max(0.1*cell->info().volume(),abs(cell->info().volume()) - flow.volumeSolidPore(cell)) ); }
	}
        if ( debug ) cout << "Volumes initialised." << endl;
}

void PeriodicFlowEngine::buildTriangulation ( double pZero, FlowSolver& flow)
{
        if (first) flow.currentTes=0;
        else {
                flow.currentTes=!flow.currentTes;
                if ( debug ) cout << "--------RETRIANGULATION-----------" << endl;}
	flow.resetNetwork();
        initSolver(flow);
        addBoundary ( flow );
        if ( debug ) cout << endl << "Added boundaries------" << endl << endl;
        triangulate (flow);
        if ( debug ) cout << endl << "Tesselating------" << endl << endl;
        flow.T[flow.currentTes].compute();
        flow.defineFictiousCells();

        //FIXME: this is already done in addBoundary(?)
        boundaryConditions ( flow );
	if ( debug ) cout << endl << "boundaryConditions------" << endl << endl;
        flow.initializePressure ( pZero );
	if ( debug ) cout << endl << "initializePressure------" << endl << endl;
        // Define the ghost cells and add indexes to the cells inside the period (the ones that will contain the pressure unknowns)
        //This must be done after boundary conditions and initialize pressure, else the indexes are not good (not accounting imposedP): FIXME
        unsigned int index=0;
	int baseIndex=-1;
        FlowSolver::Tesselation& Tes = flow.T[flow.currentTes];
	Tes.cellHandles.resize(Tes.Triangulation().number_of_finite_cells());
	const FiniteCellsIterator cellend=Tes.Triangulation().finite_cells_end();
        for ( FiniteCellsIterator cell=Tes.Triangulation().finite_cells_begin(); cell!=cellend; cell++ ){
                locateCell ( cell,index,baseIndex,flow );
		if (flow.errorCode>0) return;
		//Fill this vector than can be later used to speedup loops
		if (!cell->info().isGhost) Tes.cellHandles[cell->info().baseIndex]=cell;
	}
	Tes.cellHandles.resize(baseIndex+1);

	if ( debug ) cout << endl << "locateCell------" << endl << endl;
        flow.computePermeability ( );
        porosity = flow.vPoralPorosity/flow.vTotalPorosity;
        flow.tolerance=tolerance;flow.relax=relax;
	
        flow.displayStatistics ();
        //FIXME: check interpolate() for the periodic case, at least use the mean pressure from previous step.
	if ( !first && !multithread && (useSolver==0 || fluidBulkModulus>0 || doInterpolate)) flow.interpolate ( flow.T[!flow.currentTes], Tes );
// 	if ( !first && (useSolver==0 || fluidBulkModulus>0)) flow.interpolate ( flow.T[!flow.currentTes], flow.T[flow.currentTes] );
	
        if ( waveAction ) flow.applySinusoidalPressure ( Tes.Triangulation(), sineMagnitude, sineAverage, 30 );

        if (normalLubrication || shearLubrication || viscousShear) flow.computeEdgesSurfaces();
	if ( debug ) cout << endl << "end buildTri------" << endl << endl;
}

void PeriodicFlowEngine::preparePShifts()
{
	CellInfo::gradP = makeCgVect ( gradP );
        CellInfo::hSize[0] = makeCgVect ( scene->cell->hSize.col ( 0 ) );
        CellInfo::hSize[1] = makeCgVect ( scene->cell->hSize.col ( 1 ) );
        CellInfo::hSize[2] = makeCgVect ( scene->cell->hSize.col ( 2 ) );
        CellInfo::deltaP=CGT::CVector (
                                              CellInfo::hSize[0]*CellInfo::gradP,
                                              CellInfo::hSize[1]*CellInfo::gradP,
                                              CellInfo::hSize[2]*CellInfo::gradP );
}


YADE_PLUGIN((PeriodicFlowEngine));
#endif //FLOW_ENGINE

#endif /* YADE_CGAL */
/*************************************************************************
*  Copyright (C) 2009 by Emanuele Catalano                               *
*  emanuele.catalano@hmg.inpg.fr                                         *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef YADE_CGAL

#ifdef FLOW_ENGINE
#include "FlowEngine.hpp"
#include<yade/core/Scene.hpp>
#include<yade/lib/base/Math.hpp>
#include<yade/pkg/dem/TesselationWrapper.hpp>
#include<yade/pkg/common/Sphere.hpp>
#include<yade/pkg/common/Wall.hpp>
#include<yade/pkg/common/Box.hpp>
#include <sys/stat.h>
#include <sys/types.h>

CREATE_LOGGER ( FlowEngine );
CREATE_LOGGER ( PeriodicFlowEngine );

CGT::Vecteur makeCgVect ( const Vector3r& yv ) {return CGT::Vecteur ( yv[0],yv[1],yv[2] );}
CGT::Point makeCgPoint ( const Vector3r& yv ) {return CGT::Point ( yv[0],yv[1],yv[2] );}
Vector3r makeVector3r ( const CGT::Point& yv ) {return Vector3r ( yv[0],yv[1],yv[2] );}
Vector3r makeVector3r ( const CGT::Vecteur& yv ) {return Vector3r ( yv[0],yv[1],yv[2] );}


FlowEngine::~FlowEngine()
{
}

const int facetVertices [4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};

void FlowEngine::action()
{
        if ( !isActivated ) return;
        timingDeltas->start();

        if ( first ) Build_Triangulation ( P_zero,solver );
        timingDeltas->checkpoint ( "Triangulating" );
        UpdateVolumes ( solver );

        Eps_Vol_Cumulative += eps_vol_max;
        if ( ( Eps_Vol_Cumulative > EpsVolPercent_RTRG || retriangulationLastIter>1000 ) && retriangulationLastIter>10 ) {
                Update_Triangulation = true;
                Eps_Vol_Cumulative=0;
                retriangulationLastIter=0;
                ReTrg++;
        } else  retriangulationLastIter++;
        timingDeltas->checkpoint ( "Update_Volumes" );

        ///Compute flow and and forces here
	solver->GaussSeidel(scene->dt);
        timingDeltas->checkpoint ( "Gauss-Seidel" );
        if ( save_mgpost ) solver->MGPost();
        if ( !CachedForces ) solver->ComputeFacetForces();
        else solver->ComputeFacetForcesWithCache();
        timingDeltas->checkpoint ( "Compute_Forces" );

        ///Application of vicscous forces
        scene->forces.sync();
        if ( viscousShear ) ApplyViscousForces ( *solver );

        Finite_vertices_iterator vertices_end = solver->T[solver->currentTes].Triangulation().finite_vertices_end();
        for ( Finite_vertices_iterator V_it = solver->T[solver->currentTes].Triangulation().finite_vertices_begin(); V_it !=  vertices_end; V_it++ ) {
                if ( !viscousShear )
                        scene->forces.addForce ( V_it->info().id(), Vector3r ( ( V_it->info().forces ) [0],V_it->info().forces[1],V_it->info().forces[2] ) );
                else
                        scene->forces.addForce ( V_it->info().id(), Vector3r ( ( V_it->info().forces ) [0],V_it->info().forces[1],V_it->info().forces[2] ) +solver->viscousShearForces[V_it->info().id() ] );
        }
        ///End Compute flow and forces
        timingDeltas->checkpoint ( "Applying Forces" );

        if ( Update_Triangulation && !first ) {
                Build_Triangulation ( P_zero, solver );
                Update_Triangulation = false;
        }

        if ( velocity_profile ) /*flow->FluidVelocityProfile();*/solver->Average_Fluid_Velocity();

        first=false;
        timingDeltas->checkpoint ( "Ending" );
}

template<class Solver>

void FlowEngine::BoundaryConditions ( Solver& flow )
{
        if ( flow->y_min_id>=0 ) {
                flow->boundary ( flow->y_min_id ).flowCondition=Flow_imposed_BOTTOM_Boundary;
                flow->boundary ( flow->y_min_id ).value=Pressure_BOTTOM_Boundary;
                flow->boundary ( flow->y_min_id ).velocity = Vector3r::Zero();
        }
        if ( flow->y_max_id>=0 ) {
                flow->boundary ( flow->y_max_id ).flowCondition=Flow_imposed_TOP_Boundary;
                flow->boundary ( flow->y_max_id ).value=Pressure_TOP_Boundary;
                flow->boundary ( flow->y_max_id ).velocity = topBoundaryVelocity;
        }
        if ( flow->x_max_id>=0 ) {
                flow->boundary ( flow->x_max_id ).flowCondition=Flow_imposed_RIGHT_Boundary;
                flow->boundary ( flow->x_max_id ).value=Pressure_RIGHT_Boundary;
                flow->boundary ( flow->x_max_id ).velocity = Vector3r::Zero();
        }
        if ( flow->x_min_id>=0 ) {
                flow->boundary ( flow->x_min_id ).flowCondition=Flow_imposed_LEFT_Boundary;
                flow->boundary ( flow->x_min_id ).value=Pressure_LEFT_Boundary;
                flow->boundary ( flow->x_min_id ).velocity = Vector3r::Zero();
        }
        if ( flow->z_max_id>=0 ) {
                flow->boundary ( flow->z_max_id ).flowCondition=Flow_imposed_FRONT_Boundary;
                flow->boundary ( flow->z_max_id ).value=Pressure_FRONT_Boundary;
                flow->boundary ( flow->z_max_id ).velocity = Vector3r::Zero();
        }
        if ( flow->z_min_id>=0 ) {
                flow->boundary ( flow->z_min_id ).flowCondition=Flow_imposed_BACK_Boundary;
                flow->boundary ( flow->z_min_id ).value=Pressure_BACK_Boundary;
                flow->boundary ( flow->z_min_id ).velocity = Vector3r::Zero();
        }
}

template<class Solver>
void FlowEngine::setImposedPressure ( unsigned int cond, Real p,Solver& flow )
{
        if ( cond>=flow->imposedP.size() ) LOG_ERROR ( "Setting p with cond higher than imposedP size." );
        flow->imposedP[cond].second=p;
        //force immediate update of boundary conditions
        Update_Triangulation=true;
}

template<class Solver>
void FlowEngine::imposeFlux ( Vector3r pos, Real flux,Solver& flow )
{
        typename Solver::RTriangulation& Tri = flow.T[flow.currentTes].Triangulation();
        double flux_base=0.f;
        double perm_base=0.f;
        typename Solver::Cell_handle cell = Tri.locate ( CGT::Point ( pos[0],pos[1],pos[2] ) );
        for ( int ngb=0;ngb<4;ngb++ ) {
                if ( !cell->neighbor ( ngb )->info().Pcondition ) {
                        flux_base += cell->info().k_norm() [ngb]* ( cell->neighbor ( ngb )->info().p() );
                        perm_base += cell->info().k_norm() [ngb];
                }
        }

        flow.imposedP.push_back ( pair<CGT::Point,Real> ( CGT::Point ( pos[0],pos[1],pos[2] ), ( flux_base-flux ) /perm_base ) );
        //force immediate update of boundary conditions
        Update_Triangulation=true;
}
template<class Solver>
void FlowEngine::clearImposedPressure ( Solver& flow ) { flow->imposedP.clear(); flow->IPCells.clear();}

template<class Solver>
Real FlowEngine::getFlux ( unsigned int cond,Solver& flow )
{
        if ( cond>=flow->imposedP.size() ) LOG_ERROR ( "Getting flux with cond higher than imposedP size." );
        double flux=0;
        Cell_handle& cell= flow->IPCells[cond];
        for ( int ngb=0;ngb<4;ngb++ ) {
                /*if (!cell->neighbor(ngb)->info().Pcondition)*/
                flux+= cell->info().k_norm() [ngb]* ( cell->info().p()-cell->neighbor ( ngb )->info().p() );
        }
        return flux+cell->info().dv();
}

template<class Solver>
void FlowEngine::Build_Triangulation ( Solver& flow )
{
        Build_Triangulation ( 0.f,flow );
}
template<class Solver>
void FlowEngine::Build_Triangulation ( double P_zero, Solver& flow )
{
        flow->ResetNetwork();
        if ( first ) flow->currentTes=0;
        else {
                flow->currentTes=!flow->currentTes;
                if ( Debug ) cout << "--------RETRIANGULATION-----------" << endl;
        }

        flow->Vtotalissimo=0; flow->Vsolid_tot=0; flow->Vporale=0; flow->Ssolid_tot=0;
        flow->SLIP_ON_LATERALS=slip_boundary;
        flow->k_factor = permeability_factor;
        flow->DEBUG_OUT = Debug;
        flow->useSolver = useSolver;
        flow->VISCOSITY = viscosity;
        flow->areaR2Permeability=areaR2Permeability;
	flow->TOLERANCE=Tolerance;
	flow->RELAX=Relax;
        flow->meanK_LIMIT = meanK_correction;
        flow->meanK_STAT = meanK_opt;
        flow->permeability_map = permeability_map;
	if ( fluidBulkModulus > 0 ) flow->compressible = 1;
	flow->fluidBulkModulus = fluidBulkModulus;
        flow->T[flow->currentTes].Clear();
        flow->T[flow->currentTes].max_id=-1;
        flow->x_min = 1000.0, flow->x_max = -10000.0, flow->y_min = 1000.0, flow->y_max = -10000.0, flow->z_min = 1000.0, flow->z_max = -10000.0;

        AddBoundary ( flow );
        Triangulate ( flow );
        if ( Debug ) cout << endl << "Tesselating------" << endl << endl;
        flow->T[flow->currentTes].Compute();

        flow->Define_fictious_cells();
        flow->DisplayStatistics ();
        flow->Compute_Permeability ( );

        porosity = flow->V_porale_porosity/flow->V_totale_porosity;

	if ( compute_K ) {BoundaryConditions ( flow ); K = flow->Sample_Permeability ( flow->x_min, flow->x_max, flow->y_min, flow->y_max, flow->z_min, flow->z_max );}

        BoundaryConditions ( flow );
        flow->Initialize_pressures ( P_zero );

        if ( !first && useSolver==0 ) flow->Interpolate ( flow->T[!flow->currentTes], flow->T[flow->currentTes] );
        if ( WaveAction ) flow->ApplySinusoidalPressure ( flow->T[flow->currentTes].Triangulation(), Sinus_Amplitude, Sinus_Average, 30 );

        Initialize_volumes ( flow );
        if ( viscousShear ) flow->ComputeEdgesSurfaces();
}
template<class Solver>
void FlowEngine::AddBoundary ( Solver& flow )
{
        shared_ptr<Sphere> sph ( new Sphere );
        int Sph_Index = sph->getClassIndexStatic();

        solver->x_min = 1000.0, solver->x_max = -10000.0, solver->y_min = 1000.0, solver->y_max = -10000.0, solver->z_min = 1000.0, solver->z_max = -10000.0;
        FOREACH ( const shared_ptr<Body>& b, *scene->bodies ) {
                if ( !b ) continue;
                if ( b->shape->getClassIndex() ==  Sph_Index ) {
                        Sphere* s=YADE_CAST<Sphere*> ( b->shape.get() );
                        //const Body::id_t& id = b->getId();
                        Real& rad = s->radius;
                        Real& x = b->state->pos[0];
                        Real& y = b->state->pos[1];
                        Real& z = b->state->pos[2];
                        flow->x_min = min ( flow->x_min, x-rad );
                        flow->x_max = max ( flow->x_max, x+rad );
                        flow->y_min = min ( flow->y_min, y-rad );
                        flow->y_max = max ( flow->y_max, y+rad );
                        flow->z_min = min ( flow->z_min, z-rad );
                        flow->z_max = max ( flow->z_max, z+rad );
                }
        }

        id_offset = flow->T[flow->currentTes].max_id+1;

        flow->id_offset = id_offset;

        flow->SectionArea = ( flow->x_max - flow->x_min ) * ( flow->z_max-flow->z_min );
        flow->Vtotale = ( flow->x_max-flow->x_min ) * ( flow->y_max-flow->y_min ) * ( flow->z_max-flow->z_min );

        flow->y_min_id=wallBottomId;
        flow->y_max_id=wallTopId;
        flow->x_max_id=wallRightId;
        flow->x_min_id=wallLeftId;
        flow->z_min_id=wallBackId;
        flow->z_max_id=wallFrontId;

        if ( flow->y_min_id>=0 ) flow->boundary ( flow->y_min_id ).useMaxMin = BOTTOM_Boundary_MaxMin;
        if ( flow->y_max_id>=0 ) flow->boundary ( flow->y_max_id ).useMaxMin = TOP_Boundary_MaxMin;
        if ( flow->x_max_id>=0 ) flow->boundary ( flow->x_max_id ).useMaxMin = RIGHT_Boundary_MaxMin;
        if ( flow->x_min_id>=0 ) flow->boundary ( flow->x_min_id ).useMaxMin = LEFT_Boundary_MaxMin;
        if ( flow->z_max_id>=0 ) flow->boundary ( flow->z_max_id ).useMaxMin = FRONT_Boundary_MaxMin;
        if ( flow->z_min_id>=0 ) flow->boundary ( flow->z_min_id ).useMaxMin = BACK_Boundary_MaxMin;

        //FIXME: Id's order in boundsIds is done according to the enumeration of boundaries from TXStressController.hpp, line 31. DON'T CHANGE IT!
        flow->boundsIds[0]= &flow->x_min_id;
        flow->boundsIds[1]= &flow->x_max_id;
        flow->boundsIds[2]= &flow->y_min_id;
        flow->boundsIds[3]= &flow->y_max_id;
        flow->boundsIds[4]= &flow->z_min_id;
        flow->boundsIds[5]= &flow->z_max_id;

        flow->Corner_min = CGT::Point ( flow->x_min, flow->y_min, flow->z_min );
        flow->Corner_max = CGT::Point ( flow->x_max, flow->y_max, flow->z_max );

        if ( Debug ) {
                cout << "Section area = " << flow->SectionArea << endl;
                cout << "Vtotale = " << flow->Vtotale << endl;
// 	cout << "Rmoy " << Rmoy << endl;
                cout << "x_min = " << flow->x_min << endl;
                cout << "x_max = " << flow->x_max << endl;
                cout << "y_max = " << flow->y_max << endl;
                cout << "y_min = " << flow->y_min << endl;
                cout << "z_min = " << flow->z_min << endl;
                cout << "z_max = " << flow->z_max << endl;
                cout << endl << "Adding Boundary------" << endl;
        }

        //assign BCs types and values
        BoundaryConditions ( flow );

        double center[3];

        for ( int i=0; i<6; i++ ) {
                if ( *flow->boundsIds[i]<0 ) continue;
                CGT::Vecteur Normal ( normal[i].x(), normal[i].y(), normal[i].z() );
                if ( flow->boundary ( *flow->boundsIds[i] ).useMaxMin ) flow->AddBoundingPlane ( true, Normal, *flow->boundsIds[i] );
                else {
                        const shared_ptr<Body>& wll = Body::byId ( *flow->boundsIds[i] , scene );
                        for ( int h=0;h<3;h++ ) {center[h] = wll->state->pos[h];}
                        flow->AddBoundingPlane ( center, wall_thickness, Normal,*flow->boundsIds[i] );
                }
        }
// 	flow->AddBoundingPlanes(true);
}

template<class Solver>
void FlowEngine::Triangulate ( Solver& flow )
{
///Using Tesselation wrapper (faster)
// 	TesselationWrapper TW;
// 	if (TW.Tes) delete TW.Tes;
// 	TW.Tes = &(flow->T[flow->currentTes]);//point to the current Tes we have in Flowengine
// 	TW.insertSceneSpheres();//TW is now really inserting in FlowEngine, using the faster insert(begin,end)
// 	TW.Tes = NULL;//otherwise, Tes would be deleted by ~TesselationWrapper() at the end of the function.
///Using one-by-one insertion
        shared_ptr<Sphere> sph ( new Sphere );
        int Sph_Index = sph->getClassIndexStatic();
        FOREACH ( const shared_ptr<Body>& b, *scene->bodies ) {
                if ( !b ) continue;
                if ( b->shape->getClassIndex() ==  Sph_Index ) {
                        Sphere* s=YADE_CAST<Sphere*> ( b->shape.get() );
                        const Body::id_t& id = b->getId();
                        Real rad = s->radius;
                        Real x = b->state->pos[0];
                        Real y = b->state->pos[1];
                        Real z = b->state->pos[2];
                        flow->T[flow->currentTes].insert ( x, y, z, rad, id );
                }
        }
        flow->viscousShearForces.resize ( flow->T[flow->currentTes].max_id+1 );
}
template<class Solver>
void FlowEngine::Initialize_volumes ( Solver& flow )
{
	typedef typename Solver::element_type Flow;
	typedef typename Flow::Finite_vertices_iterator Finite_vertices_iterator;
	typedef typename Solver::element_type Flow;
	typedef typename Flow::Finite_cells_iterator Finite_cells_iterator;
	
	Finite_vertices_iterator vertices_end = flow->T[flow->currentTes].Triangulation().finite_vertices_end();
	CGT::Vecteur Zero(0,0,0);
	for (Finite_vertices_iterator V_it = flow->T[flow->currentTes].Triangulation().finite_vertices_begin(); V_it!= vertices_end; V_it++) V_it->info().forces=Zero;
	
	Finite_cells_iterator cell_end = flow->T[flow->currentTes].Triangulation().finite_cells_end();
	for ( Finite_cells_iterator cell = flow->T[flow->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ )
	{
		switch ( cell->info().fictious() )
		{
			case ( 0 ) : cell->info().volume() = Volume_cell ( cell ); break;
			case ( 1 ) : cell->info().volume() = Volume_cell_single_fictious ( cell ); break;
			case ( 2 ) : cell->info().volume() = Volume_cell_double_fictious ( cell ); break;
			case ( 3 ) : cell->info().volume() = Volume_cell_triple_fictious ( cell ); break;
			default: break; 
		}

		if (flow->compressible) { cell->info().invVoidVolume() = (1 / ( cell->info().volume() - flow->volumeSolidPore(cell) )); }
	}
	if (Debug) cout << "Volumes initialised." << endl;
}

void FlowEngine::Average_real_cell_velocity()
{
        solver->Average_Relative_Cell_Velocity();
        Vector3r Vel ( 0,0,0 );
        //AVERAGE CELL VELOCITY
        Finite_cells_iterator cell_end = solver->T[solver->currentTes].Triangulation().finite_cells_end();
        for ( Finite_cells_iterator cell = solver->T[solver->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {
                for ( int g=0;g<4;g++ ) {
                        if ( !cell->vertex ( g )->info().isFictious ) {
                                const shared_ptr<Body>& sph = Body::byId ( cell->vertex ( g )->info().id(), scene );
                                for ( int i=0;i<3;i++ ) Vel[i] = Vel[i] + sph->state->vel[i]/4;
                        }
                }
                RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
                CGT::Point pos_av_facet;
                double volume_facet_translation = 0;
                CGT::Vecteur Vel_av ( Vel[0], Vel[1], Vel[2] );
                for ( int i=0; i<4; i++ ) {
                        volume_facet_translation = 0;
                        if ( !Tri.is_infinite ( cell->neighbor ( i ) ) ) {
                                CGT::Vecteur Surfk = cell->info()-cell->neighbor ( i )->info();
                                Real area = sqrt ( Surfk.squared_length() );
                                Surfk = Surfk/area;
                                CGT::Vecteur branch = cell->vertex ( facetVertices[i][0] )->point() - cell->info();
                                pos_av_facet = ( CGT::Point ) cell->info() + ( branch*Surfk ) *Surfk;
                                volume_facet_translation += Vel_av*cell->info().facetSurfaces[i];
                                cell->info().av_vel() = cell->info().av_vel() - volume_facet_translation/cell->info().volume() * ( pos_av_facet-CGAL::ORIGIN );
                        }
                }
        }
}
template<class Solver>
void FlowEngine::UpdateVolumes ( Solver& flow )
{
        if ( Debug ) cout << "Updating volumes.............." << endl;
        Real invDeltaT = 1/scene->dt;
        Finite_cells_iterator cell_end = flow->T[flow->currentTes].Triangulation().finite_cells_end();

        double newVol, dVol;
        eps_vol_max=0;
        Real totVol=0; Real totDVol=0; Real totVol0=0; Real totVol1=0; Real totVol2=0; Real totVol3=0;

        for ( Finite_cells_iterator cell = flow->T[flow->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {
                switch ( cell->info().fictious() ) {
                case ( 3 ) : newVol = Volume_cell_triple_fictious ( cell ); totVol3+=newVol; break;
                case ( 2 ) : newVol = Volume_cell_double_fictious ( cell ); totVol2+=newVol; break;
                case ( 1 ) : newVol = Volume_cell_single_fictious ( cell ); totVol1+=newVol; break;
                case ( 0 ) : newVol = Volume_cell ( cell ); totVol0+=newVol; break;
                default: newVol = 0; break;
                }
                totVol+=newVol;
                dVol=cell->info().volumeSign* ( newVol - cell->info().volume() );
                totDVol+=dVol;
                eps_vol_max = max ( eps_vol_max, abs ( dVol/newVol ) );
                cell->info().dv() = dVol*invDeltaT;
                cell->info().volume() = newVol;
        }
        if ( Debug ) cout << "Updated volumes, total =" <<totVol<<", dVol="<<totDVol<<endl;
}
template<class Cellhandle>
Real FlowEngine::Volume_cell_single_fictious ( Cellhandle cell )
{
        Vector3r V[3];
        int b=0;
        int w=0;
        cell->info().volumeSign=1;
        Real Wall_coordinate=0;

        for ( int y=0;y<4;y++ ) {
                if ( ! ( cell->vertex ( y )->info().isFictious ) ) {
                        const shared_ptr<Body>& sph = Body::byId ( cell->vertex ( y )->info().id(), scene );
                        V[w]=sph->state->pos;
                        w++;
                } else {
                        b = cell->vertex ( y )->info().id();
                        const shared_ptr<Body>& wll = Body::byId ( b , scene );
                        if ( !solver->boundary ( b ).useMaxMin ) Wall_coordinate = wll->state->pos[solver->boundary ( b ).coordinate]+ ( solver->boundary ( b ).normal[solver->boundary ( b ).coordinate] ) *wall_thickness/2;
                        else Wall_coordinate = solver->boundary ( b ).p[solver->boundary ( b ).coordinate];
                }
        }

        double v1[3], v2[3];

        for ( int g=0;g<3;g++ ) { v1[g]=V[0][g]-V[1][g]; v2[g]=V[0][g]-V[2][g];}

        Real Volume = 0.5* ( ( V[0]-V[1] ).cross ( V[0]-V[2] ) ) [solver->boundary ( b ).coordinate] * ( 0.33333333333* ( V[0][solver->boundary ( b ).coordinate]+ V[1][solver->boundary ( b ).coordinate]+ V[2][solver->boundary ( b ).coordinate] ) - Wall_coordinate );

        return abs ( Volume );
}
template<class Cellhandle>
Real FlowEngine::Volume_cell_double_fictious ( Cellhandle cell )
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
                        const shared_ptr<Body>& wll = Body::byId ( b[j] , scene );
                        if ( !solver->boundary ( b[j] ).useMaxMin ) Wall_coordinate[j] = wll->state->pos[coord[j]] + ( solver->boundary ( b[j] ).normal[coord[j]] ) *wall_thickness/2;
                        else Wall_coordinate[j] = solver->boundary ( b[j] ).p[coord[j]];
                        j++;
                } else if ( first_sph ) {
                        const shared_ptr<Body>& sph1 = Body::byId ( cell->vertex ( g )->info().id(), scene );
                        A=AS=/*AT=*/ ( sph1->state->pos );
                        first_sph=false;
                } else {
                        const shared_ptr<Body>& sph2 = Body::byId ( cell->vertex ( g )->info().id(), scene );
                        B=BS=/*BT=*/ ( sph2->state->pos );
                }
        }
        AS[coord[0]]=BS[coord[0]] = Wall_coordinate[0];

        //first pyramid with triangular base (A,B,BS)
        Real Vol1=0.5* ( ( A-BS ).cross ( B-BS ) ) [coord[1]]* ( 0.333333333* ( 2*B[coord[1]]+A[coord[1]] )-Wall_coordinate[1] );
        //second pyramid with triangular base (A,AS,BS)
        Real Vol2=0.5* ( ( AS-BS ).cross ( A-BS ) ) [coord[1]]* ( 0.333333333* ( B[coord[1]]+2*A[coord[1]] )-Wall_coordinate[1] );
        return abs ( Vol1+Vol2 );
}
template<class Cellhandle>
Real FlowEngine::Volume_cell_triple_fictious ( Cellhandle cell )
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
                        if ( !solver->boundary ( b[j] ).useMaxMin ) Wall_coordinate[j] = wll->state->pos[coord[j]] + ( solver->boundary ( b[j] ).normal[coord[j]] ) *wall_thickness/2;
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
template<class Cellhandle>
Real FlowEngine::Volume_cell ( Cellhandle cell )
{
        Real volume = CGT::Tetraedre ( makeCgPoint ( Body::byId ( cell->vertex ( 0 )->info().id(), scene )->state->pos ),
                                       makeCgPoint ( Body::byId ( cell->vertex ( 1 )->info().id(), scene )->state->pos ),
                                       makeCgPoint ( Body::byId ( cell->vertex ( 2 )->info().id(), scene )->state->pos ),
                                       makeCgPoint ( Body::byId ( cell->vertex ( 3 )->info().id(), scene )->state->pos ) )
                      .volume();

        if ( ! ( cell->info().volumeSign ) ) cell->info().volumeSign= ( volume>0 ) ?1:-1;
        return volume;
}
template<class Solver>
void FlowEngine::ApplyViscousForces ( Solver& flow )
{
//   flow->ComputeEdgesSurfaces(); //only done in buildTriangulation
        if ( Debug ) cout << "Application of viscous forces" << endl;
        if ( Debug ) cout << "Number of edges = " << flow.Edge_ids.size() << endl;
        for ( unsigned int k=0; k<flow.viscousShearForces.size(); k++ ) flow.viscousShearForces[k]=Vector3r::Zero();

	typedef typename Solver::Tesselation Tesselation;
        const Tesselation& Tes = flow.T[flow.currentTes];
        for ( int i=0; i< ( int ) flow.Edge_ids.size(); i++ ) {
                int hasFictious= Tes.vertex ( flow.Edge_ids[i].first )->info().isFictious +  Tes.vertex ( flow.Edge_ids[i].second )->info().isFictious;
                const shared_ptr<Body>& sph1 = Body::byId ( flow.Edge_ids[i].first, scene );
                const shared_ptr<Body>& sph2 = Body::byId ( flow.Edge_ids[i].second, scene );
                Sphere* s1=YADE_CAST<Sphere*> ( sph1->shape.get() );
                Sphere* s2=YADE_CAST<Sphere*> ( sph2->shape.get() );
                Vector3r x = sph1->state->pos - sph2->state->pos;
                Vector3r n = x / sqrt(makeCgVect(x).squared_length());
                Vector3r deltaV;
                if ( !hasFictious )
                        deltaV = (sph2->state->vel + sph2->state->angVel.cross(s2->radius * n)) - (sph1->state->vel+ sph1->state->angVel.cross(s1->radius * n));
                else {
                        if ( hasFictious==1 ) {//for the fictious sphere, use velocity of the boundary, not of the body
                                Vector3r v1 = ( Tes.vertex ( flow.Edge_ids[i].first )->info().isFictious ) ? flow.boundary ( flow.Edge_ids[i].first ).velocity:sph1->state->vel + sph1->state->angVel.cross(s1->radius * n);
                                Vector3r v2 = ( Tes.vertex ( flow.Edge_ids[i].second )->info().isFictious ) ? flow.boundary ( flow.Edge_ids[i].second ).velocity:sph2->state->vel + sph2->state->angVel.cross(s2->radius * n);
                                deltaV = v2-v1;
                        } else {//both fictious, ignore
                                deltaV = Vector3r::Zero();
                        }
                }
                deltaV = deltaV - ( flow.Edge_normal[i].dot ( deltaV ) ) *flow.Edge_normal[i];
                Vector3r visc_f = flow.ComputeViscousForce ( deltaV, i );
//                 if ( Debug ) cout << "la force visqueuse entre " << flow->Edge_ids[i].first << " et " << flow->Edge_ids[i].second << "est " << visc_f << endl;
///    //(1) directement sur le body Yade...
//     scene->forces.addForce(flow->Edge_ids[i].first,visc_f);
//     scene->forces.addForce(flow->Edge_ids[i].second,-visc_f);
///   //(2) ou dans CGAL? On a le choix (on pourrait même avoir info->viscousF pour faire la différence entre les deux types de forces... mais ça prend un peu plus de mémoire et de temps de calcul)
//     Tes.vertex(flow->Edge_ids[i].first)->info().forces=Tes.vertex(flow->Edge_ids[i].first)->info().forces+makeCgVect(visc_f);
//     Tes.vertex(flow->Edge_ids[i].second)->info().forces=Tes.vertex(flow->Edge_ids[i].second)->info().forces+makeCgVect(visc_f);
/// //(3) ou dans un vecteur séparé (rapide)
                flow.viscousShearForces[flow.Edge_ids[i].first]+=visc_f;
                flow.viscousShearForces[flow.Edge_ids[i].second]-=visc_f;
        }
        if(Debug) cout<<"number of viscousShearForce"<<flow.viscousShearForces.size()<<endl;
}

YADE_PLUGIN ( ( FlowEngine ) );

//______________________________________________________________

//___________________ PERIODIC VERSION _________________________
//______________________________________________________________

PeriodicFlowEngine::~PeriodicFlowEngine()
{
}

void PeriodicFlowEngine:: action()
{
        if ( !isActivated ) return;
        CGT::PeriodicCellInfo::gradP = makeCgVect ( gradP );
        CGT::PeriodicCellInfo::hSize[0] = makeCgVect ( scene->cell->hSize.col ( 0 ) );
        CGT::PeriodicCellInfo::hSize[1] = makeCgVect ( scene->cell->hSize.col ( 1 ) );
        CGT::PeriodicCellInfo::hSize[2] = makeCgVect ( scene->cell->hSize.col ( 2 ) );
        CGT::PeriodicCellInfo::deltaP=CGT::Vecteur (
                                              CGT::PeriodicCellInfo::hSize[0]*CGT::PeriodicCellInfo::gradP,
                                              CGT::PeriodicCellInfo::hSize[1]*CGT::PeriodicCellInfo::gradP,
                                              CGT::PeriodicCellInfo::hSize[2]*CGT::PeriodicCellInfo::gradP );
// 	timingDeltas->start();

        if ( first ) {Build_Triangulation ( P_zero ); Update_Triangulation = false;}
// 	timingDeltas->checkpoint("Triangulating");
        UpdateVolumes ( );
        Eps_Vol_Cumulative += eps_vol_max;
        if ( ( Eps_Vol_Cumulative > EpsVolPercent_RTRG || retriangulationLastIter>1000 ) && retriangulationLastIter>10 ) {
                Update_Triangulation = true;
                Eps_Vol_Cumulative=0;
                retriangulationLastIter=0;
                ReTrg++;
        } else  retriangulationLastIter++;
// 		timingDeltas->checkpoint("Update_Volumes");

	///Compute flow and and forces here
	solver->GaussSeidel(scene->dt);
	solver->ComputeFacetForcesWithCache();
// 	timingDeltas->checkpoint("Gauss-Seidel");
// 	timingDeltas->checkpoint("Compute_Forces");

        ///Application of vicscous forces
        scene->forces.sync();
        if ( viscousShear ) ApplyViscousForces(*solver);

        Finite_vertices_iterator vertices_end = solver->T[solver->currentTes].Triangulation().finite_vertices_end();
        for ( Finite_vertices_iterator V_it = solver->T[solver->currentTes].Triangulation().finite_vertices_begin(); V_it !=  vertices_end; V_it++ ) {
                if ( !viscousShear )
                        scene->forces.addForce ( V_it->info().id(), Vector3r ( ( V_it->info().forces ) [0],V_it->info().forces[1],V_it->info().forces[2] ) );
                else
                        if (V_it->info().isGhost) continue;
                        else
                                scene->forces.addForce ( V_it->info().id(), Vector3r ( ( V_it->info().forces ) [0],V_it->info().forces[1],V_it->info().forces[2] ) +solver->viscousShearForces[V_it->info().id() ] );
        }
        ///End Compute flow and forces
// 		timingDeltas->checkpoint("Applying Forces");

        if ( Update_Triangulation && !first ) {
                Build_Triangulation ( P_zero );
                Update_Triangulation = false;
        }
// 	if (velocity_profile) /*flow->FluidVelocityProfile();*/solver->Average_Fluid_Velocity();
        first=false;
// 	timingDeltas->checkpoint("Ending");
}


void PeriodicFlowEngine::Triangulate()
{
        shared_ptr<Sphere> sph ( new Sphere );
        int Sph_Index = sph->getClassIndexStatic();
        FOREACH ( const shared_ptr<Body>& b, *scene->bodies ) {
                if ( !b || b->shape->getClassIndex() != Sph_Index ) continue;
                Vector3r wpos;
                Sphere* s=YADE_CAST<Sphere*> ( b->shape.get() );
                const Body::id_t& id = b->getId();
                Real rad = s->radius;
                Vector3i period;
                // FIXME: use "sheared" variant if the cell is sheared
                wpos=scene->cell->wrapPt ( b->state->pos,period );
                Real x = wpos[0];
                Real y = wpos[1];
                Real z = wpos[2];
                Vertex_handle vh0=solver->T[solver->currentTes].insert ( x, y, z, rad, id );
                // FIXME: it can be out of the period!
                for ( int k=0;k<3;k++ ) vh0->info().period[k]=-period[k];
                const Vector3r& cellSize ( scene->cell->getSize() );
//                 wpos=scene->cell->unshearPt ( wpos );
                // traverse all periodic cells around the body, to see if any of them touches
                Vector3r halfSize= ( rad+duplicateThreshold ) *Vector3r ( 1,1,1 );
                Vector3r pmin,pmax;
                Vector3i i;
                for ( i[0]=-1; i[0]<=1; i[0]++ )
                        for ( i[1]=-1;i[1]<=1; i[1]++ )
                                for ( i[2]=-1; i[2]<=1; i[2]++ ) {
                                        if ( i[0]==0 && i[1]==0 && i[2]==0 ) continue; // middle; already rendered above
                                        Vector3r pos2=wpos+Vector3r ( cellSize[0]*i[0],cellSize[1]*i[1],cellSize[2]*i[2] ); // shift, but without shear!
                                        pmin=pos2-halfSize;
                                        pmax=pos2+halfSize;
                                        if ( pmin[0]<=cellSize[0] && pmax[0]>=0 && pmin[1]<=cellSize[1] && pmax[1]>=0 && pmin[2]<=cellSize[2] && pmax[2]>=0 ) {
                                                //with shear:
                                                //Vector3r pt=scene->cell->shearPt ( pos2 );
                                                //without shear:
                                                const Vector3r& pt= pos2;
                                                Vertex_handle vh=solver->T[solver->currentTes].insert ( pt[0],pt[1],pt[2],rad,id,false,id );
                                                for ( int k=0;k<3;k++ ) vh->info().period[k]=i[k]-period[k];}}
        }
        solver -> viscousShearForces.resize ( solver -> T[solver -> currentTes].max_id+1 );
}


Real PeriodicFlowEngine::Volume_cell ( Cell_handle cell )
{
        Real volume = CGT::Tetraedre (
		makeCgPoint ( Body::byId ( cell->vertex ( 0 )->info().id(), scene )->state->pos ) + cell->vertex ( 0 )->info().ghostShift(),
		makeCgPoint ( Body::byId ( cell->vertex ( 1 )->info().id(), scene )->state->pos ) + cell->vertex ( 1 )->info().ghostShift(),
		makeCgPoint ( Body::byId ( cell->vertex ( 2 )->info().id(), scene )->state->pos ) + cell->vertex ( 2 )->info().ghostShift(),
		makeCgPoint ( Body::byId ( cell->vertex ( 3 )->info().id(), scene )->state->pos ) + cell->vertex ( 3 )->info().ghostShift() )
		.volume();
        if ( ! ( cell->info().volumeSign ) ) cell->info().volumeSign= ( volume>0 ) ?1:-1;
        return volume;
}

Real PeriodicFlowEngine::Volume_cell_single_fictious ( Cell_handle cell )
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
                        if ( !solver->boundary ( b ).useMaxMin ) Wall_coordinate = wll->state->pos[solver->boundary ( b ).coordinate]+ ( solver->boundary ( b ).normal[solver->boundary ( b ).coordinate] ) *wall_thickness/2;
                        else Wall_coordinate = solver->boundary ( b ).p[solver->boundary ( b ).coordinate];
                }
        }
        double v1[3], v2[3];
        for ( int g=0;g<3;g++ ) { v1[g]=V[0][g]-V[1][g]; v2[g]=V[0][g]-V[2][g];}
        Real Volume = 0.5* ( ( V[0]-V[1] ).cross ( V[0]-V[2] ) ) [solver->boundary ( b ).coordinate] * ( 0.33333333333* ( V[0][solver->boundary ( b ).coordinate]+ V[1][solver->boundary ( b ).coordinate]+ V[2][solver->boundary ( b ).coordinate] ) - Wall_coordinate );
        return abs ( Volume );
}


void PeriodicFlowEngine::locateCell ( Cell_handle baseCell, unsigned int& index, unsigned int count)
{
        if (count>10) LOG_ERROR("More than 10 attempts to locate a cell, duplicateThreshold may be too small, resulting in periodicity inconsistencies.");
	PeriFlowTesselation::Cell_Info& base_info = baseCell->info();
        //already located, return FIXME: is inline working correctly? else move this test outside the function, just before the calls
	if ( base_info.index>0 || base_info.isGhost ) return;
	RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	Vector3r center ( 0,0,0 );
	Vector3i period;
	if (baseCell->info().fictious()==0)
		for ( int k=0;k<4;k++ ) center+= 0.25*makeVector3r (baseCell->vertex(k)->point());
	else {
		Real boundPos=0; int coord=0;
		for ( int k=0;k<4;k++ ) {
			if ( !baseCell->vertex ( k )->info().isFictious ) center+= 0.3333333333*makeVector3r ( baseCell->vertex ( k )->point() );
			else {
				coord=solver->boundary ( baseCell->vertex ( k )->info().id() ).coordinate;
				boundPos=solver->boundary ( baseCell->vertex ( k )->info().id() ).p[coord];}
		}
		center[coord]=boundPos;
	}
	Vector3r wdCenter= scene->cell->wrapPt ( center,period );
	if ( period[0]!=0 || period[1]!=0 || period[2]!=0 ) {
		if ( baseCell->info().index>0 ) {
			cout<<"indexed cell is found ghost!"<<base_info.index <<endl;
			base_info.isGhost=false;
			return;
		}
		Cell_handle ch= Tri.locate ( CGT::Point ( wdCenter[0],wdCenter[1],wdCenter[2] ) );//T[currentTes].vertexHandles[id]
		base_info.period[0]=period[0];
		base_info.period[1]=period[1];
		base_info.period[2]=period[2];
		//call recursively, since the returned cell could be also a ghost (especially if baseCell is a non-periodic type from the external contour
		locateCell ( ch,index,++count );
		if ( ch==baseCell ) cerr<<"WTF!!"<<endl;
		base_info.isGhost=true;
		base_info._pression=& ( ch->info().p() );
		base_info.index=ch->info().index;
		base_info.Pcondition=ch->info().Pcondition;
	} else {
		base_info.isGhost=false;
		//index is 1-based, if it is zero it is not initialized, we define it here
		if ( !base_info.Pcondition && base_info.index==0 ) base_info.index=++index;
	}
}

Vector3r PeriodicFlowEngine::meanVelocity()
{
        solver->Average_Relative_Cell_Velocity();
        Vector3r meanVel ( 0,0,0 );
        Real volume=0;
        Finite_cells_iterator cell_end = solver->T[solver->currentTes].Triangulation().finite_cells_end();
        for ( Finite_cells_iterator cell = solver->T[solver->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {
                if ( cell->info().isGhost ) continue;
                for ( int i=0;i<3;i++ )
                        meanVel[i]=meanVel[i]+ ( ( cell->info().av_vel() ) [i] * abs ( cell->info().volume() ) );
                volume+=abs ( cell->info().volume() );
        }
        return ( meanVel/volume );
}

void PeriodicFlowEngine::UpdateVolumes ()
{
        if ( Debug ) cout << "Updating volumes.............." << endl;
        Real invDeltaT = 1/scene->dt;
        Finite_cells_iterator cell_end = solver->T[solver->currentTes].Triangulation().finite_cells_end();

        double newVol, dVol;
        eps_vol_max=0;
        Real totVol=0;
        Real totDVol=0;
        Real totVol0=0;
        Real totVol1=0;

	for ( Finite_cells_iterator cell = solver->T[solver->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {
                if ( cell->info().isGhost ) continue;
                switch ( cell->info().fictious() ) {
                case ( 1 ) :
                        newVol = Volume_cell_single_fictious ( cell );
                        totVol1+=newVol;
                        break;
                case ( 0 ) :
                        newVol = Volume_cell ( cell );
                        totVol0+=newVol;
                        break;
                default:
                        newVol = 0;
                        break;
                }
                totVol+=newVol;
                dVol=cell->info().volumeSign * ( newVol - cell->info().volume() );
                totDVol+=dVol;
                eps_vol_max = max ( eps_vol_max, abs ( dVol/newVol ) );
                cell->info().dv() = dVol * invDeltaT;
                cell->info().volume() = newVol;
        }
        if ( Debug ) cout << "Updated volumes, total =" <<totVol<<", dVol="<<totDVol<<" "<< totVol0<<" "<< totVol1<<endl;
}

void PeriodicFlowEngine::Initialize_volumes ()
{
        Finite_vertices_iterator vertices_end = solver->T[solver->currentTes].Triangulation().finite_vertices_end();
        CGT::Vecteur Zero ( 0,0,0 );
        for ( Finite_vertices_iterator V_it = solver->T[solver->currentTes].Triangulation().finite_vertices_begin(); V_it!= vertices_end; V_it++ ) V_it->info().forces=Zero;

        Finite_cells_iterator cell_end = solver->T[solver->currentTes].Triangulation().finite_cells_end();
        for ( Finite_cells_iterator cell = solver->T[solver->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {
                if ( cell->info().isGhost ) continue;
		switch ( cell->info().fictious() )
		{
			case ( 0 ) : cell->info().volume() = Volume_cell ( cell ); break;
			case ( 1 ) : cell->info().volume() = Volume_cell_single_fictious ( cell ); break;
// 			case ( 2 ) : cell->info().volume() = Volume_cell_double_fictious ( cell ); break;
// 			case ( 3 ) : cell->info().volume() = Volume_cell_triple_fictious ( cell ); break;

			if (solver->compressible) { cell->info().invVoidVolume() = 1 / ( Volume_cell ( cell ) - solver->volumeSolidPore(cell) ); }

			default:  cell->info().volume() = 0; break;
		}
        }
        if ( Debug ) cout << "Volumes initialised." << endl;
}

void PeriodicFlowEngine::Build_Triangulation ( double P_zero )
{
        solver->ResetNetwork();
        if ( first ) solver->currentTes=0;
        else {
                solver->currentTes=!solver->currentTes;
                if ( Debug ) cout << "--------RETRIANGULATION-----------" << endl;
        }

        solver->Vtotalissimo=0; solver->Vsolid_tot=0; solver->Vporale=0; solver->Ssolid_tot=0;
        solver->SLIP_ON_LATERALS=slip_boundary;
        solver->k_factor = permeability_factor;
        solver->DEBUG_OUT = Debug;
        solver->useSolver = useSolver;
        solver->VISCOSITY = viscosity;
        solver->areaR2Permeability=areaR2Permeability;
	if ( fluidBulkModulus > 0 ) solver->compressible = 1;
	solver->fluidBulkModulus = fluidBulkModulus;
        solver->T[solver->currentTes].Clear();
        solver->T[solver->currentTes].max_id=-1;
        AddBoundary ( solver );
        if ( Debug ) cout << endl << "Added boundaries------" << endl << endl;
        Triangulate ();
        if ( Debug ) cout << endl << "Tesselating------" << endl << endl;
        solver->T[solver->currentTes].Compute();
        solver->Define_fictious_cells();
        solver->meanK_LIMIT = meanK_correction;
        solver->meanK_STAT = meanK_opt;
        solver->permeability_map = permeability_map;
        solver->Compute_Permeability ( );
        porosity = solver->V_porale_porosity/solver->V_totale_porosity;
        solver->TOLERANCE=Tolerance;solver->RELAX=Relax;
        Finite_cells_iterator cell_end = solver->T[solver->currentTes].Triangulation().finite_cells_end();
        for ( Finite_cells_iterator cell = solver->T[solver->currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ) {cell->info().dv() = 0; cell->info().p() = P_zero;}
        BoundaryConditions ( solver );

        solver->Initialize_pressures ( P_zero );
        // Define the ghost cells and add indexes to the cells inside the period (the ones that will contain the pressure unknowns)
        //This must be done after boundary conditions and before initialize pressure, else the indexes are not good (not accounting imposedP):
        Finite_cells_iterator cellend=solver->T[solver->currentTes].Triangulation().finite_cells_end();
        unsigned int index=0;
        for ( Finite_cells_iterator cell=solver -> T[solver -> currentTes].Triangulation().finite_cells_begin(); cell!=cellend; cell++ )
                locateCell ( cell,index );
        solver->DisplayStatistics ();
        //FIXME: check interpolate() for the periodic case, at least use the mean pressure from previous step.
// 	if (!first && useSolver==0) solver->Interpolate (solver->T[!solver->currentTes], solver->T[solver->currentTes]);

        if ( WaveAction ) solver->ApplySinusoidalPressure ( solver->T[solver->currentTes].Triangulation(), Sinus_Amplitude, Sinus_Average, 30 );
        Initialize_volumes();
        if ( viscousShear ) solver->ComputeEdgesSurfaces();
}

YADE_PLUGIN ( ( PeriodicFlowEngine ) );

#endif //FLOW_ENGINE

#endif /* YADE_CGAL */


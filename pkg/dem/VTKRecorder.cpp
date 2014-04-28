#ifdef YADE_VTK

#include"VTKRecorder.hpp"

#include<vtkCellArray.h>
#include<vtkPoints.h>
#include<vtkPointData.h>
#include<vtkCellData.h>
#include<vtkSmartPointer.h>
#include<vtkFloatArray.h>
#include<vtkUnstructuredGrid.h>
#include<vtkPolyData.h>
#include<vtkXMLUnstructuredGridWriter.h>
#include<vtkXMLPolyDataWriter.h>
#include<vtkZLibDataCompressor.h>
#include<vtkTriangle.h>
#include<vtkLine.h>
#include<vtkQuad.h>
#include<vtkHexahedron.h>
#ifdef YADE_VTK_MULTIBLOCK
  #include<vtkXMLMultiBlockDataWriter.h>
  #include<vtkMultiBlockDataSet.h>
#endif

#include<yade/core/Scene.hpp>
#include<yade/pkg/common/Sphere.hpp>
#include<yade/pkg/common/Facet.hpp>
#include<yade/pkg/common/Box.hpp>
#include<yade/pkg/dem/ConcretePM.hpp>
#include<yade/pkg/dem/WirePM.hpp>
#include<yade/pkg/dem/JointedCohesiveFrictionalPM.hpp>
#include<yade/pkg/dem/Shop.hpp>
#ifdef YADE_LIQCONTROL
	#include<yade/pkg/dem/ViscoelasticCapillarPM.hpp>
#endif

YADE_PLUGIN((VTKRecorder));
CREATE_LOGGER(VTKRecorder);

void VTKRecorder::action(){
	vector<bool> recActive(REC_SENTINEL,false);
	FOREACH(string& rec, recorders){
		if(rec=="all"){
			recActive[REC_SPHERES]=true;
			recActive[REC_VELOCITY]=true;
			recActive[REC_FACETS]=true;
			recActive[REC_BOXES]=true;
			recActive[REC_COLORS]=true;
			recActive[REC_MASS]=true;
			recActive[REC_INTR]=true;
			recActive[REC_ID]=true;
			recActive[REC_MASK]=true;
			recActive[REC_CLUMPID]=true;
			recActive[REC_MATERIALID]=true;
			recActive[REC_STRESS]=true;
			if (scene->isPeriodic) { recActive[REC_PERICELL]=true; }
		}
		else if(rec=="spheres") recActive[REC_SPHERES]=true;
		else if(rec=="velocity") recActive[REC_VELOCITY]=true;
		else if(rec=="facets") recActive[REC_FACETS]=true;
		else if(rec=="boxes") recActive[REC_BOXES]=true;
		else if(rec=="mass") recActive[REC_MASS]=true;
		else if((rec=="colors") || (rec=="color"))recActive[REC_COLORS]=true;
		else if(rec=="cpm") recActive[REC_CPM]=true;
		else if(rec=="wpm") recActive[REC_WPM]=true;
		else if(rec=="intr") recActive[REC_INTR]=true;
		else if((rec=="ids") || (rec=="id")) recActive[REC_ID]=true;
		else if(rec=="mask") recActive[REC_MASK]=true;
		else if((rec=="clumpids") || (rec=="clumpId")) recActive[REC_CLUMPID]=true;
		else if(rec=="materialId") recActive[REC_MATERIALID]=true;
		else if(rec=="stress") recActive[REC_STRESS]=true;
		else if(rec=="jcfpm") recActive[REC_JCFPM]=true;
		else if(rec=="cracks") recActive[REC_CRACKS]=true;
		else if(rec=="pericell" && scene->isPeriodic) recActive[REC_PERICELL]=true;
		else if(rec=="liquidcontrol") recActive[REC_LIQ]=true;
		else LOG_ERROR("Unknown recorder named `"<<rec<<"' (supported are: all, spheres, velocity, facets, boxes, color, stress, cpm, wpm, intr, id, clumpId, materialId, jcfpm, cracks, pericell). Ignored.");
	}
	// cpm needs interactions
	if(recActive[REC_CPM]) recActive[REC_INTR]=true;
	
	// jcfpm needs interactions
	if(recActive[REC_JCFPM]) recActive[REC_INTR]=true;

	// wpm needs interactions
	if(recActive[REC_WPM]) recActive[REC_INTR]=true;

	// liquid control needs interactions
	if(recActive[REC_LIQ]) recActive[REC_INTR]=true;


	// spheres
	vtkSmartPointer<vtkPoints> spheresPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> spheresCells = vtkSmartPointer<vtkCellArray>::New();
	
	vtkSmartPointer<vtkFloatArray> radii = vtkSmartPointer<vtkFloatArray>::New();
	radii->SetNumberOfComponents(1);
	radii->SetName("radii");
	
	vtkSmartPointer<vtkFloatArray> spheresMass = vtkSmartPointer<vtkFloatArray>::New();
	spheresMass->SetNumberOfComponents(1);
	spheresMass->SetName("mass");
	
	vtkSmartPointer<vtkFloatArray> spheresId = vtkSmartPointer<vtkFloatArray>::New();
	spheresId->SetNumberOfComponents(1);
	spheresId->SetName("id");

#ifdef YADE_SPH
	vtkSmartPointer<vtkFloatArray> spheresCsSPH = vtkSmartPointer<vtkFloatArray>::New();
	spheresCsSPH->SetNumberOfComponents(1);
	spheresCsSPH->SetName("SPH_Cs");
	
	vtkSmartPointer<vtkFloatArray> spheresRhoSPH = vtkSmartPointer<vtkFloatArray>::New();
	spheresRhoSPH->SetNumberOfComponents(1);
	spheresRhoSPH->SetName("SPH_Rho");
	
	vtkSmartPointer<vtkFloatArray> spheresPressSPH = vtkSmartPointer<vtkFloatArray>::New();
	spheresPressSPH->SetNumberOfComponents(1);
	spheresPressSPH->SetName("SPH_Press");
	
	vtkSmartPointer<vtkFloatArray> spheresCoordNumbSPH = vtkSmartPointer<vtkFloatArray>::New();
	spheresCoordNumbSPH->SetNumberOfComponents(1);
	spheresCoordNumbSPH->SetName("SPH_Neigh");
#endif

#ifdef YADE_LIQCONTROL
	vtkSmartPointer<vtkFloatArray> spheresLiqVol = vtkSmartPointer<vtkFloatArray>::New();
	spheresLiqVol->SetNumberOfComponents(1);
	spheresLiqVol->SetName("Liq_Vol");
	
	vtkSmartPointer<vtkFloatArray> spheresLiqVolIter = vtkSmartPointer<vtkFloatArray>::New();
	spheresLiqVolIter->SetNumberOfComponents(1);
	spheresLiqVolIter->SetName("Liq_VolIter");
	
	vtkSmartPointer<vtkFloatArray> spheresLiqVolTotal = vtkSmartPointer<vtkFloatArray>::New();
	spheresLiqVolTotal->SetNumberOfComponents(1);
	spheresLiqVolTotal->SetName("Liq_VolTotal");
#endif

	vtkSmartPointer<vtkFloatArray> spheresMask = vtkSmartPointer<vtkFloatArray>::New();
	spheresMask->SetNumberOfComponents(1);
	spheresMask->SetName("mask");
	
	vtkSmartPointer<vtkFloatArray> clumpId = vtkSmartPointer<vtkFloatArray>::New();
	clumpId->SetNumberOfComponents(1);
	clumpId->SetName("clumpId");
	
	vtkSmartPointer<vtkFloatArray> spheresColors = vtkSmartPointer<vtkFloatArray>::New();
	spheresColors->SetNumberOfComponents(3);
	spheresColors->SetName("color");
	
	vtkSmartPointer<vtkFloatArray> spheresLinVelVec = vtkSmartPointer<vtkFloatArray>::New();
	spheresLinVelVec->SetNumberOfComponents(3);
	spheresLinVelVec->SetName("linVelVec");		//Linear velocity in Vector3 form
	
	vtkSmartPointer<vtkFloatArray> spheresLinVelLen = vtkSmartPointer<vtkFloatArray>::New();
	spheresLinVelLen->SetNumberOfComponents(1);
	spheresLinVelLen->SetName("linVelLen");		//Length (magnitude) of linear velocity
	
	vtkSmartPointer<vtkFloatArray> spheresAngVelVec = vtkSmartPointer<vtkFloatArray>::New();
	spheresAngVelVec->SetNumberOfComponents(3);
	spheresAngVelVec->SetName("angVelVec");		//Angular velocity in Vector3 form
	
	vtkSmartPointer<vtkFloatArray> spheresAngVelLen = vtkSmartPointer<vtkFloatArray>::New();
	spheresAngVelLen->SetNumberOfComponents(1);
	spheresAngVelLen->SetName("angVelLen");		//Length (magnitude) of angular velocity
	
	vtkSmartPointer<vtkFloatArray> spheresNormalStressVec = vtkSmartPointer<vtkFloatArray>::New();
	spheresNormalStressVec->SetNumberOfComponents(3);
	spheresNormalStressVec->SetName("normalStress");
	
	vtkSmartPointer<vtkFloatArray> spheresShearStressVec = vtkSmartPointer<vtkFloatArray>::New();
	spheresShearStressVec->SetNumberOfComponents(3);
	spheresShearStressVec->SetName("shearStress");
	
	vtkSmartPointer<vtkFloatArray> spheresNormalStressNorm = vtkSmartPointer<vtkFloatArray>::New();
	spheresNormalStressNorm->SetNumberOfComponents(1);
	spheresNormalStressNorm->SetName("normalStressNorm");
	
	vtkSmartPointer<vtkFloatArray> spheresMaterialId = vtkSmartPointer<vtkFloatArray>::New();
	spheresMaterialId->SetNumberOfComponents(1);
	spheresMaterialId->SetName("materialId");

	// facets
	vtkSmartPointer<vtkPoints> facetsPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> facetsCells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkFloatArray> facetsColors = vtkSmartPointer<vtkFloatArray>::New();
	facetsColors->SetNumberOfComponents(3);
	facetsColors->SetName("color");
	
	vtkSmartPointer<vtkFloatArray> facetsForceVec = vtkSmartPointer<vtkFloatArray>::New();
	facetsForceVec->SetNumberOfComponents(3);
	facetsForceVec->SetName("stressVec");
	
	vtkSmartPointer<vtkFloatArray> facetsForceLen = vtkSmartPointer<vtkFloatArray>::New();
	facetsForceLen->SetNumberOfComponents(1);
	facetsForceLen->SetName("stressLen");
	
	vtkSmartPointer<vtkFloatArray> facetsMaterialId = vtkSmartPointer<vtkFloatArray>::New();
	facetsMaterialId->SetNumberOfComponents(1);
	facetsMaterialId->SetName("materialId");
	
	vtkSmartPointer<vtkFloatArray> facetsMask = vtkSmartPointer<vtkFloatArray>::New();
	facetsMask->SetNumberOfComponents(1);
	facetsMask->SetName("mask");

	// boxes
	vtkSmartPointer<vtkPoints> boxesPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> boxesCells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkFloatArray> boxesColors = vtkSmartPointer<vtkFloatArray>::New();
	boxesColors->SetNumberOfComponents(3);
	boxesColors->SetName("color");
	
	vtkSmartPointer<vtkFloatArray> boxesForceVec = vtkSmartPointer<vtkFloatArray>::New();
	boxesForceVec->SetNumberOfComponents(3);
	boxesForceVec->SetName("stressVec");
	
	vtkSmartPointer<vtkFloatArray> boxesForceLen = vtkSmartPointer<vtkFloatArray>::New();
	boxesForceLen->SetNumberOfComponents(1);
	boxesForceLen->SetName("stressLen");
	
	vtkSmartPointer<vtkFloatArray> boxesMaterialId = vtkSmartPointer<vtkFloatArray>::New();
	boxesMaterialId->SetNumberOfComponents(1);
	boxesMaterialId->SetName("materialId");
	
	vtkSmartPointer<vtkFloatArray> boxesMask = vtkSmartPointer<vtkFloatArray>::New();
	boxesMask->SetNumberOfComponents(1);
	boxesMask->SetName("mask");

	// interactions
	vtkSmartPointer<vtkPoints> intrBodyPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> intrCells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkFloatArray> intrForceN = vtkSmartPointer<vtkFloatArray>::New();
	intrForceN->SetNumberOfComponents(1);
	intrForceN->SetName("forceN");
	vtkSmartPointer<vtkFloatArray> intrAbsForceT = vtkSmartPointer<vtkFloatArray>::New();
	intrAbsForceT->SetNumberOfComponents(3);
	intrAbsForceT->SetName("absForceT");

	// pericell
	vtkSmartPointer<vtkPoints> pericellPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> pericellHexa = vtkSmartPointer<vtkCellArray>::New();

	// extras for CPM
	if(recActive[REC_CPM]){ CpmStateUpdater csu; csu.update(scene); }
	vtkSmartPointer<vtkFloatArray> cpmDamage = vtkSmartPointer<vtkFloatArray>::New();
	cpmDamage->SetNumberOfComponents(1);
	cpmDamage->SetName("cpmDamage");
	vtkSmartPointer<vtkFloatArray> cpmStress = vtkSmartPointer<vtkFloatArray>::New();
	cpmStress->SetNumberOfComponents(9);
	cpmStress->SetName("cpmStress");

	// extras for JCFpm
	vtkSmartPointer<vtkFloatArray> damage = vtkSmartPointer<vtkFloatArray>::New();
	damage->SetNumberOfComponents(1);;
	damage->SetName("damage");
	vtkSmartPointer<vtkFloatArray> damageRel = vtkSmartPointer<vtkFloatArray>::New();
	damageRel->SetNumberOfComponents(1);;
	damageRel->SetName("damageRel");
	vtkSmartPointer<vtkFloatArray> intrIsCohesive = vtkSmartPointer<vtkFloatArray>::New();
	intrIsCohesive->SetNumberOfComponents(1);
	intrIsCohesive->SetName("isCohesive");
	vtkSmartPointer<vtkFloatArray> intrIsOnJoint = vtkSmartPointer<vtkFloatArray>::New();
	intrIsOnJoint->SetNumberOfComponents(1);
	intrIsOnJoint->SetName("isOnJoint");
	
	// extras for cracks
	vtkSmartPointer<vtkPoints> crackPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> crackCells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkFloatArray> iter = vtkSmartPointer<vtkFloatArray>::New();
	iter->SetNumberOfComponents(1);
	iter->SetName("iter");
	vtkSmartPointer<vtkFloatArray> crackType = vtkSmartPointer<vtkFloatArray>::New();
	crackType->SetNumberOfComponents(1);
	crackType->SetName("crackType");
	vtkSmartPointer<vtkFloatArray> crackSize = vtkSmartPointer<vtkFloatArray>::New();
	crackSize->SetNumberOfComponents(1);
	crackSize->SetName("crackSize");
	vtkSmartPointer<vtkFloatArray> crackNorm = vtkSmartPointer<vtkFloatArray>::New();
	crackNorm->SetNumberOfComponents(3);
	crackNorm->SetName("crackNorm");
	
#ifdef YADE_LIQCONTROL
	vtkSmartPointer<vtkFloatArray> liqVol = vtkSmartPointer<vtkFloatArray>::New();
	liqVol->SetNumberOfComponents(1);
	liqVol->SetName("liqVol");
	
	vtkSmartPointer<vtkFloatArray> liqVolNorm = vtkSmartPointer<vtkFloatArray>::New();
	liqVolNorm->SetNumberOfComponents(1);
	liqVolNorm->SetName("liqVolNorm");
#endif
	// the same for newly created cracks
// 	vtkSmartPointer<vtkPoints> crackPosNew = vtkSmartPointer<vtkPoints>::New();
// 	vtkSmartPointer<vtkCellArray> crackCellsNew = vtkSmartPointer<vtkCellArray>::New();
// 	vtkSmartPointer<vtkFloatArray> iterNew = vtkSmartPointer<vtkFloatArray>::New();
// 	iterNew->SetNumberOfComponents(1);
// 	iterNew->SetName("iter");
// 	vtkSmartPointer<vtkFloatArray> crackTypeNew = vtkSmartPointer<vtkFloatArray>::New();
// 	crackTypeNew->SetNumberOfComponents(1);
// 	crackTypeNew->SetName("crackType");
// 	vtkSmartPointer<vtkFloatArray> crackSizeNew = vtkSmartPointer<vtkFloatArray>::New();
// 	crackSizeNew->SetNumberOfComponents(1);
// 	crackSizeNew->SetName("crackSize");
// 	vtkSmartPointer<vtkFloatArray> crackNormNew = vtkSmartPointer<vtkFloatArray>::New();
// 	crackNormNew->SetNumberOfComponents(3);
// 	crackNormNew->SetName("crackNorm");
	
	// extras for WireMatPM
	vtkSmartPointer<vtkFloatArray> wpmNormalForce = vtkSmartPointer<vtkFloatArray>::New();
	wpmNormalForce->SetNumberOfComponents(1);
	wpmNormalForce->SetName("wpmNormalForce");
	vtkSmartPointer<vtkFloatArray> wpmLimitFactor = vtkSmartPointer<vtkFloatArray>::New();
	wpmLimitFactor->SetNumberOfComponents(1);
	wpmLimitFactor->SetName("wpmLimitFactor");

	if(recActive[REC_INTR]){
		// holds information about cell distance between spatial and displayed position of each particle
		vector<Vector3i> wrapCellDist; if (scene->isPeriodic){ wrapCellDist.resize(scene->bodies->size()); }
		// save body positions, referenced by ids by vtkLine
		FOREACH(const shared_ptr<Body>& b, *scene->bodies){
			if (!b) {
				/* must keep ids contiguous, so that position in the array corresponds to Body::id */
				intrBodyPos->InsertNextPoint(NaN,NaN,NaN);
				continue;
			}
			if(!scene->isPeriodic){ intrBodyPos->InsertNextPoint(b->state->pos[0],b->state->pos[1],b->state->pos[2]); }
			else {
				Vector3r pos=scene->cell->wrapShearedPt(b->state->pos,wrapCellDist[b->id]);
				intrBodyPos->InsertNextPoint(pos[0],pos[1],pos[2]);
			}
			assert(intrBodyPos->GetNumberOfPoints()==b->id+1);
		}
		FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
			if(!I->isReal()) continue;
			if(skipFacetIntr){
				if(!(Body::byId(I->getId1()))) continue;
				if(!(Body::byId(I->getId2()))) continue;
				if(!(dynamic_cast<Sphere*>(Body::byId(I->getId1())->shape.get()))) continue;
				if(!(dynamic_cast<Sphere*>(Body::byId(I->getId2())->shape.get()))) continue;
			}
			/* For the periodic boundary conditions,
				find out whether the interaction crosses the boundary of the periodic cell;
				if it does, display the interaction on both sides of the cell, with one of the
				points sticking out in each case.
				Since vtkLines must connect points with an ID assigned, we will create a new additional
				point for each point outside the cell. It might create some data redundancy, but
				let us suppose that the number of interactions crossing the cell boundary is low compared
				to total numer of interactions
			*/
			// how many times to add values defined on interactions, depending on how many times the interaction is saved
			int numAddValues=1;
			// aperiodic boundary, or interaction is inside the cell
			if(!scene->isPeriodic || (scene->isPeriodic && (I->cellDist==wrapCellDist[I->getId2()]-wrapCellDist[I->getId1()]))){
				vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0,I->getId1());
				line->GetPointIds()->SetId(1,I->getId2());
				intrCells->InsertNextCell(line);
			} else {
				assert(scene->isPeriodic);
				// spatial positions of particles
				const Vector3r& p01(Body::byId(I->getId1())->state->pos); const Vector3r& p02(Body::byId(I->getId2())->state->pos);
				// create two line objects; each of them has one endpoint inside the cell and the other one sticks outside
				// A,B are the "fake" bodies outside the cell for id1 and id2 respectively, p1,p2 are the displayed points
				// distance in cell units for shifting A away from p1; negated value is shift of B away from p2
				Vector3r ptA(p01+scene->cell->hSize*(wrapCellDist[I->getId2()]-I->cellDist).cast<Real>());
				Vector3r ptB(p02+scene->cell->hSize*(wrapCellDist[I->getId1()]-I->cellDist).cast<Real>());
				vtkIdType idPtA=intrBodyPos->InsertNextPoint(ptA[0],ptA[1],ptA[2]), idPtB=intrBodyPos->InsertNextPoint(ptB[0],ptB[1],ptB[2]);
				vtkSmartPointer<vtkLine> line1B(vtkSmartPointer<vtkLine>::New()); line1B->GetPointIds()->SetId(0,I->getId1()); line1B->GetPointIds()->SetId(1,idPtB);
				vtkSmartPointer<vtkLine> lineA2(vtkSmartPointer<vtkLine>::New()); lineA2->GetPointIds()->SetId(0,idPtA); lineA2->GetPointIds()->SetId(1,I->getId2());
				numAddValues=2;
			}
			const NormShearPhys* phys = YADE_CAST<NormShearPhys*>(I->phys.get());
			const GenericSpheresContact* geom = YADE_CAST<GenericSpheresContact*>(I->geom.get());
			// gives _signed_ scalar of normal force, following the convention used in the respective constitutive law
			float fn=phys->normalForce.dot(geom->normal); 
			float fs[3]={ (float) abs(phys->shearForce[0]), (float) abs(phys->shearForce[1]), (float) abs(phys->shearForce[2])};
			// add the value once for each interaction object that we created (might be 2 for the periodic boundary)
			for(int i=0; i<numAddValues; i++){
				intrAbsForceT->InsertNextTupleValue(fs);
				if(recActive[REC_WPM]) {
					const WirePhys* wirephys = dynamic_cast<WirePhys*>(I->phys.get());
					if (wirephys!=NULL && wirephys->isLinked) {
						wpmLimitFactor->InsertNextValue(wirephys->limitFactor);
						wpmNormalForce->InsertNextValue(fn);
						intrForceN->InsertNextValue(NaN);
					}
					else {
						intrForceN->InsertNextValue(fn);
						wpmNormalForce->InsertNextValue(NaN);
						wpmLimitFactor->InsertNextValue(NaN);
					}
				}
				else if (recActive[REC_JCFPM]){
					const JCFpmPhys* jcfpmphys = YADE_CAST<JCFpmPhys*>(I->phys.get());
					intrIsCohesive->InsertNextValue(jcfpmphys->isCohesive);
					intrIsOnJoint->InsertNextValue(jcfpmphys->isOnJoint);
					intrForceN->InsertNextValue(fn);
				} else {
					intrForceN->InsertNextValue(fn);
				}
#ifdef YADE_LIQCONTROL
				if (recActive[REC_LIQ]) {
					const ViscElCapPhys* capphys = YADE_CAST<ViscElCapPhys*>(I->phys.get());
					liqVol->InsertNextValue(capphys->Vb);
					liqVolNorm->InsertNextValue(capphys->Vb/capphys->Vmax);
				}
#endif
			}
		}
	}

	//Additional Vector for storing forces
	vector<Shop::bodyState> bodyStates;
	if(recActive[REC_STRESS]) Shop::getStressForEachBody(bodyStates);
	
	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
		if (!b) continue;
		if(mask!=0 && (b->groupMask & mask)==0) continue;
		if (recActive[REC_SPHERES]){
			const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get()); 
			if (sphere){
				if(skipNondynamic && b->state->blockedDOFs==State::DOF_ALL) continue;
				vtkIdType pid[1];
				Vector3r pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				pid[0] = spheresPos->InsertNextPoint(pos[0], pos[1], pos[2]);
				spheresCells->InsertNextCell(1,pid);
				radii->InsertNextValue(sphere->radius);
				if (recActive[REC_ID]) spheresId->InsertNextValue(b->getId()); 
				if (recActive[REC_MASK]) spheresMask->InsertNextValue(b->groupMask);
				if (recActive[REC_MASS]) spheresMass->InsertNextValue(b->state->mass);
				if (recActive[REC_CLUMPID]) clumpId->InsertNextValue(b->clumpId);
				if (recActive[REC_COLORS]){
					const Vector3r& color = sphere->color;
					float c[3] = { (float) color[0], (float) color[1], (float) color[2]};
					spheresColors->InsertNextTupleValue(c);
				}
				if(recActive[REC_VELOCITY]){
					const Vector3r& vel = b->state->vel;
					float v[3] = { (float) vel[0], (float) vel[1], (float) vel[2] };
					spheresLinVelVec->InsertNextTupleValue(v);
					spheresLinVelLen->InsertNextValue(vel.norm());
					
					const Vector3r& angVel = b->state->angVel;
					float av[3] = { (float) angVel[0], (float) angVel[1], (float) angVel[2] };
					spheresAngVelVec->InsertNextTupleValue(av);
					spheresAngVelLen->InsertNextValue(angVel.norm());
				}
				if(recActive[REC_STRESS]){
					const Vector3r& stress = bodyStates[b->getId()].normStress;
					const Vector3r& shear = bodyStates[b->getId()].shearStress;
					float n[3] = { (float)  stress[0], (float) stress[1], (float) stress[2] };
					float s[3] = { (float)  shear [0], (float) shear [1], (float) shear [2] };
					spheresNormalStressVec->InsertNextTupleValue(n);
					spheresShearStressVec->InsertNextTupleValue(s);
					spheresNormalStressNorm->InsertNextValue(stress.norm());
				}
				
				if (recActive[REC_CPM]){
					cpmDamage->InsertNextValue(YADE_PTR_CAST<CpmState>(b->state)->normDmg);
					const Matrix3r& ss=YADE_PTR_CAST<CpmState>(b->state)->stress;
					//float s[3]={ss[0],ss[1],ss[2]};
					float s[9]={ (float) ss(0,0), (float) ss(0,1), (float) ss(0,2), (float) ss(1,0), (float) ss(1,1), (float) ss(1,2), (float) ss(2,0), (float) ss(2,1), (float) ss(2,2)};
					cpmStress->InsertNextTupleValue(s);
				}
				
				if (recActive[REC_JCFPM]){
					damage->InsertNextValue(YADE_PTR_CAST<JCFpmState>(b->state)->tensBreak + YADE_PTR_CAST<JCFpmState>(b->state)->shearBreak);
					damageRel->InsertNextValue(YADE_PTR_CAST<JCFpmState>(b->state)->tensBreakRel + YADE_PTR_CAST<JCFpmState>(b->state)->shearBreakRel);
				}
#ifdef YADE_SPH
				spheresCsSPH->InsertNextValue(b->Cs); 
				spheresRhoSPH->InsertNextValue(b->rho); 
				spheresPressSPH->InsertNextValue(b->press); 
				spheresCoordNumbSPH->InsertNextValue(b->coordNumber()); 
#endif
#ifdef YADE_LIQCONTROL
				spheresLiqVol->InsertNextValue(b->Vf);
				const Real tmpVolIter = liqVolIterBody(b);
				spheresLiqVolIter->InsertNextValue(tmpVolIter);
				spheresLiqVolTotal->InsertNextValue(tmpVolIter + b->Vf);
#endif
				if (recActive[REC_MATERIALID]) spheresMaterialId->InsertNextValue(b->material->id);
				continue;
			}
		}
		if (recActive[REC_FACETS]){
			const Facet* facet = dynamic_cast<Facet*>(b->shape.get()); 
			if (facet){
				Vector3r pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				const vector<Vector3r>& localPos = facet->vertices;
				Matrix3r facetAxisT=b->state->ori.toRotationMatrix();
				vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();
				vtkIdType nbPoints=facetsPos->GetNumberOfPoints();
				for (int i=0;i<3;++i){
					Vector3r globalPos = pos + facetAxisT * localPos[i];
					facetsPos->InsertNextPoint(globalPos[0], globalPos[1], globalPos[2]);
					tri->GetPointIds()->SetId(i,nbPoints+i);
				}
				facetsCells->InsertNextCell(tri);
				if (recActive[REC_COLORS]){
					const Vector3r& color = facet->color;
					float c[3] = { (float) color[0], (float) color[1], (float) color[2]};
					facetsColors->InsertNextTupleValue(c);
				}
				if(recActive[REC_STRESS]){
					const Vector3r& stress = bodyStates[b->getId()].normStress+bodyStates[b->getId()].shearStress;
					float s[3] = { (float) stress[0], (float) stress[1], (float) stress[2] };
					facetsForceVec->InsertNextTupleValue(s);
					facetsForceLen->InsertNextValue(stress.norm());
				}
				if (recActive[REC_MATERIALID]) facetsMaterialId->InsertNextValue(b->material->id);
				if (recActive[REC_MASK]) facetsMask->InsertNextValue(b->groupMask);
				continue;
			}
		}
		if (recActive[REC_BOXES]){
			const Box* box = dynamic_cast<Box*>(b->shape.get()); 
			if (box){
				Vector3r pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				Vector3r ext(box->extents);
				vtkSmartPointer<vtkQuad> boxes = vtkSmartPointer<vtkQuad>::New();
				
				Vector3r A = Vector3r(pos[0]-ext[0], pos[1]-ext[1], pos[2]-ext[2]);
				Vector3r B = Vector3r(pos[0]-ext[0], pos[1]+ext[1], pos[2]-ext[2]);
				Vector3r C = Vector3r(pos[0]+ext[0], pos[1]+ext[1], pos[2]-ext[2]);
				Vector3r D = Vector3r(pos[0]+ext[0], pos[1]-ext[1], pos[2]-ext[2]);
				
				Vector3r E = Vector3r(pos[0]-ext[0], pos[1]-ext[1], pos[2]+ext[2]);
				Vector3r F = Vector3r(pos[0]-ext[0], pos[1]+ext[1], pos[2]+ext[2]);
				Vector3r G = Vector3r(pos[0]+ext[0], pos[1]+ext[1], pos[2]+ext[2]);
				Vector3r H = Vector3r(pos[0]+ext[0], pos[1]-ext[1], pos[2]+ext[2]);
				
				addWallVTK(boxes, boxesPos, A, B, C, D);
				boxesCells->InsertNextCell(boxes);
				
				addWallVTK(boxes, boxesPos, E, H, G, F);
				boxesCells->InsertNextCell(boxes);
				
				addWallVTK(boxes, boxesPos, A, E, F, B);
				boxesCells->InsertNextCell(boxes);
				
				addWallVTK(boxes, boxesPos, G, H, D, C);
				boxesCells->InsertNextCell(boxes);
				
				addWallVTK(boxes, boxesPos, F, G, C, B);
				boxesCells->InsertNextCell(boxes);
				
				addWallVTK(boxes, boxesPos, D, H, E, A);
				boxesCells->InsertNextCell(boxes);
				
				for(int i=0; i<6; i++){
					if (recActive[REC_COLORS]){
						const Vector3r& color = box->color;
						float c[3] = { (float) color[0], (float) color[1], (float) color[2]};
						boxesColors->InsertNextTupleValue(c);
					}
					if(recActive[REC_STRESS]){
						const Vector3r& stress = bodyStates[b->getId()].normStress+bodyStates[b->getId()].shearStress;
						float s[3] = { (float) stress[0], (float) stress[1], (float) stress[2] };
						boxesForceVec->InsertNextTupleValue(s);
						boxesForceLen->InsertNextValue(stress.norm());
					}
					if (recActive[REC_MATERIALID]) boxesMaterialId->InsertNextValue(b->material->id);
					if (recActive[REC_MASK]) boxesMask->InsertNextValue(b->groupMask);
				}
				continue;
			}
		}
	}

	if (recActive[REC_PERICELL]) {
		const Matrix3r& hSize = scene->cell->hSize;
		Vector3r v0 = hSize*Vector3r(0,0,1);
		Vector3r v1 = hSize*Vector3r(0,1,1);
		Vector3r v2 = hSize*Vector3r(1,1,1);
		Vector3r v3 = hSize*Vector3r(1,0,1);
		Vector3r v4 = hSize*Vector3r(0,0,0);
		Vector3r v5 = hSize*Vector3r(0,1,0);
		Vector3r v6 = hSize*Vector3r(1,1,0);
		Vector3r v7 = hSize*Vector3r(1,0,0);
		pericellPoints->InsertNextPoint(v0[0],v0[1],v0[2]);
		pericellPoints->InsertNextPoint(v1[0],v1[1],v1[2]);
		pericellPoints->InsertNextPoint(v2[0],v2[1],v2[2]);
		pericellPoints->InsertNextPoint(v3[0],v3[1],v3[2]);
		pericellPoints->InsertNextPoint(v4[0],v4[1],v4[2]);
		pericellPoints->InsertNextPoint(v5[0],v5[1],v5[2]);
		pericellPoints->InsertNextPoint(v6[0],v6[1],v6[2]);
		pericellPoints->InsertNextPoint(v7[0],v7[1],v7[2]);
		vtkSmartPointer<vtkHexahedron> h = vtkSmartPointer<vtkHexahedron>::New();
		vtkIdList* l = h->GetPointIds();
		for (int i=0; i<8; i++) {
			l->SetId(i,i);
		}
		pericellHexa->InsertNextCell(h);
	}

	
	vtkSmartPointer<vtkDataCompressor> compressor;
	if(compress) compressor=vtkSmartPointer<vtkZLibDataCompressor>::New();

	vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
	if (recActive[REC_SPHERES]){
		spheresUg->SetPoints(spheresPos);
		spheresUg->SetCells(VTK_VERTEX, spheresCells);
		spheresUg->GetPointData()->AddArray(radii);
		if (recActive[REC_ID]) spheresUg->GetPointData()->AddArray(spheresId);
		if (recActive[REC_MASK]) spheresUg->GetPointData()->AddArray(spheresMask);
		if (recActive[REC_MASS]) spheresUg->GetPointData()->AddArray(spheresMass);
		if (recActive[REC_CLUMPID]) spheresUg->GetPointData()->AddArray(clumpId);
		if (recActive[REC_COLORS]) spheresUg->GetPointData()->AddArray(spheresColors);
		if (recActive[REC_VELOCITY]){
			spheresUg->GetPointData()->AddArray(spheresLinVelVec);
			spheresUg->GetPointData()->AddArray(spheresAngVelVec);
			spheresUg->GetPointData()->AddArray(spheresLinVelLen);
			spheresUg->GetPointData()->AddArray(spheresAngVelLen);
		}
#ifdef YADE_SPH
		spheresUg->GetPointData()->AddArray(spheresCsSPH);
		spheresUg->GetPointData()->AddArray(spheresRhoSPH);
		spheresUg->GetPointData()->AddArray(spheresPressSPH);
		spheresUg->GetPointData()->AddArray(spheresCoordNumbSPH);
#endif
#ifdef YADE_LIQCONTROL
		spheresUg->GetPointData()->AddArray(spheresLiqVol);
		spheresUg->GetPointData()->AddArray(spheresLiqVolIter);
		spheresUg->GetPointData()->AddArray(spheresLiqVolTotal);
#endif
		if (recActive[REC_STRESS]){
			spheresUg->GetPointData()->AddArray(spheresNormalStressVec);
			spheresUg->GetPointData()->AddArray(spheresShearStressVec);
			spheresUg->GetPointData()->AddArray(spheresNormalStressNorm);
		}
		if (recActive[REC_CPM]){
			spheresUg->GetPointData()->AddArray(cpmDamage);
			spheresUg->GetPointData()->AddArray(cpmStress);
		}

		if (recActive[REC_JCFPM]) {
			spheresUg->GetPointData()->AddArray(damage);
		}

		if (recActive[REC_MATERIALID]) spheresUg->GetPointData()->AddArray(spheresMaterialId);

		#ifdef YADE_VTK_MULTIBLOCK
		if(!multiblock)
		#endif
			{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if(compress) writer->SetCompressor(compressor);
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+"spheres."+lexical_cast<string>(scene->iter)+".vtu";
			writer->SetFileName(fn.c_str());
			#ifdef YADE_VTK6
				writer->SetInputData(spheresUg);
			#else
				writer->SetInput(spheresUg);
			#endif
			writer->Write();
		}
	}
	vtkSmartPointer<vtkUnstructuredGrid> facetsUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
	if (recActive[REC_FACETS]){
		facetsUg->SetPoints(facetsPos);
		facetsUg->SetCells(VTK_TRIANGLE, facetsCells);
		if (recActive[REC_COLORS]) facetsUg->GetCellData()->AddArray(facetsColors);
		if (recActive[REC_STRESS]){
			facetsUg->GetCellData()->AddArray(facetsForceVec);
			facetsUg->GetCellData()->AddArray(facetsForceLen);
		}
		if (recActive[REC_MATERIALID]) facetsUg->GetCellData()->AddArray(facetsMaterialId);
		if (recActive[REC_MASK]) facetsUg->GetCellData()->AddArray(facetsMask);
		#ifdef YADE_VTK_MULTIBLOCK
			if(!multiblock)
		#endif
			{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if(compress) writer->SetCompressor(compressor);
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+"facets."+lexical_cast<string>(scene->iter)+".vtu";
			writer->SetFileName(fn.c_str());
			#ifdef YADE_VTK6
				writer->SetInputData(facetsUg);
			#else
				writer->SetInput(facetsUg);
			#endif
			writer->Write();	
		}
	}
	vtkSmartPointer<vtkUnstructuredGrid> boxesUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
	if (recActive[REC_BOXES]){
		boxesUg->SetPoints(boxesPos);
		boxesUg->SetCells(VTK_QUAD, boxesCells);
		if (recActive[REC_COLORS]) boxesUg->GetCellData()->AddArray(boxesColors);
		if (recActive[REC_STRESS]){
			boxesUg->GetCellData()->AddArray(boxesForceVec);
			boxesUg->GetCellData()->AddArray(boxesForceLen);
		}
		if (recActive[REC_MATERIALID]) boxesUg->GetCellData()->AddArray(boxesMaterialId);
		if (recActive[REC_MASK]) boxesUg->GetCellData()->AddArray(boxesMask);
		#ifdef YADE_VTK_MULTIBLOCK
			if(!multiblock)
		#endif
			{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if(compress) writer->SetCompressor(compressor);
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+"boxes."+lexical_cast<string>(scene->iter)+".vtu";
			writer->SetFileName(fn.c_str());
			#ifdef YADE_VTK6
				writer->SetInputData(boxesUg);
			#else
				writer->SetInput(boxesUg);
			#endif
			writer->Write();	
		}
	}
	vtkSmartPointer<vtkPolyData> intrPd = vtkSmartPointer<vtkPolyData>::New();
	if (recActive[REC_INTR]){
		intrPd->SetPoints(intrBodyPos);
		intrPd->SetLines(intrCells);
		intrPd->GetCellData()->AddArray(intrForceN);
		intrPd->GetCellData()->AddArray(intrAbsForceT);
#ifdef YADE_LIQCONTROL
		if (recActive[REC_LIQ]) { 
			intrPd->GetCellData()->AddArray(liqVol);
			intrPd->GetCellData()->AddArray(liqVolNorm);
		}
#endif
		if (recActive[REC_JCFPM]) { 
			intrPd->GetCellData()->AddArray(intrIsCohesive);
			intrPd->GetCellData()->AddArray(intrIsOnJoint);
		}
		if (recActive[REC_WPM]){
			intrPd->GetCellData()->AddArray(wpmNormalForce);
			intrPd->GetCellData()->AddArray(wpmLimitFactor);
		}
		#ifdef YADE_VTK_MULTIBLOCK
			if(!multiblock)
		#endif
			{
			vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
			if(compress) writer->SetCompressor(compressor);
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+"intrs."+lexical_cast<string>(scene->iter)+".vtp";
			writer->SetFileName(fn.c_str());
			#ifdef YADE_VTK6
				writer->SetInputData(intrPd);
			#else
				writer->SetInput(intrPd);
			#endif
			writer->Write();
		}
	}
	vtkSmartPointer<vtkUnstructuredGrid> pericellUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
	if (recActive[REC_PERICELL]){
		pericellUg->SetPoints(pericellPoints);
		pericellUg->SetCells(12,pericellHexa);
		#ifdef YADE_VTK_MULTIBLOCK
			if(!multiblock)
		#endif
			{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if(compress) writer->SetCompressor(compressor);
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+"pericell."+lexical_cast<string>(scene->iter)+".vtu";
			writer->SetFileName(fn.c_str());
			#ifdef YADE_VTK6
				writer->SetInputData(pericellUg);
			#else
				writer->SetInput(pericellUg);
			#endif
			writer->Write();
		}
	}

	if (recActive[REC_CRACKS]) {
		string fileCracks = "cracks_"+Key+".txt";
		std::ifstream file (fileCracks.c_str(),std::ios::in);
		vtkSmartPointer<vtkUnstructuredGrid> crackUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkSmartPointer<vtkUnstructuredGrid> crackUgNew = vtkSmartPointer<vtkUnstructuredGrid>::New();
		
		 if(file){
			 while ( !file.eof() ){
				std::string line;
				Real i,p0,p1,p2,t,s,n0,n1,n2;
				while ( std::getline(file, line)/* writes into string "line", a line of file "file". To go along diff. lines*/ ) 
				{
					file >> i >> p0 >> p1 >> p2 >> t >> s >> n0 >> n1 >> n2;
					vtkIdType pid[1];
					pid[0] = crackPos->InsertNextPoint(p0, p1, p2);
					crackCells->InsertNextCell(1,pid);
					crackType->InsertNextValue(t);
					crackSize->InsertNextValue(s);
					iter->InsertNextValue(i);
					float n[3] = { n0, n1, n2 };
					crackNorm->InsertNextTupleValue(n);
					// Then, taking care only of newly created cracks :
// 					if (i > scene->iter - iterPeriod)
// 					{
// 					  pid[0] = crackPosNew->InsertNextPoint(p0, p1, p2);
// 					  crackCellsNew->InsertNextCell(1,pid);
// 					  crackTypeNew->InsertNextValue(t);
// 					  crackSizeNew->InsertNextValue(s);
// 					  iterNew->InsertNextValue(i);
// 					  crackNormNew->InsertNextTupleValue(n);
// 					}  
				}
			 }
			 file.close();
		} 

		crackUg->SetPoints(crackPos);
		crackUg->SetCells(VTK_VERTEX, crackCells);
		crackUg->GetPointData()->AddArray(iter);
		crackUg->GetPointData()->AddArray(crackType);
		crackUg->GetPointData()->AddArray(crackSize);
		crackUg->GetPointData()->AddArray(crackNorm); //orientation of 2D glyphs does not match this direction (some work to do in order to have the good orientation) 
		
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		if(compress) writer->SetCompressor(compressor);
		if(ascii) writer->SetDataModeToAscii();
		string fn=fileName+"cracks."+lexical_cast<string>(scene->iter)+".vtu";
		writer->SetFileName(fn.c_str());
		#ifdef YADE_VTK6
			writer->SetInputData(crackUg);
		#else
			writer->SetInput(crackUg);
		#endif
		writer->Write();
		
		// Same operations, for newly created cracks :
// 		crackUgNew->SetPoints(crackPosNew);
// 		crackUgNew->SetCells(VTK_VERTEX, crackCellsNew);
// 		crackUgNew->GetPointData()->AddArray(iterNew);
// 		crackUgNew->GetPointData()->AddArray(crackTypeNew);
// 		crackUgNew->GetPointData()->AddArray(crackSizeNew);
// 		crackUgNew->GetPointData()->AddArray(crackNormNew); //same remark about the orientation...
// 	
// 		fn=fileName+"newcracks."+lexical_cast<string>(scene->iter)+".vtu";
// 		writer->SetFileName(fn.c_str());
// 		#ifdef YADE_VTK6
// 			writer->SetInputData(crackUgNew);
// 		#else
// 			writer->SetInput(crackUgNew);
// 		#endif
// 		writer->Write();
	}

	#ifdef YADE_VTK_MULTIBLOCK
		if(multiblock){
			vtkSmartPointer<vtkMultiBlockDataSet> multiblockDataset = vtkSmartPointer<vtkMultiBlockDataSet>::New();
			int i=0;
			if(recActive[REC_SPHERES]) multiblockDataset->SetBlock(i++,spheresUg);
			if(recActive[REC_FACETS]) multiblockDataset->SetBlock(i++,facetsUg);
			if(recActive[REC_INTR]) multiblockDataset->SetBlock(i++,intrPd);
			if(recActive[REC_PERICELL]) multiblockDataset->SetBlock(i++,pericellUg);
			vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+lexical_cast<string>(scene->iter)+".vtm";
			writer->SetFileName(fn.c_str());
			#ifdef YADE_VTK6
				writer->SetInputData(multiblockDataset);
			#else
				writer->SetInput(multiblockDataset);
			#endif
			writer->Write();	
		}
	#endif
};

void VTKRecorder::addWallVTK (vtkSmartPointer<vtkQuad>& boxes, vtkSmartPointer<vtkPoints>& boxesPos, Vector3r& W1, Vector3r& W2, Vector3r& W3, Vector3r& W4) {
	//Function for exporting walls of boxes
	vtkIdType nbPoints=boxesPos->GetNumberOfPoints();
	
	boxesPos->InsertNextPoint(W1[0], W1[1], W1[2]);
	boxes->GetPointIds()->SetId(0,nbPoints+0);
	
	boxesPos->InsertNextPoint(W2[0], W2[1], W2[2]);
	boxes->GetPointIds()->SetId(1,nbPoints+1);
	
	boxesPos->InsertNextPoint(W3[0], W3[1], W3[2]);
	boxes->GetPointIds()->SetId(2,nbPoints+2);
	
	boxesPos->InsertNextPoint(W4[0], W4[1], W4[2]);
	boxes->GetPointIds()->SetId(3,nbPoints+3);
};

#endif /* YADE_VTK */

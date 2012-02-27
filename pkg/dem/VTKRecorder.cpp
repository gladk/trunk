#include"VTKRecorder.hpp"

// avoid warnings in VTK code (does not work?!)
#pragma GCC diagnostic ignored "-Wdeprecated"
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
	#ifdef YADE_VTK_MULTIBLOCK
		#include<vtkXMLMultiBlockDataWriter.h>
		#include<vtkMultiBlockDataSet.h>
	#endif
#pragma GCC diagnostic warning "-Wdeprecated"

#include<yade/core/Scene.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Facet.hpp>
#include<yade/pkg-dem/ConcretePM.hpp>
#include<yade/pkg-dem/RockPM.hpp>
#include<yade/pkg-dem/Shop.hpp>


YADE_PLUGIN((VTKRecorder));
YADE_REQUIRE_FEATURE(VTK)
CREATE_LOGGER(VTKRecorder);

void VTKRecorder::action(){
	vector<bool> recActive(REC_SENTINEL,false);
	FOREACH(string& rec, recorders){
		if(rec=="all"){
			recActive[REC_SPHERES]=true;
			recActive[REC_VELOCITY]=true;
			recActive[REC_FACETS]=true;
			recActive[REC_COLORS]=true;
			recActive[REC_INTR]=true;
			recActive[REC_ID]=true;
			recActive[REC_MASK]=true;
			recActive[REC_CLUMPID]=true;
			recActive[REC_MATERIALID]=true;
			recActive[REC_STRESS]=true;
		}
		else if(rec=="spheres") recActive[REC_SPHERES]=true;
		else if(rec=="velocity") recActive[REC_VELOCITY]=true;
		else if(rec=="facets") recActive[REC_FACETS]=true;
		else if((rec=="colors") || (rec=="color"))recActive[REC_COLORS]=true;
		else if(rec=="cpm") recActive[REC_CPM]=true;
		else if(rec=="rpm") recActive[REC_RPM]=true;
		else if(rec=="intr") recActive[REC_INTR]=true;
		else if((rec=="ids") || (rec=="id")) recActive[REC_ID]=true;
		else if(rec=="mask") recActive[REC_MASK]=true;
		else if((rec=="clumpids") || (rec=="clumpId")) recActive[REC_CLUMPID]=true;
		else if(rec=="materialId") recActive[REC_MATERIALID]=true;
		else if(rec=="stress") recActive[REC_STRESS]=true;
		else LOG_ERROR("Unknown recorder named `"<<rec<<"' (supported are: all, spheres, velocity, facets, color, stress, cpm, rpm, intr, id, clumpId, materialId). Ignored.");
	}
	// cpm needs interactions
	if(recActive[REC_CPM]) recActive[REC_INTR]=true;

	// spheres
	vtkSmartPointer<vtkPoints> spheresPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> spheresCells = vtkSmartPointer<vtkCellArray>::New();
	
	vtkSmartPointer<vtkFloatArray> radii = vtkSmartPointer<vtkFloatArray>::New();
	radii->SetNumberOfComponents(1);
	radii->SetName("radii");
	
	vtkSmartPointer<vtkFloatArray> spheresId = vtkSmartPointer<vtkFloatArray>::New();
	spheresId->SetNumberOfComponents(1);
	spheresId->SetName("id");
	
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

	// interactions
	vtkSmartPointer<vtkPoints> intrBodyPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> intrCells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkFloatArray> intrForceN = vtkSmartPointer<vtkFloatArray>::New();
	intrForceN->SetNumberOfComponents(1);
	intrForceN->SetName("forceN");
	vtkSmartPointer<vtkFloatArray> intrAbsForceT = vtkSmartPointer<vtkFloatArray>::New();
	intrAbsForceT->SetNumberOfComponents(3);
	intrAbsForceT->SetName("absForceT");

	// extras for CPM
	if(recActive[REC_CPM]){ CpmStateUpdater csu; csu.update(scene); }
	vtkSmartPointer<vtkFloatArray> cpmDamage = vtkSmartPointer<vtkFloatArray>::New();
	cpmDamage->SetNumberOfComponents(1);
	cpmDamage->SetName("cpmDamage");
	vtkSmartPointer<vtkFloatArray> cpmSigma = vtkSmartPointer<vtkFloatArray>::New();
	cpmSigma->SetNumberOfComponents(3);
	cpmSigma->SetName("cpmSigma");
	vtkSmartPointer<vtkFloatArray> cpmSigmaM = vtkSmartPointer<vtkFloatArray>::New();
	cpmSigmaM->SetNumberOfComponents(1);
	cpmSigmaM->SetName("cpmSigmaM");
	vtkSmartPointer<vtkFloatArray> cpmTau = vtkSmartPointer<vtkFloatArray>::New();
	cpmTau->SetNumberOfComponents(3);
	cpmTau->SetName("cpmTau");
	
	// extras for RPM
	vtkSmartPointer<vtkFloatArray> rpmSpecNum = vtkSmartPointer<vtkFloatArray>::New();
	rpmSpecNum->SetNumberOfComponents(1);
	rpmSpecNum->SetName("rpmSpecNum");
	vtkSmartPointer<vtkFloatArray> rpmSpecMass = vtkSmartPointer<vtkFloatArray>::New();
	rpmSpecMass->SetNumberOfComponents(1);
	rpmSpecMass->SetName("rpmSpecMass");
	vtkSmartPointer<vtkFloatArray> rpmSpecDiam = vtkSmartPointer<vtkFloatArray>::New();
	rpmSpecDiam->SetNumberOfComponents(1);
	rpmSpecDiam->SetName("rpmSpecDiam");

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
				Vector3r ptA(p01+scene->cell->Hsize*(wrapCellDist[I->getId2()]-I->cellDist).cast<Real>());
				Vector3r ptB(p02+scene->cell->Hsize*(wrapCellDist[I->getId1()]-I->cellDist).cast<Real>());
				vtkIdType idPtA=intrBodyPos->InsertNextPoint(ptA[0],ptA[1],ptA[2]), idPtB=intrBodyPos->InsertNextPoint(ptB[0],ptB[1],ptB[2]);
				vtkSmartPointer<vtkLine> line1B(vtkSmartPointer<vtkLine>::New()); line1B->GetPointIds()->SetId(0,I->getId1()); line1B->GetPointIds()->SetId(1,idPtB);
				vtkSmartPointer<vtkLine> lineA2(vtkSmartPointer<vtkLine>::New()); lineA2->GetPointIds()->SetId(0,idPtA); line1B->GetPointIds()->SetId(1,I->getId2());
				numAddValues=2;
			}
			const NormShearPhys* phys = YADE_CAST<NormShearPhys*>(I->phys.get());
			const GenericSpheresContact* geom = YADE_CAST<GenericSpheresContact*>(I->geom.get());
			// gives _signed_ scalar of normal force, following the convention used in the respective constitutive law
			float fn=phys->normalForce.dot(geom->normal); 
			float fs[3]={abs(phys->shearForce[0]),abs(phys->shearForce[1]),abs(phys->shearForce[2])};
			// add the value once for each interaction object that we created (might be 2 for the periodic boundary)
			for(int i=0; i<numAddValues; i++){
				intrForceN->InsertNextValue(fn);
				intrAbsForceT->InsertNextTupleValue(fs);
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
				if(skipNondynamic && !b->isDynamic()) continue;
				vtkIdType pid[1];
				Vector3r pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				pid[0] = spheresPos->InsertNextPoint(pos[0], pos[1], pos[2]);
				spheresCells->InsertNextCell(1,pid);
				radii->InsertNextValue(sphere->radius);
				if (recActive[REC_ID]) spheresId->InsertNextValue(b->getId()); 
				if (recActive[REC_MASK]) spheresMask->InsertNextValue(b->groupMask);
				if (recActive[REC_CLUMPID]) clumpId->InsertNextValue(b->clumpId);
				if (recActive[REC_COLORS]){
					const Vector3r& color = sphere->color;
					float c[3] = {color[0],color[1],color[2]};
					spheresColors->InsertNextTupleValue(c);
				}
				if(recActive[REC_VELOCITY]){
					const Vector3r& vel = b->state->vel;
					float v[3] = { vel[0],vel[1],vel[2] };
					spheresLinVelVec->InsertNextTupleValue(v);
					spheresLinVelLen->InsertNextValue(vel.norm());
					
					const Vector3r& angVel = b->state->angVel;
					float av[3] = { angVel[0],angVel[1],angVel[2] };
					spheresAngVelVec->InsertNextTupleValue(av);
					spheresAngVelLen->InsertNextValue(angVel.norm());
				}
				if(recActive[REC_STRESS]){
					const Vector3r& stress = bodyStates[b->getId()].normStress;
					const Vector3r& shear = bodyStates[b->getId()].shearStress;
					float n[3] = { stress[0],stress[1],stress[2] };
					float s[3] = { shear [0],shear [1],shear [2] };
					spheresNormalStressVec->InsertNextTupleValue(n);
					spheresShearStressVec->InsertNextTupleValue(s);
					spheresNormalStressNorm->InsertNextValue(stress.norm());
				}
				
				if (recActive[REC_CPM]){
					cpmDamage->InsertNextValue(YADE_PTR_CAST<CpmState>(b->state)->normDmg);
					const Vector3r& ss=YADE_PTR_CAST<CpmState>(b->state)->sigma;
					const Vector3r& tt=YADE_PTR_CAST<CpmState>(b->state)->tau;
					float s[3]={ss[0],ss[1],ss[2]};
					float t[3]={tt[0],tt[1],tt[2]};
					cpmSigma->InsertNextTupleValue(s);
					cpmSigmaM->InsertNextValue((ss[0]+ss[1]+ss[2])/3.);
					cpmTau->InsertNextTupleValue(t);
				}
				if (recActive[REC_RPM]){
					rpmSpecNum->InsertNextValue(YADE_PTR_CAST<RpmState>(b->state)->specimenNumber);
					rpmSpecMass->InsertNextValue(YADE_PTR_CAST<RpmState>(b->state)->specimenMass);
					rpmSpecDiam->InsertNextValue(YADE_PTR_CAST<RpmState>(b->state)->specimenMaxDiam);
				}
				
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
					float c[3] = {color[0],color[1],color[2]};
					facetsColors->InsertNextTupleValue(c);
				}
				if(recActive[REC_STRESS]){
					const Vector3r& stress = bodyStates[b->getId()].normStress+bodyStates[b->getId()].shearStress;
					float s[3] = { stress[0],stress[1],stress[2] };
					facetsForceVec->InsertNextTupleValue(s);
					facetsForceLen->InsertNextValue(stress.norm());
				}
				if (recActive[REC_MATERIALID]) facetsMaterialId->InsertNextValue(b->material->id);
				if (recActive[REC_MASK]) facetsMask->InsertNextValue(b->groupMask);
				continue;
			}
		}
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
		if (recActive[REC_CLUMPID]) spheresUg->GetPointData()->AddArray(clumpId);
		if (recActive[REC_COLORS]) spheresUg->GetPointData()->AddArray(spheresColors);
		if (recActive[REC_VELOCITY]){
			spheresUg->GetPointData()->AddArray(spheresLinVelVec);
			spheresUg->GetPointData()->AddArray(spheresAngVelVec);
			spheresUg->GetPointData()->AddArray(spheresLinVelLen);
			spheresUg->GetPointData()->AddArray(spheresAngVelLen);
		}
		if (recActive[REC_STRESS]){
			spheresUg->GetPointData()->AddArray(spheresNormalStressVec);
			spheresUg->GetPointData()->AddArray(spheresShearStressVec);
			spheresUg->GetPointData()->AddArray(spheresNormalStressNorm);
		}
		if (recActive[REC_CPM]){
			spheresUg->GetPointData()->AddArray(cpmDamage);
			spheresUg->GetPointData()->AddArray(cpmSigma);
			spheresUg->GetPointData()->AddArray(cpmSigmaM);
			spheresUg->GetPointData()->AddArray(cpmTau);
		}
		if (recActive[REC_RPM]){
			spheresUg->GetPointData()->AddArray(rpmSpecNum);
			spheresUg->GetPointData()->AddArray(rpmSpecMass);
			spheresUg->GetPointData()->AddArray(rpmSpecDiam);
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
			writer->SetInput(spheresUg);
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
			writer->SetInput(facetsUg);
			writer->Write();	
		}
	}
	vtkSmartPointer<vtkPolyData> intrPd = vtkSmartPointer<vtkPolyData>::New();
	if (recActive[REC_INTR]){
		intrPd->SetPoints(intrBodyPos);
		intrPd->SetLines(intrCells);
		intrPd->GetCellData()->AddArray(intrForceN);
		intrPd->GetCellData()->AddArray(intrAbsForceT);
		#ifdef YADE_VTK_MULTIBLOCK
			if(!multiblock)
		#endif
			{
			vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
			if(compress) writer->SetCompressor(compressor);
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+"intrs."+lexical_cast<string>(scene->iter)+".vtp";
			writer->SetFileName(fn.c_str());
			writer->SetInput(intrPd);
			writer->Write();
		}
	}
	#ifdef YADE_VTK_MULTIBLOCK
		if(multiblock){
			vtkSmartPointer<vtkMultiBlockDataSet> multiblockDataset = vtkSmartPointer<vtkMultiBlockDataSet>::New();
			int i=0;
			if(recActive[REC_SPHERES]) multiblockDataset->SetBlock(i++,spheresUg);
			if(recActive[REC_FACETS]) multiblockDataset->SetBlock(i++,facetsUg);
			if(recActive[REC_INTR]) multiblockDataset->SetBlock(i++,intrPd);
			vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
			if(ascii) writer->SetDataModeToAscii();
			string fn=fileName+lexical_cast<string>(scene->iter)+".vtm";
			writer->SetFileName(fn.c_str());
			writer->SetInput(multiblockDataset);
			writer->Write();	
		}
	#endif
}

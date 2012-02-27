/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeLaw.hpp"
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LatticeBeamAngularSpring.hpp>
#include<yade/pkg-lattice/NonLocalDependency.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
#include<yade/pkg-lattice/LatticeSetParameters.hpp>
#include<yade/core/BodyContainer.hpp>
#include<yade/core/Scene.hpp>
#include<yade/lib-base/yadeWm3Extra.hpp>

LatticeLaw::LatticeLaw() : InteractionSolver()
{
//	respect_non_destroy = 500.0;
	backward_compatible = true;
}


LatticeLaw::~LatticeLaw()
{

}

bool LatticeLaw::deleteBeam(Scene* metaBody , LatticeBeamParameters* beam, Body* b)
{
	Real strain;

//	if ( nonlocal) // calculate strain with non-local law
//		strain = beam->nonLocalStrain / beam->nonLocalDivisor;
//	else

	
	strain = beam->strain();
	return     strain < -beam->criticalCompressiveStrain
		|| strain >  beam->criticalTensileStrain;
}


void LatticeLaw::calcBeamPositionOrientationNewLength(Body* body, BodyContainer* bodies)
{
// FIXME - verify that this updating of length, position, orientation and color is in correct place/plugin
////	LatticeBeamParameters* beam 	  = static_cast<LatticeBeamParameters*>(body->physicalParameters.get());

	LatticeBeamParameters* beam 	  = YADE_CAST<LatticeBeamParameters*>(body->physicalParameters.get());
	if(beam==0 || (!(*(bodies)).exists(beam->id1)) || (!(*(bodies)).exists(beam->id2)))
	{
		std::string error(boost::lexical_cast<std::string>(body->getId()) + " is not a beam, but it SHOULD really be a beam, ERROR. You should check if bodies are not renumbered somehow after reloading a simulation. That's probably a container error, that is changing the IDs when not asked to do so.\n");
		std::cerr << "================\nERROR: ";
		std::cerr << error;
		std::cerr << "================\n";
		throw error;
	}
//std::cerr << body->getId() << " " << beam->id1 << " " << beam->id2 << "\n";

	Body* bodyA 			  = (*(bodies))[beam->id1].get();
	Body* bodyB 			  = (*(bodies))[beam->id2].get();
	Se3r& se3A 			  = bodyA->physicalParameters->se3;
	Se3r& se3B 			  = bodyB->physicalParameters->se3;
	
	Se3r se3Beam;
	se3Beam.position 		  = (se3A.position + se3B.position)*0.5;
	Vector3r dist 			  =  se3A.position - se3B.position;
	
	Real length 			  = dist.Normalize();

	beam->direction 		  = dist;
	beam->length 			  = length;

/////////////////////////////////////////
//	// rotate otherDirection
//	{
		beam->otherDirection -= beam->direction.Dot(beam->otherDirection)*beam->direction;
		Quaternionr tor(beam->direction,beam->torsionAngle);
		beam->otherDirection = tor * beam->otherDirection;
		beam->otherDirection.Normalize();
//	}
/////////////////////////////////////////
	
	Vector3r CP			  = beam->direction.Cross(beam->otherDirection);
	se3Beam.orientation.FromRotationMatrix( Matrix3r( beam->direction , beam->otherDirection , CP , true ) );
	
	beam->se3Displacement.position    = se3Beam.position - beam->se3.position;
	beam->se3Displacement.orientation = se3Beam.orientation * beam->se3.orientation.Conjugate();
	
	beam->se3 			  = se3Beam;

	// PhysicalActionResetter
	beam->bendingRotation		= Quaternionr(1.0,0.0,0.0,0.0);
	beam->torsionalRotation		= Quaternionr(1.0,0.0,0.0,0.0);
	beam->torsionAngle		= 0.0;
	beam->count			= 0.0;

//	double zzzzz=beam->direction.dot(beam->otherDirection);
//	if(zzzzz < -3e-16 || zzzzz > 3e-16 )
//		std::cerr << "============== " << zzzzz << "\n";
//		std::cerr << " (" << beam->otherDirection << ")";

}

void LatticeLaw::action(Scene* lattice)
{
	futureDeletes.clear();


	int nodeGroupMask  = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->nodeGroupMask;
	int beamGroupMask  = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->beamGroupMask;

	BodyContainer* bodies = lattice->bodies.get();

	InteractionContainer::iterator angles     = lattice->interactions->begin();
	InteractionContainer::iterator angles_end = lattice->interactions->end();

	{ // 'B' calculate how many neighbors each beam currently has
		for(  ; angles != angles_end; ++angles )
		{
			if(        bodies->exists( (*angles)->getId1() )        // FIXME - remove this test ....
				&& bodies->exists( (*angles)->getId2() ) ) 	//
			{
				++(static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId1()])->physicalParameters.get())->count);
                                ++(static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId2()])->physicalParameters.get())->count);
			}
			// else FIXME - delete unused angularSpring
		}
	}
	{ // calculate new angles between beams
		angles     = lattice->interactions->begin();
		angles_end = lattice->interactions->end();
		for(  ; angles != angles_end; ++angles )
		{
			if( 	   bodies->exists( (*angles)->getId1() )
				&& bodies->exists( (*angles)->getId2() ) )
			{
				LatticeBeamParameters* beam1 = static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId1()])->physicalParameters.get());
				LatticeBeamParameters* beam2 = static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId2()])->physicalParameters.get());
				
				LatticeBeamAngularSpring* an = static_cast<LatticeBeamAngularSpring*>((*angles)->interactionPhysics.get());
				
				Vector3r	newCP;
				Vector3r&	lastCP		= an->lastCrossProduct;
				bool&		planeSwap180	= an->planeSwap180;
				bool		oldPS180	= planeSwap180;
				Real		planeAngle;

				newCP		= beam1->direction.Cross(beam2->direction);
				if(newCP.Dot(lastCP)< 0.0)
					planeSwap180=!planeSwap180;
				
				// beam1->direction and beam2->direction are unit vectors
				planeAngle	= ( planeSwap180      ? -1.0 : 1.0 ) * unitVectorsAngle(beam1->direction,beam2->direction);

				Real sinAngleSquared = newCP.SquaredLength();
				//Real cosAngleSquared = std::pow(beam1->direction.Dot(beam2->direction),2);
				if(sinAngleSquared > 0.0001)
					lastCP = newCP;
				else
					planeSwap180 = oldPS180;	
				//sinAngleSquared = newCP.Length();
				
		 // 'B' calculate needed beam rotations 
				Real angleDifference = planeAngle - an->initialPlaneAngle;
			
				// FIXME - is it really necessary?
				if( angleDifference > 0 )
					while(angleDifference > Mathr::PI) angleDifference -= Mathr::TWO_PI;
				else
					while(angleDifference < -Mathr::PI) angleDifference += Mathr::TWO_PI;
				// FIXME END
				
				// axis of bending rotation is orthogonal to the plane between two beams and must have unit length
				lastCP.Normalize();
				Quaternionr rotationDifference1,rotationDifference2;

				rotationDifference1.FromAxisAngle(lastCP, ( planeSwap180 ? -1.0 : 1.0 ) * angleDifference/beam1->count);
				rotationDifference2.FromAxisAngle(lastCP, ( planeSwap180 ? 1.0 : -1.0 ) * angleDifference/beam2->count);
				
				beam1->bendingRotation = beam1->bendingRotation * rotationDifference1;
				beam2->bendingRotation = beam2->bendingRotation * rotationDifference2;

		// 'T' calculate torsional rotation
			if(calcTorsion)
			{

				bool sameFlow			= true;
				if( (beam1->id1 == beam2->id1) || (beam1->id2 == beam2->id2) )
					sameFlow		= false;

				Real offPlaneAngle1		= unitVectorsAngle(beam1->otherDirection,lastCP);
				Real offPlaneAngle2		= unitVectorsAngle(beam2->otherDirection,lastCP);
		
				/////////////////////////////////////////////////////////////
		// FIXME [1] look below
				Quaternionr	aligner1,aligner2;
				aligner1.FromAxisAngle(beam1->direction , offPlaneAngle1);
				aligner2.FromAxisAngle(beam2->direction , offPlaneAngle2);
						
				Vector3r	dir1			= aligner1 * lastCP;
				Vector3r	dir2			= aligner2 * lastCP;
		
				// FIXME
				if( dir1.Dot(beam1->otherDirection) < 0.999999 )
					offPlaneAngle1   *= -1.0;
				if( dir2.Dot(beam2->otherDirection) < 0.999999 )
					offPlaneAngle2   *= -1.0;
	//			if( (lastCP.Cross(beam1->direction)).Dot(beam1->otherDirection) > 0.0 )
	//				offPlaneAngle1   *= -1.0;
	//			if( (lastCP.Cross(beam2->direction)).Dot(beam2->otherDirection) > 0.0 )
	//				offPlaneAngle2   *= -1.0;
	//	// FIXME [1] END
				/////////////////////////////////////////////////////////////

				Real offPlaneAngleDifference1	= an->initialOffPlaneAngle1 - offPlaneAngle1;
				Real offPlaneAngleDifference2	= an->initialOffPlaneAngle2 - offPlaneAngle2;
				
				if( offPlaneAngleDifference1 > 0 )
					while(offPlaneAngleDifference1 > Mathr::PI) offPlaneAngleDifference1 -= Mathr::TWO_PI;
				else
					while(offPlaneAngleDifference1 < -Mathr::PI) offPlaneAngleDifference1 += Mathr::TWO_PI;
				
				if( offPlaneAngleDifference2 > 0 )
					while(offPlaneAngleDifference2 > Mathr::PI) offPlaneAngleDifference2 -= Mathr::TWO_PI;
				else
					while(offPlaneAngleDifference2 < -Mathr::PI) offPlaneAngleDifference2 += Mathr::TWO_PI;

				Real torsionAngle1, torsionAngle2;
//				if(sameFlow)
//				{
//					torsionAngle1 = offPlaneAngleDifference1 - offPlaneAngleDifference2;
//					torsionAngle2 = -torsionAngle1;
//				}
//				else
//				{
//					torsionAngle1 = offPlaneAngleDifference1 + offPlaneAngleDifference2;
//					torsionAngle2 = torsionAngle1;
//				}
				torsionAngle1 = offPlaneAngleDifference1 + (sameFlow?-1:1)*offPlaneAngleDifference2;
		//		torsionAngle1 *= cosAngleSquared;
				torsionAngle2 = (sameFlow?-1:1)*torsionAngle1;

				
				if( torsionAngle1 > 0 )
					while(torsionAngle1 > Mathr::PI) torsionAngle1 -= Mathr::TWO_PI;
				else
					while(torsionAngle1 < -Mathr::PI) torsionAngle1 += Mathr::TWO_PI;
				
				if( torsionAngle2 > 0 )
					while(torsionAngle2 > Mathr::PI) torsionAngle2 -= Mathr::TWO_PI;
				else
					while(torsionAngle2 < -Mathr::PI) torsionAngle2 += Mathr::TWO_PI;

				// calculate beam swirl
		// FIXME [1] - this can solve this double checking of angle +/-180
				if( an->lastOffPlaneAngleDifference1 >  3.0 && offPlaneAngleDifference1 < -3.0 ) ++(an->swirl1);
				if( an->lastOffPlaneAngleDifference1 < -3.0 && offPlaneAngleDifference1 >  3.0 ) --(an->swirl1);
				if( an->lastOffPlaneAngleDifference2 >  3.0 && offPlaneAngleDifference2 < -3.0 ) ++(an->swirl2);
				if( an->lastOffPlaneAngleDifference2 < -3.0 && offPlaneAngleDifference2 >  3.0 ) --(an->swirl2);

				an->lastOffPlaneAngleDifference1   = offPlaneAngleDifference1;
				an->lastOffPlaneAngleDifference2   = offPlaneAngleDifference2;

				Real swirlAngle1 = (offPlaneAngleDifference1 + Mathr::TWO_PI * an->swirl1)*sinAngleSquared;
				Real swirlAngle2 = (offPlaneAngleDifference2 + Mathr::TWO_PI * an->swirl2)*sinAngleSquared;
				
				beam1->torsionAngle += 0.5*( swirlAngle1 + torsionAngle1)/beam1->count;
				beam2->torsionAngle += 0.5*( swirlAngle2 + torsionAngle2)/beam2->count;
				
				//axis of torsional rotation is the other's beam direction, the beam tries to rotate its neighbours to fit its orientation
				rotationDifference1.FromAxisAngle(beam2->direction, ((sameFlow ? -1.0 : 1.0)*swirlAngle1-torsionAngle2)/beam1->count);
				rotationDifference2.FromAxisAngle(beam1->direction, ((sameFlow ? -1.0 : 1.0)*swirlAngle2-torsionAngle1)/beam2->count);
			
				beam1->torsionalRotation = beam1->torsionalRotation * rotationDifference1;
				beam2->torsionalRotation = beam2->torsionalRotation * rotationDifference2;
			}
			}
		}
	}

	BodyContainer::iterator bi    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();

	for(  ; bi!=biEnd ; ++bi )  // loop over all beams
	{
		Body* body = (*bi).get();
		if( ! ( body->getGroupMask() & beamGroupMask ) )
			continue; // skip non-beams

		// next beam
		LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>(body->physicalParameters.get() );

		Real      stretch      = beam->length - beam->initialLength;

		// 'D' from picture. how much beam wants to change length at each node, to bounce back through original length to mirror position.
		Vector3r  displacementLongitudal = beam->direction * stretch;

		{ // 'E_min' 'E_max' criterion
			if( deleteBeam(lattice , beam, body) ) // calculates strain
			{
//				std::cerr << "one beam deleted\n";
			
	if(respect_non_destroy == beam->longitudalStiffness)
	{
		std::cerr << "Tried to delete a non-destroy beam, skipping.\n";
	}
	else
	{
				futureDeletes.push_back(body->getId());
				continue;
	}
			}
		}
		
		LatticeNodeParameters* node1 = static_cast<LatticeNodeParameters*>(((*(bodies))[beam->id1])->physicalParameters.get());
		LatticeNodeParameters* node2 = static_cast<LatticeNodeParameters*>(((*(bodies))[beam->id2])->physicalParameters.get());

		++(node1->countIncremental);
		++(node2->countIncremental);

		{ // 'W' from picture - previous displacement of the beam. try to do it again.
			node1->displacementIncremental += beam->se3Displacement.position;
			node2->displacementIncremental += beam->se3Displacement.position;
		}

		{ // 'R' from picture - previous rotation of the beam. try to do it again.
			Vector3r halfLength = beam->length * beam->direction * 0.5;
			Vector3r R = beam->se3Displacement.orientation * ( halfLength) - halfLength;
			node1->displacementIncremental += R;
			node2->displacementIncremental -= R;
			//node1->displacementIncremental += beam->se3Displacement.orientation * ( halfLength) - halfLength;
			//node2->displacementIncremental += beam->se3Displacement.orientation * (-halfLength) + halfLength;

		//	beam->otherDirection		= beam->se3Displacement.orientation * beam->otherDirection;
		}

		{ // 'D' to nodes
			node1->countStiffness += (beam->longitudalStiffness)/(beam->initialLength);
			node2->countStiffness += (beam->longitudalStiffness)/(beam->initialLength);
			node1->displacementAlignmental -= (displacementLongitudal * beam->longitudalStiffness)/(beam->initialLength);
			node2->displacementAlignmental += (displacementLongitudal * beam->longitudalStiffness)/(beam->initialLength);
		}

		if( beam->count != 0 )
		{
		Vector3r beam_vec = beam->length * beam->direction;
		{ // 'B' from picture - rotate to align with neighbouring beams

/*   
 *   Bending stiffness different for compression and tension.
 */
			static bool first=true;
			if(first)
			{
if(tension_compression_different_stiffness){
				std::cerr << "\nusing k.b tension=0.6, k.b compression=0.2 !\n/beam->initialLength !\n";
}else{
				std::cerr << "\nNOT! using k.b tension=0.6, k.b compression=0.2 ! (just a classical formula)\n\n";
}
				first=false;
			}

			Real kb = beam->bendingStiffness;

if(tension_compression_different_stiffness)
{

			Real Em = beam->criticalTensileStrain/3.0;
			Real x  = beam->strain();
			//const Real howmuch = 0.5;
			const Real howmuch = 0.3333333333333333333333333333;
		
		//	x<-E ? kb*0.2 : ((kb*0.2+kb)/2+(kb*0.2-kb)/(-2*E)*x)

			if( x < -Em )
				kb *= howmuch;
			else if( x < Em )
				kb = (kb*howmuch+kb)/2.0+x*(kb*howmuch-kb)/(-2.0*Em);
			// if strain > criticalTensileStrain/2.0 then kb is not changed
}

			node1->countStiffness += kb/beam->initialLength;
			node2->countStiffness += kb/beam->initialLength;

			node1->displacementAlignmental += ((beam->bendingRotation * ( beam_vec) - beam_vec) * kb)/(beam->initialLength);
			node2->displacementAlignmental += ((beam->bendingRotation * (-beam_vec) + beam_vec) * kb)/(beam->initialLength);


/*			node1->countStiffness += (beam->bendingStiffness)/(beam->initialLength);
			node2->countStiffness += (beam->bendingStiffness)/(beam->initialLength);

			node1->displacementAlignmental += ((beam->bendingRotation * ( beam_vec) - beam_vec) * beam->bendingStiffness)/(beam->initialLength);
			node2->displacementAlignmental += ((beam->bendingRotation * (-beam_vec) + beam_vec) * beam->bendingStiffness)/(beam->initialLength);
*/
			beam->shearing_strain = ((beam->bendingRotation * (-beam_vec) + beam_vec) * beam->bendingStiffness)/(beam->initialLength);

		}
		
		if(calcTorsion)
		{ // 'T' from picture - torsion

			node1->countStiffness += (beam->torsionalStiffness)/(beam->initialLength);
			node2->countStiffness += (beam->torsionalStiffness)/(beam->initialLength);

			node1->displacementAlignmental += ((beam->torsionalRotation * ( beam_vec) - beam_vec) * beam->torsionalStiffness)/(beam->initialLength);
			node2->displacementAlignmental += ((beam->torsionalRotation * (-beam_vec) + beam_vec) * beam->torsionalStiffness)/(beam->initialLength);
		}
		}
	}

	//std::cout << "iter:" << Omega::instance().getCurrentIteration() << " beams_deleted: " << futureDeletes.size() << std::endl;
	
	{ // move nodes
		bi    = bodies->begin();
		biEnd = bodies->end();
		for(  ; bi!=biEnd ; ++bi )  // loop over all nodes -- move them according to displacements from beams
		{
			Body* body = (*bi).get();
	
			if( ! ( body->getGroupMask() & nodeGroupMask ) )
				continue; // skip non-nodes
			
			LatticeNodeParameters* node = static_cast<LatticeNodeParameters*>(body->physicalParameters.get() );
			
			{ // check nodes if they have any beams left
				if(node->countIncremental == 0) 
				{ // node not moving (marked for deletion)
					futureDeletes.push_back(body->getId());
					continue; 
				}
			}
			Vector3r straight_line_movement		= node->displacementIncremental / node->countIncremental;
			Vector3r elastic_deformation_movement	= node->displacementAlignmental / node->countStiffness;

/////////// No damping
			Vector3r displacementTotal		= straight_line_movement + elastic_deformation_movement;
				//  node->displacementIncremental / node->countIncremental 
				//+ node->displacementAlignmental / node->countStiffness;
/////////// No damping END

///////////// ZZ
//			Real DAMPING=0.3;
//			{
//				Vector3r sign =  elastic_deformation_movement - node->previousDisplacement;
//				for(int i=0; i<3; i++)
//					if(sign[i]<0)
//						elastic_deformation_movement[i] *= (1.0-DAMPING);
//			}
//			Vector3r displacementTotal	= straight_line_movement + elastic_deformation_movement;
//			node->previousDisplacement	= elastic_deformation_movement;
///////////// ZZ END
//
///////////// ZZ2
//			Vector3r displacementTotal	= straight_line_movement + elastic_deformation_movement;
//			Real DAMPING=0.3;
//			{
//				for(int i=0; i<3; i++)
//					if(displacementTotal[i] * node->previousDisplacement[i] < 0)
//						displacementTotal[i] *= (1.0-DAMPING);
//			}
//
//			node->previousDisplacement	= displacementTotal;
///////////// ZZ2 END
//
///////////// ZZ3
//			Real DAMPING=0.3;
//			{
//				for(int i=0; i<3; i++)
//					if(elastic_deformation_movement[i] * node->previousDisplacement[i] < 0)
//						elastic_deformation_movement[i] *= (1.0-DAMPING);
//			}
//			Vector3r displacementTotal	= straight_line_movement + elastic_deformation_movement;
//			node->previousDisplacement	= elastic_deformation_movement;
///////////// ZZ3 END

			// FIXME - ponder changing names:
			// 
			// uniform straight line movement:
			//  displacementIncremental	(numerator) 
			//  countIncremental		(denominator)
			//
			// elastic deformation:
			//  displacementAlignmental	(numerator)
			//  countStiffness		(denominator)

			node->countIncremental		= 0;
			node->countStiffness            = 0;
			node->displacementIncremental 	= Vector3r(0.0,0.0,0.0);
			node->displacementAlignmental   = Vector3r(0.0,0.0,0.0);


			{
				static bool first1(true);
				if(first1)
				{
					first1=false;
					std::cerr << "====== backward_compatible is:" << (int)(backward_compatible) << "\n";
				}
			}

			if(backward_compatible)
			{
						if(body->isDynamic)
						{
							node->se3.position      += displacementTotal;
													// DAMPING: *(1.0 - damping_with_energy_loss_0_to_1);
						//	node->se3.orientation	+= ;
						}
						// FIXME FIXME FIXME FIXME FIXME FIXME FIXME [2]
						else
						{
							if(!roughEdges)// FIXME - else move only in x and z directions
							{
								node->se3.position[0]   += displacementTotal[0];
								node->se3.position[2]   += displacementTotal[2];
								//std::cerr << "zz!\n";
							}
						}
						if(ensure2D)
							node->se3.position[2] = 0; // ensure 2D
			}
			else
			{
				static bool first(true);
				if(first)
				{
					first=false;
					std::cerr << "Using DOFs !\n";
				}
				
				{// FIXED (was FIXME [2]) here is good.
					if(node->blockedDOFs == PhysicalParameters::DOF_NONE)
					{
						node->se3.position      += displacementTotal;
					}
					else
					{
						if((node->blockedDOFs & PhysicalParameters::DOF_X) == 0) node->se3.position[0] += displacementTotal[0];
						if((node->blockedDOFs & PhysicalParameters::DOF_Y) == 0) node->se3.position[1] += displacementTotal[1];
						if((node->blockedDOFs & PhysicalParameters::DOF_Z) == 0) node->se3.position[2] += displacementTotal[2];
					}
				}
			}
		}
	}
	
	{ // calc new beam position:    X_b = ( X_n1 + X_n2 ) / 2 ,
	  // calc new beam orientation: ------------||-----------
	  // store position/orientation displacement, for 'W' and 'R' in se3Displacement
		bi    = bodies->begin();
		biEnd = bodies->end(); 
		for(  ; bi!=biEnd ; ++bi )  // loop over all beams
		{
			Body* body = (*bi).get();
			if( body->getGroupMask() & beamGroupMask )
				calcBeamPositionOrientationNewLength(body,bodies);
		}
	} 
	
	{ // delete all beams and nodes marked for deletion 
	  // - beams that exceeded critical tensile/compressive strain
	  // - nodes that have no more beams
		std::list<unsigned int>::iterator vend = futureDeletes.end();
		for( std::list<unsigned int>::iterator vsta = futureDeletes.begin() ; vsta != vend ; ++vsta)
			bodies->erase(*vsta); 
	}
}









/* 
 * this method should be placed in a class derived from this LatticeLaw:
 * it will be a LatticeLaw with softening.
 *
 *
 *
 *
	//bool FIXME_useBendTensileSoftening = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->useBendTensileSoftening;
	//bool FIXME_useStiffnessSoftening   = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->useStiffnessSoftening;
	//
	//void* nonl = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->nonl;
	//unsigned long int& total = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->total;
	//nonlocal = ( total != 0 );
	//Real FIXME_range                   = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->range;

 *
 * this non-local calculations stuff was right after 'B' calculate needed beam rotations
 * it should be placed in a derived class.
 * Or some other place, where it could be turned on/off easily
 *
	for(  ; bi!=biEnd ; ++bi )  // clear all non-local vaules in all beams.
	{
		Body* body = (*bi).get();
		if( ! ( body->getGroupMask() & beamGroupMask ) )
			continue; // skip non-beams

		LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>(body->physicalParameters.get() );

		beam->nonLocalStrain = 0.0;
		beam->nonLocalDivisor = 0.0;
	}

        //if (lattice->interactions->size() != 0) // calculate non-local strain of beams
        if ( nonlocal ) // calculate non-local strain of beams
	{
		unsigned long int counter = 0;
		LatticeSetParameters::NonLocalInteraction* nb = reinterpret_cast<LatticeSetParameters::NonLocalInteraction*>(nonl);
		for( ; counter != total ; ++nb , ++counter )
		{
			if(	   bodies->exists( nb->id1 )
					&& bodies->exists( nb->id2 ) ) // both beams exist ..
			{
				LatticeBeamParameters* beam1 = static_cast<LatticeBeamParameters*>(((*(bodies))[nb->id1])->physicalParameters.get());
				LatticeBeamParameters* beam2 = static_cast<LatticeBeamParameters*>(((*(bodies))[nb->id2])->physicalParameters.get());

				Real cosAngle = std::abs(beam1->direction.dot(beam2->direction)); // direction is a unit vector, abs() - because I want a number between 0 .. 1

				Real sqDist = (beam1->se3.position - beam2->se3.position).squaredLength();
				Real len = FIXME_range;
				static Real sqPi = std::sqrt(Mathr::PI);
				Real gaussValue  = std::exp( - (sqDist/(len*len)) ) / ( len * sqPi );
				//	Real gaussValue  = std::exp( - std::pow( dist / len , 2) ) / ( len * sqPi );

				Real b2 = beam2->length * gaussValue * cosAngle;
				Real b1 = beam1->length * gaussValue * cosAngle;
				beam1->nonLocalStrain += beam2->strain() * b2;// * cosAngle;
				beam2->nonLocalStrain += beam1->strain() * b1;// * cosAngle;
				beam1->nonLocalDivisor += b2;
				beam2->nonLocalDivisor += b1;
			}
			// else FIXME - delete unused NonLocalDependency
		}
	}

	bi    = bodies->begin();
	biEnd = bodies->end();

 *
 *
 * this stiffness softening must go somewhere else, too
 *
 *


		Real FIXME_useStiffnessSoftening_Factor=1.0;
		if(FIXME_useStiffnessSoftening)
		{ // give 'D' to nodes
			Real sTrAiN = beam->strain();                                                  // only with over half of tension. (compression is ignored)
			FIXME_useStiffnessSoftening_Factor = sTrAiN > beam->criticalTensileStrain*0.5 
				? (2.0-sTrAiN*2.0/beam->criticalTensileStrain) : 1.0 ;
			if(FIXME_useStiffnessSoftening_Factor<0.01) FIXME_useStiffnessSoftening_Factor=0.01;

			assert(FIXME_useStiffnessSoftening_Factor<=1.0);

			Real FIXME_longitudalStiffness = beam->longitudalStiffness*FIXME_useStiffnessSoftening_Factor;
			node1->countStiffness += FIXME_longitudalStiffness;
			node2->countStiffness += FIXME_longitudalStiffness;
			node1->displacementAlignmental -= displacementLongitudal * FIXME_longitudalStiffness;
			node2->displacementAlignmental += displacementLongitudal * FIXME_longitudalStiffness;
		}
 
 
 *
 *
 *
 *
		if(FIXME_useBendTensileSoftening)
		{ // 'B' from picture - rotate to align with neighbouring beams
			Real FIXME_bendingStiffness = beam->bendingStiffness * FIXME_useStiffnessSoftening_Factor;
			Real strain                     = beam->strain();

			Real factor;

			if( strain > 0 )
			{
				//                              (   0 ... 0.5 ... 1 ... inf          ) 
				factor                  = 1.0 - (strain / beam->criticalTensileStrain); // positive
				if(factor< 0.0) factor  = 0.0;
				FIXME_bendingStiffness *=factor;
			}
			else
			{
				//                        (  0 ... -0.5 ... -1 ... -inf            )
				factor                  = (strain / beam->criticalCompressiveStrain); // negative
				if(factor<-1.0) factor  =-1.0;
				FIXME_bendingStiffness -= (1.0-FIXME_bendingStiffness)*factor;
			}
			node1->countStiffness += FIXME_bendingStiffness;
			node2->countStiffness += FIXME_bendingStiffness;
			Vector3r axis           = Vector3r( 0.0 , 0.0 , 1.0 ); // 2D - always rotation around z 
			Real angle              = beam->bendingRotation;
			angle /= beam->count;
			Quaternionr rotation;
			rotation.fromAxisAngle(axis,angle);
			Vector3r length = beam->length * beam->direction;// * 0.5;// beam->bendingStiffness / beam->count;// * 0.5; // FIXME - duplicate line
			node1->displacementAlignmental += (rotation * ( length) - length) * FIXME_bendingStiffness;
			node2->displacementAlignmental += (rotation * (-length) + length) * FIXME_bendingStiffness;
			beam->bendingRotation	= 0.0;
			beam->torsionalRotation	= 0.0;
			beam->count             = 0.0;
		}

*/

YADE_PLUGIN((LatticeLaw));

YADE_REQUIRE_FEATURE(PHYSPAR);


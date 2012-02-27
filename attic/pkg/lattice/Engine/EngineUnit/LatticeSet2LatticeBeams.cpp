/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeSet2LatticeBeams.hpp"
#include<yade/pkg-lattice/LatticeSetParameters.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LineSegment.hpp>
#include<yade/core/Scene.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

void LatticeSet2LatticeBeams::go(	  const shared_ptr<PhysicalParameters>& ph
					, shared_ptr<GeometricalModel>& 
					, const Body* body)
{
	int beamGroupMask = static_cast<const LatticeSetParameters*>(ph.get())->beamGroupMask;
	const Scene * ncb = static_cast<const Scene*>(body);
	const shared_ptr<BodyContainer>& bodies = ncb->bodies;

	BodyContainer::iterator bi    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();
	for(  ; bi!=biEnd ; ++bi )
	{
		Body* b = (*bi).get();
		if( b->getGroupMask() & beamGroupMask )
		{
			LineSegment* line 		= static_cast<LineSegment*>	     (b->geometricalModel.get());
			LatticeBeamParameters* beam 	= static_cast<LatticeBeamParameters*>(b->physicalParameters.get());

// FIXME - this copying of length between latticeBeam geometry and physics, inside Scene could
//         be done just once, if length was inside shared_ptr. This can be improved once we make
//         indexable Parameters: Velocity, Position, Orientation, ....
			line->length                    = beam->length;

			// FIXME - display aggregates as brown, bonds as dark brown.
			                                             line->diffuseColor = Vector3r(0.6,0.6,0.6); // DEFAULT color

			      if(beam->longitudalStiffness == 14.0 ) line->diffuseColor = Vector3r(0.06,0.06,0.06); // BOND  
			      if(beam->longitudalStiffness == 20.0 ) line->diffuseColor = Vector3r(0.46,0.46,0.46); // CEMENT   
			      if(beam->longitudalStiffness == 60.0 ) line->diffuseColor = Vector3r(0.0,0.0,0.0); // AGGREGATE
			      if(beam->longitudalStiffness == 160.0) line->diffuseColor = Vector3r(1.0,1.0,0.35); // STEEL

			      if(beam->longitudalStiffness == 500.0 ) line->diffuseColor = Vector3r(0.0,1.0,0.0); // NON-DESTROY
			      
			      
//			      if(beam->longitudalStiffness == 0.69999999999999996 ) line->diffuseColor = Vector3r(0.06,0.06,0.06); // BOND  
//			      if(beam->longitudalStiffness == 1.0 ) line->diffuseColor = Vector3r(0.46,0.46,0.46); // CEMENT   
//			      if(beam->longitudalStiffness == 3.0 ) line->diffuseColor = Vector3r(0.0,0.0,0.0); // AGGREGATE
//			      if(beam->longitudalStiffness == 8.0) line->diffuseColor = Vector3r(1.0,1.0,0.35); // STEEL
//
//			      if(beam->longitudalStiffness == 10.0 ) line->diffuseColor = Vector3r(0.0,1.0,0.0); // NON-DESTROY
		}
	}
}

YADE_PLUGIN((LatticeSet2LatticeBeams));

YADE_REQUIRE_FEATURE(PHYSPAR);


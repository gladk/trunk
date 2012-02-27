/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once


#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/core/InteractionSolver.hpp>
#include<yade/core/BodyContainer.hpp>
#include <list>



class LatticeLaw : public InteractionSolver
{
	private:
		std::list<unsigned int> futureDeletes;
		//bool nonlocal;
                
                bool deleteBeam(Scene* lattice , LatticeBeamParameters* beam, Body*);
                void calcBeamPositionOrientationNewLength(Body* body, BodyContainer* bodies);
	public :
		bool roughEdges,ensure2D,calcTorsion,tension_compression_different_stiffness;// FIXME, FIXME, FIXME 
		bool backward_compatible;
		Real respect_non_destroy;
		LatticeLaw();
		virtual ~LatticeLaw();
		void action();
		REGISTER_ATTRIBUTES(InteractionSolver,
			(roughEdges)
			(ensure2D)
			(calcTorsion)
			(backward_compatible)
			(tension_compression_different_stiffness)
			(respect_non_destroy)
			// must go to derived class
			//(nonlocal)
		);
/// Serializtion
	REGISTER_CLASS_NAME(LatticeLaw);
	REGISTER_BASE_CLASS_NAME(InteractionSolver);
};


REGISTER_SERIALIZABLE(LatticeLaw);





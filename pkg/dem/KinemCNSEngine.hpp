/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-dem/KinemSimpleShearBox.hpp>



class KinemCNSEngine : public KinemSimpleShearBox
{
	private :
		int	temoin,it_stop
			;

	public :
		void 	action()
			;



	YADE_CLASS_BASE_DOC_ATTRS_CTOR(KinemCNSEngine,KinemSimpleShearBox,
		"To apply a Constant Normal Stifness (CNS) shear for a parallelogram box (simple shear)\n\nThis engine, useable in simulations implying one deformable parallelepipedic box, allows to translate horizontally the upper plate while the lateral ones rotate so that they always keep contact with the lower and upper walls. The upper plate can move not only horizontally but also vertically, so that the normal rigidity defined by DeltaF(upper plate)/DeltaU(upper plate) = constant (= *KnC* defined by the user).\n\nThe movement is moreover controlled by the user via a *shearSpeed* which is the horizontal speed of the upper wall, and by a maximum value of horizontal displacement *gammalim* (of the upper plate), after which the shear stops.\n\n.. note::\n\t not only the positions of walls are updated but also their speeds, which is all but useless considering the fact that in the contact laws these velocities of bodies are used to compute values of tangential relative displacements.\n\n.. warning::\n\tBut, because of this last point, if you want to use later saves of simulations executed with this Engine, but without that stopMovement was executed, your boxes will keep their speeds => you will have to cancel them by hand in the .xml",
		((Real,shearSpeed,0.0,,"the speed at wich the shearing is performed : speed of the upper plate [m/s]"))
		((Real,gammalim,0.0,,"the value of tangential displacement (of upper plate) at wich the shearing is stopped [m]"))
		((Real,gamma,0.0,,"current value of tangential displacement [m]"))
		((Real,KnC,10.0e6,,"the normal rigidity chosen by the user [MPa/mm] - the conversion in Pa/m will be made")),
		temoin=0;
		it_stop=0;
	);

};

REGISTER_SERIALIZABLE(KinemCNSEngine);



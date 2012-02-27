/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-dem/KinemSimpleShearBox.hpp>




class KinemCNLEngine : public KinemSimpleShearBox
{
	private :

		int	temoin,// utile pour savoir ou on en est
			it_stop
			;

	public :
		void	action()
			;

	
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(KinemCNLEngine,KinemSimpleShearBox,
				 "To apply a constant normal stress shear (i.e. Constant Normal Load : CNL) for a parallelogram box (simple shear box : :yref:`SimpleShear` Preprocessor or scripts/simpleShear.py)\n\nThis engine allows one to translate horizontally the upper plate while the lateral ones rotate so that they always keep contact with the lower and upper walls.\n\nIn fact the upper plate can move not only horizontally but also vertically, so that the normal stress acting on it remains constant (this constant value is not chosen by the user but is the one that exists at the beginning of the simulation)\n\nThe right vertical displacements which will be allowed are computed from the rigidity Kn of the sample over the wall (so to cancel a deltaSigma, a normal dplt deltaSigma*S/(Kn) is set)\n\nThe movement is moreover controlled by the user via a *shearSpeed* which will be the speed of the upper wall, and by a maximum value of horizontal displacement *gammalim*, after which the shear stops.\n\n.. note::\n\tNot only the positions of walls are updated but also their speeds, which is all but useless considering the fact that in the contact laws these velocities of bodies are used to compute values of tangential relative displacements.\n\n.. warning::\n\tBecause of this last point, if you want to use later saves of simulations executed with this Engine, but without that stopMovement was executed, your boxes will keep their speeds => you will have to cancel them 'by hand' in the .xml.\n",
				 ((Real,shearSpeed,0.0,,"the speed at wich the shearing is performed : speed of the upper plate [m/s]"))
				 ((Real,gammalim,0.0,,"the value of tangential displacement (of upper plate) at wich the shearing is stopped [m]"))
				 ((Real,gamma,0.0,,"current value of tangential displacement [m]"))
				 ((std::vector<Real>,gamma_save,,,"vector with the values of gamma at which a save of the simulation is performed [m]")),
// 				 ((Real,coeff_dech,1.0,,"in the case of the use of 'Law2_ScGeom6D_NormalInelasticityPhys_NormalInelasticity' for ex, where kn(unload)#kn(load). The engine cares to find the value at the first run BROKEN actually")),
				temoin=0;
				it_stop=0;
				 );
};

REGISTER_SERIALIZABLE(KinemCNLEngine);



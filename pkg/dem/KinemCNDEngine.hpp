/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-dem/KinemSimpleShearBox.hpp>


class KinemCNDEngine : public KinemSimpleShearBox
{
	private :
		int temoinfin;



	public :
		void action();

	protected :

		YADE_CLASS_BASE_DOC_ATTRS_CTOR(KinemCNDEngine,KinemSimpleShearBox,
			"To apply a Constant Normal Displacement (CND) shear for a parallelogram box\n\n \tThis engine, designed for simulations implying a simple shear box (:yref:`SimpleShear` Preprocessor or scripts/simpleShear.py), allows one to perform a constant normal displacement shear, by translating horizontally the upper plate, while the lateral ones rotate so that they always keep contact with the lower and upper walls.",
			((Real,shearSpeed,0.0,,"the speed at which the shear is performed : speed of the upper plate [m/s]"))
			((Real,gammalim,0.0,,"the value of the tangential displacement at wich the displacement is stopped [m]"))
			((Real,gamma,0.0,,"the current value of the tangential displacement"))
			((std::vector<Real>,gamma_save,,,"vector with the values of gamma at which a save of the simulation is performed [m]")),
			temoinfin=0;
						);
};

REGISTER_SERIALIZABLE(KinemCNDEngine);



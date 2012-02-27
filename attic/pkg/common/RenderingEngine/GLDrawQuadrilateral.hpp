/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>

class GLDrawQuadrilateral : public GlGeometricalModelFunctor
{
	private:
		// FIXME - this does not belong here
		void calculateStrainQuadrilateral(
			double ox1 ,double oy1 ,double ox2 ,double oy2 ,double ox3 ,double oy3 ,double ox4 ,double oy4 
			,double nx1 ,double ny1 ,double nx2 ,double ny2 ,double nx3 ,double ny3 ,double nx4 ,double ny4 
			,double r ,double s
			,double& e11 ,double& e22 ,double& e12 ,double& e21,double& x,double& y);
	public :
		virtual void go(const shared_ptr<GeometricalModel>&, const shared_ptr<PhysicalParameters>&,bool);

	RENDERS(Quadrilateral);
	REGISTER_CLASS_NAME(GLDrawQuadrilateral);
	REGISTER_BASE_CLASS_NAME(GlGeometricalModelFunctor);
};

REGISTER_SERIALIZABLE(GLDrawQuadrilateral);



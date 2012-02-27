/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>
#include<yade/lib-computational-geometry/MarchingCube.hpp>
#include <vector>

class GLDrawLatticeSetGeometry : public GlGeometricalModelFunctor
{
	public :
		virtual void go(const shared_ptr<GeometricalModel>&, const shared_ptr<PhysicalParameters>&,bool);

	RENDERS(LatticeSetGeometry);
	REGISTER_CLASS_NAME(GLDrawLatticeSetGeometry);
	REGISTER_BASE_CLASS_NAME(GlGeometricalModelFunctor);

	private :
		MarchingCube mc;
		Vector3r min,max;
		int sizeX,sizeY,sizeZ;
		vector<vector<vector<float > > > 	scalarField,weights;
		void generateScalarField();
		void calcMinMax();
		float oldIsoValue,oldIsoSec,oldIsoThick;
		Vector3r isoStep;
		void drawLine(Vector3r a,Vector3r b,float color,float thickness);
		void drawPoint(Vector3r a,float color,float thickness);


};

REGISTER_SERIALIZABLE(GLDrawLatticeSetGeometry);



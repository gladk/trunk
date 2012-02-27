/*************************************************************************
*  Copyright (C) 2008 by Jerome Duriez                                   *
*  duriez@geo.hmg.inpg.fr                                                *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/FileGenerator.hpp>


typedef pair<Vector3r, Real> BasicSphere;
//! make a list of spheres non-overlapping sphere

using namespace std;


class SimpleShear : public FileGenerator
{
		void createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents);
		void createSphere(shared_ptr<Body>& body, Vector3r position, Real radius);
		void createActors(shared_ptr<Scene>& scene);
		//method to create a list (containing the positions of centers and radii) of n non interpenetrating spheres, occupying a rectangle with a given (rather high) porosity (issued from TriaxialTest) :
		string GenerateCloud(vector<BasicSphere>& sphere_list,Vector3r lowerCorner,Vector3r upperCorner,long number,Real rad_std_dev, Real porosity);
// 		to create the same list but by reading a text file containing the information :
		std::pair<string,bool> ImportCloud(vector<BasicSphere>& sphere_list,string importFilename);


	public :
		~SimpleShear ();
		bool generate(std::string& message);


	YADE_CLASS_BASE_DOC_ATTRS(SimpleShear,FileGenerator,"Preprocessor for creating a numerical model of a simple shear box.\n\n\t- Boxes (6) constitute the different sides of the box itself\n \t- Spheres are contained in the box. The sample could be generated via the same method used in TriaxialTest Preprocesor (=> see GenerateCloud) or by reading a text file containing positions and radii of a sample (=> see ImportCloud). This last one is the one by default used by this PreProcessor as it is written here => you need to have such a file.\n \n \t Thanks to the Engines (in pkg/common/Engine/PartialEngine) KinemCNDEngine, KinemCNSEngine and KinemCNLEngine, respectively constant normal displacement, constant normal rigidity and constant normal stress are possible to execute over such samples.\n \n \tNB about micro-parameters : their values correspond to those used in [Duriez2009a]_.",
				  ((string,filename,"../porosite0_44.txt",,"file with the list of spheres centers and radii"))
				  ((Vector3r,gravity,Vector3r(0,-9.81,0),,"vector corresponding to used gravity [$m/s^2$]"))
				  ((Real,thickness,0.001,,"thickness of the boxes constituting the shear box [$m$]"))
				  ((Real,length,0.1,,"initial length (along x-axis) of the shear box [$m$]"))
				  ((Real,height,0.02,,"initial height (along y-axis) of the shear box [$m$]"))
				  ((Real,width,0.04,,"initial width (along z-axis) of the shear box [$m$]"))
				  ((Real,density,2600,,"density of the spheres [$kg/m^3$]"))
				  ((Real,sphereYoungModulus,4.0e9,,"value of :yref:`ElastMat::young` for the spheres [$Pa$]"))
				  ((Real,spherePoissonRatio,0.04,,"value of :yref:`ElastMat::poisson` for the spheres [-]"))
				  ((Real,sphereFrictionDeg,37,,"value of :yref:`ElastMat::poisson` for the spheres [$^\\circ$] (the necessary conversion in $rad$ is done automatically)"))
				  ((Real,boxYoungModulus,4.0e9,,"value of :yref:`ElastMat::young` for the boxes [$Pa$]"))
				  ((Real,boxPoissonRatio,0.04,,"value of :yref:`ElastMat::poisson` for the spheres [-]"))
				  ((bool,gravApplied,false,,"depending on this, :yref:`GravityEngine` is added or not to the scene to take into account the weight of particles"))
				  ((int,timeStepUpdateInterval,50,,"value of :yref:`TimeStepper::timeStepUpdateInterval` for the :yref:`TimeStepper` used here"))
				  );
};

REGISTER_SERIALIZABLE(SimpleShear);




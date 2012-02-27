/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/FileGenerator.hpp>
#include<Wm3Vector3.h>
#include<yade/lib-base/yadeWm3.hpp>

#include<set>
#include<string>
#include<boost/foreach.hpp>
#include<boost/archive/binary_oarchive.hpp>
#include<boost/archive/binary_iarchive.hpp>
#include<boost/serialization/vector.hpp>
#include<boost/serialization/set.hpp>
#include<boost/serialization/list.hpp>
#include<boost/serialization/utility.hpp>

#include<boost/iostreams/filtering_stream.hpp>
#include<boost/iostreams/filter/gzip.hpp>
#include<boost/iostreams/filter/bzip2.hpp>
#include<boost/iostreams/device/file.hpp>
#include<boost/shared_ptr.hpp>

#include"Voxel/VoxelEnvelope.hpp"
#include<yade/pkg-snow/BshSnowGrain.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>

class VoxelDocument
{
	public:
		VoxelEnvelope	m_data;
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			//boost::mutex::scoped_lock scoped_lock(m_document_mutex);
			ar & m_data;
		}
};

class SnowVoxelsLoader : public FileGenerator
{
	private:
		Vector3r	 gravity
				,lowerCorner
				,upperCorner;

		Vector3r	 spheresColor;

		Real		 sphereYoungModulus
				,spherePoissonRatio
				,sphereFrictionDeg
				,boxYoungModulus
				,boxPoissonRatio
				,boxFrictionDeg
				,density
				
				,dampingForce
				,dampingMomentum

				,defaultDt

				,normalCohesion
				,shearCohesion
				
				,creep_viscosity
				,sigma_iso
				,thickness
				
				,strainRate
				,StabilityCriterion
				
				,maxMultiplier ///max multiplier of diameters during internal compaction
				;
		bool		 setCohesionOnNewContacts
				,autoCompressionActivation
				,internalCompaction
				,use_grain_shear_creep
				,use_grain_twist_creep
				,enable_layers_creep
				,use_gravity_engine
				;
		Real		 one_voxel_in_meters_is
				,layer_distance_voxels
				,angle_increment_radians
				,layers_creep_viscosity;

		int		 timeStepUpdateInterval
				,radiusControlInterval
				,wallStiffnessUpdateInterval

				,recordIntervalIter

				,skip_small_grains
				;
				
		std::string	WallStressRecordFile;

		std::string	voxel_binary_data_file;
		std::string	voxel_txt_dir;
		std::string	voxel_caxis_file;
		std::string	voxel_colors_file;
		std::string	grain_binary_data_file;
		
		shared_ptr<TriaxialCompressionEngine> triaxialcompressionEngine;
		shared_ptr<TriaxialStressController> triaxialstressController;
		shared_ptr<TriaxialStateRecorder> triaxialStateRecorder;

		VoxelDocument	m_voxel;
		std::vector<boost::shared_ptr<BshSnowGrain> > m_grains;

	public : 
		SnowVoxelsLoader();
		~SnowVoxelsLoader();
		virtual bool generate();
		bool load_voxels();
		void createActors(shared_ptr<Scene>& rootBody);
		void positionRootBody(shared_ptr<Scene>& rootBody);
		void create_grain(shared_ptr<Body>& body, Vector3r position, bool dynamic , boost::shared_ptr<BshSnowGrain> grain);
		void create_box(shared_ptr<Body>& body, Vector3r position, Vector3r extents, bool wire);
	
	REGISTER_ATTRIBUTES(FileGenerator,
		(voxel_binary_data_file)
		// that was for integrating snow-read with this preprocessor, but for now it's not integrated. snow-read does the conversion separately
		//	(voxel_txt_dir)
		//	(voxel_caxis_file)
		//	(voxel_colors_file)
	//	(grain_binary_data_file)
		(one_voxel_in_meters_is)
		(layer_distance_voxels)
		(angle_increment_radians)
		(skip_small_grains)

		(shearCohesion)
		(normalCohesion)
			  (creep_viscosity)
		(use_grain_shear_creep)
		(use_grain_twist_creep)
		(enable_layers_creep)
		(layers_creep_viscosity)
		(sigma_iso)
		(setCohesionOnNewContacts)

		(sphereYoungModulus)
		(spherePoissonRatio)
		(sphereFrictionDeg)
		(boxYoungModulus)
		(boxPoissonRatio)
		(boxFrictionDeg)
		(density)
		(use_gravity_engine)
		(gravity)
		(dampingForce)
		(dampingMomentum)

		(WallStressRecordFile)
	);
	REGISTER_CLASS_NAME(SnowVoxelsLoader);
	REGISTER_BASE_CLASS_NAME(FileGenerator);
};

REGISTER_SERIALIZABLE(SnowVoxelsLoader);


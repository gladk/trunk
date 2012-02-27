#pragma once

#include "Config.hpp"
#include "DataVoxel.hpp"
#include "DataSurface.hpp"
#include "SafeVectors3.hpp"

class VoxelEnvelope
{
	private:
		boost::mutex		m_voxel_envelope_mutex;
		bool			m_centers_calculated;
		bool			m_grain_voxel_surfaces_calculated;
		std::string		m_bad_grain_count;
		std::string		m_myname;

		DataVoxel		m_data_voxel;
		SafeVectors3		m_axes;
		SafeVectors3		m_colors;
		std::vector<Vector3r>	m_centers;
		DataSurface		m_data_surface;
	private:
		void			calculate_surfaces_loop(Config& c);
		void			calculate_surfaces_loop_spawn(Config& c);
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			boost::mutex::scoped_lock scoped_lock(m_voxel_envelope_mutex);
			ar & m_centers_calculated;
			ar & m_grain_voxel_surfaces_calculated;
			ar & m_bad_grain_count;
			ar & m_myname;
			ar & m_data_voxel;
			ar & m_axes;
			ar & m_colors;
			ar & m_centers;
			ar & m_data_surface;
		}
	public:
		VoxelEnvelope() : m_myname("unknown") {};
		VoxelEnvelope(Config& c,std::string name);
		void			load_experimental_data(Config& c);
		void			save_txt(Config& c);
		void			create_singled_particles(VoxelEnvelope& other,Config& config,std::string name);
		void			grow_them_a_little(VoxelEnvelope& other,VoxelEnvelope& source,Config& config,std::string name);
		void			calculate_difference(VoxelEnvelope& a,VoxelEnvelope& b,Config& config,std::string name);
		void			calculate_final_version(VoxelEnvelope& a,VoxelEnvelope& b,Config& config,std::string name);
		void			cut_off_all_tentacles(VoxelEnvelope& other,std::string name,Config& config);
		void			draw(Config& c);

		void			calculate_centers(Config& c);
		void			calculate_grain_voxel_surfaces(Config& c);
		Vector3r		center(Config& c);
		std::string		message();
		bool			centers_calculated();
		bool			grain_voxel_surfaces_calculated();
		bool			ready();
		const DataVoxel&	get_data_voxel_const_ref(){return m_data_voxel;};
		const SafeVectors3&	get_axes_const_ref()const{return m_axes;};
		const SafeVectors3&	get_colors_const_ref()const{return m_colors;};
};


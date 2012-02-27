#pragma once

#include "Config.hpp"
#include "SafeVectors3.hpp"
#include <vector>
#include <string>
#include <boost/thread.hpp>

typedef unsigned char				t_voxel;
typedef std::vector<std::vector<t_voxel> >	t_voxel_slice;
typedef std::vector<t_voxel_slice>		t_voxel_data;

class DataVoxel
{
	private:
		boost::mutex	m_voxel_mutex;
		bool		m_ready;
		t_voxel_data	m_voxel_data;
	
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			boost::mutex::scoped_lock lock(m_voxel_mutex);
			ar & m_ready;
			ar & m_voxel_data;
		}

	public:
		DataVoxel();
		void draw(SafeVectors3& axes,SafeVectors3& colors,Config& c);

		void load(Config& c);
		void save_txt_result(Config& c);
		void simple_volume_finder(const t_voxel_data& other,Config& config,std::set<int> grains);
		void simple_volume_grow(const t_voxel_data& other,const t_voxel_data& source);
		void diff(const t_voxel_data& a,const t_voxel_data& b);
		void final(const t_voxel_data& a,const t_voxel_data& b,Config& config);
		void cut_off_junk(const t_voxel_data& other,std::set<int> grains,const std::vector<Vector3r>& centers);
		void grow_balloon(Vector3r center, int grain_id,const t_voxel_data& other, std::vector<std::vector<std::vector<unsigned char> > >& box);
		void find_volumes(const t_voxel_data& data, std::vector<std::vector<std::vector<unsigned char> > >& box, int& how_many, int GRAIN_ID);
		int  volume_of_grain_id(const t_voxel_data& data, int GRAIN_ID);
		void read_txt(std::string file,t_voxel_slice& dat);

		bool ready();
		//void ready(bool t);
		t_voxel_data get_a_voxel_copy();
		const t_voxel_data& get_a_voxel_const_ref() const { return m_voxel_data;};

		std::pair<std::set<int>,std::string> how_many_grains();
};


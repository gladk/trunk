#pragma once

#include "Config.hpp"
#include "GrainSurface.hpp"
#include <vector>
#include <boost/thread.hpp>

class VoxelEnvelope;

//struct Contact
//{
//	Vector3r normal;
//	Vector3r center;
//	Vector3r color;
//	float count;
//};

class DataSurface
{
	private:
		mutable boost::mutex		m_surface_mutex;
		std::list<GrainSurface>		m_voxel_surfaces;
	
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			boost::mutex::scoped_lock scoped_lock(m_surface_mutex);
			ar & m_voxel_surfaces;
		}

	public:
		DataSurface();
		void create_surfaces(t_voxel_data dat, std::set<int>& grains,VoxelEnvelope* voxel_envelope);
		void draw(SafeVectors3& axes,SafeVectors3& colors,Config& c,const std::vector<Vector3r>& centers);
		DataSurface& operator= (const DataSurface & other)
		{
			boost::mutex::scoped_lock scoped_lock1(this->m_surface_mutex);
			boost::mutex::scoped_lock scoped_lock2(other.m_surface_mutex);
			this->m_voxel_surfaces = other.m_voxel_surfaces;
			return *this;
		}

};



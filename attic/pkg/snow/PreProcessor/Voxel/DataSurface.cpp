YADE_REQUIRE_FEATURE(COMPILED_BY_INCLUDE_FROM_SnowVoxelsLoader.cpp)
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "DataSurface.hpp"
#include<yade/pkg-snow/Config.hpp>
#include "VoxelEnvelope.hpp"
#include <unistd.h>

DataSurface::DataSurface()
{
	m_voxel_surfaces.clear();
}
		
void DataSurface::create_surfaces(t_voxel_data dat, std::set<int>& grains, VoxelEnvelope* voxel_envelope)
{
	sleep(5);
	if(grains.size() < 1)
	{
		std::cerr << "Not enough grains, yet.\n";
		return;
	}
	BOOST_FOREACH(int i, grains)
	{
		GrainSurface v(dat,i);
		{
			boost::mutex::scoped_lock scoped_lock(m_surface_mutex);
			m_voxel_surfaces.push_back(v);

			if(voxel_envelope->ready() && !voxel_envelope->centers_calculated())
			{
				std::cerr << "\nLoading was just finished, now we can proceed with final calculation\n";
				return;
			}
		}
	}
}

	const std::vector<Vector3r>*   s_local_centers;
	Vector3r                 s_local_camera;

	float point_plane_dist(Vector3r& point, Vector3r& plane_normal)
	{
		float a = plane_normal[0];
		float b = plane_normal[1];
		float c = plane_normal[2];

		float x0 = point[0];
		float y0 = point[1];
		float z0 = point[2];

		return ( a*x0 + b*y0 + c*z0 ); // (is already normalized) / std:sqrt(a*a + b*b + c*c)
	}

	bool depth_check(GrainSurface& first, GrainSurface& second)
	{
		t_voxel id1(first.id());
		t_voxel id2(second.id());

		Vector3r p1((*s_local_centers)[id1]);
		Vector3r p2((*s_local_centers)[id2]);

		return point_plane_dist(p1,s_local_camera) > point_plane_dist(p2,s_local_camera);
	}

inline Vector3r toVec(qglviewer::Vec v){return Vector3r(v[0],v[1],v[2]);};

void DataSurface::draw(SafeVectors3& axes,SafeVectors3& colors,Config& c,const std::vector<Vector3r>& centers)
{
	if(!c.draw_using_surfaces())
		return;

	boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	bool selected=c.selected();
	bool names=c.draw_with_names();
	boost::mutex::scoped_lock lock(m_surface_mutex);

	if(!names)
		c.apply_transparency();
	if(c.transparency() && (c.selection_size() > 1 || !c.display_selection()))
	{
		s_local_centers = &centers;
		s_local_camera  = toVec(c.view_direction());
		m_voxel_surfaces.sort(depth_check);
	}
	BOOST_FOREACH(GrainSurface& v , m_voxel_surfaces)
	{
		t_voxel gr = v.id();
		if( (selected && c.display_selection()) ? c.selected(gr) : gr!=0 )
		{
			if(names) glPushName(gr);
			v.draw(colors,c);
			if(names) glPopName();
		}
		if(c.time_draw_limit() &&
		   boost::posix_time::time_duration(boost::posix_time::microsec_clock::local_time() - start) > boost::posix_time::milliseconds(400))
			return;
	}
	if(!names)
		c.remove_transparency();
}


YADE_REQUIRE_FEATURE(PHYSPAR);


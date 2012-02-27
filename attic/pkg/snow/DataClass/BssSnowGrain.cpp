#include"BshSnowGrain.hpp"
#include"BssSnowGrain.hpp"
#include<Wm3Quaternion.h>

YADE_REQUIRE_FEATURE(geometricalmodel);

BssSnowGrain::BssSnowGrain():Sphere()
{
	createIndex();
};

BssSnowGrain::BssSnowGrain(BshSnowGrain* grain, Real one_voxel_in_meters_is)
{
	createIndex();
	
	m_copy.center   =grain->center;
	m_copy.c_axis   =grain->c_axis;
	m_copy.start    =grain->start;
	m_copy.end      =grain->end;
	m_copy.color    =grain->color;
	m_copy.selection=grain->selection;
	m_copy.slices   =grain->slices;
	m_copy.layer_distance = grain->layer_distance;

//	Real LEN=(m_copy.start - m_copy.end).Length();
//	std::cerr << LEN/one_voxel_in_meters_is << " ---------------\n";

	// calculate volume
	Real area=0;
	Real vol=0;
	Vector3r prev(0,0,0);
	Vector3r pos(m_copy.start);
	BOOST_FOREACH(std::vector<Vector3r>& g,m_copy.slices)
	{
		BOOST_FOREACH(Vector3r& v,g)
		{
			if(prev != Vector3r(0,0,0))
			{
				Real a = (prev - v).Length();
				Real h = (((prev + v)*0.5) - pos).Length();
				area += a*h*0.5;
			}
			prev = v;
		}
		vol += area * one_voxel_in_meters_is * m_copy.layer_distance;
		area = 0;
//		std::cerr << one_voxel_in_meters_is << " " << m_copy.layer_distance << " " << m_copy.c_axis << "\n";
//		std::cerr << ((pos-m_copy.end).Length()/one_voxel_in_meters_is)  << " - " << ((one_voxel_in_meters_is * m_copy.layer_distance * m_copy.c_axis)/one_voxel_in_meters_is)  << " ?\n";
		pos = pos - one_voxel_in_meters_is * m_copy.layer_distance * m_copy.c_axis;
//		std::cerr << ((pos-m_copy.end).Length()/one_voxel_in_meters_is) << "\n";
	}

	// now we can get radius from volume
	//vol/(4.0/3.0*Mathr::PI)    = radius*radius*radius;
	radius = std::pow(vol/(4.0/3.0*Mathr::PI) , 0.333333333333333333333 );

	std::cerr << "radius: " << radius << "\n";

}


YADE_PLUGIN((BssSnowGrain));

YADE_REQUIRE_FEATURE(PHYSPAR);


#pragma once

#include<yade/core/GeometricalModel.hpp>
#include<Wm3Vector3.h>
#include<yade/lib-base/yadeWm3.hpp>
#include<vector>
#include<boost/serialization/vector.hpp>
#include<boost/serialization/shared_ptr.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-snow/BshSnowGrain.hpp>
#include<boost/tuple/tuple.hpp>

struct depth_one
{
	mutable Real current_depth;
	Real original_depth;

	int i;
	int j;

	depth_one(Real D,int I,int J):current_depth(D),original_depth(D),i(I),j(J){};
	bool operator< (const depth_one& b) const
	{
		/* BAD, WRONG
		if(i < b.i);
			return true;
		if(i == b.i && j < b.j);
			return true;
		return false;
		*/
		// GOOD, CORRECT

		if (i < b.i) return true;
		if (i > b.i) return false;

		if (j < b.j) return true;
		if (j > b.j) return false;

		return false;
	};
	void set_current_depth(Real d) const
	{
		current_depth=d;
	};
};

class BssSnowGrain : public Sphere
{
	public: 
		std::map<int,std::set<depth_one> >	depths; // body_id , < depth, i, j (indexes of point on a slice) >
		std::map<int,Real>	sphere_depth;
		BshSnowGrain	m_copy;

		BssSnowGrain();
		BssSnowGrain(BshSnowGrain* grain, Real one_voxel_in_meters_is);
	
	// FIXME: BssSnowGrain has no copy-constructor; will be inaccessible from python
	// REGISTER_ATTRIBUTES(Sphere,(m_copy));
	virtual void registerAttributes(){ Sphere::registerAttributes(); REGISTER_ATTRIBUTE_(m_copy); }
	REGISTER_CLASS_NAME(BssSnowGrain);
	REGISTER_BASE_CLASS_NAME(Sphere);
	REGISTER_CLASS_INDEX(BssSnowGrain,Sphere);
};

REGISTER_SERIALIZABLE(BssSnowGrain);


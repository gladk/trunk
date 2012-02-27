YADE_REQUIRE_FEATURE(COMPILED_BY_INCLUDE_FROM_SnowVoxelsLoader.cpp)
#include "SafeVectors3.hpp"

SafeVectors3::SafeVectors3(int size,std::string myname,std::string fname,bool has_id) : m_myname(myname)
{
	m_data.resize(size+1);
	std::ifstream file(fname.c_str());
	if(!file)
		std::cerr << "\n" << myname << " ERROR: file " << fname <<" could not be opened for reading\n", exit(0);
	float id,x,y,z;
	for(int j=has_id?1:0;j<size+1;++j)
	{
		if(has_id) 
			file >> id;
		file >> x >> y >> z;
		if(has_id && id!=(float)j)
			std::cerr << "\n" << myname << " ERROR: bad grain number: " << id << " should be " << j << "\n", exit(0);
		m_data[j]=Vector3r(x,y,z);
	}
	std::cerr << myname << " file loaded\n";
}

Vector3r SafeVectors3::operator[](size_t i)const
{
	if(i>=0 && i<m_data.size()) {
		return m_data[i];
	} else {
		std::cerr << m_myname << " WARNING: trying to get Vector3r outside the container: [" << i << "], returning [0]\n";
		return m_data[0];
	}
}


YADE_REQUIRE_FEATURE(PHYSPAR);


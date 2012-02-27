#pragma once

#include<yade/lib-base/yadeWm3Extra.hpp>
#include<Wm3Vector3.h>
#include<Wm3Quaternion.h>
#include<vector>
#include<boost/archive/binary_oarchive.hpp>
#include<boost/archive/binary_iarchive.hpp>
//#include<boost/serialization/vector.hpp>

class SafeVectors3
{
	private:
		std::vector<Vector3r>	m_data;
		std::string		m_myname;
	private:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & m_data;
			ar & m_myname;
		}
	public:
		SafeVectors3() : m_myname("undefined") {};
		SafeVectors3(int size,std::string myname,std::string fname,bool has_id);
		SafeVectors3(int size,std::string myname):m_myname(myname){m_data.resize(size+1,Vector3r(0,0,0));};
		Vector3r operator[](size_t i)const;
};


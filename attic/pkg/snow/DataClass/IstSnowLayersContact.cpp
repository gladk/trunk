#include<yade/pkg-snow/IstSnowLayersContact.hpp>

IstSnowLayersContact::IstSnowLayersContact()
{
	createIndex();
}

IstSnowLayersContact::~IstSnowLayersContact()
{}
		

YADE_REQUIRE_FEATURE(PHYSPAR);


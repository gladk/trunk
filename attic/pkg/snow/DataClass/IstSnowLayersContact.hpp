// 2009 Janek Kozicki <cosurgi@berlios.de>

#pragma once

#include<yade/core/InteractionGeometry.hpp>
#include<yade/pkg-dem/ScGeom.hpp>

class IstSnowLayersContact : public ScGeom
{
	public:
		IstSnowLayersContact();
		virtual ~IstSnowLayersContact();

	REGISTER_ATTRIBUTES(ScGeom,/* */);
	REGISTER_CLASS_NAME(IstSnowLayersContact);
	REGISTER_BASE_CLASS_NAME(ScGeom);
	REGISTER_CLASS_INDEX(IstSnowLayersContact,ScGeom);
};

REGISTER_SERIALIZABLE(IstSnowLayersContact);


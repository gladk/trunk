// © 2009 Janek Kozicki <cosurgi@mail.berlios.de>

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>

class Ef1_IstSnowLayersContact_glDraw : public GlInteractionGeometryFunctor
{
	public :
		Ef1_IstSnowLayersContact_glDraw(){}
		virtual void go(const shared_ptr<InteractionGeometry>&,const shared_ptr<Interaction>&,const shared_ptr<Body>&,const shared_ptr<Body>&,bool wireFrame);

	RENDERS(IstSnowLayersContact);
	REGISTER_CLASS_NAME(Ef1_IstSnowLayersContact_glDraw);
	REGISTER_BASE_CLASS_NAME(GlInteractionGeometryFunctor);
};

REGISTER_SERIALIZABLE(Ef1_IstSnowLayersContact_glDraw);


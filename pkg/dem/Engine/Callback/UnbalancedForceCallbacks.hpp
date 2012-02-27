// 2010 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include<yade/pkg-common/Callbacks.hpp>
#include<yade/lib-base/openmp-accu.hpp>

class SumIntrForcesCb: public IntrCallback{
	public:
		// zero values, so that we can pass pointers to them to OpenMPAccumulator
		static int int0;
		static Real Real0;
		OpenMPAccumulator<int,&SumIntrForcesCb::int0> numIntr;
		OpenMPAccumulator<Real,&SumIntrForcesCb::Real0> force;
		static void go(IntrCallback*,Interaction*);
		virtual IntrCallback::FuncPtr stepInit();
	YADE_CLASS_BASE_DOC(SumIntrForcesCb,IntrCallback,"Callback summing magnitudes of forces over all interactions. :yref:`InteractionPhysics` of interactions must derive from :yref:`NormShearPhys` (responsability fo the user).");
};
REGISTER_SERIALIZABLE(SumIntrForcesCb);

class SumBodyForcesCb: public BodyCallback{
	Scene* scene;
	public:
		OpenMPAccumulator<int,&SumIntrForcesCb::int0> numBodies;
		OpenMPAccumulator<Real,&SumIntrForcesCb::Real0> force;
		static void go(BodyCallback*,Body*);
		virtual BodyCallback::FuncPtr stepInit();
	YADE_CLASS_BASE_DOC(SumBodyForcesCb,BodyCallback,"Callback summing magnitudes of resultant forces over :yref:`dynamic<Body::isDynamic>` bodies.");
};
REGISTER_SERIALIZABLE(SumBodyForcesCb);

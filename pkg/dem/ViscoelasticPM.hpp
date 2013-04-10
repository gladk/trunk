// 2009 © Sergei Dorofeenko <sega@users.berlios.de>
// This file contains a set of classes for modelling of viscoelastic
// particles.

#pragma once

#include<yade/core/Material.hpp>
#include<yade/pkg/dem/FrictPhys.hpp>
#include<yade/pkg/common/Dispatching.hpp>
#include<yade/pkg/dem/ScGeom.hpp>

/* Simple viscoelastic model */

/// Material
/// Note: Shop::getViscoelasticFromSpheresInteraction can get kn,cn,ks,cs from a analytical solution of a pair spheres interaction problem.
class ViscElMat : public Material {
	public:
		virtual ~ViscElMat();
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(ViscElMat,Material,"Material for simple viscoelastic model of contact.\n\n.. note::\n\t ``Shop::getViscoelasticFromSpheresInteraction`` (and :yref:`yade.utils.getViscoelasticFromSpheresInteraction` in python) compute :yref:`kn<ViscElMat::kn>`, :yref:`cn<ViscElMat::cn>`,  :yref:`ks<ViscElMat::ks>`,  :yref:`cs<ViscElMat::cs>` from analytical solution of a pair spheres interaction problem.",
		((Real,kn,NaN,,"Normal elastic stiffness"))
		((Real,cn,NaN,,"Normal viscous constant"))
		((Real,ks,NaN,,"Shear elastic stiffness"))
		((Real,cs,NaN,,"Shear viscous constant"))
		((Real,frictionAngle,NaN,,"Friction angle [rad]"))
		((bool,Capillar,false,,"True, if capillar forces need to be added."))
		((Real,Vb,NaN,,"Liquid bridge volume [m^3]"))
		((Real,gamma,NaN,,"Surface tension [N/m]"))
		((Real,theta,NaN,,"Contact angle [°]"))
		((std::string,CapillarType,"",,"Different types of capillar interaction: Willett, Weigert, Herminghaus")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(ViscElMat,Material);
};
REGISTER_SERIALIZABLE(ViscElMat);

/// Interaction physics
class ViscElPhys : public FrictPhys{
	public:
		virtual ~ViscElPhys();
		Real R;
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(ViscElPhys,FrictPhys,"IPhys created from :yref:`ViscElMat`, for use with :yref:`Law2_ScGeom_ViscElPhys_Basic`.",
		((Real,cn,NaN,,"Normal viscous constant"))
		((Real,cs,NaN,,"Shear viscous constant"))
		((bool,Capillar,false,,"True, if capillar forces need to be added."))
		((bool,liqBridgeCreated,false,,"Whether liquid bridge was created, only after a normal contact of spheres"))
		((Real,sCrit,false,,"Critical bridge length [m]"))
		((Real,Vb,NaN,,"Liquid bridge volume [m^3]"))
		((Real,gamma,NaN,,"Surface tension [N/m]"))
		((Real,theta,NaN,,"Contact angle [rad]"))
		((std::string,CapillarType,"",,"Different types of capillar interaction: Willett, Weigert, Herminghaus")),
		createIndex();
	)
};
REGISTER_SERIALIZABLE(ViscElPhys);

/// Convert material to interaction physics.
// Uses the rule of consecutively connection.
class Ip2_ViscElMat_ViscElMat_ViscElPhys: public IPhysFunctor {
	public :
		virtual void go(const shared_ptr<Material>& b1,
					const shared_ptr<Material>& b2,
					const shared_ptr<Interaction>& interaction);
	YADE_CLASS_BASE_DOC(Ip2_ViscElMat_ViscElMat_ViscElPhys,IPhysFunctor,"Convert 2 instances of :yref:`ViscElMat` to :yref:`ViscElPhys` using the rule of consecutive connection.");
	FUNCTOR2D(ViscElMat,ViscElMat);

};
REGISTER_SERIALIZABLE(Ip2_ViscElMat_ViscElMat_ViscElPhys);

/// Constitutive law
/// This class provides linear viscoelastic contact model
class Law2_ScGeom_ViscElPhys_Basic: public LawFunctor {
	public :
		virtual void go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
	FUNCTOR2D(ScGeom,ViscElPhys);
	YADE_CLASS_BASE_DOC(Law2_ScGeom_ViscElPhys_Basic,LawFunctor,"Linear viscoelastic model operating on :yref:`ScGeom` and :yref:`ViscElPhys`.");
};
REGISTER_SERIALIZABLE(Law2_ScGeom_ViscElPhys_Basic);

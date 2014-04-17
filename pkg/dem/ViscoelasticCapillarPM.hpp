#include"ViscoelasticPM.hpp"
#include <boost/unordered_map.hpp>

class ViscElCapMat : public ViscElMat {
	public:
		virtual ~ViscElCapMat();
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(ViscElCapMat,ViscElMat,"Material for extended viscoelastic model of contact with capillary parameters.",
		((bool,Capillar,false,,"True, if capillar forces need to be added."))
		((Real,Vb,NaN,,"Liquid bridge volume [m^3]"))
		((Real,gamma,NaN,,"Surface tension [N/m]"))
		((Real,theta,NaN,,"Contact angle [°]"))
		((std::string,CapillarType,"",,"Different types of capillar interaction: Willett_numeric, Willett_analytic [Willett2000]_ , Weigert [Weigert1999]_ , Rabinovich [Rabinov2005]_ , Lambert (simplified, corrected Rabinovich model) [Lambert2008]_ ")),
		createIndex();
	);
	REGISTER_CLASS_INDEX(ViscElCapMat,ViscElMat);
};
REGISTER_SERIALIZABLE(ViscElCapMat);

/// Interaction physics
enum CapType {None_Capillar, Willett_numeric, Willett_analytic, Weigert, Rabinovich, Lambert, Soulie};
class ViscElCapPhys : public ViscElPhys{
	typedef Real (* CapillarFunction)(const ScGeom& geom, ViscElCapPhys& phys);
	public:
		virtual ~ViscElCapPhys();
		Real R;
		CapillarFunction CapFunct;
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(ViscElCapPhys,ViscElPhys,"IPhys created from :yref:`ViscElCapMat`, for use with :yref:`Law2_ScGeom_ViscElCapPhys_Basic`.",
		((bool,Capillar,false,,"True, if capillar forces need to be added."))
		((bool,liqBridgeCreated,false,,"Whether liquid bridge was created, only after a normal contact of spheres"))
		((bool,liqBridgeActive,false,, "Whether liquid bridge is active at the moment"))
		((Real,sCrit,false,,"Critical bridge length [m]"))
		((Real,Vb,NaN,,"Liquid bridge volume [m^3]"))
		((Real,gamma,NaN,,"Surface tension [N/m]"))
		((Real,theta,NaN,,"Contact angle [rad]"))
		((CapType,CapillarType,None_Capillar,,"Different types of capillar interaction: Willett_numeric, Willett_analytic, Weigert, Rabinovich, Lambert, Soulie"))
#ifdef YADE_LIQCONTROL
		((Real,Vmax,NaN,,"Maximal liquid bridge volume [m^3]"))
#endif
		,
		createIndex();
	)
	REGISTER_CLASS_INDEX(ViscElCapPhys,ViscElPhys);
};
REGISTER_SERIALIZABLE(ViscElCapPhys);

/// Convert material to interaction physics.
class Ip2_ViscElCapMat_ViscElCapMat_ViscElCapPhys: public Ip2_ViscElMat_ViscElMat_ViscElPhys {
	public :
		virtual void go(const shared_ptr<Material>& b1,
					const shared_ptr<Material>& b2,
					const shared_ptr<Interaction>& interaction);
	YADE_CLASS_BASE_DOC(Ip2_ViscElCapMat_ViscElCapMat_ViscElCapPhys,Ip2_ViscElMat_ViscElMat_ViscElPhys,"Convert 2 instances of :yref:`ViscElCapMat` to :yref:`ViscElCapPhys` using the rule of consecutive connection.");
	FUNCTOR2D(ViscElCapMat,ViscElCapMat);
};
REGISTER_SERIALIZABLE(Ip2_ViscElCapMat_ViscElCapMat_ViscElCapPhys);

/// Constitutive law
class Law2_ScGeom_ViscElCapPhys_Basic: public LawFunctor {
	public :
		virtual void go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*);
		static Real Willett_numeric_f     (const ScGeom& geom, ViscElCapPhys& phys);
		static Real Willett_analytic_f    (const ScGeom& geom, ViscElCapPhys& phys);
		static Real Weigert_f             (const ScGeom& geom, ViscElCapPhys& phys);
		static Real Rabinovich_f          (const ScGeom& geom, ViscElCapPhys& phys);
		static Real Lambert_f             (const ScGeom& geom, ViscElCapPhys& phys);
		static Real Soulie_f              (const ScGeom& geom, ViscElCapPhys& phys);
		static Real None_f                (const ScGeom& geom, ViscElCapPhys& phys);
		Real critDist(const Real& Vb, const Real& R, const Real& Theta);
	FUNCTOR2D(ScGeom,ViscElCapPhys);
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom_ViscElCapPhys_Basic,LawFunctor,"Extended version of Linear viscoelastic model with capillary parameters.",
		((OpenMPAccumulator<Real>,VLiqBridg,,Attr::noSave,"The total volume of liquid bridges"))
		((OpenMPAccumulator<int>, NLiqBridg,,Attr::noSave,"The total number of liquid bridges"))
		,/* ctor */
		,/* py */
		;
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom_ViscElCapPhys_Basic);

#ifdef YADE_LIQCONTROL
typedef boost::unordered_map<Body::id_t, int> mapBodyInt;
class LiqControl: public PartialEngine{
	public:
		virtual void action();
		void addBodyMapInt( mapBodyInt & m, Body::id_t b );
		Real vMax(shared_ptr<Body> b1, shared_ptr<Body> b2);
		void updateLiquid(shared_ptr<Body> b);
	YADE_CLASS_BASE_DOC_ATTRS(LiqControl,PartialEngine,"Apply given torque (momentum) value at every subscribed particle, at every step. ",
		((int, mask,-1,, "Bitmask for particles."))
  );
};

REGISTER_SERIALIZABLE(LiqControl);
#endif

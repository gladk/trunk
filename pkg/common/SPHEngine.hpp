#ifdef YADE_SPH
#pragma once

#include<yade/core/PartialEngine.hpp>
#include<yade/pkg/dem/ScGeom.hpp>

typedef Real (* KernelFunction)(const double & r, const double & h);

enum KernFunctions {Poly6=1, Spiky=2, Visco=3, Lucy=4, Monaghan=5};
#define KERNELFUNCDESCR throw runtime_error("Type of kernel function undefined! The following kernel functions are available: Poly6=1, Spiky=2, Visco=3, Lucy=4, Monaghan=5.");

enum typeKernFunctions {Norm, Grad, Lapl};
class SPHEngine: public PartialEngine{
  public:
    void calculateSPHRho(const shared_ptr<Body>& b);
    virtual void action();
  YADE_CLASS_BASE_DOC_ATTRS(SPHEngine,PartialEngine,"Apply given torque (momentum) value at every subscribed particle, at every step. ",
    ((int, mask,-1,, "Bitmask for SPH-particles."))
    ((Real,k,-1,,    "Gas constant for SPH-interactions (only for SPH-model). See Mueller [Mueller2003]_ .")) // [Mueller2003], (11)
    ((Real,rho0,-1,, "Rest density. See Mueller [Mueller2003]_ ."))                                           // [Mueller2003], (1)
    ((Real,h,-1,,    "Core radius. See Mueller [Mueller2003]_ ."))                                            // [Mueller2003], (1)
    ((int,KernFunctionDensity, Poly6,, "Kernel function for density calculation (by default - Poly6). The following kernel functions are available: Poly6=1, Spiky=2, Visco=3, Lucy=4, Monaghan=5."))
  );
};
REGISTER_SERIALIZABLE(SPHEngine);
Real smoothkernelPoly6(const double & r, const double & h);       // [Mueller2003] (20)
Real smoothkernelPoly6Grad(const double & r, const double & h);
Real smoothkernelPoly6Lapl(const double & r, const double & h);
Real smoothkernelSpiky(const double & r, const double & h);       // [Mueller2003] (21)
Real smoothkernelSpikyGrad(const double & r, const double & h);
Real smoothkernelSpikyLapl(const double & r, const double & h);
Real smoothkernelVisco(const double & r, const double & h);       // [Mueller2003] (22)
Real smoothkernelViscoGrad(const double & r, const double & h);   // [Mueller2003] (22)
Real smoothkernelViscoLapl(const double & r, const double & h);
Real smoothkernelLucy(const double & r, const double & h);         
Real smoothkernelLucyGrad(const double & r, const double & h);     
Real smoothkernelLucyLapl(const double & r, const double & h);
Real smoothkernelMonaghan(const double & r, const double & h);         
Real smoothkernelMonaghanGrad(const double & r, const double & h);     
Real smoothkernelMonaghanLapl(const double & r, const double & h);

KernelFunction returnKernelFunction(const int a, const int b, const typeKernFunctions typeF);
KernelFunction returnKernelFunction(const int a, const typeKernFunctions typeF);

bool computeForceSPH(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I, Vector3r & force);
#endif


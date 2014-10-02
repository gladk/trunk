#ifdef YADE_SPH
#include"SPHEngine.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg/dem/ViscoelasticPM.hpp>
#include<yade/pkg/common/Sphere.hpp>

#include<yade/core/State.hpp>
#include<yade/core/Omega.hpp>

void SPHEngine::action(){
  {
    YADE_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){
      if(mask>0 && (b->groupMask & mask)==0) continue;
      this->calculateSPHRho(b);
      b->press=std::max(0.0, k*(b->rho - b->rho0));
    } YADE_PARALLEL_FOREACH_BODY_END();
  }
  /*
  {
    YADE_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies){
      if(mask>0 && (b->groupMask & mask)==0) continue;
      this->calculateSPHCs(b);
    } YADE_PARALLEL_FOREACH_BODY_END();
  }
  */
}

void SPHEngine::calculateSPHRho(const shared_ptr<Body>& b) {
  if (b->rho0<0) {
    b->rho0 = rho0;
  }
  Real rho = 0;
  
  // Pointer to kernel function
  KernelFunction kernelFunctionCurDensity = returnKernelFunction (KernFunctionDensity, KernFunctionDensity, Norm);
  
  // Calculate rho for every particle
  for(Body::MapId2IntrT::iterator it=b->intrs.begin(),end=b->intrs.end(); it!=end; ++it) {
    const shared_ptr<Body> b2 = Body::byId((*it).first,scene);
    Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
    if(!s) continue;
    
    if (((*it).second)->geom and ((*it).second)->phys) {
      const ScGeom geom = *(YADE_PTR_CAST<ScGeom>(((*it).second)->geom));
      const ViscElPhys phys=*(YADE_PTR_CAST<ViscElPhys>(((*it).second)->phys));
      
      if((b2->groupMask & mask)==0)  continue;
      
      Real Mass = b2->state->mass;
      if (Mass == 0) Mass = b->state->mass;
      
      const Real SmoothDist = -geom.penetrationDepth + phys.h;
     
      // [Mueller2003], (3)
      rho += b2->state->mass*kernelFunctionCurDensity(SmoothDist, phys.h);
    }
    // [Mueller2003], (3), we need to consider the density of the current body (?)
    // rho += b->state->mass*kernelFunctionCurDensity(0.0, s->radius);
  }
  b->rho = rho;
}

void SPHEngine::calculateSPHCs(const shared_ptr<Body>& b) {
  if (b->Cs<0) {
    b->Cs = 0.0;
  }
  Real Cs = 0;
  
  // Pointer to kernel function
  KernelFunction kernelFunctionCurDensity = returnKernelFunction (KernFunctionDensity, KernFunctionDensity, Norm);
  
  // Calculate Cs for every particle
  for(Body::MapId2IntrT::iterator it=b->intrs.begin(),end=b->intrs.end(); it!=end; ++it) {
    const shared_ptr<Body> b2 = Body::byId((*it).first,scene);
    Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
    if(!s) continue;
    
    if (((*it).second)->geom and ((*it).second)->phys) {
      const ScGeom geom = *(YADE_PTR_CAST<ScGeom>(((*it).second)->geom));
      const ViscElPhys phys=*(YADE_PTR_CAST<ViscElPhys>(((*it).second)->phys));
      
      if((b2->groupMask & mask)==0)  continue;
      
      Real Mass = b2->state->mass;
      if (Mass == 0) Mass = b->state->mass;
      
      const Real SmoothDist = -geom.penetrationDepth + phys.h;
     
      // [Mueller2003], (15)
      Cs += b2->state->mass/b2->rho*kernelFunctionCurDensity(SmoothDist, phys.h);
    }
  }
  b->Cs = Cs;
}

Real smoothkernelPoly6(const double & r, const double & h) {
  if (r<=h) {
    return 315/(64*M_PI*pow(h,9))*pow((h*h - r*r), 3);
  } else {
    return 0;
  }
}

Real smoothkernelPoly6Grad(const double & r, const double & h) {
  if (r<=h) {
    return -315/(64*M_PI*pow(h,9))*(-6*r*pow((h*h-r*r), 2));
  } else {
    return 0;
  }
}

Real smoothkernelPoly6Lapl(const double & r, const double & h) {
  if (r<=h) {
    return 315/(64*M_PI*pow(h,9))*(-6*(h*h-r*r)*(h*h - 5*r*r));
  } else {
    return 0;
  }
}

Real smoothkernelSpiky(const double & r, const double & h) {
  if (r<=h) {
    return 15/(M_PI*pow(h,6))*(pow((h-r), 3));             // [Mueller2003], (21)
  } else {
    return 0;
  }
}

Real smoothkernelSpikyGrad(const double & r, const double & h) {
  if (r<=h) {
    return -15/(M_PI*pow(h,6))*(-3*pow((h-r),2));
  } else {
    return 0;
  }
}

Real smoothkernelSpikyLapl(const double & r, const double & h) {
  if (r<=h) {
    return 15/(M_PI*pow(h,6))*(6*(h-r));
  } else {
    return 0;
  }
}

Real smoothkernelVisco(const double & r, const double & h) {
  if (r<=h and r!=0 and h!=0) {
    return 15/(2*M_PI*pow(h,3))*(-r*r*r/(2*h*h*h) + r*r/(h*h) + h/(2*r) -1);   // [Mueller2003], (21)
  } else {
    return 0;
  }
}

Real smoothkernelViscoGrad(const double & r, const double & h) {
  if (r<=h and r!=0 and h!=0) {
    return -15/(2*M_PI*pow(h,3))*(-3*r*r/(2*h*h*h) + 2*r/(h*h) - h/(2*r*r));
  } else {
    return 0;
  }
}

Real smoothkernelViscoLapl(const double & r, const double & h) {
  if (r<=h and r!=0 and h!=0) {
    return 45/(M_PI*pow(h,6))*(h - r);                     // [Mueller2003], (22+)
  } else {
    return 0;
  }
}

Real smoothkernelLucy(const double & r, const double & h) {
  if (r<=h and r!=0 and h!=0) {
    return 5/(9*M_PI)*(1+3*r)*pow((1-r),3);
  } else {
    return 0;
  }
}

Real smoothkernelLucyGrad(const double & r, const double & h) {
  if (r<=h and r!=0 and h!=0) {
    return -5/(9*M_PI*pow(h,2))*(-12*r/(h*h))*pow((1-r/h),2);
  } else {
    return 0;
  }
}

Real smoothkernelLucyLapl(const double & r, const double & h) {
  if (r<=h and r!=0 and h!=0) {
    return  5/(9*M_PI*pow(h,2))*(-12/(h*h))*(1-r/h)*(1-3*r/h);
  } else {
    return 0;
  }
}

Real smoothkernelMonaghan(const double & r, const double & h) {
  Real ret = 0.0;
  if (h!=0) {
    if (r/h<0.5) {
      ret = 40/(7*M_PI*h*h)*(1 - 6*pow((r/h),2) + 6*pow((r/h),3));
    } else {
      ret = 80/(7*M_PI*h*h)*pow((1 - (r/h)), 3);
    }
  }
  return ret;
}

Real smoothkernelMonaghanGrad(const double & r, const double & h) {
  Real ret = 0.0;
  if (h!=0) {
    if (r/h<0.5) {
      ret = -40/(7*M_PI*h*h)*( - 6*r/(h*h))*(2 - 3 * r/(h*h*h));
    } else {
      ret = -80/(7*M_PI*h*h)*( -3/h)*pow((1 -r/h), 2);
    }
  }
  return ret;
}

Real smoothkernelMonaghanLapl(const double & r, const double & h) {
  Real ret = 0.0;
  if (h!=0) {
    if (r/h<0.5) {
      ret = 40/(7*M_PI*h*h)*( - 12/(h*h))*(1 - 3 * r/(h*h*h*h*h));
    } else {
      ret = 80/(7*M_PI*h*h)*( 6/(h*h))*(1 -r/h);
    }
  }
  return ret;
}

KernelFunction returnKernelFunction(const int a, const int b, const typeKernFunctions typeF) {
  if (a != b) {
    throw runtime_error("Kernel types should be equal!");
  }
  if (a==Poly6) {
    if (typeF==Norm) {
      return smoothkernelPoly6;
    } else if (typeF==Grad) {
      return smoothkernelPoly6Grad;
    } else if (typeF==Lapl) {
      return smoothkernelPoly6Lapl;
    } else {
      KERNELFUNCDESCR
    }
  } else if (a==Spiky) {
    if (typeF==Norm) {
      return smoothkernelSpiky;
    } else if (typeF==Grad) {
      return smoothkernelSpikyGrad;
    } else if (typeF==Lapl) {
      return smoothkernelSpikyLapl;
    } else {
      KERNELFUNCDESCR
    }
  } else if (a==Visco) {
    if (typeF==Norm) {
      return smoothkernelVisco;
    } else if (typeF==Grad) {
      return smoothkernelViscoGrad;
    } else if (typeF==Lapl) {
      return smoothkernelViscoLapl;
    } else {
    }
  } else if (a==Lucy) {
    if (typeF==Norm) {
      return smoothkernelLucy;
    } else if (typeF==Grad) {
      return smoothkernelLucyGrad;
    } else if (typeF==Lapl) {
      return smoothkernelLucyLapl;
    } else {
      KERNELFUNCDESCR
    }
  } else if (a==Monaghan) {
    if (typeF==Norm) {
      return smoothkernelMonaghan;
    } else if (typeF==Grad) {
      return smoothkernelMonaghanGrad;
    } else if (typeF==Lapl) {
      return smoothkernelMonaghanLapl;
    } else {
      KERNELFUNCDESCR
    }
  } else {
    KERNELFUNCDESCR
  }
}

bool computeForceSPH(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I, Vector3r & force) {
  
  const ScGeom& geom=*static_cast<ScGeom*>(_geom.get());
  Scene* scene=Omega::instance().getScene().get();
  ViscElPhys& phys=*static_cast<ViscElPhys*>(_phys.get());
  
  const int id1 = I->getId1();
  const int id2 = I->getId2();
  
  const BodyContainer& bodies = *scene->bodies;
  
  //////////////////////////////////////////////////////////////////
  // Copy-paste
  
  const State& de1 = *static_cast<State*>(bodies[id1]->state.get());
  const State& de2 = *static_cast<State*>(bodies[id2]->state.get());
  
    // Handle periodicity.
  const Vector3r shift2 = scene->isPeriodic ? scene->cell->intrShiftPos(I->cellDist): Vector3r::Zero(); 
  const Vector3r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(I->cellDist): Vector3r::Zero(); 

  const Vector3r c1x = (geom.contactPoint - de1.pos);
  const Vector3r c2x = (geom.contactPoint - de2.pos - shift2);
  
  const Vector3r relativeVelocity = (de1.vel+de1.angVel.cross(c1x)) - (de2.vel+de2.angVel.cross(c2x)) + shiftVel;
  const Real normalVelocity	= geom.normal.dot(relativeVelocity);
  //const Vector3r shearVelocity	= relativeVelocity-normalVelocity*geom.normal;
  
  // Copy-paste
  //////////////////////////////////////////////////////////////////
  
  if (phys.h<0) {
    Sphere* s1=dynamic_cast<Sphere*>(bodies[id1]->shape.get());
    Sphere* s2=dynamic_cast<Sphere*>(bodies[id2]->shape.get());
    if (s1 and s2) {
      phys.h = s1->radius + s2->radius;
    } else if (s1 and not(s2)) {
      phys.h = s1->radius;
    } else {
      phys.h = s2->radius;
    }
  }
  
  Real Mass = bodies[id2]->state->mass;
  if (Mass==0.0 and bodies[id1]->state->mass!= 0.0) {
    Mass = bodies[id1]->state->mass;
  }
  
  Real Rho = bodies[id2]->rho;
  if (Rho==0.0 and bodies[id1]->rho!=0.0) {
    Rho = bodies[id1]->rho;
  }
  
  const Vector3r xixj = (-geom.penetrationDepth + phys.h)*geom.normal;
  
  if (xixj.norm() < phys.h) {
    Real fpressure = 0.0;
    if (Rho!=0.0) {
      // [Mueller2003], (10)
      fpressure = -Mass * 
                  (bodies[id1]->press + bodies[id2]->press)/(2.0*Rho) *
                  phys.kernelFunctionCurrentPressure(xixj.norm(), phys.h);
    }
    
    Vector3r fvisc = Vector3r::Zero();
    if (Rho!=0.0) {
      fvisc = phys.mu * Mass * 
                  normalVelocity*geom.normal/Rho *
                  phys.kernelFunctionCurrentVisco(xixj.norm(), phys.h);
    }
    force = fpressure*geom.normal + fvisc;
    return true;
  } else {
    return false;
  }
}
YADE_PLUGIN((SPHEngine));
#endif


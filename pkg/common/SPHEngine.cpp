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
      //std::cerr<<"(b->rho - b->rho0) = " << (b->rho - b->rho0)<<"; b->press = " << b->press<<std::endl;
    } YADE_PARALLEL_FOREACH_BODY_END();
  }
}

void SPHEngine::calculateSPHRho(const shared_ptr<Body>& b) {
  if (b->rho0<0) {
    b->rho0 = rho0;
  }
  Real rho = 0;
  
  // Pointer to kernel function
  KernelFunction kernelFunctionCurDensity = returnKernelFunction (KernFunctionDensity, Norm);
  
  // Calculate rho for every particle
  for(Body::MapId2IntrT::iterator it=b->intrs.begin(),end=b->intrs.end(); it!=end; ++it) {
    const shared_ptr<Body> b2 = Body::byId((*it).first,scene);
    Sphere* s=dynamic_cast<Sphere*>(b->shape.get());
    if(!s) continue;
    
    // std::cerr<<"Density, b->id"<<b->id<<"; b2->id"<<b2->id<<"   ";
    if (((*it).second)->geom and ((*it).second)->phys) {
      const ScGeom geom = *(YADE_PTR_CAST<ScGeom>(((*it).second)->geom));
      const ViscElPhys phys=*(YADE_PTR_CAST<ViscElPhys>(((*it).second)->phys));
      
      if((b2->groupMask & mask)==0)  continue;
      
      Real Mass = b2->state->mass;
      if (Mass == 0) Mass = b->state->mass;
      
      const Real SmoothDist = (b2->state->pos - b->state->pos).norm();
     
      // [Mueller2003], (3)
      rho += b2->state->mass*kernelFunctionCurDensity(SmoothDist, h);
      // std::cerr<<"SmoothDist = "<<SmoothDist<<"; h = "<<h<<"; kern(Smooth, h)="<<kernelFunctionCurDensity(SmoothDist, h)<<"; b2->state->mass*kernelFunctionCurDensity(SmoothDist, h)="<<b2->state->mass*kernelFunctionCurDensity(SmoothDist, h);
    }
    // std::cerr<<std::endl;
  }
  // Self mass contribution
  rho += b->state->mass*kernelFunctionCurDensity(0.0, h);
  b->rho = rho;
  //std::cerr<<"New rho = "<<b->rho<<std::endl;
}

Real smoothkernelPoly6(const double & rr, const double & hh) {
  if (rr<=hh) {
    return (2.1541870227086614782e0 * pow((1.0 - rr*rr / (hh*hh)), 4) /(hh*hh*hh));
  } else {
    return 0;
  }
}

Real smoothkernelPoly6Grad(const double & rr, const double & hh) {
  if (rr<=hh and hh != 0 ) {
    
    return -15/(M_PI)*(-3*pow((1-rr/hh),2));
    //return -25.066903536973515383e0 * rr * pow((hh - rr), 2) * pow(hh, -7);
  } else {
    return 0;
  }
}

Real smoothkernelPoly6Lapl(const double & rr, const double & hh) {
  if (rr<=hh) {
    
    return 45/(M_PI)*(1 - rr/hh);
    //return -25.066903536973515383e0 * rr * pow((hh - rr), 2) * pow(hh, -7);
  } else {
    return 0;
  }
}

Real smoothkernelSpiky(const double & rr, const double & hh) {
  if (rr<=hh) {
    const Real h = 1; const Real r = rr/hh;
    //return 15/(M_PI*pow(h,6))*(pow((h-r), 3));             // [Mueller2003], (21)
    return 15/(M_PI)*(pow((h-r), 3));             // [Mueller2003], (21)
  } else {
    return 0;
  }
}

Real smoothkernelSpikyGrad(const double & rr, const double & hh) {
  if (rr<=hh) {
    const Real h = 1; const Real r = rr/hh;
    //return -15/(M_PI*pow(h,6))*(-3*pow((h-r),2));
    return -15/(M_PI)*(-3*pow((h-r),2));
  } else {
    return 0;
  }
}

Real smoothkernelSpikyLapl(const double & rr, const double & hh) {
  if (rr<=hh) {
    const Real h = 1; const Real r = rr/hh;
    //return 15/(M_PI*pow(h,6))*(6*(h-r));
    return 15/(M_PI)*(6*(h-r));
  } else {
    return 0;
  }
}

Real smoothkernelVisco(const double & rr, const double & hh) {
  if (rr<=hh and rr!=0 and hh!=0) {
    const Real h = 1; const Real r = rr/hh;
    //return 15/(2*M_PI*pow(h,3))*(-r*r*r/(2*h*h*h) + r*r/(h*h) + h/(2*r) -1);   // [Mueller2003], (21)
    return 15/(2*M_PI)*(-r*r*r/(2) + r*r + h/(2*r) -1);   // [Mueller2003], (21)
  } else {
    return 0;
  }
}

Real smoothkernelViscoGrad(const double & rr, const double & hh) {
  if (rr<=hh and rr!=0 and hh!=0) {
    const Real h = 1; const Real r = rr/hh;
    //return -15/(2*M_PI*pow(h,3))*(-3*r*r/(2*h*h*h) + 2*r/(h*h) - h/(2*r*r));
    return -15/(2*M_PI)*(-3*r*r/(2) + 2*r - h/(2*r*r));
  } else {
    return 0;
  }
}

Real smoothkernelViscoLapl(const double & rr, const double & hh) {
  if (rr<=hh and rr!=0 and hh!=0) {
    const Real h = 1; const Real r = rr/hh;
    //return 45/(M_PI*pow(h,6))*(h - rrj);                     // [Mueller2003], (22+)
    //return 15/(2*M_PI*pow(h,3))*(-3*r*r/(2*h*h*h) + 2*r/(h*h) - h/(2*r*r));
    return 15/(2*M_PI)*(-3*r*r/(2) + 2*r - h/(2*r*r));
  } else {
    return 0;
  }
}

Real smoothkernelLucy(const double & rr, const double & hh) {
  if (rr<=hh and rr!=0 and hh!=0) {
    //const Real h = 1; 
    //return 5/(9*M_PI*pow(h,2))*(1+3*r/h)*pow((1-r/h),3);
    const Real r = rr/hh;
    return 5/(9*M_PI)*(1+3*r)*pow((1-r),3);
  } else {
    return 0;
  }
}

Real smoothkernelLucyGrad(const double & rr, const double & hh) {
  if (rr<=hh and rr!=0 and hh!=0) {
    //const Real h = 1; 
    //return -5/(9*M_PI*pow(h,2))*(-12*r/(h*h))*pow((1-r/h),2);
    const Real r = rr/hh;
    return -5/(9*M_PI)*(-12*r)*pow((1-r),2);
  } else {
    return 0;
  }
}

Real smoothkernelLucyLapl(const double & rr, const double & hh) {
  if (rr<=hh and rr!=0 and hh!=0) {
    //const Real h = 1; 
    //return  5/(9*M_PI*pow(h,2))*(-12/(h*h))*(1-r/h)*(1-3*r/h);
    const Real r = rr/hh;
    return  5/(9*M_PI)*(-12)*(1-r)*(1-3*r);
  } else {
    return 0;
  }
}

Real smoothkernelMonaghan(const double & rr, const double & hh) {
  Real ret = 0.0;
  if (hh!=0) {
    //const Real h = 1; 
    const Real r = rr/hh;
    if (rr/hh<0.5) {
      //ret = 40/(7*M_PI*h*h)*(1 - 6*pow((r/h),2) + 6*pow((r/h),3));
      ret = 40/(7*M_PI)*(1 - 6*pow((r),2) + 6*pow((r),3));
    } else {
      //ret = 80/(7*M_PI*h*h)*pow((1 - (r/h)), 3);
      ret = 80/(7*M_PI)*pow((1 - (r)), 3);
    }
  }
  return ret;
}

Real smoothkernelMonaghanGrad(const double & rr, const double & hh) {
  Real ret = 0.0;
  if (hh!=0) {
    const Real h = 1; const Real r = rr/hh;
    if (rr/hh<0.5) {
      //ret = -40/(7*M_PI*h*h)*( - 6*r/(h*h))*(2 - 3 * r/(h*h*h));
      ret = -40/(7*M_PI)*( - 6*r)*(2 - 3 * r);
    } else {
      //ret = -80/(7*M_PI*h*h)*( -3/h)*pow((1 -r/h), 2);
      ret = -80/(7*M_PI)*( -3/h)*pow((1 -r), 2);
    }
  }
  return ret;
}

Real smoothkernelMonaghanLapl(const double & rr, const double & hh) {
  Real ret = 0.0;
  if (hh!=0) {
    //const Real h = 1; 
    const Real r = rr/hh;
    if (rr/hh<0.5) {
      //ret = 40/(7*M_PI*h*h)*( - 12/(h*h))*(1 - 3 * r/(h*h*h*h*h));
      ret = 40/(7*M_PI)*( - 12)*(1 - 3 * r);
    } else {
      //ret = 80/(7*M_PI*h*h)*( 6/(h*h))*(1 -r/h);
      ret = 80/(7*M_PI)*( 6)*(1 -r);
    }
  }
  return ret;
}

KernelFunction returnKernelFunction(const int a, const typeKernFunctions typeF) {
  return returnKernelFunction(a, a, typeF);
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
  
  Real Mass = bodies[id2]->state->mass;
  if (Mass==0.0 and bodies[id1]->state->mass!= 0.0) {
    Mass = bodies[id1]->state->mass;
  }
  
  Real Rho = bodies[id2]->rho;
  if (Rho==0.0 and bodies[id1]->rho!=0.0) {
    Rho = bodies[id1]->rho;
  }
  
  const Vector3r xixj = de2.pos - de1.pos;
  //std::cerr<<"Rho="<<Rho<<std::endl;
  if (xixj.norm() < phys.h) {
    Real fpressure = 0.0;
    if (Rho!=0.0) {
      // [Mueller2003], (10)
      fpressure = Mass * 
                  (bodies[id1]->press + bodies[id2]->press)/(2.0*Rho) *
                  phys.kernelFunctionCurrentPressure(xixj.norm(), phys.h);
      /*            
      std::cerr<<"bodies[id1]->press = "<<bodies[id1]->press<<"; bodies[id2]->press = "<<bodies[id2]->press<<"; 2.0*Rho = "<<2.0*Rho
               <<"; phys.kern = "<<phys.kernelFunctionCurrentPressure(xixj.norm(), phys.h)
               <<"; Mass = "<<Mass
               <<"; fpressure = "<<fpressure<<std::endl;
      */
    }
    
    Vector3r fvisc = Vector3r::Zero();
    if (Rho!=0.0 and phys.mu !=0.0 ) {
      fvisc = -phys.mu * Mass * 
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


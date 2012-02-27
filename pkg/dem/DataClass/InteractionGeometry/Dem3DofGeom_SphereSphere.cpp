#include "Dem3DofGeom_SphereSphere.hpp"

#include<yade/pkg-common/Sphere.hpp>
#include<yade/core/Omega.hpp>
YADE_PLUGIN((Dem3DofGeom_SphereSphere)(Dem6DofGeom_SphereSphere)
	#ifdef YADE_OPENGL
		(Gl1_Dem3DofGeom_SphereSphere)
	#endif
	(Ig2_Sphere_Sphere_Dem3DofGeom)(Ig2_Sphere_Sphere_Dem6DofGeom));


Dem3DofGeom_SphereSphere::~Dem3DofGeom_SphereSphere(){}

/*! Project point from sphere surface to tangent plane,
 * such that the angle of shortest arc from (1,0,0) pt on the sphere to the point itself is the same
 * as the angle of segment of the same length on the tangent plane.
 *
 * This function is (or should be) inverse of ScGeom::rollPlanePtToSphere.
 * 
 * @param fromXtoPtOri gives orientation of the vector from sphere center to the sphere point from the global +x axis.
 * @param radius the distance from sphere center to the contact plane
 * @param planeNormal unit vector pointing away from the sphere center, determining plane orientation on which the projected point lies.
 * @returns The projected point coordinates (with origin at the contact point).
 */
Vector3r Dem3DofGeom_SphereSphere::unrollSpherePtToPlane(const Quaternionr& fromXtoPtOri, const Real& radius, const Vector3r& planeNormal){
	Quaternionr normal2pt; normal2pt.setFromTwoVectors(planeNormal,fromXtoPtOri*Vector3r::UnitX());
	AngleAxisr aa(angleAxisFromQuat(normal2pt));
	return (aa.angle()*radius) /* length */ *(aa.axis().cross(planeNormal)) /* direction: both are unit vectors */;
}

/*! Project point from tangent plane to the sphere.
 *
 * This function is (or should be) inverse of ScGeom::unrollSpherePtToPlane.
 *
 * @param planePt point on the tangent plane, with origin at the contact point (i.e. at sphere center + normal*radius)
 * @param radius sphere radius
 * @param planeNormal _unit_ vector pointing away from sphere center
 * @returns orientation that transforms +x axis to the vector between sphere center and point on the sphere that corresponds to planePt.
 *
 * @note It is not checked whether planePt relly lies on the tangent plane. If not, result will be incorrect.
 */
Quaternionr Dem3DofGeom_SphereSphere::rollPlanePtToSphere(const Vector3r& planePt, const Real& radius, const Vector3r& planeNormal){
	if (planePt!=Vector3r::Zero()) {
		Quaternionr normal2pt;
		Vector3r axis=planeNormal.cross(planePt); axis.normalize();
		Real angle=planePt.norm()/radius;
		normal2pt=Quaternionr(AngleAxisr(angle,axis));
		Quaternionr ret; return ret.setFromTwoVectors(Vector3r::UnitX(),normal2pt*planeNormal);
	} else {
		Quaternionr ret; return ret.setFromTwoVectors(Vector3r::UnitX(),planeNormal);
	}
}



/* Set contact points on both spheres such that their projection is the one given
 * (should be on the plane passing through origin and oriented with normal; not checked!)
 */
void Dem3DofGeom_SphereSphere::setTgPlanePts(Vector3r p1new, Vector3r p2new){
	cp1rel=ori1.conjugate()*rollPlanePtToSphere(p1new,effR1,normal);
	cp2rel=ori2.conjugate()*rollPlanePtToSphere(p2new,effR2,-normal);
}



/*! Perform slip of the projected contact points so that their distance becomes equal (or remains smaller) than the given one.
 * The slipped distance is returned.
 */
Real Dem3DofGeom_SphereSphere::slipToDisplacementTMax(Real displacementTMax){
	// very close, reset shear
	if(displacementTMax<=0.){ setTgPlanePts(Vector3r(0,0,0),Vector3r(0,0,0)); return displacementTMax;}
	// otherwise
	Vector3r p1=contPtInTgPlane1(), p2=contPtInTgPlane2();
	Real currDistSq=(p2-p1).squaredNorm();
	if(currDistSq<pow(displacementTMax,2)) return 0; // close enough, no slip needed
	Vector3r diff=.5*(displacementTMax/sqrt(currDistSq)-1)*(p2-p1);
	setTgPlanePts(p1-diff,p2+diff);
	return 2*diff.norm();
}

/*! As above : perform slip of the projected contact points. Here, we directly give the multiplier applied on the distance for faster results.
 * The plastic displacement (vector) is returned.
 */
Vector3r Dem3DofGeom_SphereSphere::scaleDisplacementT(Real multiplier){
	assert(multiplier>=0 && multiplier<=1);
	Vector3r p1=contPtInTgPlane1(), p2=contPtInTgPlane2();
	Vector3r diff=.5*(multiplier-1)*(p2-p1);
	setTgPlanePts(p1-diff,p2+diff);
	return diff*2.0;
}


/* Move contact point on both spheres in such way that their relative position (displacementT) is the same;
 * this should be done regularly to ensure that the angle doesn't go over π, since then quaternion would
 * flip axis and the point would project on other side of the tangent plane piece. */
void Dem3DofGeom_SphereSphere::relocateContactPoints(){
	relocateContactPoints(contPtInTgPlane1(),contPtInTgPlane2());
}

/*! Like Dem3DofGeom_SphereSphere::relocateContactPoints(), but use already computed tangent plane points. */
void Dem3DofGeom_SphereSphere::relocateContactPoints(const Vector3r& p1, const Vector3r& p2){
	Vector3r midPt=(effR1/(effR1+effR2))*(p1+p2); // proportionally to radii, so that angle would be the same
	if((p1.squaredNorm()>pow(effR1,2) || p2.squaredNorm()>pow(effR2,2)) && midPt.squaredNorm()>pow(min(effR1,effR2),2)){
		//cerr<<"RELOCATION with displacementT="<<displacementT(); // should be the same before and after relocation
		setTgPlanePts(p1-midPt,p2-midPt);
		//cerr<<" → "<<displacementT()<<endl;
	}
}


Dem6DofGeom_SphereSphere::~Dem6DofGeom_SphereSphere(){}

Vector3r Dem6DofGeom_SphereSphere::relRotVector() const{
	// FIXME: this is not correct, as it assumes normal will not change (?)
	Quaternionr relOri12=ori1.conjugate()*ori2;
	Quaternionr oriDiff=initRelOri12.conjugate()*relOri12;
	AngleAxisr aa(angleAxisFromQuat(oriDiff));
	if(aa.angle()>Mathr::PI)aa.angle()-=Mathr::TWO_PI;
	// cerr<<axis<<";"<<angle<<";"<<ori1<<";"<<ori2<<";"<<oriDiff<<endl;
	return aa.angle()*aa.axis();
}

void Dem6DofGeom_SphereSphere::bendTwistAbs(Vector3r& bend, Real& twist){
	const Vector3r& relRot=relRotVector();
	twist=relRot.dot(normal);
	bend=relRot-twist*normal;
}


#ifdef YADE_OPENGL
	#include<yade/lib-opengl/OpenGLWrapper.hpp>
	#include<yade/lib-opengl/GLUtils.hpp>
	bool Gl1_Dem3DofGeom_SphereSphere::normal=false;
	bool Gl1_Dem3DofGeom_SphereSphere::rolledPoints=false;
	bool Gl1_Dem3DofGeom_SphereSphere::unrolledPoints=false;
	bool Gl1_Dem3DofGeom_SphereSphere::shear=false;
	bool Gl1_Dem3DofGeom_SphereSphere::shearLabel=false;

	void Gl1_Dem3DofGeom_SphereSphere::go(const shared_ptr<InteractionGeometry>& ig, const shared_ptr<Interaction>& ip, const shared_ptr<Body>& b1, const shared_ptr<Body>& b2, bool wireFrame){
		Dem3DofGeom_SphereSphere* ss = static_cast<Dem3DofGeom_SphereSphere*>(ig.get());
		//const Se3r& se31=b1->physicalParameters->dispSe3,se32=b2->physicalParameters->dispSe3;
		const Se3r& se31=b1->state->se3,se32=b2->state->se3;
		const Vector3r& pos1=se31.position,pos2=se32.position;
		Vector3r& contPt=ss->contactPoint;
		
		if(normal){
			GLUtils::GLDrawArrow(contPt,contPt+ss->normal*.5*ss->refLength); // normal of the contact
		}
		#if 0
			// never used, since bending/torsion not used
			//Vector3r contPt=se31.position+(ss->effR1/ss->refLength)*(se32.position-se31.position); // must be recalculated to not be unscaled if scaling displacements ...
			GLUtils::GLDrawLine(pos1,pos2,Vector3r(.5,.5,.5));
			Vector3r bend; Real tors;
			ss->bendingTorsionRel(bend,tors);
			GLUtils::GLDrawLine(contPt,contPt+10*ss->radius1*(bend+ss->normal*tors),Vector3r(1,0,0));
			#if 0
				GLUtils::GLDrawNum(bend[0],contPt-.2*ss->normal*ss->radius1,Vector3r(1,0,0));
				GLUtils::GLDrawNum(bend[1],contPt,Vector3r(0,1,0));
				GLUtils::GLDrawNum(bend[2],contPt+.2*ss->normal*ss->radius1,Vector3r(0,0,1));
				GLUtils::GLDrawNum(tors,contPt+.5*ss->normal*ss->radius2,Vector3r(1,1,0));
			#endif
		#endif
		// sphere center to point on the sphere
		if(rolledPoints){
			GLUtils::GLDrawLine(pos1,pos1+(ss->ori1*ss->cp1rel*Vector3r::UnitX()*ss->effR1),Vector3r(0,.5,1));
			GLUtils::GLDrawLine(pos2,pos2+(ss->ori2*ss->cp2rel*Vector3r::UnitX()*ss->effR2),Vector3r(0,1,.5));
		}
		//TRVAR4(pos1,ss->ori1,pos2,ss->ori2);
		//TRVAR2(ss->cp2rel,pos2+(ss->ori2*ss->cp2rel*Vector3r::UnitX()*ss->effR2));
		// contact point to projected points
		if(unrolledPoints||shear){
			Vector3r ptTg1=ss->contPtInTgPlane1(), ptTg2=ss->contPtInTgPlane2();
			if(unrolledPoints){
				//TRVAR3(ptTg1,ptTg2,ss->normal)
				GLUtils::GLDrawLine(contPt,contPt+ptTg1,Vector3r(0,.5,1)); GLUtils::GLDrawLine(pos1,contPt+ptTg1,Vector3r(0,.5,1));
				GLUtils::GLDrawLine(contPt,contPt+ptTg2,Vector3r(0,1,.5)); GLUtils::GLDrawLine(pos2,contPt+ptTg2,Vector3r(0,1,.5));
			}
			if(shear){
				GLUtils::GLDrawLine(contPt+ptTg1,contPt+ptTg2,Vector3r(1,1,1));
				if(shearLabel) GLUtils::GLDrawNum(ss->displacementT().norm(),contPt,Vector3r(1,1,1));
			}
		}
	}
#endif

CREATE_LOGGER(Ig2_Sphere_Sphere_Dem3DofGeom);

bool Ig2_Sphere_Sphere_Dem3DofGeom::go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c){
	Sphere *s1=static_cast<Sphere*>(cm1.get()), *s2=static_cast<Sphere*>(cm2.get());
	Vector3r normal=(state2.pos+shift2)-state1.pos;
	Real penetrationDepthSq=pow((distFactor>0?distFactor:1.)*(s1->radius+s2->radius),2)-normal.squaredNorm();
	if (penetrationDepthSq<0 && !c->isReal() && !force){
		return false;
	}

	Real dist=normal.norm(); normal/=dist; /* normal is unit vector now */
	shared_ptr<Dem3DofGeom_SphereSphere> ss;
	if(c->interactionGeometry) ss=YADE_PTR_CAST<Dem3DofGeom_SphereSphere>(c->interactionGeometry);
	else {
		ss=shared_ptr<Dem3DofGeom_SphereSphere>(new Dem3DofGeom_SphereSphere());
		c->interactionGeometry=ss;
		
		// constants
		if(distFactor>0) ss->refLength=dist;
		else ss->refLength=s1->radius+s2->radius;
		ss->refR1=s1->radius; ss->refR2=s2->radius;
		
		Real penetrationDepth=s1->radius+s2->radius-ss->refLength;
		
		if(Omega::instance().getCurrentIteration()<=10){
			ss->effR1=s1->radius-.5*penetrationDepth; ss->effR2=s2->radius-.5*penetrationDepth;
		} else {ss->effR1=s1->radius; ss->effR2=s2->radius;}
		
		// for bending only: ss->initRelOri12=state1.ori.Conjugate()*state2.ori;
		// quasi-constants
		ss->cp1rel.setFromTwoVectors(Vector3r::UnitX(),state1.ori.conjugate()*normal);
		ss->cp2rel.setFromTwoVectors(Vector3r::UnitX(),state2.ori.conjugate()*(-normal));
		ss->cp1rel.normalize(); ss->cp2rel.normalize();
	}
	ss->normal=normal;
	ss->contactPoint=state1.pos+(ss->effR1-.5*(ss->refLength-dist))*ss->normal;
	ss->se31=state1.se3; ss->se32=state2.se3; ss->se32.position+=shift2;
	return true;
}

CREATE_LOGGER(Ig2_Sphere_Sphere_Dem6DofGeom);
bool Ig2_Sphere_Sphere_Dem6DofGeom::go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c){
	bool hadIntrGeom=c->interactionGeometry;
	if(!Ig2_Sphere_Sphere_Dem3DofGeom::go(cm1,cm2,state1,state2,shift2,force,c)) return false;
	// HACK: dem3dof functor creates a dem3dof instance; we need to copy-construct dem6dof from it instead
	// proper solution would be to factor out the computation part from the dem3dof to separate functions and call those from here, or make the dem3dof functor templated on the dem3dof/dem6dof
	if(!hadIntrGeom){
		assert(c->interactionGeometry);
		assert(c->interactionGeometry->getClassName()=="Dem3DofGeom_SphereSphere");
		const shared_ptr<Dem6DofGeom_SphereSphere> geom(new Dem6DofGeom_SphereSphere(*YADE_CAST<Dem3DofGeom_SphereSphere*>(c->interactionGeometry.get())));
		geom->initRelOri12=state1.ori.conjugate()*state2.ori;
		c->interactionGeometry=geom;
		//TRVAR3(geom->refLength,geom->contactPoint,geom->initRelOri12)
	}
	return true;
}

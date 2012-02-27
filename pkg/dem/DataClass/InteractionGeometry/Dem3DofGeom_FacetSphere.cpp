// © Václav Šmilauer <eudoxos@arcig.cz>
#include "Dem3DofGeom_FacetSphere.hpp"
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Facet.hpp>
YADE_PLUGIN((Dem3DofGeom_FacetSphere)
	#ifdef YADE_OPENGL
		(Gl1_Dem3DofGeom_FacetSphere)
	#endif	
		(Ig2_Facet_Sphere_Dem3DofGeom));

CREATE_LOGGER(Dem3DofGeom_FacetSphere);
Dem3DofGeom_FacetSphere::~Dem3DofGeom_FacetSphere(){}

void Dem3DofGeom_FacetSphere::setTgPlanePts(const Vector3r& p1new, const Vector3r& p2new){
	TRVAR3(cp1pt,cp2rel,contPtInTgPlane2()-contPtInTgPlane1());	
	cp1pt=se31.orientation.conjugate()*(turnPlanePt(p1new,normal,se31.orientation*localFacetNormal)+contactPoint-se31.position);
	cp2rel=se32.orientation.conjugate()*Dem3DofGeom_SphereSphere::rollPlanePtToSphere(p2new,effR2,-normal);
	TRVAR3(cp1pt,cp2rel,contPtInTgPlane2()-contPtInTgPlane1());	
}

void Dem3DofGeom_FacetSphere::relocateContactPoints(const Vector3r& p1, const Vector3r& p2){
	//TRVAR2(p2.norm(),effR2);
	if(p2.squaredNorm()>pow(effR2,2)){
		setTgPlanePts(Vector3r::Zero(),p2-p1);
	}
}

Real Dem3DofGeom_FacetSphere::slipToDisplacementTMax(Real displacementTMax){
	//FIXME: not yet tested
	// negative or zero: reset shear
	if(displacementTMax<=0.){ setTgPlanePts(Vector3r(0,0,0),Vector3r(0,0,0)); return displacementTMax;}
	// otherwise
	Vector3r p1=contPtInTgPlane1(), p2=contPtInTgPlane2();
	Real currDistSq=(p2-p1).squaredNorm();
	if(currDistSq<pow(displacementTMax,2)) return 0; // close enough, no slip needed
	//Vector3r diff=.5*(sqrt(currDistSq)/displacementTMax-1)*(p2-p1); setTgPlanePts(p1+diff,p2-diff);
	Real scale=displacementTMax/sqrt(currDistSq); setTgPlanePts(scale*p1,scale*p2);
	return (displacementTMax/scale)*(1-scale);
}

Vector3r Dem3DofGeom_FacetSphere::scaleDisplacementT(Real multiplier){
	assert(multiplier>=0 && multiplier<=1);
	Vector3r p1=contPtInTgPlane1(), p2=contPtInTgPlane2();
	Vector3r diff=.5*(multiplier-1)*(p2-p1);
	setTgPlanePts(p1-diff,p2+diff);
	return diff*2.0;
}


CREATE_LOGGER(Ig2_Facet_Sphere_Dem3DofGeom);
bool Ig2_Facet_Sphere_Dem3DofGeom::go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c){
	Facet* facet=static_cast<Facet*>(cm1.get());
	Real sphereRadius=static_cast<Sphere*>(cm2.get())->radius;

	// InteractionGeometryFunctor::go(cm1,cm2,state1,state2,shift2,force,c);

	#if 1
		/* new code written from scratch, to make sure the algorithm is correct; it is about the same speed 
			as sega's algo below, but seems more readable to me.
			The FACET_TOPO thing is still missing here but can be copied literally once it is tested */
		// begin facet-local coordinates
			Vector3r cogLine=state1.ori.conjugate()*(state2.pos+shift2-state1.pos); // connect centers of gravity
			//TRVAR4(state1.pos,state1.ori,state2.pos,cogLine);
			Vector3r normal=facet->nf;
			Real planeDist=normal.dot(cogLine);
			if(planeDist<0){normal*=-1; planeDist*=-1; }
			if(planeDist>sphereRadius && !c->isReal()) { /* LOG_TRACE("Sphere too far ("<<planeDist<<") from plane"); */ return false;  }
			Vector3r planarPt=cogLine-planeDist*normal; // project sphere center to the facet plane
			Real normDotPt[3];
			Vector3r contactPt(Vector3r::Zero());
			for(int i=0; i<3; i++) normDotPt[i]=facet->ne[i].dot(planarPt-facet->vertices[i]);
			short w=(normDotPt[0]>0?1:0)+(normDotPt[1]>0?2:0)+(normDotPt[2]>0?4:0);
			//TRVAR4(planarPt,normDotPt[0],normDotPt[1],normDotPt[2]);
			//TRVAR2(normal,cogLine);
			//TRVAR3(facet->vertices[0],facet->vertices[1],facet->vertices[2]);
			switch(w){
				case 0: contactPt=planarPt; break; // inside triangle
				case 1: contactPt=getClosestSegmentPt(planarPt,facet->vertices[0],facet->vertices[1]); break; // +-- (n1)
				case 2: contactPt=getClosestSegmentPt(planarPt,facet->vertices[1],facet->vertices[2]); break; // -+- (n2)
				case 4: contactPt=getClosestSegmentPt(planarPt,facet->vertices[2],facet->vertices[0]); break; // --+ (n3)
				case 3: contactPt=facet->vertices[1]; break; // ++- (v1)
				case 5: contactPt=facet->vertices[0]; break; // +-+ (v0)
				case 6: contactPt=facet->vertices[2]; break; // -++ (v2)
				case 7: throw logic_error("Impossible triangle intersection?"); // +++ (impossible)
				default: throw logic_error("Nonsense intersection value!");
			}
			normal=cogLine-contactPt; // called normal, but it is no longer the facet's normal (for compat)
			//TRVAR3(normal,contactPt,sphereRadius);
			if(!c->isReal() && normal.squaredNorm()>sphereRadius*sphereRadius && !force) { /* LOG_TRACE("Sphere too far from closest point"); */ return false; } // fast test before sqrt
			Real norm=normal.norm(); normal/=norm; // normal is unit vector now
			Real penetrationDepth=sphereRadius-norm;
	#else
		/* This code was mostly copied from InteractingFacet2InteractinSphere4SpheresContactGeometry */
		// begin facet-local coordinates 
			Vector3r contactLine=state1.ori.Conjugate()*(state2.pos+shift2-state1.pos);
			Vector3r normal=facet->nf;
			Real L=normal.Dot(contactLine); // height/depth of sphere's center from facet's plane
			if(L<0){normal*=-1; L*=-1;}
			if(L>sphereRadius && !c->isReal()) return false; // sphere too far away from the plane

			Vector3r contactPt=contactLine-L*normal; // projection of sphere's center to facet's plane (preliminary contact point)
			const Vector3r* edgeNormals=facet->ne; // array[3] of edge normals (in facet plane)
			int edgeMax=0; Real distMax=edgeNormals[0].Dot(contactPt);
			for(int i=1; i<3; i++){
				Real dist=edgeNormals[i].Dot(contactPt);
				if(distMax<dist){edgeMax=i; distMax=dist;}
			}
			//TRVAR2(distMax,edgeMax);
			// OK, what's the logic here? Copying from IF2IS4SCG…
			Real sphereRReduced=shrinkFactor*sphereRadius;
			Real inCircleR=facet->icr-sphereRReduced;
			Real penetrationDepth;
			if(inCircleR<0){inCircleR=facet->icr; sphereRReduced=0;}
			if(distMax<inCircleR){// contact with facet's surface
				penetrationDepth=sphereRadius-L;	
				normal.normalize();
			} else { // contact with the edge
				contactPt+=edgeNormals[edgeMax]*(inCircleR-distMax);
				bool noVertexContact=false;
				//TRVAR3(edgeNormals[edgeMax],inCircleR,distMax);
				// contact with vertex no. edgeMax
				// FIXME: this is the original version, but why (edgeMax-1)%3? IN that case, edgeNormal to edgeMax would never be tried
				//    if     (contactPt.Dot(edgeNormals[        (edgeMax-1)%3])>inCircleR) contactPt=facet->vu[edgeMax]*(facet->vl[edgeMax]-sphereRReduced);
				if     (contactPt.Dot(edgeNormals[        edgeMax      ])>inCircleR) contactPt=facet->vu[edgeMax]*(facet->vl[edgeMax]-sphereRReduced);
				// contact with vertex no. edgeMax+1
				else if(contactPt.Dot(edgeNormals[edgeMax=(edgeMax+1)%3])>inCircleR) contactPt=facet->vu[edgeMax]*(facet->vl[edgeMax]-sphereRReduced);
				// contact with edge no. edgeMax
				else noVertexContact=true;
				normal=contactLine-contactPt;
				#ifdef FACET_TOPO
					if(noVertexContact && facet->edgeAdjIds[edgeMax]!=Body::ID_NONE){
						// find angle between our normal and the facet's normal (still local coords)
						Quaternionr q; q.Align(facet->nf,normal); AngleAxisr aa(angleAxisFromQuat(q));
						assert(aa.angle()>=0 && aa.angle()<=Mathr::PI);
						if(edgeNormals[edgeMax].Dot(aa.axis())<0) aa.angle()*=-1.;
						bool negFace=normal.Dot(facet->nf)<0; // contact in on the negative facet's face
						Real halfAngle=(negFace?-1.:1.)*facet->edgeAdjHalfAngle[edgeMax]; 
						if(halfAngle<0 && aa.angle()>halfAngle) return false; // on concave boundary, and if in the other facet's sector, no contact
						// otherwise the contact will be created
					}
				#endif
				//TRVAR4(contactLine,contactPt,normal,normal.norm());
				//TRVAR3(se31.orientation*contactLine,se31.position+se31.orientation*contactPt,se31.orientation*normal);
				Real norm=normal.norm(); normal/=norm;
				penetrationDepth=sphereRadius-norm;
				//TRVAR1(penetrationDepth);
			}
		// end facet-local coordinates
	#endif

	if(penetrationDepth<0 && !c->isReal()) return false;


	shared_ptr<Dem3DofGeom_FacetSphere> fs;
	Vector3r normalGlob=state1.ori*normal;
	bool isNew=false;
	if(c->interactionGeometry) fs=YADE_PTR_CAST<Dem3DofGeom_FacetSphere>(c->interactionGeometry);
	else {
		// LOG_TRACE("Creating new Dem3DofGeom_FacetSphere");
		fs=shared_ptr<Dem3DofGeom_FacetSphere>(new Dem3DofGeom_FacetSphere());
		c->interactionGeometry=fs;
		isNew=true;
		fs->effR2=sphereRadius-penetrationDepth;
		fs->refR1=-1; fs->refR2=sphereRadius;
		// postponed till below, to avoid numeric issues
		// see https://lists.launchpad.net/yade-dev/msg02794.html
		// since displacementN() is computed from fs->contactPoint,
		// it was returning something like +1e-16 at the very first step
		// when it was created ⇒ the constitutive law was erasing the
		// contact as soon as it was created.
		// fs->refLength=…
		fs->cp1pt=contactPt; // facet-local intial contact point
		fs->localFacetNormal=facet->nf;
		fs->cp2rel.setFromTwoVectors(Vector3r::UnitX(),state2.ori.conjugate()*(-normalGlob)); // initial sphere-local center-contactPt orientation WRT +x
		fs->cp2rel.normalize();
	}
	fs->se31=state1.se3; fs->se32=state2.se3; fs->se32.position+=shift2;
	fs->normal=normalGlob;
	fs->contactPoint=state2.pos+shift2+(-normalGlob)*(sphereRadius-penetrationDepth);
	// this refLength computation mimics what displacementN() does inside
	// displcementN will therefore return exactly zero at the step the contact
	// was created, which is what we want
	if(isNew) fs->refLength=(state2.pos+shift2-fs->contactPoint).norm();
	return true;
}

#ifdef YADE_OPENGL

	#include<yade/lib-opengl/OpenGLWrapper.hpp>
	#include<yade/lib-opengl/GLUtils.hpp>

	bool Gl1_Dem3DofGeom_FacetSphere::normal=false;
	bool Gl1_Dem3DofGeom_FacetSphere::rolledPoints=false;
	bool Gl1_Dem3DofGeom_FacetSphere::unrolledPoints=false;
	bool Gl1_Dem3DofGeom_FacetSphere::shear=false;
	bool Gl1_Dem3DofGeom_FacetSphere::shearLabel=false;

	void Gl1_Dem3DofGeom_FacetSphere::go(const shared_ptr<InteractionGeometry>& ig, const shared_ptr<Interaction>& ip, const shared_ptr<Body>& b1, const shared_ptr<Body>& b2, bool wireFrame){
		Dem3DofGeom_FacetSphere* fs = static_cast<Dem3DofGeom_FacetSphere*>(ig.get());
		const Se3r& se31=b1->state->se3,se32=b2->state->se3;
		const Vector3r& pos1=se31.position; const Vector3r& pos2=se32.position;
		const Quaternionr& ori1=se31.orientation; const Quaternionr& ori2=se32.orientation;
		const Vector3r& contPt=fs->contactPoint;
		
		if(normal){
			GLUtils::GLDrawArrow(contPt,contPt+fs->refLength*fs->normal); // normal of the contact
		}
		// sphere center to point on the sphere
		if(rolledPoints){
			//cerr<<pos1<<" "<<pos1+ori1*fs->cp1pt<<" "<<contPt<<endl;
			GLUtils::GLDrawLine(pos1+ori1*fs->cp1pt,contPt,Vector3r(0,.5,1));
			GLUtils::GLDrawLine(pos2,pos2+(ori2*fs->cp2rel*Vector3r::UnitX()*fs->effR2),Vector3r(0,1,.5));
		}
		// contact point to projected points
		if(unrolledPoints||shear){
			Vector3r ptTg1=fs->contPtInTgPlane1(), ptTg2=fs->contPtInTgPlane2();
			if(unrolledPoints){
				//TRVAR3(ptTg1,ptTg2,ss->normal)
				GLUtils::GLDrawLine(contPt,contPt+ptTg1,Vector3r(0,.5,1));
				GLUtils::GLDrawLine(contPt,contPt+ptTg2,Vector3r(0,1,.5)); GLUtils::GLDrawLine(pos2,contPt+ptTg2,Vector3r(0,1,.5));
			}
			if(shear){
				GLUtils::GLDrawLine(contPt+ptTg1,contPt+ptTg2,Vector3r(1,1,1));
				if(shearLabel) GLUtils::GLDrawNum(fs->displacementT().norm(),contPt,Vector3r(1,1,1));
			}
		}
	}

#endif


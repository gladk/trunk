# -*- coding: utf-8 -*-

# ---- Script to detect spheres which are "onJoint", according to JCFpm. -----
# To be called directly within an other script, for example, with execfile('identifBis.py')
# The sample (spheres + facets) has to exist already, with their JCFpmMat


### engines definition, according to our only goal that is to detect spheres concerned by joint surfaces

O.engines=[
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()],verletDist=0), #verletDist=0 to avoid introducing NewtonIntegrator in engines list
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True,label='interactionLaw')]
	)
]

O.dt=1 # whatever value is ok


############################ Identification spheres on joint
#### color set for particles on joint
jointcolor1=(0,1,0)
jointcolor2=(1,0,0)
jointcolor3=(0,0,1)
jointcolor4=(1,1,1)
jointcolor5=(0,0,0)


O.step(); # one step to know interactions

#### first step-> find spheres on facet
for i in O.interactions:
    ##if not i.isReal : continue
    ### Rk: facet are only stored in id1 
    if isinstance(O.bodies[i.id1].shape,Facet) and isinstance(O.bodies[i.id2].shape,Sphere): 
	vertices=O.bodies[i.id1].shape.vertices
	normalRef=vertices[0].cross(vertices[1]) # defines the normal to the facet normalRef
	nRef=normalRef/(normalRef.norm()) ## normalizes normalRef
	normalFacetSphere=i.geom.normal # geom.normal is oriented from id1 to id2 -> normalFacetSphere from facet (i.id1) to sphere (i.id2)

	if O.bodies[i.id2].state.onJoint==False : ## particles has not yet been identified as belonging to a joint plane
	    O.bodies[i.id2].state.onJoint=True
	    O.bodies[i.id2].state.joint=1
	    O.bodies[i.id2].shape.color=jointcolor1
	    if nRef.dot(normalFacetSphere)>=0 :
		O.bodies[i.id2].state.jointNormal1=nRef
	    elif nRef.dot(normalFacetSphere)<0 :
		O.bodies[i.id2].state.jointNormal1=-nRef
	elif O.bodies[i.id2].state.onJoint==True :  ## particles has already been identified as belonging to, at least, 1 facet
	    if O.bodies[i.id2].state.joint==1 and ((O.bodies[i.id2].state.jointNormal1.cross(nRef)).norm()>0.05) : ## particles has already been identified as belonging to only 1 facet
		O.bodies[i.id2].state.joint=2
		O.bodies[i.id2].shape.color=jointcolor2
		if nRef.dot(normalFacetSphere)>=0 :
		    O.bodies[i.id2].state.jointNormal2=nRef
		elif nRef.dot(normalFacetSphere)<0 :
		    O.bodies[i.id2].state.jointNormal2=-nRef
	    elif O.bodies[i.id2].state.joint==2 and ((O.bodies[i.id2].state.jointNormal1.cross(nRef)).norm()>0.05) and ((O.bodies[i.id2].state.jointNormal2.cross(nRef)).norm()>0.05): ## particles has already been identified as belonging to more than 1 facet
		O.bodies[i.id2].state.joint=3
		O.bodies[i.id2].shape.color=jointcolor3
		if nRef.dot(normalFacetSphere)>=0 :
		    O.bodies[i.id2].state.jointNormal3=nRef
		elif nRef.dot(normalFacetSphere)<0 :
		    O.bodies[i.id2].state.jointNormal3=-nRef
	    elif O.bodies[i.id2].state.joint==3 and ((O.bodies[i.id2].state.jointNormal1.cross(nRef)).norm()>0.05) and ((O.bodies[i.id2].state.jointNormal2.cross(nRef)).norm()>0.05) and ((O.bodies[i.id2].state.jointNormal3.cross(nRef)).norm()>0.05):
		O.bodies[i.id2].state.joint=4
		O.bodies[i.id2].shape.color=jointcolor5

#### second step -> find spheres interacting with spheres on facet (could be executed in the same timestep as step 1?)
for j in O.interactions:
    #if not i.isReal : continue
    ## Rk: facet are only stored in id1 
    if isinstance(O.bodies[j.id1].shape,Facet) and isinstance(O.bodies[j.id2].shape,Sphere): 
	vertices=O.bodies[j.id1].shape.vertices
	normalRef=vertices[0].cross(vertices[1]) # defines the normal to the facet normalRef
	nRef=normalRef/(normalRef.norm()) ## normalizes normalRef
	if ((O.bodies[j.id2].state.jointNormal1.cross(nRef)).norm()<0.05) :
	    jointNormalRef=O.bodies[j.id2].state.jointNormal1
	elif ((O.bodies[j.id2].state.jointNormal2.cross(nRef)).norm()<0.05) :
	    jointNormalRef=O.bodies[j.id2].state.jointNormal2
	elif ((O.bodies[j.id2].state.jointNormal3.cross(nRef)).norm()<0.05) :
	    jointNormalRef=O.bodies[j.id2].state.jointNormal3
	else : continue
	facetCenter=O.bodies[j.id1].state.pos
	#### seek for each sphere interacting with the identified sphere j.id2
	for n in O.interactions.withBody(j.id2) :
            if isinstance(O.bodies[n.id1].shape,Sphere) and isinstance(O.bodies[n.id2].shape,Sphere):
                if j.id2==n.id1: # the sphere that was detected on facet (that is, j.id2) is id1 of interaction n
                    sphOnF=n.id1
                    othSph=n.id2
                elif j.id2==n.id2: # here, this sphere that was detected on facet (that is, j.id2) is id2 of interaction n
                    sphOnF=n.id2
                    othSph=n.id1
		facetSphereDir=(O.bodies[othSph].state.pos-facetCenter)
		if O.bodies[othSph].state.onJoint==True :
		    if O.bodies[othSph].state.joint==3 and ((O.bodies[othSph].state.jointNormal1.cross(jointNormalRef)).norm()>0.05) and ((O.bodies[othSph].state.jointNormal2.cross(jointNormalRef)).norm()>0.05) and ((O.bodies[othSph].state.jointNormal3.cross(jointNormalRef)).norm()>0.05):
			O.bodies[othSph].state.joint=4
			O.bodies[othSph].shape.color=jointcolor5
		    elif O.bodies[othSph].state.joint==2 and ((O.bodies[othSph].state.jointNormal1.cross(jointNormalRef)).norm()>0.05) and ((O.bodies[othSph].state.jointNormal2.cross(jointNormalRef)).norm()>0.05):
			O.bodies[othSph].state.joint=3
			if facetSphereDir.dot(jointNormalRef)>=0:
			    O.bodies[othSph].state.jointNormal3=jointNormalRef
			elif facetSphereDir.dot(jointNormalRef)<0:
			    O.bodies[othSph].state.jointNormal3=-jointNormalRef 
		    elif O.bodies[othSph].state.joint==1 and ((O.bodies[othSph].state.jointNormal1.cross(jointNormalRef)).norm()>0.05) :
			O.bodies[othSph].state.joint=2
			if facetSphereDir.dot(jointNormalRef)>=0:
			    O.bodies[othSph].state.jointNormal2=jointNormalRef
			elif facetSphereDir.dot(jointNormalRef)<0:
			    O.bodies[othSph].state.jointNormal2=-jointNormalRef
		elif  O.bodies[othSph].state.onJoint==False :
		    O.bodies[othSph].state.onJoint=True
		    O.bodies[othSph].state.joint=1
		    O.bodies[othSph].shape.color=jointcolor4
		    if facetSphereDir.dot(jointNormalRef)>=0:
			O.bodies[othSph].state.jointNormal1=jointNormalRef
		    elif facetSphereDir.dot(jointNormalRef)<0:
			O.bodies[othSph].state.jointNormal1=-jointNormalRef



##### to delete facets
for b in O.bodies:
    if isinstance(b.shape,Facet):
	O.bodies.erase(b.id)

O.resetTime()
O.interactions.clear()

print '\nIdentificationSpheresOnJoint executed ! Spheres onJoint (and so on...) detected, facets deleted, simulation may go on.\n\n'




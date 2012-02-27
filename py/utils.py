# encoding: utf-8
#
# utility functions for yade
#
# 2008-2009 © Václav Šmilauer <eudoxos@arcig.cz>

"""Heap of functions that don't (yet) fit anywhere else.

Devs: please DO NOT ADD more functions here, it is getting too crowded!
"""

import math,random
from yade import *
from yade.wrapper import *
from miniEigen import *
try: # use psyco if available
	import psyco
	psyco.full()
except ImportError: pass
import doctest

# c++ implementations for performance reasons
from yade._utils import *

def saveVars(mark='',loadNow=True,**kw):
	"""Save passed variables into the simulation so that it can be recovered when the simulation is loaded again.

	For example, variables *a*, *b* and *c* are defined. To save them, use::

		>>> from yade import utils
		>>> utils.saveVars('mark',a=1,b=2,c=3)
		>>> from yade.params.mark import *
		>>> a,b,c
		(1, 2, 3)

	those variables will be save in the .xml file, when the simulation itself is saved. To recover those variables once the .xml is loaded again, use

		>>> utils.loadVars('mark')

	and they will be defined in the yade.params.\ *mark* module

	m*==True, variables will be loaded immediately after saving. That effectively makes *\*\*kw* available in builtin namespace.
	"""
	import cPickle
	Omega().tags['pickledPythonVariablesDictionary'+mark]=cPickle.dumps(kw)
	if loadNow: loadVars(mark)

def loadVars(mark=None):
	"""Load variables from saveVars, which are saved inside the simulation.
	If mark==None, all save variables are loaded. Otherwise only those with
	the mark passed."""
	import cPickle, types, sys, warnings
	def loadOne(d,mark=None):
		"""Load given dictionary into a synthesized module yade.params.name (or yade.params if *name* is not given). Update yade.params.__all__ as well."""
		import yade.params
		if mark:
			if mark in yade.params.__dict__: warnings.warn('Overwriting yade.params.%s which already exists.'%mark)
			modName='yade.params.'+mark
			mod=types.ModuleType(modName)
			mod.__dict__.update(d)
			mod.__all__=list(d.keys()) # otherwise params starting with underscore would not be imported
			sys.modules[modName]=mod
			yade.params.__all__.append(mark)
			yade.params.__dict__[mark]=mod
		else:
			yade.params.__all__+=list(d.keys())
			yade.params.__dict__.update(d)
	if mark!=None:
		d=cPickle.loads(Omega().tags['pickledPythonVariablesDictionary'+mark])
		loadOne(d,mark)
	else: # load everything one by one
		for m in Omega().tags.keys():
			if m.startswith('pickledPythonVariablesDictionary'):
				loadVars(m[len('pickledPythonVariableDictionary')+1:])


def SpherePWaveTimeStep(radius,density,young):
	r"""Compute P-wave critical timestep for a single (presumably representative) sphere, using formula for P-Wave propagation speed $\Delta t_{c}=\frac{r}{\sqrt{E/\rho}}$.
	If you want to compute minimum critical timestep for all spheres in the simulation, use :yref:`yade.utils.PWaveTimeStep` instead.
	
	>>> SpherePWaveTimeStep(1e-3,2400,30e9)
	2.8284271247461903e-07
	"""
	from math import sqrt
	return radius/sqrt(young/density)

def randomColor():
	"Return random Vector3 with each component in interval 0…1 (uniform distribution)"
	return Vector3(random.random(),random.random(),random.random())

def typedEngine(name):
	"""Return first engine from current O.engines, identified by its type (as string). For example:

	>>> from yade import utils
	>>> O.engines=[InsertionSortCollider(),NewtonIntegrator(),GravityEngine()]
	>>> utils.typedEngine("NewtonIntegrator") == O.engines[1]
	True
	"""
	return [e for e in Omega().engines if e.__class__.__name__==name][0]

def downCast(obj,newClassName):
	"""Cast given object to class deriving from the same yade root class and copy all parameters from given object.
	Obj should be up in the inheritance tree, otherwise some attributes may not be defined in the new class."""
	return obj.__class__(newClassName,dict([ (key,obj[key]) for key in obj.keys() ]))


def defaultMaterial():
	"""Return default material, when creating bodies with :yref:`yade.utils.sphere` and friends, material is unspecified and there is no shared material defined yet. By default, this function returns::

		FrictMat(density=1e3,young=1e7,poisson=.3,frictionAngle=.5,label='defaultMat')
	"""
	return FrictMat(density=1e3,young=1e7,poisson=.3,frictionAngle=.5,label='defaultMat')

def _commonBodySetup(b,volume,geomInertia,material,noBound=False,resetState=True):
	"""Assign common body parameters."""
	if isinstance(material,int):
		if material<0 and len(O.materials)==0: O.materials.append(defaultMaterial());
		b.mat=O.materials[material]
	elif isinstance(material,str): b.mat=O.materials[material]
	elif isinstance(material,Material): b.mat=material
	elif callable(material): b.mat=material()
	else: raise TypeError("The 'material' argument must be None (for defaultMaterial), string (for shared material label), int (for shared material id) or Material instance.");
	## resets state (!!)
	if resetState: b.state=b.mat.newAssocState()
	mass=volume*b.mat.density
	b.state.mass,b.state.inertia=mass,geomInertia*b.mat.density
	b.bounded=(not noBound)

def sphere(center,radius,dynamic=True,wire=False,color=None,highlight=False,material=-1,mask=1):
	"""Create sphere with given parameters; mass and inertia computed automatically.

	Last assigned material is used by default (*material*=-1), and utils.defaultMaterial() will be used if no material is defined at all.

	:param Vector3 center: center
	:param float radius: radius
	:param Vector3-or-None: body's color, as normalized RGB; random color will be assigned if ``None`.
	:param material:
		specify :yref:`Body.material`; different types are accepted:
			* int: O.materials[material] will be used; as a special case, if material==-1 and there is no shared materials defined, utils.defaultMaterial() will be assigned to O.materials[0]
			* string: label of an existing material that will be used
			* :yref:`Material` instance: this instance will be used
			* callable: will be called without arguments; returned Material value will be used (Material factory object, if you like)
	:param int mask: :yref:`Body.mask` for the body

	:return:
		A Body instance with desired characteristics.


	Creating default shared material if none exists neither is given::

		>>> O.reset()
		>>> from yade import utils
		>>> len(O.materials)
		0
		>>> s0=utils.sphere([2,0,0],1)
		>>> len(O.materials)
		1

	Instance of material can be given::

		>>> s1=utils.sphere([0,0,0],1,wire=False,color=(0,1,0),material=ElastMat(young=30e9,density=2e3))
		>>> s1.shape.wire
		False
		>>> s1.shape.color
		Vector3(0,1,0)
		>>> s1.mat.density
		2000.0

	Material can be given by label::

		>>> O.materials.append(FrictMat(young=10e9,poisson=.11,label='myMaterial'))
		1
		>>> s2=utils.sphere([0,0,2],1,material='myMaterial')
		>>> s2.mat.label
		'myMaterial'
		>>> s2.mat.poisson
		0.11

	Finally, material can be a callable object (taking no arguments), which returns a Material instance.
	Use this if you don't call this function directly (for instance, through yade.pack.randomDensePack), passing
	only 1 *material* parameter, but you don't want material to be shared.

	For instance, randomized material properties can be created like this:

		>>> import random
		>>> def matFactory(): return ElastMat(young=1e10*random.random(),density=1e3+1e3*random.random())
		... 
		>>> s3=utils.sphere([0,2,0],1,material=matFactory)
		>>> s4=utils.sphere([1,2,0],1,material=matFactory)

	"""
	b=Body()
	b.shape=Sphere(radius=radius,color=color if color else randomColor(),wire=wire,highlight=highlight)
	V=(4./3)*math.pi*radius**3
	geomInert=(2./5.)*V*radius**2
	_commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material)
	b.state.pos=b.state.refPos=center
	b.dynamic=dynamic
	b.aspherical=False
	b.mask=mask
	return b

def box(center,extents,orientation=[1,0,0,0],dynamic=True,wire=False,color=None,highlight=False,material=-1,mask=1):
	"""Create box (cuboid) with given parameters.

	:param Vector3 extents: half-sizes along x,y,z axes
	
	See :yref:`yade.utils.sphere`'s documentation for meaning of other parameters."""
	b=Body()
	b.shape=Box(extents=extents,color=color if color else randomColor(),wire=wire,highlight=highlight)
	V=8*extents[0]*extents[1]*extents[2]
	geomInert=Vector3(4*(extents[1]**2+extents[2]**2),4*(extents[0]**2+extents[2]**2),4*(extents[0]**2+extents[1]**2))
	_commonBodySetup(b,V,geomInert,material)
	b.state.pos=b.state.refPos=center
	b.dynamic=dynamic
	b.mask=mask
	b.aspherical=True
	return b


def chainedCylinder(begin=Vector3(0,0,0),end=Vector3(1.,0.,0.),radius=0.2,dynamic=True,wire=False,color=None,highlight=False,material=-1,mask=1):
	"""Create and chain a MinkCylinder with given parameters. This shape is the Minkowski sum of line and sphere.

	:param Real radius: radius of sphere in the Minkowski sum.
	:param Vector3 begin: first point positioning the line in the Minkowski sum
	:param Vector3 last: last point positioning the line in the Minkowski sum

	In order to build a correct chain, last point of element of rank N must correspond to first point of element of rank N+1 in the same chain (with some tolerance, since bounding boxes will be used to create connections.
	
	:return: Body object with the :yref:`ChainedCylinder` :yref:`shape<Body.shape>`.
	"""
	segment=end-begin
	b=Body()
	b.shape=ChainedCylinder(radius=radius,length=segment.norm(),color=color if color else randomColor(),wire=wire,highlight=highlight)
	b.shape.segment=segment
	V=2*(4./3)*math.pi*radius**3
	geomInert=(2./5.)*V*radius**2+b.shape.length*b.shape.length*2*(4./3)*math.pi*radius**3
	b.state=ChainedState()
	b.state.addToChain(b.id)
	#b.state.blockedDOFs=['rx','ry','rz']
	_commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,resetState=False)
	b.state.pos=b.state.refPos=begin
	b.dynamic=dynamic
	b.mask=mask
	b.bound=Aabb(color=[0,1,0])
	b.state.ori=b.state.ori.setFromTwoVectors(Vector3(0.,0.,1.),segment)
	return b

def wall(position,axis,sense=0,color=None,material=-1,mask=1):
	"""Return ready-made wall body.

	:Parameters:
	:param float-or-Vector3 position: center of the wall. If float, it is the position along given axis, the other 2 components being zero
	:param ∈{0,1,2} axis: orientation of the wall normal (0,1,2) for x,y,z (sc. planes yz, xz, xy)
	:param ∈{-1,0,1} sense: sense in which to interact (0: both, -1: negative, +1: positive; see :yref:`Wall`)

	See :yref:`yade.utils.sphere`'s documentation for meaning of other parameters."""
	b=Body()
	b.shape=Wall(sense=sense,axis=axis,color=color if color else randomColor())
	_commonBodySetup(b,0,Vector3(0,0,0),material)
	if isinstance(position,(int,long,float)):
		pos2=Vector3(0,0,0); pos2[axis]=position
	else: pos2=position
	b.state.pos=b.state.refPos=pos2
	b.dynamic=False
	b.aspherical=True # never used, since the wall is not dynamic
	b.mask=mask
	return b

def facet(vertices,dynamic=False,wire=True,color=None,highlight=False,noBound=False,material=-1,mask=1):
	"""Create facet with given parameters.

	:Parameters:
	:param [Vector3,Vector3,Vector3] vertices: coordinates of vertices in the global coordinate system.
	:param bool wire: if ``True``, facets are shown as skeleton; otherwise facets are filled
	:param bool noBound: set :yref:`Body.bounded`
	:param Vector3-or-None color: color of the facet; random color will be assigned if ``None``.
	
	See :yref:`yade.utils.sphere`'s documentation for meaning of other parameters."""
	b=Body()
	center=inscribedCircleCenter(vertices[0],vertices[1],vertices[2])
	vertices=Vector3(vertices[0])-center,Vector3(vertices[1])-center,Vector3(vertices[2])-center
	b.shape=Facet(color=color if color else randomColor(),wire=wire,highlight=highlight,vertices=vertices)
	_commonBodySetup(b,0,Vector3(0,0,0),material,noBound=noBound)
	b.state.pos=b.state.refPos=center
	b.dynamic=dynamic
	b.aspherical=True
	b.mask=mask
	return b

def facetBox(center,extents,orientation=Quaternion.Identity,wallMask=63,**kw):
	"""
	Create arbitrarily-aligned box composed of facets, with given center, extents and orientation.
	If any of the box dimensions is zero, corresponding facets will not be created. The facets are oriented outwards from the box.

	
	:Parameters:
			`center`: Vector3
				center of the created box
			`extents`: (eX,eY,eZ)
				lengths of the box sides
			`orientation`: quaternion
				orientation of the box
			`wallMask`: bitmask
				determines which walls will be created, in the order -x (1), +x (2), -y (4), +y (8), -z (16), +z (32). The numbers are ANDed; the default 63 means to create all walls;
			`**kw`: (unused keyword arguments)
				passed to utils.facet

	:Returns: list of facets forming the box.
	"""
	
	
	#Defense from zero dimensions
	if (wallMask>63):
		print "wallMask must be 63 or less"
		wallMask=63
	if (extents[0]==0):
		wallMask=1
	elif (extents[1]==0):
		wallMask=4
	elif (extents[2]==0):
		wallMask=16
	if (((extents[0]==0) and (extents[1]==0)) or ((extents[0]==0) and (extents[2]==0)) or ((extents[1]==0) and (extents[2]==0))):
		raise RuntimeError("Please, specify at least 2 none-zero dimensions in extents!");
	# ___________________________
	
	mn,mx=[-extents[i] for i in 0,1,2],[extents[i] for i in 0,1,2]
	def doWall(a,b,c,d):
		return [facet((a,b,c),**kw),facet((a,c,d),**kw)]
	ret=[]

	A=orientation*Vector3(mn[0],mn[1],mn[2])+center
	B=orientation*Vector3(mx[0],mn[1],mn[2])+center
	C=orientation*Vector3(mx[0],mx[1],mn[2])+center
	D=orientation*Vector3(mn[0],mx[1],mn[2])+center
	E=orientation*Vector3(mn[0],mn[1],mx[2])+center
	F=orientation*Vector3(mx[0],mn[1],mx[2])+center
	G=orientation*Vector3(mx[0],mx[1],mx[2])+center
	H=orientation*Vector3(mn[0],mx[1],mx[2])+center
	if wallMask&1:  ret+=doWall(A,D,H,E)
	if wallMask&2:  ret+=doWall(B,F,G,C)
	if wallMask&4:  ret+=doWall(A,E,F,B)
	if wallMask&8:  ret+=doWall(D,C,G,H)
	if wallMask&16: ret+=doWall(A,B,C,D)
	if wallMask&32: ret+=doWall(E,H,G,F)
	return ret
	
def facetCylinder(center,radius,height,orientation=Quaternion.Identity,segmentsNumber=10,wallMask=7,angleRange=None,closeGap=False,**kw):
	"""
	Create arbitrarily-aligned cylinder composed of facets, with given center, radius, height and orientation.
	Return List of facets forming the cylinder;
	
	:param Vector3 center: center of the created cylinder
	:param float radius:  cylinder radius
	:param float height: cylinder height
	:param Quaternion orientation: orientation of the cylinder; the reference orientation has axis along the $+x$ axis.
	:param int segmentsNumber: number of edges on the cylinder surface (>=5)
	:param bitmask wallMask: determines which walls will be created, in the order up (1), down (2), side (4). The numbers are ANDed; the default 7 means to create all walls
	:param (θmin,Θmax) angleRange: allows one to create only part of cylinder by specifying range of angles; if ``None``, (0,2*pi) is assumed.
	:param bool closeGap: close range skipped in angleRange with triangular facets at cylinder bases.
	:param **kw: (unused keyword arguments) passed to utils.facet;
	"""
	
	# check zero dimentions
	if (segmentsNumber<3): raise RuntimeError("The segmentsNumber should be at least 3");
	if (height<=0): raise RuntimeError("The height should have the positive value");
	if (radius<=0): raise RuntimeError("The radius should have the positive value");

	import numpy
	if angleRange==None: angleRange=(0,2*math.pi)
	if isinstance(angleRange,float):
		print u'WARNING: utils.facetCylinder,angleRange should be (Θmin,Θmax), not just Θmax (one number), update your code.'
		angleRange=(0,angleRange)
	anglesInRad = numpy.linspace(angleRange[0], angleRange[1], segmentsNumber+1, endpoint=True)
	P1=[]; P2=[]
	P1.append(Vector3(0,0,-height/2))
	P2.append(Vector3(0,0,+height/2))
	
	for i in anglesInRad:
		X=radius*math.cos(i)
		Y=radius*math.sin(i)
		P1.append(Vector3(X,Y,-height/2))
		P2.append(Vector3(X,Y,+height/2))
	for i in range(0,len(P1)):
		P1[i]=orientation*P1[i]+center
		P2[i]=orientation*P2[i]+center
		
	ret=[]
	for i in range(2,len(P1)):
		if wallMask&2:
			ret.append(facet((P1[0],P1[i],P1[i-1]),**kw))
		if wallMask&1:
			ret.append(facet((P2[0],P2[i-1],P2[i]),**kw))
		if wallMask&4:
			ret.append(facet((P1[i],P2[i],P2[i-1]),**kw))
			ret.append(facet((P2[i-1],P1[i-1],P1[i]),**kw))
	if closeGap and (angleRange[0]%(2*math.pi))!=(angleRange[1]%(2*math.pi)): # some part is skipped
		pts=[(radius*math.cos(angleRange[i]),radius*math.sin(angleRange[i])) for i in (0,1)]
		for Z in -height/2,height/2:
			#print (pts[0][0],pts[0][1],Z),(pts[1][0],pts[1][1],Z),(0,0,Z)
			pp=[(pts[0][0],pts[0][1],Z),(pts[1][0],pts[1][1],Z),(0,0,Z)]
			pp=[orientation*p+center for p in pp]
			ret.append(facet(pp,**kw))
		
	return ret
	
	
def aabbWalls(extrema=None,thickness=None,oversizeFactor=1.5,**kw):
	"""Return 6 boxes that will wrap existing packing as walls from all sides;
	extrema are extremal points of the Aabb of the packing (will be calculated if not specified)
	thickness is wall thickness (will be 1/10 of the X-dimension if not specified)
	Walls will be enlarged in their plane by oversizeFactor.
	returns list of 6 wall Bodies enclosing the packing, in the order minX,maxX,minY,maxY,minZ,maxZ.
	"""
	walls=[]
	if not extrema: extrema=aabbExtrema()
	if not thickness: thickness=(extrema[1][0]-extrema[0][0])/10.
	for axis in [0,1,2]:
		mi,ma=extrema
		center=[(mi[i]+ma[i])/2. for i in range(3)]
		extents=[.5*oversizeFactor*(ma[i]-mi[i]) for i in range(3)]
		extents[axis]=thickness/2.
		for j in [0,1]:
			center[axis]=extrema[j][axis]+(j-.5)*thickness
			walls.append(box(center=center,extents=extents,dynamic=False,**kw))
			walls[-1].shape.wire=True
	return walls


def aabbDim(cutoff=0.,centers=False):
	"""Return dimensions of the axis-aligned bounding box, optionally with relative part *cutoff* cut away."""
	a=aabbExtrema(cutoff,centers)
	return (a[1][0]-a[0][0],a[1][1]-a[0][1],a[1][2]-a[0][2])

def aabbExtrema2d(pts):
	"""Return 2d bounding box for a sequence of 2-tuples."""
	inf=float('inf')
	min,max=[inf,inf],[-inf,-inf]
	for pt in pts:
		if pt[0]<min[0]: min[0]=pt[0]
		elif pt[0]>max[0]: max[0]=pt[0]
		if pt[1]<min[1]: min[1]=pt[1]
		elif pt[1]>max[1]: max[1]=pt[1]
	return tuple(min),tuple(max)

def perpendicularArea(axis):
	"""Return area perpendicular to given axis (0=x,1=y,2=z) generated by bodies
	for which the function consider returns True (defaults to returning True always)
	and which is of the type :yref:`Sphere`.
	"""
	ext=aabbExtrema()
	other=((axis+1)%3,(axis+2)%3)
	return (ext[1][other[0]]-ext[0][other[0]])*(ext[1][other[1]]-ext[0][other[1]])

def fractionalBox(fraction=1.,minMax=None):
	"""retrurn (min,max) that is the original minMax box (or aabb of the whole simulation if not specified)
	linearly scaled around its center to the fraction factor"""
	if not minMax: minMax=aabbExtrema()
	half=[.5*(minMax[1][i]-minMax[0][i]) for i in [0,1,2]]
	return (tuple([minMax[0][i]+(1-fraction)*half[i] for i in [0,1,2]]),tuple([minMax[1][i]-(1-fraction)*half[i] for i in [0,1,2]]))


def randomizeColors(onlyDynamic=False):
	"""Assign random colors to :yref:`Shape::color`.

	If onlyDynamic is true, only dynamic bodies will have the color changed.
	"""
	for b in O.bodies:
		color=(random.random(),random.random(),random.random())
		if b.dynamic or not onlyDynamic: b.shape.color=color

def avgNumInteractions(cutoff=0.,skipFree=False):
	r"""Return average numer of interactions per particle, also known as *coordination number* $Z$. This number is defined as 

	.. math:: Z=2C/N

	where $C$ is number of contacts and $N$ is number of particles.

	With *skipFree*, particles not contributing to stable state of the packing are skipped, following equation (8) given in [Thornton2000]_:

	.. math:: Z_m=\frac{2C-N_1}{N-N_0-N_1}

	:param cutoff: cut some relative part of the sample's bounding box away.
	:param skipFree: 	"""
	if cutoff==0 and not skipFree: return 2*O.interactions.countReal()*1./len(O.bodies)
	else:
		nums,counts=bodyNumInteractionsHistogram(aabbExtrema(cutoff))
		## CC is 2*C
		CC=sum([nums[i]*counts[i] for i in range(len(nums))]); N=sum(counts)
		if not skipFree: return CC*1./N
		## find bins with 0 and 1 spheres
		N0=0 if (0 not in nums) else counts[nums.index(0)]
		N1=0 if (1 not in nums) else counts[nums.index(1)]
		NN=N-N0-N1
		return (CC-N1)*1./NN if NN>0 else float('nan')


def plotNumInteractionsHistogram(cutoff=0.):
	"Plot histogram with number of interactions per body, optionally cutting away *cutoff* relative axis-aligned box from specimen margin."
	nums,counts=bodyNumInteractionsHistogram(aabbExtrema(cutoff))
	import pylab
	pylab.bar(nums,counts)
	pylab.title('Number of interactions histogram, average %g (cutoff=%g)'%(avgNumInteractions(cutoff),cutoff))
	pylab.xlabel('Number of interactions')
	pylab.ylabel('Body count')
	pylab.show()

def plotDirections(aabb=(),mask=0,bins=20,numHist=True,noShow=False):
	"""Plot 3 histograms for distribution of interaction directions, in yz,xz and xy planes and
	(optional but default) histogram of number of interactions per body.
	
	:returns: If *noShow* is ``False``, displays the figure and returns nothing. If *noShow*, the figure object is returned without being displayed (works the same way as :yref:`yade.plot.plot`).
	"""
	import pylab,math
	from yade import utils
	for axis in [0,1,2]:
		d=utils.interactionAnglesHistogram(axis,mask=mask,bins=bins,aabb=aabb)
		fc=[0,0,0]; fc[axis]=1.
		subp=pylab.subplot(220+axis+1,polar=True);
		# 1.1 makes small gaps between values (but the column is a bit decentered)
		pylab.bar(d[0],d[1],width=math.pi/(1.1*bins),fc=fc,alpha=.7,label=['yz','xz','xy'][axis])
		#pylab.title(['yz','xz','xy'][axis]+' plane')
		pylab.text(.5,.25,['yz','xz','xy'][axis],horizontalalignment='center',verticalalignment='center',transform=subp.transAxes,fontsize='xx-large')
	if numHist:
		pylab.subplot(224,polar=False)
		nums,counts=utils.bodyNumInteractionsHistogram(aabb if len(aabb)>0 else utils.aabbExtrema())
		avg=sum([nums[i]*counts[i] for i in range(len(nums))])/(1.*sum(counts))
		pylab.bar(nums,counts,fc=[1,1,0],alpha=.7,align='center')
		pylab.xlabel('Interactions per body (avg. %g)'%avg)
		pylab.axvline(x=avg,linewidth=3,color='r')
		pylab.ylabel('Body count')
	if noShow: return pylab.gcf()
	else: pylab.show()


def encodeVideoFromFrames(*args,**kw):
	"|ydeprecated|"
	_deprecatedUtilsFunction('utils.encodeVideoFromFrames','utils.makeVideo')
	return makeVideo(*args,**kw)

def makeVideo(frameSpec,out,renameNotOverwrite=True,fps=24,bps=2400):
	"""Create a video from external image files using `mencoder <http://www.mplayerhq.hu>`__. Two-pass encoding using the default mencoder codec (mpeg4) is performed, running multi-threaded with number of threads equal to number of OpenMP threads allocated for Yade.

	:param frameSpec: wildcard | sequence of filenames. If list or tuple, filenames to be encoded in given order; otherwise wildcard understood by mencoder's mf:// URI option (shell wildcards such as ``/tmp/snap-*.png`` or and printf-style pattern like ``/tmp/snap-%05d.png``)
	:param str out: file to save video into
	:param bool renameNotOverwrite: if True, existing same-named video file will have -*number* appended; will be overwritten otherwise.
	:param int fps: Frames per second (``-mf fps=…``)
	:param int bps: Bitrate (``-lavcopts vbitrate=…``)
	"""
	import os,os.path,subprocess
	if renameNotOverwrite and os.path.exists(out):
		i=0
		while(os.path.exists(out+"~%d"%i)): i+=1
		os.rename(out,out+"~%d"%i); print "Output file `%s' already existed, old file renamed to `%s'"%(out,out+"~%d"%i)
	if isinstance(frameSpec,list) or isinstance(frameSpec,tuple): frameSpec=','.join(frameSpec)
	for passNo in (1,2):
		cmd=['mencoder','mf://%s'%frameSpec,'-mf','fps=%d'%int(fps),'-ovc','lavc','-lavcopts','vbitrate=%d:vpass=%d:threads=%d:%s'%(int(bps),passNo,O.numThreads,'turbo' if passNo==1 else ''),'-o',('/dev/null' if passNo==1 else out)]
		print 'Pass %d:'%passNo,' '.join(cmd)
		ret=subprocess.call(cmd)
		if ret!=0: raise RuntimeError("Error when running mencoder.")


def replaceCollider(colliderEngine):
	"""Replaces collider (Collider) engine with the engine supplied. Raises error if no collider is in engines."""
	colliderIdx=-1
	for i,e in enumerate(O.engines):
		if O.isChildClassOf(e.__class__.__name__,"Collider"):
			colliderIdx=i
			break
	if colliderIdx<0: raise RuntimeError("No Collider found within O.engines.")
	O.engines=O.engines[:colliderIdx]+[colliderEngine]+O.engines[colliderIdx+1:]

def _procStatus(name):
	import os
	for l in open('/proc/%d/status'%os.getpid()):
		if l.split(':')[0]==name: return l
	raise "No such line in /proc/[pid]/status: "+name
def vmData():
	"Return memory usage data from Linux's /proc/[pid]/status, line VmData."
	l=_procStatus('VmData'); ll=l.split(); assert(ll[2]=='kB')
	return int(ll[1])

def uniaxialTestFeatures(filename=None,areaSections=10,axis=-1,**kw):
	"""Get some data about the current packing useful for uniaxial test:

#. Find the dimensions that is the longest (uniaxial loading axis)

#. Find the minimum cross-section area of the specimen by examining several (areaSections) sections perpendicular to axis, computing area of the convex hull for each one. This will work also for non-prismatic specimen.

#. Find the bodies that are on the negative/positive boundary, to which the straining condition should be applied.

:parameters:
	`filename`:
		if given, spheres will be loaded from this file (ASCII format); if not, current simulation will be used.
	`areaSection`:
		number of section that will be used to estimate cross-section
	`axis`:
		if given, force strained axis, rather than computing it from predominant length

:return: dictionary with keys 'negIds', 'posIds', 'axis', 'area'.

.. warning::
	The function :yref:`yade.utils.approxSectionArea` uses convex hull algorithm to find the area, but the implementation is reported to be *buggy* (bot works in some cases). Always check this number, or fix the convex hull algorithm (it is documented in the source, see :ysrc:`py/_utils.cpp`).

	"""
	if filename: ids=spheresFromFile(filename,**kw)
	else: ids=[b.id for b in O.bodies]
	mm,mx=aabbExtrema()
	dim=aabbDim();
	if axis<0: axis=list(dim).index(max(dim)) # list(dim) for compat with python 2.5 which didn't have index defined for tuples yet (appeared in 2.6 first)
	assert(axis in (0,1,2))
	import numpy
	areas=[approxSectionArea(coord,axis) for coord in numpy.linspace(mm[axis],mx[axis],num=10)[1:-1]]
	negIds,posIds=negPosExtremeIds(axis=axis,distFactor=2.2)
	return {'negIds':negIds,'posIds':posIds,'axis':axis,'area':min(areas)}

def NormalRestitution2DampingRate(en):
	r"""Compute the normal damping rate as a function of the normal coefficient of restitution $e_n$. For $e_n\in\langle0,1\rangle$ damping rate equals
	
	.. math:: -\frac{\log e_n}{\sqrt{e_n^2+\pi^2}}

	"""
	if en == 0.0: return 0.999999999
	if en == 1.0: return 0.0
	from math import sqrt,log,pi
	ln_en = math.log(en)
	return (-ln_en/math.sqrt((math.pow(ln_en,2) + math.pi*math.pi)))

def xMirror(half):
	"""Mirror a sequence of 2d points around the x axis (changing sign on the y coord).
The sequence should start up and then it will wrap from y downwards (or vice versa).
If the last point's x coord is zero, it will not be duplicated."""
	return list(half)+[(x,-y) for x,y in reversed(half[:-1] if half[-1][0]==0 else half)]

#############################
##### deprecated functions


def _deprecatedUtilsFunction(old,new):
	"Wrapper for deprecated functions, example below."
	import warnings
	warnings.warn('Function utils.%s is deprecated, use %s instead.'%(old,new),stacklevel=2,category=DeprecationWarning)

# example of _deprecatedUtilsFunction usage:
#
# def import_mesh_geometry(*args,**kw):
#    "|ydeprecated|"
#    _deprecatedUtilsFunction('import_mesh_geometry','yade.import.gmsh')
#    import yade.ymport
#    return yade.ymport.stl(*args,**kw)


class TableParamReader():
	"""Class for reading simulation parameters from text file.
	
Each parameter is represented by one column, each parameter set by one line. Colums are separated by blanks (no quoting).

First non-empty line contains column titles (without quotes).
You may use special column named 'description' to describe this parameter set;
if such colum is absent, description will be built by concatenating column names and corresponding values (``param1=34,param2=12.22,param4=foo``)

* from columns ending in ``!`` (the ``!`` is not included in the column name)
* from all columns, if no columns end in ``!``.

Empty lines within the file are ignored (although counted); ``#`` starts comment till the end of line. Number of blank-separated columns must be the same for all non-empty lines.

A special value ``=`` can be used instead of parameter value; value from the previous non-empty line will be used instead (works recursively).

This class is used by :yref:`yade.utils.readParamsFromTable`.
	"""
	def __init__(self,file):
		"Setup the reader class, read data into memory."
		import re
		# read file in memory, remove newlines and comments; the [''] makes lines 1-indexed
		ll=[re.sub('\s*#.*','',l[:-1]) for l in ['']+open(file,'r').readlines()]
		# usable lines are those that contain something else than just spaces
		usableLines=[i for i in range(len(ll)) if not re.match(r'^\s*(#.*)?$',ll[i])]
		headings=ll[usableLines[0]].split()
		# use all values of which heading has ! after its name to build up the description string
		# if there are none, use all columns
		if not 'description' in headings:
			bangHeads=[h[:-1] for h in headings if h[-1]=='!'] or headings
			headings=[(h[:-1] if h[-1]=='!' else h) for h in headings]
		usableLines=usableLines[1:] # and remove headinds from usableLines
		values={}
		for l in usableLines:
			val={}
			for i in range(len(headings)):
				val[headings[i]]=ll[l].split()[i]
			values[l]=val
		lines=values.keys(); lines.sort()
		# replace '=' by previous value of the parameter
		for i,l in enumerate(lines):
			for j in values[l].keys():
				if values[l][j]=='=':
					try:
						values[l][j]=values[lines[i-1]][j]
					except IndexError,KeyError:
						raise RuntimeError("The = specifier on line %d refers to nonexistent value on previous line?"%l)
		#import pprint; pprint.pprint(headings); pprint.pprint(values)
		# add descriptions, but if they repeat, append line number as well
		if not 'description' in headings:
			descs=set()
			for l in lines:
				dd=','.join(head.replace('!','')+'='+('%g'%values[head] if isinstance(values[l][head],float) else str(values[l][head])) for head in bangHeads).replace("'",'').replace('"','')
				if dd in descs: dd+='__line=%d__'%l
				values[l]['description']=dd
				descs.add(dd)
		self.values=values
				
	def paramDict(self):
		"""Return dictionary containing data from file given to constructor. Keys are line numbers (which might be non-contiguous and refer to real line numbers that one can see in text editors), values are dictionaries mapping parameter names to their values given in the file. The special value '=' has already been interpreted, ``!`` (bangs) (if any) were already removed from column titles, ``description`` column has already been added (if absent)."""
		return self.values
		
if __name__=="__main__":
	tryTable="""head1 important2! !OMP_NUM_THREADS! abcd
	1 1.1 1.2 1.3
	'a' 'b' 'c' 'd'  ### comment
	
	# empty line
	1 = = g
"""
	file='/tmp/try-tbl.txt'
	f=open(file,'w')
	f.write(tryTable)
	f.close()
	from pprint import *
	pprint(TableParamReader(file).paramDict())
	
def runningInBatch():
	'Tell whether we are running inside the batch or separately.'
	import os
	return os.environ.has_key('PARAM_TABLE')

def waitIfBatch():
	'Block the simulation if running inside a batch. Typically used at the end of script so that it does not finish prematurely in batch mode (the execution would be ended in such a case).'
	if runningInBatch(): O.wait()

def readParamsFromTable(tableFileLine=None,noTableOk=False,unknownOk=False,**kw):
	"""
	Read parameters from a file and assign them to __builtin__ variables.

	The format of the file is as follows (commens starting with # and empty lines allowed)::
		
		# commented lines allowed anywhere
		name1 name2 … # first non-blank line are column headings
					# empty line is OK, with or without comment
		val1  val2  … # 1st parameter set
		val2  val2  … # 2nd 
		…

	Assigned tags:

	* *description* column is assigned to Omega().tags['description']; this column is synthesized if absent (see :yref:`yade.utils.TableParamReader`);
	* ``Omega().tags['params']="name1=val1,name2=val2,…"``
	* ``Omega().tags['defaultParams']="unassignedName1=defaultValue1,…"``
	* ``Omega().tags['d.id']=O.tags['id']+'.'+O.tags['description']``
	* ``Omega().tags['id.d']=O.tags['description']+'.'+O.tags['id']``

	All parameters (default as well as settable) are saved using :yref:`yade.utils.saveVars`\ ``('table')``.

	:parameters:
		`tableFile`:
			text file (with one value per blank-separated columns)
		`tableLine`:
			number of line where to get the values from.
		`noTableOk`: bool
			do not raise exception if the file cannot be open; use default values
		`unknownOk`: bool
			do not raise exception if unknown column name is found in the file; assign it as well

	:return:
		number of assigned parameters.

	"""
	tagsParams=[]
	# dictParams is what eventually ends up in yade.params.table (default+specified values)
	dictDefaults,dictParams,dictAssign={},{},{}
	import os, __builtin__,re,math
	if not tableFileLine and not os.environ.has_key('PARAM_TABLE'):
		if not noTableOk: raise EnvironmentError("PARAM_TABLE is not defined in the environment")
		O.tags['line']='l!'
	else:
		if not tableFileLine: tableFileLine=os.environ['PARAM_TABLE']
		env=tableFileLine.split(':')
		tableFile,tableLine=env[0],int(env[1])
		allTab=TableParamReader(tableFile).paramDict()
		if not allTab.has_key(tableLine): raise RuntimeError("Table %s doesn't contain valid line number %d"%(tableFile,tableLine))
		vv=allTab[tableLine]
		O.tags['line']='l%d'%tableLine
		O.tags['description']=vv['description']
		O.tags['id.d']=O.tags['id']+'.'+O.tags['description']; O.tags['d.id']=O.tags['description']+'.'+O.tags['id']
		# assign values specified in the table to python vars
		# !something cols are skipped, those are env vars we don't treat at all (they are contained in description, though)
		for col in vv.keys():
			if col=='description' or col[0]=='!': continue
			if col not in kw.keys() and (not unknownOk): raise NameError("Parameter `%s' has no default value assigned"%names[i])
			if vv[col]=='*': vv[col]=kw[col] # use default value for * in the table
			elif col in kw.keys(): kw.pop(col) # remove the var from kw, so that it contains only those that were default at the end of this loop
			#print 'ASSIGN',col,vv[col]
			tagsParams+=['%s=%s'%(col,vv[col])];
			dictParams[col]=eval(vv[col],math.__dict__)
	# assign remaining (default) keys to python vars
	defaults=[]
	for k in kw.keys():
		dictDefaults[k]=kw[k]
		defaults+=["%s=%s"%(k,kw[k])]; 
	O.tags['defaultParams']=",".join(defaults)
	O.tags['params']=",".join(tagsParams)
	dictParams.update(dictDefaults)
	saveVars('table',loadNow=True,**dictParams)
	return len(tagsParams)


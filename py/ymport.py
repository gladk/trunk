"""
Import geometry from various formats ('import' is python keyword, hence the name 'ymport').
"""

from yade.wrapper import *
from miniEigen import *
from yade import utils

def textExt(fileName,format='x_y_z_r',shift=[0.0,0.0,0.0],scale=1.0,**kw):
	"""Load sphere coordinates from file in specific format, create spheres, insert them to the simulation.
	
	:Parameters:
		`filename`: string
		`format`:
			the name of output format. Supported `x_y_z_r`(default), `x_y_z_r_matId`
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.sphere`
	:Returns: list of spheres.
	Lines starting with # are skipped
	"""
	infile = open(fileName,"r")
	lines = infile.readlines()
	infile.close()
	ret=[]
	for line in lines:
		data = line.split()
		if (data[0] == "#format"):
			format=data[1]
			continue
		elif (data[0][0] == "#"): continue
		
		if (format=='x_y_z_r'):
			ret.append(utils.sphere([shift[0]+scale*float(data[0]),shift[1]+scale*float(data[1]),shift[2]+scale*float(data[2])],scale*float(data[3]),**kw))
			
		elif (format=='x_y_z_r_matId'):
			ret.append(utils.sphere([shift[0]+scale*float(data[0]),shift[1]+scale*float(data[1]),shift[2]+scale*float(data[2])],scale*float(data[3]),material=int(data[4]),**kw))
			
		else:
			raise RuntimeError("Please, specify a correct format output!");
	return ret

def text(fileName,shift=[0.0,0.0,0.0],scale=1.0,**kw):
	"""Load sphere coordinates from file, create spheres, insert them to the simulation.

	:Parameters:
		`filename`: string
			file which has 4 colums [x, y, z, radius].
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.sphere`
	:Returns: list of spheres.
	Lines starting with # are skipped
	"""
	
	return textExt(fileName=fileName,format='x_y_z_r',shift=shift,scale=scale,**kw)
	
def stl(file, dynamic=False,wire=True,color=None,highlight=False,noBound=False,material=-1):
	""" Import geometry from stl file, return list of created facets."""
	imp = STLImporter()
	facets=imp.ymport(file)
	for b in facets:
		b.dynamic=dynamic
		b.shape.postProcessAttributes(True)
		b.shape.Color=color if color else utils.randomColor()
		b.shape.wire=wire
		b.shape.highlight=highlight
		pos,ori=b.state.pos,b.state.ori
		utils._commonBodySetup(b,0,Vector3(0,0,0),noBound=noBound,material=material)
		b.state.pos,b.state.ori=pos,ori
	return facets

def gts(meshfile,shift=(0,0,0),scale=1.0,**kw):
	""" Read given meshfile in gts format.

	:Parameters:
		`meshfile`: string
			name of the input file.
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.facet`
	:Returns: list of facets.
	"""
	import gts,yade.pack
	surf=gts.read(open(meshfile))
	surf.scale(scale)
	surf.translate(shift) 
	yade.pack.gtsSurface2Facets(surf,**kw)

def gmsh(meshfile="file.mesh",shift=[0.0,0.0,0.0],scale=1.0,orientation=Quaternion().Identity,**kw):
	""" Imports geometry from mesh file and creates facets.

	:Parameters:
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`orientation`: quaternion
			orientation of the imported mesh
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.facet`
	:Returns: list of facets forming the specimen.
	
	mesh files can be easily created with `GMSH <http://www.geuz.org/gmsh/>`_.
	Example added to :ysrc:`examples/regular-sphere-pack/regular-sphere-pack.py`
	
	Additional examples of mesh-files can be downloaded from 
	http://www-roc.inria.fr/gamma/download/download.php
	"""
	infile = open(meshfile,"r")
	lines = infile.readlines()
	infile.close()

	nodelistVector3=[]
	findVerticesString=0
	
	while (lines[findVerticesString].split()[0]<>'Vertices'): #Find the string with the number of Vertices
		findVerticesString+=1
	findVerticesString+=1
	numNodes = int(lines[findVerticesString].split()[0])
	
	for i in range(numNodes):
		nodelistVector3.append(Vector3(0.0,0.0,0.0))
	id = 0
	
	qTemp = Quaternion(Vector3(orientation[0],orientation[1],orientation[2]),orientation[3])
	for line in lines[findVerticesString+1:numNodes+findVerticesString+1]:
		data = line.split()
		nodelistVector3[id] = qTemp.Rotate(Vector3(float(data[0])*scale,float(data[1])*scale,float(data[2])*scale))+Vector3(shift[0],shift[1],shift[2])
		id += 1

	
	findTriangleString=findVerticesString+numNodes
	while (lines[findTriangleString].split()[0]<>'Triangles'): #Find the string with the number of Triangles
		findTriangleString+=1
	findTriangleString+=1
	numTriangles = int(lines[findTriangleString].split()[0])

	triList = []
	for i in range(numTriangles):
		triList.append([0,0,0,0])
	
	tid = 0
	for line in lines[findTriangleString+1:findTriangleString+numTriangles+1]:
		data = line.split()
		id1 = int(data[0])-1
		id2 = int(data[1])-1
		id3 = int(data[2])-1
		triList[tid][0] = tid
		triList[tid][1] = id1
		triList[tid][2] = id2
		triList[tid][3] = id3
		tid += 1
		ret=[]
	for i in triList:
		a=nodelistVector3[i[1]]
		b=nodelistVector3[i[2]]
		c=nodelistVector3[i[3]]
		ret.append(utils.facet((nodelistVector3[i[1]],nodelistVector3[i[2]],nodelistVector3[i[3]]),**kw))
	return ret

def gengeoFile(fileName="file.geo",shift=[0.0,0.0,0.0],scale=1.0,**kw):
	""" Imports geometry from LSMGenGeo .geo file and creates spheres.
	
	:Parameters:
		`filename`: string
			file which has 4 colums [x, y, z, radius].
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.sphere`
	:Returns: list of spheres.
	
	LSMGenGeo library allows to create pack of spheres
	with given [Rmin:Rmax] with null stress inside the specimen.
	Can be useful for Mining Rock simulation.
	
	Example: :ysrc:`examples/regular-sphere-pack/regular-sphere-pack.py`, usage of LSMGenGeo library in :ysrc:`scripts/test/genCylLSM.py`.
	
	* https://answers.launchpad.net/esys-particle/+faq/877
	* http://www.access.edu.au/lsmgengeo_python_doc/current/pythonapi/html/GenGeo-module.html
	* https://svn.esscc.uq.edu.au/svn/esys3/lsm/contrib/LSMGenGeo/"""
	from yade.utils import sphere

	infile = open(fileName,"r")
	lines = infile.readlines()
	infile.close()

	numSpheres = int(lines[6].split()[0])
	ret=[]
	for line in lines[7:numSpheres+7]:
		data = line.split()
		ret.append(utils.sphere([shift[0]+scale*float(data[0]),shift[1]+scale*float(data[1]),shift[2]+scale*float(data[2])],scale*float(data[3]),**kw))
	return ret

def gengeo(mntable,shift=Vector3().Zero,scale=1.0,**kw):
	""" Imports geometry from LSMGenGeo library and creates spheres.

	:Parameters:
		`mntable`: mntable
			object, which creates by LSMGenGeo library, see example
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.sphere`
	
	LSMGenGeo library allows to create pack of spheres
	with given [Rmin:Rmax] with null stress inside the specimen.
	Can be useful for Mining Rock simulation.
	
	Example: :ysrc:`examples/regular-sphere-pack/regular-sphere-pack.py`, usage of LSMGenGeo library in :ysrc:`scripts/test/genCylLSM.py`.
	
	* https://answers.launchpad.net/esys-particle/+faq/877
	* http://www.access.edu.au/lsmgengeo_python_doc/current/pythonapi/html/GenGeo-module.html
	* https://svn.esscc.uq.edu.au/svn/esys3/lsm/contrib/LSMGenGeo/"""
	try:
		from GenGeo import MNTable3D,Sphere
	except ImportError:
		from gengeo import MNTable3D,Sphere
	ret=[]
	sphereList=mntable.getSphereListFromGroup(0)
	for i in range(0, len(sphereList)):
		r=sphereList[i].Radius()
		c=sphereList[i].Centre()
		ret.append(utils.sphere([shift[0]+scale*float(c.X()),shift[1]+scale*float(c.Y()),shift[2]+scale*float(c.Z())],scale*float(r),**kw))
	return ret

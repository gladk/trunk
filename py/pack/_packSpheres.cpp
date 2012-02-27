// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

#include<yade/pkg-dem/SpherePack.hpp>
#include<yade/lib-pyutil/doc_opts.hpp>

BOOST_PYTHON_MODULE(_packSpheres){
	python::scope().attr("__doc__")="Creation, manipulation, IO for generic sphere packings.";
	YADE_SET_DOCSTRING_OPTS;
	python::class_<SpherePack>("SpherePack","Set of spheres represented as centers and radii. This class is returned by :yref:`yade.pack.randomDensePack`, :yref:`yade.pack.randomPeriPack` and others. The object supports iteration over spheres, as in \n\n\t>>> sp=SpherePack()\n\t>>> for center,radius in sp: print center,radius\n\n\t>>> for sphere in sp: print sphere[0],sphere[1]   ## same, but without unpacking the tuple automatically\n\n\t>>> for i in range(0,len(sp)): print sp[i][0], sp[i][1]   ## same, but accessing spheres by index\n\n\n.. admonition:: Special constructors\n\n\tConstruct from list of ``[(c1,r1),(c2,r2),…]``. To convert two same-length lists of ``centers`` and ``radii``, construct with ``zip(centers,radii)``.\n",python::init<python::optional<python::list> >(python::args("list"),"Empty constructor, optionally taking list [ ((cx,cy,cz),r), … ] for initial data." ))
		.def("add",&SpherePack::add,"Add single sphere to packing, given center as 3-tuple and radius")
		.def("toList",&SpherePack::toList,"Return packing data as python list.")
		.def("fromList",&SpherePack::fromList,"Make packing from given list, same format as for constructor. Discards current data.")
		.def("fromList",&SpherePack::fromLists,(python::arg("centers"),python::arg("radii")),"Make packing from given list, same format as for constructor. Discards current data.")
		.def("load",&SpherePack::fromFile,(python::arg("fileName")),"Load packing from external text file (current data will be discarded).")
		.def("save",&SpherePack::toFile,(python::arg("fileName")),"Save packing to external text file (will be overwritten).")
		.def("fromSimulation",&SpherePack::fromSimulation,"Make packing corresponding to the current simulation. Discards current data.")
		.def("makeCloud",&SpherePack::makeCloud,(python::arg("minCorner"),python::arg("maxCorner"),python::arg("rMean")=-1,python::arg("rRelFuzz")=0,python::arg("num")=-1,python::arg("periodic")=false,python::arg("porosity")=-1,python::arg("psdSizes")=vector<Real>(),python::arg("psdCumm")=vector<Real>(),python::arg("distributeMass")=false),"Create random loose packing enclosed in box."
		"\nSphere radius distribution can be specified using one of the following ways (they are mutually exclusive):\n\n#. *rMean* and *rRelFuzz* gives uniform radius distribution between *rMean*(1 ± *rRelFuzz*).\n#. *porosity*, *num* and *rRelFuzz* which estimates mean radius so that *porosity* is attained at the end\n#. *psdSizes* and *psdCumm*, two arrays specifying points of the `particle size distribution <http://en.wikipedia.org/wiki/Particle_size_distribution>`__ function.\n\nBy default (with ``distributeMass==False``), the distribution is applied to particle radii. The usual sense of \"particle size distribution\" is the distribution of *mass fraction* (rather than particle count); this can be achieved with ``distributeMass=True``."
		"\n\n:param Vector3 minCorner: lower corner of the box\n:param Vector3 maxCorner: upper corner of the box\n:param float rMean: mean radius or spheres\n:param float rRelFuzz: dispersion of radius relative to rMean\n:param int num: number of spheres to be generated (if negavite, generate as many as possible, ending after a fixed number of tries to place the sphere in space)\n:param bool periodic: whether the packing to be generated should be periodic\n:param float porosity: if specified, estimate mean radius $r_m$ (*rMean* must not be given) using *rRelFuzz* ($z$) and *num* ($N$) so that the porosity given ($\\rho$) is approximately achieved at the end of generation, $r_m=\\sqrt[3]{\\frac{V(1-\\rho)}{\\frac{4}{3}\\pi(1+z^2)N}}$. The value of $\\rho$=0.65 is recommended.\n:param psdSizes: sieve sizes (particle diameters) when particle size distribution (PSD) is specified\n:param psdCumm: cummulative fractions of particle sizes given by *psdSizes*; must be the same length as *psdSizes* and should be non-decreasing\n:param bool distributeMass: if ``True``, given distribution will be used to distribute sphere's mass rather than radius of them.\n:returns: number of created spheres, which can be lower than *num* is the packing is too tight.\n")
		.def("psd",&SpherePack::psd,(python::arg("bins")=10,python::arg("mass")=false),"Return `particle size distribution <http://en.wikipedia.org/wiki/Particle_size_distribution>`__ of the packing.\n:param int bins: number of bins between minimum and maximum diameter\n:param mass: Compute relative mass rather than relative particle count for each bin. Corresponds to :yref:`distributeMass parameter for makeCloud<yade.pack.SpherePack.makeCloud>`.\n:returns: tuple of ``(cumm,edges)``, where ``cumm`` are cummulative fractions for respective diameters  and ``edges`` are those diameter values. Dimension of both arrays is equal to ``bins+1``.")
		// new psd
		.def("particleSD",&SpherePack::particleSD,(python::arg("minCorner"),python::arg("maxCorner"),python::arg("rMean"),python::arg("periodic")=false,python::arg("name"),python::arg("numSph"),python::arg("radii")=vector<Real>(),python::arg("passing")=vector<Real>(),python::arg("passingIsNotPercentageButCount")=false),"Create random packing enclosed in box given by minCorner and maxCorner, containing numSph spheres. Returns number of created spheres, which can be < num if the packing is too tight. The computation is done according to the given psd.")
		.def("makeClumpCloud",&SpherePack::makeClumpCloud,(python::arg("minCorner"),python::arg("maxCorner"),python::arg("clumps"),python::arg("periodic")=false,python::arg("num")=-1),"Create random loose packing of clumps within box given by *minCorner* and *maxCorner*. Clumps are selected with equal probability. At most *num* clumps will be positioned if *num* is positive; otherwise, as many clumps as possible will be put in space, until maximum number of attemps to place a new clump randomly is attained.")
		//
		.def("aabb",&SpherePack::aabb_py,"Get axis-aligned bounding box coordinates, as 2 3-tuples.")
		.def("dim",&SpherePack::dim,"Return dimensions of the packing in terms of aabb(), as a 3-tuple.")
		.def("center",&SpherePack::midPt,"Return coordinates of the bounding box center.")
		.def_readwrite("cellSize",&SpherePack::cellSize,"Size of periodic cell; is Vector3(0,0,0) if not periodic. (Change this property only if you know what you're doing).")
		.def("cellFill",&SpherePack::cellFill,"Repeat the packing (if periodic) so that the results has dim() >= given size. The packing retains periodicity, but changes cellSize. Raises exception for non-periodic packing.")
		.def("cellRepeat",&SpherePack::cellRepeat,"Repeat the packing given number of times in each dimension. Periodicity is retained, cellSize changes. Raises exception for non-periodic packing.")
		.def("relDensity",&SpherePack::relDensity,"Relative packing density, measured as sum of spheres' volumes / aabb volume.\n(Sphere overlaps are ignored.)")
		.def("translate",&SpherePack::translate,"Translate all spheres by given vector.")
		.def("rotate",&SpherePack::rotate,(python::arg("axis"),python::arg("angle")),"Rotate all spheres around packing center (in terms of aabb()), given axis and angle of the rotation.")
		.def("scale",&SpherePack::scale,"Scale the packing around its center (in terms of aabb()) by given factor (may be negative).")
		.def("hasClumps",&SpherePack::hasClumps,"Whether this object contains clumps.")
		.def("getClumps",&SpherePack::getClumps,"Return lists of sphere ids sorted by clumps they belong to. The return value is (standalones,[clump1,clump2,…]), where each item is list of id's of spheres.")
		.def("__len__",&SpherePack::len,"Get number of spheres in the packing")
		.def("__getitem__",&SpherePack::getitem,"Get entry at given index, as tuple of center and radius.")
		.def("__iter__",&SpherePack::getIterator,"Return iterator over spheres.")
		.def_readwrite("psdScaleExponent",&SpherePack::psdScaleExponent,"Exponent used to scale cummulative distribution function values so that standard uniform distribution is mapped to uniform mass distribution. Theoretically, this value is 3, but that usually overfavors small particles. The dafault value is 2.5.")
		;
	python::class_<SpherePack::_iterator>("SpherePackIterator",python::init<SpherePack::_iterator&>())
		.def("__iter__",&SpherePack::_iterator::iter)
		.def("next",&SpherePack::_iterator::next)
	;
}


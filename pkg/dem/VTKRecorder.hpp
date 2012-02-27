#pragma once
#include<yade/pkg-common/PeriodicEngines.hpp>

// multiblock features don't seem to exist prioor to 5.2 
#if (VTK_MAJOR_VERSION==5 && VTK_MINOR_VERSION>=2) || (VTK_MAJOR_VERSION > 5)
	#define YADE_VTK_MULTIBLOCK
#endif


class VTKRecorder: public PeriodicEngine {
	public:
		enum {REC_SPHERES=0,REC_FACETS,REC_COLORS,REC_CPM,REC_INTR,REC_VELOCITY,REC_ID,REC_CLUMPID,REC_SENTINEL,REC_MATERIALID,REC_STRESS,REC_MASK,REC_RPM};
		virtual void action();
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(VTKRecorder,PeriodicEngine,"Engine recording snapshots of simulation into series of *.vtu files, readable by VTK-based postprocessing programs such as Paraview. Both bodies (spheres and facets) and interactions can be recorded, with various vector/scalar quantities that are defined on them.\n\n:yref:`PeriodicEngine.initRun` is initialized to ``True`` automatically.",
		((bool,compress,false,,"Compress output XML files [experimental]."))
		((bool,ascii,false,,"Store data as readable text in the XML file (sets `vtkXMLWriter <http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html>`__ data mode to ``vtkXMLWriter::Ascii``, while the default is ``Appended``"))
		((bool,skipFacetIntr,true,,"Skip interactions with facets, when saving interactions"))
		((bool,skipNondynamic,false,,"Skip non-dynamic spheres (but not facets)."))
		#ifdef YADE_VTK_MULTIBLOCK
			((bool,multiblock,false,,"Use multi-block (``.vtm``) files to store data, rather than separate ``.vtu`` files."))
		#endif
		((string,fileName,"",,"Base file name; it will be appended with {spheres,intrs,facets}-243100.vtu (unless *multiblock* is ``True``) depending on active recorders and step number (243100 in this case). It can contain slashes, but the directory must exist already."))
		((vector<string>,recorders,vector<string>(1,string("all")),,"List of active recorders (as strings). ``all`` (the default value) enables all base and generic recorders.\n\n.. admonition:: Base recorders\n\n\tBase recorders save the geometry (unstructured grids) on which other data is defined. They are implicitly activated by many of the other recorders. Each of them creates a new file (or a block, if :yref:`multiblock <VTKRecorder.multiblock>` is set).\n\n\t``spheres``\n\t\tSaves positions and radii (``radii``) of :yref:`spherical<Sphere>` particles.\n\t``facets``\n\t\tSave :yref:`facets<Facet>` positions (vertices).\n\t``intr``\n\t\tStore interactions as lines between nodes at respective particles positions. Additionally stores magnitude of normal (``forceN``) and shear (``absForceT``) forces on interactions (the :yref:`geom<Interaction.geom> must be of type :yref:`NormShearPhys`). \n\n.. admonition:: Generic recorders\n\n\tGeneric recorders do not depend on specific model being used and save commonly useful data.\n\n\t``id``\n\t\tSaves id's (field ``id``) of spheres; active only if ``spheres`` is active.\n\t``clumpId``\n\t\tSaves id's of clumps to which each sphere belongs (field ``clumpId``); active only if ``spheres`` is active.\n\t``colors``\n\t\tSaves colors of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>` (field ``color``); only active if ``spheres`` or ``facets`` are activated.\n\t``mask``\n\t\tSaves groupMasks of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>` (field ``mask``); only active if ``spheres`` or ``facets`` are activated.\n\t``materialId``\n\t\tSaves materialID of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>`; only active if ``spheres`` or ``facets`` are activated.\n\t``velocity``\n\t\tSaves linear and angular velocities of spherical particles as Vector3 and length(fields ``linVelVec``, ``linVelLen`` and ``angVelVec``, ``angVelLen`` respectively``); only effective with ``spheres``.\n\t``stress``\n\t\tSaves stresses of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>`  as Vector3 and length; only active if ``spheres`` or ``facets`` are activated.\n\n.. admonition:: Specific recorders\n\n\tThe following should only be activated in appropriate cases, otherwise crashes can occur due to violation of type presuppositions.\n\n\t``cpm``\n\t\tSaves data pertaining to the :yref:`concrete model<Law2_Dem3DofGeom_CpmPhys_Cpm>`: ``cpmDamage`` (normalized residual strength averaged on particle), ``cpmSigma`` (stress on particle, normal components), ``cpmTau`` (shear components of stress on particle), ``cpmSigmaM`` (mean stress around particle); ``intr`` is activated automatically by ``cpm``\n\t``rpm``\n\t\tSaves data pertaining to the :yref:`rock particle model<Law2_Dem3DofGeom_RockPMPhys_Rpm>`: ``rpmSpecNum`` shows different pieces of separated stones, only ids. ``rpmSpecMass`` shows masses of separated stones.\n\n"))
		((int,mask,0,,"If mask defined, only bodies with corresponding groupMask will be exported. If 0, all bodies will be exported.")),
		/*ctor*/
		initRun=true;
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(VTKRecorder);

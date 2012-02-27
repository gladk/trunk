.. _chapter-formulation:

****************
DEM Background
****************
In this chapter, we mathematically describe general features of explicit DEM simulations, with some reference to Yade implementation of these algorithms. They are given roughly in the order as they appear in simulation; first, two particles might establish a new interaction, which consists in

#. detecting collision between particles;
#. creating new interaction and determining its properties (such as stiffness); they are either precomputed or derived from properties of both particles;

Then, for already existing interactions, the following is performed:

#. strain evaluation;
#. stress computation based on strains;
#. force application to particles in interaction.

This simplified description serves only to give meaning to the ordering of sections within this chapter. A more detailed description of this *simulation loop* is given later.

Collision detection
=====================

Generalities
-------------

Exact computation of collision configuration between two particles can be relatively expensive (for instance between :yref:`Sphere` and :yref:`Facet`). Taking a general pair of bodies $i$ and $j$ and their ``exact'' (In the sense of precision admissible by numerical implementation.) spatial predicates (called :yref:`Shape` in Yade) represented by point sets $P_i$, $P_j$ the detection generally proceeds in 2 passes:

#. fast collision detection using approximate predicate $\tilde P_i$ and $\tilde P_j$; they are pre-constructed in such a way as to abstract away individual features of $P_i$ and $P_j$ and satisfy the condition

   .. math:: \forall {\bf x}\in R^3: x\in P_i\Rightarrow x\in \tilde P_i
     :label: eq-bounding-predicate

   (likewise for $P_j$). The approximate predicate is called ``bounding volume'' (:yref:`Bound` in Yade) since it bounds any particle's volume from outside (by virtue of the implication). It follows that $(P_i \cap  P_j)\neq\emptyset \Rightarrow (\tilde P_i \cap  \tilde P_j)\neq\emptyset$ and, by applying *modus tollens*,

   .. math:: \bigl(\tilde P_i \cap \tilde P_j\bigr)=\emptyset\Rightarrow\bigl( P_i \cap  P_j \bigr)=\emptyset
     :label: eq-collide-exclude

   which is a candidate exclusion rule in the proper sense.

#. By filtering away impossible collisions in :eq:`eq-collide-exclude`, a more expensive, exact collision detection algorithms can be run on possible interactions, filtering out remaining spurious couples $(\tilde P_i \cap  \tilde P_j)\neq\emptyset \wedge \bigl(P_i \cap  P_j\bigr)=\emptyset$. These algorithms operate on $P_i$ and $P_j$ and have to be able to handle all possible combinations of shape types.

It is only the first step we are concerned with here.

Algorithms
-----------

Collision evaluation algorithms have been the subject of extensive research in fields such as robotics, computer graphics and simulations. They can be roughly divided in two groups:

Hierarchical algorithms
   which recursively subdivide space and restrict the number of approximate checks in the first pass, knowing that lower-level bounding volumes can intersect only if they are part of the same higher-level bounding volume. Hierarchy elements are bounding volumes of different kinds: octrees [Jung1997]_, bounding spheres [Hubbard1996]_, k-DOP's [Klosowski1998]_.

Flat algorithms
   work directly with bounding volumes without grouping them in hierarchies first; let us only mention two kinds commonly used in particle simulations:

	Sweep and prune
	   algorithm operates on axis-aligned bounding boxes, which overlap if and only if they overlap along all axes. These algorithms have roughly $\bigO{n\log n}$ complexity, where $n$ is number of particles as long as they exploit `temporal coherence <sect-temp-coherence>`_ of the simulation.
	
	Grid algorithms
	   represent continuous $R^3$ space by a finite set of regularly spaced points, leading to very fast neighbor search; they can reach the $\bigO{n}$ complexity [Munjiza1998]_ and recent research  suggests ways to overcome one of the major drawbacks of this method, which is the necessity to adjust grid cell size to the largest particle in the simulation ([Munjiza2006]_, the ``multistep'' extension).

.. _sect-temp-coherence:
Temporal coherence
   expresses the fact that motion of particles in simulation is not arbitrary but governed by physical laws. This knowledge can be exploited to optimize performance.
				
Numerical stability of integrating motion equations dictates an upper limit on $\Delta t$ (sect. :ref:`sect-formulation-dt`) and, by consequence, on displacement of particles during one step. This consideration is taken into account in [Munjiza2006]_, implying that any particle may not move further than to a neighboring grid cell during one step allowing the $\bigO{n}$ complexity; it is also explored in the periodic variant of the sweep and prune algorithm described below.

On a finer level, it is common to enlarge $\tilde P_i$ predicates in such a way that they satisfy the :eq:`eq-bounding-predicate` condition during *several* timesteps; the first collision detection pass might then be run with stride, speeding up the simulation considerably. The original publication of this optimization by Verlet [Verlet1967]_ used enlarged list of neighbors, giving this technique the name *Verlet list*. In general cases, however, where neighbor lists are not necessarily used, the term *Verlet distance* is employed.

.. _sect-sweep-and-prune:

Sweep and prune
----------------
Let us describe in detail the sweep and prune algorithm used for collision detection in Yade (class :yref:`InsertionSortCollider`). Axis-aligned bounding boxes (:yref:`Aabb`) are used as $\tilde P_i$; each :yref:`Aabb` is given by lower and upper corner $\in R^3$ (in the following, $\tilde P_i^{x0}$, $\tilde P_i^{x1}$ are minimum/maximum coordinates of $\tilde P_i$ along the $x$-axis and so on). Construction of :yref:`Aabb` from various particle :yref:`Shape`'s (such as :yref:`Sphere`, :yref:`Facet`, :yref:`Wall`) is straightforward, handled by appropriate classes deriving form :yref:`BoundFunctor` (:yref:`Bo1_Sphere_Aabb`, :yref:`Bo1_Facet_Aabb`, …).
			
Presence of overlap of two :yref:`Aabb`'s can be determined from conjunction of separate overlaps of intervals along each axis (`fig-sweep-and-prune`_):

.. math:: \left(\tilde P_i \cap  \tilde P_j\right)\neq\emptyset \Leftrightarrow \bigwedge_{w\in\{x,y,z\}}\left[\left(\left(\tilde P_i^{w0},\tilde P_i^{w1}\right) \cap \left(\tilde P_j^{w0},\tilde P_j^{w1}\right)\right)\neq\emptyset\right]

where $(a,b)$ denotes interval in $R$.

.. _fig-sweep-and-prune:
.. figure:: fig/sweep-and-prune.*

   Sweep and prune algorithm (shown in 2D), where :yref:`Aabb` of each sphere is represented by minimum and maximum value along each axis. Spatial overlap of :yref:`Aabb`'s is present if they overlap along all axes. In this case, $\tilde P_1\cap\tilde P_2\neq\emptyset$ (but note that $P_1\cap P_2=\emptyset$) and $\tilde P_2 \cap\tilde P_3\neq\emptyset$.}

The collider keeps 3 separate lists (arrays) $L_w$ for each axis $w\in\{x,y,z\}$

.. math:: L_w=\bigcup_{i} \left\{\tilde P_i^{w0}, \tilde P_i^{w1} \right\}

where $i$ traverses all particles. $L_w$ arrays (sorted sets) contain respective coordinates of minimum and maximum corners for each :yref:`Aabb` (we call these coordinates *bound* in the following); besides bound, each of list elements further carries ``id`` referring to particle it belongs to, and a flag whether it is lower or upper bound.

In the initial step, all lists are sorted (using quicksort, average $\bigO{n\log n}$) and one axis is used to create initial interactions: the range between lower and upper bound for each body is traversed, while bounds in-between indicate potential :yref:`Aabb` overlaps which must be checked on the remaining axes as well.

At each successive step, lists are already pre-sorted. Inversions occur where a particle's coordinate has just crossed another particle's coordinate; this number is limited by numerical stability of simulation and its physical meaning (giving spatio-temporal coherence to the algorithm). The insertion sort algorithm swaps neighboring elements if they are inverted, and has complexity between \bigO{n} and \bigO{n^2}, for pre-sorted and unsorted lists respectively. For our purposes, we need only to handle inversions, which by nature of the sort algorithm are detected inside the sort loop. An inversion might signify: 

* overlap along the current axis, if an upper bound inverts (swaps) with a lower bound (i.e. that the upper bound with a higher coordinate was out of order in coming before the lower bound with a lower coordinate). Overlap along the other 2 axes is checked and if there is overlap along all axes, a new potential interaction is created.
* End of overlap along the current axis, if lower bound inverts (swaps) with an upper bound. If there is only potential interaction between the two particles in question, it is deleted.
* Nothing if both bounds are upper or both lower.

Aperiodic insertion sort
^^^^^^^^^^^^^^^^^^^^^^^^^^


Let us show the sort algorithm on a sample sequence of numbers:

.. math:: \sortlines{\sortSep& 3 &\sortInv& 7 &\sortInv& 2 &\sortInv& 4 &\sortSep}

Elements are traversed from left to right; each of them keeps inverting (swapping) with neighbors to the left, moving left itself, until any of the following conditions is satisfied:

.. FIXME table

=================== ========================================================
($\leq$)            the sorting order with the left neighbor is correct, or
($||$)              the element is at the beginning of the sequence.
=================== ========================================================

We start at the leftmost element (the current element is marked $\currelem{i}$)

.. math:: \sortlines{\sortSep& \currelem{3} &\sortInv& 7 &\sortInv& 2 &\sortInv& 4 &\sortSep.}

It obviously immediately satisfies ($||$), and we move to the next element:

.. math:: \sortlines{\sortSep& 3 &\sortInv& \isleq{7} &\sortInv& 2 &\sortInv& 4 &\sortSep.}

Condition ($\leq$) holds, therefore we move to the right. The $\currelem{2}$ is not in order (violating ($\leq$)) and two inversions take place; after that, ($||$) holds:

.. math::

	\sortlines{
		\sortSep& 3 &\sortInv& 7 &\sortInv& \isnleq{2} &\sortInv& 4 &\sortSep,  \\
		\sortSep& 3 &\sortInv& \isnleq{2} &\sortInv& 7 &\sortInv& 4 &\sortSep, \\
		\sortSep& \currelem{2} &\sortInv& 3 &\sortInv& 7 &\sortInv& 4 &\sortSep.
	}

The last element $\currelem{4}$ first violates ($\leq$), but satisfies it after one inversion

.. math::

	\sortlines{ 
		\sortSep& 2 &\sortInv& 3 &\sortInv& 7 &\sortInv& \isnleq{4} &\sortSep,  \\
		\sortSep& 2 &\sortInv& 3 &\sortInv& \isleq{4} &\sortInv& 7 &\sortSep.  
	}

All elements having been traversed, the sequence is now sorted.
				
It is obvious that if the initial sequence were sorted, elements only would have to be traversed without any inversion to handle (that happens in $\mathcal{O}(n)$ time).

For each inversion during the sort in simulation, the function that investigates change in :yref:`Aabb` overlap is invoked, creating or deleting interactions.
				
The periodic variant of the sort algorithm is described in :ref:`sect-periodic-insertion-sort`, along with other periodic-boundary related topics.

Optimization with Verlet distances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As noted above, [Verlet1967]_ explored the possibility of running the collision detection only sparsely by enlarging predicates $\tilde P_i$.
				
In Yade, this is achieved by enlarging :yref:`Aabb` of particles by fixed relative length in all dimensions $\Delta L$ (:yref:`InsertionSortCollider.sweepLength`). Suppose the collider run last time at step $m$ and the current step is $n$. :yref:`NewtonIntegrator` tracks maximum distance traversed by particles (via maximum velocity magnitudes $\curr{v_{\rm max}}=\max |\curr{\dot{u}_i}|$ in each step, with the initial cummulative distance $L_{\rm max}=0$,

.. math:: \curr{L_{\rm max}}&=\prev{L_{\rm max}}+\curr{v_{\rm max}}\curr{\Dt}
  :label: eq-verlet-l0

triggering the collider re-run as soon as

.. math:: \curr{L_{\rm max}}&>\Delta L.
  :label: eq-verlet-trigger

The disadvantage of this approach is that even one fast particle determines $\curr{v_{\rm max}}$.

A solution is to track maxima per particle groups. The possibility of tracking each particle separately (that is what `ESyS-Particle <http://www.launchpad.net/esys-particle>`_ does) seemed to us too fine-grained. Instead, we assign particles to $b_n$ (:yref:`InsertionSortCollider.nBins`) *velocity bins* based on their current velocity magnitude. The bins' limit values are geometrical with the coefficient $b_c>1$ (:yref:`InsertionSortCollider.binCoeff`), the maximum velocity being the current global velocity maximum $\curr{v_{\rm max}}$ (with some constraints on its change rate, to avoid large oscillations); for bin $i\in\{0,\dots,b_n\}$ and particle $j$:

.. math:: \curr{v_{\rm max}} b_c^{-(i+1)}\leq |\curr{\dot{u}_j}| < v_{\rm max} b_c^{-i}.

(note that in this case, superscripts of $b_c$ mean exponentiation). Equations :eq:`eq-verlet-l0`--:eq:`eq-verlet-trigger` are used for each bin separately; however, when :eq:`eq-verlet-trigger` is satisfied, full collider re-run is necessary and all bins' distances are reset.

Particles in high-speed oscillatory motion could be put into a slow bin if they happen to be at the point where their instantaneous speed is low, causing the necessity of early collider re-run. This is avoided by allowing particles to only go slower by one bin rather than several at once.

Results of using Verlet distance depend highly on the nature of simulation and choice of parameters :yref:`InsertionSortCollider.nBins` and :yref:`InsertionSortColldier.binCoeff`. The binning algorithm was specifically designed for simulating local fracture of larger concrete specimen; in that way, only particles in the fracturing zone, with greater velocities, had the :yref:`Aabb`'s enlarged, without affecting quasi-still particles outside of this zone. In such cases, up to 50% overall computation time savings were observed, collider being run every ≈100 steps in average.


Creating interaction between particles
================================================
Collision detection described above is only approximate. Exact collision detection depends on the geometry of individual particles and is handled separately. In Yade terminology, the :yref:`Collider` creates only *potential* interactions; potential interactions are evaluated exactly using specialized algorithms for collision of two spheres or other combinations. Exact collision detection must be run at every timestep since it is at every step that particles can change their mutual position (the collider is only run sometimes if the Verlet distance optimization is in use). Some exact collision detection algorithms are described in :ref:`sect-strain-evaluation`; in Yade, they are implemented in classes deriving from :yref:`IGeomFunctor` (prefixed with ``Ig2``).
		
Besides detection of geometrical overlap (which corresponds to :yref:`IGeom` in Yade), there are also non-geometrical properties of the interaction to be determined (:yref:`IPhys`). In Yade, they are computed for every new interaction by calling a functor deriving from :yref:`IPhysFunctor` (prefixed with ``Ip2``) which accepts the given combination of :yref:`Material` types of both particles.

Stiffnesses
-----------
Basic DEM interaction defines two stiffnesses: normal stiffness $K_N$ and shear (tangent) stiffness $K_T$. It is desirable that $K_N$ be related to fictitious Young's modulus of the particles' material, while $K_T$ is typically determined as a given fraction of computed $K_N$. The $K_T/K_N$ ratio determines macroscopic Poisson's ratio of the arrangement, which can be shown by dimensional analysis: elastic continuum has two parameters ($E$ and $\nu$) and basic DEM model also has 2 parameters with the same dimensions $K_N$ and $K_T/K_N$; macroscopic Poisson's ratio is therefore determined solely by $K_T/K_N$ and macroscopic Young's modulus is then proportional to $K_N$ and affected by $K_T/K_N$.

Naturally, such analysis is highly simplifying and does not account for particle radius distribution, packing configuration and other possible parameters such as the interaction radius introduced later.

Normal stiffness
^^^^^^^^^^^^^^^^
The algorithm commonly used in Yade computes normal interaction stiffness as stiffness of two springs in serial configuration with lengths equal to the sphere radii (`fig-spheres-contact-stiffness`_).
				
.. _fig-spheres-contact-stiffness:
.. figure:: fig/spheres-contact-stiffness.*
   
	Series of 2 springs representing normal stiffness of contact between 2 spheres.
				
Let us define distance $l=l_1+l_2$, where $l_i$ are distances between contact point and sphere centers, which are initially (roughly speaking) equal to sphere radii. Change of distance between the spehre centers $\Delta l$ is distributed onto deformations of both spheres $\Delta l=\Delta l_1+\Delta l_2$ proportionally to their compliances. Displacement change $\Delta l_i$ generates force $F_i=K_i \Delta l_i$, where $K_i$ assures proportionality and has physical meaning and dimension of stiffness; $K_i$ is related to the sphere material modulus $E_i$ and some length $\tilde l_i$ proportional to $r_i$.

.. math::
	:nowrap:

	\begin{align*}
		\Delta l&=\Delta l_1+\Delta l_2\\
		K_i&=E_i \tilde l_i\\
		K_N\Delta l&=F=F_1=F_2\\
		K_N\left(\Delta l_1+\Delta l_2\right)&=F\\
		K_N\left(\frac{F}{K_1}+\frac{F}{K_2}\right)&=F\\
		K_1^{-1}+K_2^{-1}&=K_N^{-1} \\
		K_N&=\frac{K_1 K_2}{K_1 + K_2} \\
		K_N&=\frac{E_1 \tilde l_1 E_2 \tilde l_2}{E_1\tilde l_1+E_2\tilde l_2}
	\end{align*}

The most used class computing interaction properties :yref:`Ip2_FrictMat_FrictMat_FrictPhys` uses $\tilde l_i=2r_i$.
				
Some formulations define an equivalent cross-section $A_{\rm eq}$, which in that case appears in the $\tilde l_i$ term as $K_i=E_i\tilde l_i=E_i\frac{A_{\rm eq}}{l_i}$. Such is the case for the concrete model (:yref:`Ip2_CpmMat_CpmMat_CpmPhys`), where $A_{\rm eq}=\min(r_1,r_2)$.
				
For reasons given above, no pretense about equality of particle-level $E_i$ and macroscopic modulus $E$ should be made. Some formulations, such as [Hentz2003]_, introduce parameters to match them numerically. This is not appropriate, in our opinion, since it binds those values to particular features of the sphere arrangement that was used for calibration.

Other parameters
----------------
Non-elastic parameters differ for various material models. Usually, though, they are averaged from the particles' material properties, if it makes sense. For instance, :yref:`Ip2_CpmMat_CpmMat_CpmPhys` averages most quantities, while :yref:`Ip2_FrictMat_FrictMat_FrictPhys` computes internal friction angle as $\phi=\min(\phi_1,\phi_2)$ to avoid friction with bodies that are frictionless.

.. _sect-strain-evaluation:

Strain evaluation
=================
In the general case, mutual configuration of two particles has 6 degrees of freedom (DoFs) just like a beam in 3D space: both particles have 6 DoFs each, but the interaction itself is free to move and rotate in space (with both spheres) having 6 DoFs itself; then $12-6=6$. They are shown at `fig-spheres-dofs`_.

.. _fig-spheres-dofs:
.. figure:: fig/spheres-dofs.*

	Degrees of freedom of configuration of two spheres. Normal strain appears if there is a difference of linear velocity along the interaction axis ($n$); shearing originates from the difference of linear velocities perpendicular to $n$ *and* from the part of $\vec{\omega}_1+\vec{\omega}_2$ perpendicular to $n$; twisting is caused by the part of $\vec{\omega}_1-\vec{\omega}_2$ parallel with $n$; bending comes from the part of $\vec{\omega}_1-\vec{\omega}_2$ perpendicular to $n$.

We will only describe normal and shear components of strain in the following, leaving torsion and bending aside. The reason is that most constitutive laws for contacts do not use the latter two.

Normal strain
-------------

.. _sect-normal-strain-constants:
Constants
^^^^^^^^^^

Let us consider two spheres with *initial* centers $\bar{\vec{C}_1}$, $\bar{\vec{C}}_2$ and radii $r_1$, $r_2$ that enter into contact. The order of spheres within the contact is arbitrary and has no influence on the behavior. Then we define lengths

.. math::
	:nowrap:

	\begin{align*}
		d_0&=|\bar{\vec{C}}_2-\bar{\vec{C}}_1| \\
		d_1&=r_1+\frac{d_0-r_1-r_2}{2}, & d_2&=d_0-d_1.
	\end{align*}

These quantities are *constant* throughout the life of the interaction and are computed only once when the interaction is established. The distance $d_0$ is the *reference distance* and is used for the conversion of absolute displacements to dimensionless strain, for instance. It is also the distance where (for usual contact laws) there is neither repulsive nor attractive force between the spheres, whence the name *equilibrium distance*.

.. _fig-sphere-sphere:
.. figure:: fig/sphere-sphere.*

	Geometry of the initial contact of 2 spheres; this case pictures spheres which already overlap when the contact is created (which can be the case at the beginning of a simulation) for the sake of generality. The initial contact point $\bar{\vec{C}}$ is in the middle of the overlap zone.
				
Distances $d_1$ and $d_2$ define reduced (or expanded) radii of spheres; geometrical radii $r_1$ and $r_2$ are used only for collision detection and may not be the same as $d_1$ and $d_2$, as shown in fig. `fig-sphere-sphere`_. This difference is exploited in cases where the average number of contacts between spheres should be increased, e.g. to influence the response in compression or to stabilize the packing. In such case, interactions will be created also for spheres that do not geometrically overlap based on the *interaction radius* $R_I$, a dimensionless parameter determining „non-locality“ of contact detection. For $R_I=1$, only spheres that touch are considered in contact; the general condition reads

.. math:: d_0&\leq R_I(r_1+r_2).
	:label: eq-strain-interaction-radius

The value of $R_I$ directly influences the average number of interactions per sphere (percolation), which for some models is necessary in order to achieve realistic results. In such cases, :yref:`Aabb` (or $\tilde P_i$ predicates in general) must be enlarged accordingly (:yref:`Bo1_Sphere_Aabb.aabbEnlargeFactor`).

Contact cross-section
"""""""""""""""""""""
Some constitutive laws are formulated with strains and stresses (:yref:`Law2_Dem3DofGeom_CpmPhys_Cpm`, the concrete model described later, for instance); in that case, equivalent cross-section of the contact must be introduced for the sake of dimensionality. The exact definition is rather arbitrary; the CPM model (:yref:`Ip2_CpmMat_CpmMat_CpmPhys`) uses the relation

.. math:: A_{\rm eq}=\pi\min(r_1,r_2)^2
	:label: eq-strain-crosssection

which will be used to convert stresses to forces, if the constitutive law used is formulated in terms of stresses and strains. Note that other values than $\pi$ can be used; it will merely scale macroscopic packing stiffness; it is only for the intuitive notion of a truss-like element between the particle centers that we choose $A_{\rm eq}$ representing the circle area. Besides that, another function than $\min(r_1,r_2)$ can be used, although the result should depend linearly on $r_1$ and $r_2$ so that the equation gives consistent results if the particle dimensions are scaled.

Variables
^^^^^^^^^^^
The following state variables are updated as spheres undergo motion during the simulation (as $\currC_1$ and $\currC_2$ change):

.. math::
	:label: eq-contact-normal

	\currn=\frac{\currC_2-\currC_1}{|\currC_2-\currC_1|}\equiv \normalized{\currC_2-\currC_1}

and

.. math::
	:label: eq-contact-point

	\currC=\currC_1+\left(d_1-\frac{d_0-|\currC_2-\currC_1|}{2}\right)\vec{n}.

The contact point $\currC$ is always in the middle of the spheres' overlap zone (even if the overlap is negative, when it is in the middle of the empty space between the spheres). The *contact plane* is always perpendicular to the contact plane normal $\currn$ and passes through $\currC$.

Normal displacement and strain can be defined as

.. math:: 
	:nowrap:
	
	\begin{align*}
		u_N&=|\currC_2-\currC_1|-d_0, \\
		\eps_N&=\frac{u_N}{d_0}=\frac{|\currC_2-\currC_1|}{d_0}-1.
	\end{align*}

Since $u_N$ is always aligned with $\vec{n}$, it can be stored as a scalar value multiplied by $\vec{n}$ if necessary.

For massively compressive simulations, it might be beneficial to use the logarithmic strain, such that the strain tends to $-\infty$ (rather than $-1$) as centers of both spheres approach. Otherwise, repulsive force would remain finite and the spheres could penetrate through each other. Therefore, we can adjust the definition of normal strain as follows:

.. math::

	\eps_N = \begin{cases}
		\log \left(\frac{|\currC_2-\currC_1|}{d_0}\right)& \hbox{if } |\currC_2-\currC_1|<d_0 \\
		\frac{|\currC_2-\currC_1|}{d_0}-1 & \hbox{otherwise.}
	\end{cases}

Such definition, however, has the disadvantage of effectively increasing rigidity (up to infinity) of contacts, requiring $\Dt$ to be adjusted, lest the simulation becomes unstable. Such dynamic adjustment is possible using a stiffness-based time-stepper (:yref:`GlobalStiffnessTimeStepper` in Yade).

Shear strain
-------------
In order to keep $\vec{u}_T$ consistent (e.g. that $\vec{u}_T$ must be constant if two spheres retain mutually constant configuration but move arbitrarily in space), then either $\vec{u}_T$ must track spheres' spatial motion or must (somehow) rely on sphere-local data exclusively.

These two possibilities lead to two algorithms of computing shear strains. They should give the same results (disregarding numerical imprecision), but there is a trade-off between computational cost of the incremental method and robustness of the total one.

Geometrical meaning of shear strain is shown in `fig-shear-2d`_.

.. _fig-shear-2d:
.. figure:: fig/shear-2d.*
	
	Evolution of shear displacement $\vec{u}_T$ due to mutual motion of spheres, both linear and rotational. Left configuration is the initial contact, right configuration is after displacement and rotation of one particle.

Incremental algorithm
^^^^^^^^^^^^^^^^^^^^^
The incremental algorithm is widely used in DEM codes and is described frequently ([Luding2008]_, [Alonso2004]_). Yade implements this algorithm in the :yref:`ScGeom` class. At each step, shear displacement $\uT$ is updated; the update increment can be decomposed in 2 parts: motion of the interaction (i.e. $\vec{C}$ and $\vec{n}$) in global space and mutual motion of spheres.

#. Contact moves dues to changes of the spheres' positions $\vec{C}_1$ and $\vec{C}_2$, which updates current $\currC$ and $\currn$ as per :eq:`eq-contact-point` and :eq:`eq-contact-normal`. $\prevuT$ is perpendicular to the contact plane at the previous step $\prevn$ and must be updated so that $\prevuT+(\Delta\uT)=\curruT\perp\currn$; this is done by perpendicular projection to the plane first (which might decrease $|\uT|$) and adding what corresponds to spatial rotation of the interaction instead:
  
   .. math::
      :nowrap:

      \begin{align*}
         (\Delta \uT)_1&=-\prevuT\times(\prevn \times \currn) \\
         (\Delta \uT)_2&=-\prevuT\times\left(\frac{\Delta t}{2} \currn \cdot (\pprev{\vec{\omega}}_1+\pprev{\vec{\omega}}_2)\right) \currn
      \end{align*}

#. Mutual movement of spheres, using only its part perpendicular to $\currn$; $\vec{v}_{12}$ denotes mutual velocity of spheres at the contact point:

   .. math::
      :nowrap:

      \begin{align*}
          \vec{v}_{12}&=\left(\pprev{\vec{v}}_2+\prev{\vec{\omega}}_2\times(-d_2 \currn)\right)-\left(\pprev{\vec{v}}_1+\pprev{\vec{\omega}}_1\times(d_1 \currn)\right) \\
         \vec{v}_{12}^{\perp}&=\vec{v}_{12}-(\curr{\vec{n}} \cdot \vec{v}_{12})\currn \\
         (\Delta \uT)_3&=-\Delta t \vec{v}_{12}^{\perp}
      \end{align*}


Finally, we compute

.. math:: \curruT=\prevuT+(\Delta\uT)_1 + (\Delta\uT)_2 + (\Delta\uT)_3.

.. _sect-formulation-total-shear:
Total alogithm
^^^^^^^^^^^^^^^
The following algorithm, aiming at stabilization of response even with large rotation speeds or $\Delta t$ approaching stability limit, was designed by the author of this thesis. (A similar algorithm based on total formulation, which covers additionally bending and torsion, was proposed in [Wang2009]_.) It is based on tracking original contact points (with zero shear) in the particle-local frame.

In this section, variable symbols implicitly denote their current values unless explicitly stated otherwise.

Shear strain may have two sources: mutual rotation of spheres or transversal displacement of one sphere with respect to the other. Shear strain does not change if both spheres move or rotate but are not in linear or angular motion mutually. To accurately and reliably model this situation, for every new contact the initial contact point $\bar{\vec{C}}$ is mapped into local sphere coordinates ($\vec{p}_{01}$, $\vec{p}_{02}$). As we want to determine the distance between both points (i.e. how long the trajectory in on both spheres' surfaces together), the shortest path from current $\vec{C}$ to the initial locally mapped point on the sphere's surface is „unrolled“ to the contact plane ($\vec{p}'_{01}$, $\vec{p}'_{02}$); then we can measure their linear distance $\uT$ and define shear strain $\vec{\eps}_T=\uT/d_0$ (fig. `fig-shear-displacement`_).

More formally, taking $\bar{\vec{C}}_i$, $\bar{q}_i$ for the sphere initial positions and orientations (as quaterions) in global coordinates, the initial sphere-local contact point *orientation* (relative to sphere-local axis $\hat{x}$) is remembered:

.. math::
	:nowrap:
	
	\begin{align*}
		\bar{\vec{n}}&=\normalized{{\vec{C}}_1-{\vec{C}}_2}, \\
		\bar{q}_{01}&=\Align(\hat x,\bar{q}_1^*\bar{\vec{n}} \bar{q}_1^{**}), \\ 
		\bar{q}_{02}&=\Align(\hat x,\bar{q}_2^* (-\bar{\vec{n}}) \bar{q}_2^{**}).
	\end{align*}
	
.. (See \autoref{sect-quaternions} for definition of $\Align$.)

				
After some spheres motion, the original point can be "unrolled" to the current contact plane:

.. math::
	:nowrap:
	
	\begin{align*}
		q&=\Align(\vec{n},q_1 \bar{q}_{01} \hat x (q_1 \bar{q}_{01})^*) \quad\hbox{(auxiliary)} \\
		\vec{p}'_{01}&=q_{\theta}d_1(q_{\vec{u}} \times \vec{n})
	\end{align*}

where $q_{\vec{u}}$, $q_{\theta}$ are axis and angle components of $q$ and $p_{01}'$ is the unrolled point. Similarly,

.. math::
	:nowrap:

	\begin{align*}
		q&=\Align(\vec{n},q_2 \bar{q}_{02} \hat x (q_2 \bar{q}_{02})^*) \\
		\vec{p}'_{02}&=q_{\theta}d_1(q_{\vec{u}} \times (-\vec{n})).
	\end{align*}

Shear displacement and strain are then computed easily:

.. math::
	:nowrap:

	\begin{align*}
		\uT&=\vec{p}'_{02}-\vec{p}'_{01} \\
		\vec{\eps}_T&=\frac{\uT}{d_0}
	\end{align*}

When using material law with plasticity in shear, it may be necessary to limit maximum shear strain, in which case the mapped points are moved closer together to the requested distance (without changing $\hat{\vec{u}}_T$). This allows us to remember the previous strain direction and also avoids summation of increments of plastic strain at every step (`fig-shear-slip`_).

.. _fig-shear-displacement:
.. figure:: fig/shear-displacement.*

	Shear displacement computation for two spheres in relative motion.


.. _fig-shear-slip:
.. figure:: fig/shear-slip.*

	Shear plastic slip for two spheres.

This algorithm is straightforwardly modified to facet-sphere interactions. In Yade, it is implemented by :yref:`Dem3DofGeom` and related classes.

.. _sect-formulation-stress-cundall:
Stress evaluation (example)
===========================
Once strain on a contact is computed, it can be used to compute stresses/forces acting on both spheres.

The constitutive law presented here is the most usual DEM formulation, originally proposed by Cundall. While the strain evaluation will be similar to algorithms described in the previous section regardless of stress evaluation, stress evaluation itself depends on the nature of the material being modeled. The constitutive law presented here is the most simple non-cohesive elastic case with dry friction, which Yade implements in :yref:`Law2_Dem3DofGeom_FrictPhys_Basic` (all constitutive laws derive from base class :yref:`LawFunctor`).
		
In DEM generally, some constitutive laws are expressed using strains and stresses while others prefer displacement/force formulation. The law described here falls in the latter category.

When new contact is established (discussed in :ref:`sect-simulation-loop`) it has its properties (:yref:`IPhys`) computed from :yref:`Materials<Material>` associated with both particles. In the simple case of frictional material :yref:`FrictMat`, :yref:`Ip2_FrictMat_FrictMat_FrictPhys` creates a new :yref:`FrictPhys` instance, which defines normal stiffness $K_N$, shear stiffness $K_T$ and friction angle $\phi$.

At each step, given normal and shear displacements $u_N$, $\uT$, normal and shear forces are computed (if $u_N>0$, the contact is deleted without generating any forces):

.. math::
	:nowrap:

	\begin{align*}
		\vec{F}_N&=K_N u_N \vec{n},\\
		\vec{F}_T^t&=K_T \vec{u}_T
	\end{align*}

where $\vec{F}_N$ is normal force and $\vec{F}_T$ is trial shear force.  A simple non-associated stress return algorithm is applied to compute final shear force

.. math::
	:nowrap:
	
	\begin{align*}
		\vec{F}_T=\begin{cases}
			\vec{F}_T^t\frac{|\vec{F}_N|\tan\phi}{\vec{F}_T^t} & \hbox{if } |\vec{F}_T|>|\vec{F}_N|\tan\phi, \\
			\vec{F}_T^t & \hbox{otherwise.}
		\end{cases}
	\end{align*}
S
ummary force $\vec{F}=\vec{F}_N+\vec{F}_T$ is then applied to both particles -- each particle accumulates forces and torques acting on it in the course of each step. Because the force computed acts at contact point $\vec{C}$, which is difference from sphes' centers, torque generated by $\vec{F}$ must also be considered.

.. math:: 
	:nowrap:

	\begin{align*}
		\vec{F}_1&+=\vec{F} & \vec{F}_2&+=-\vec{F} \\
		\vec{T}_1&+=d_1(-\vec{n})\times \vec{F} & \vec{T}_2&+=d_2 \vec{n} \times \vec{F}.
	\end{align*}

Motion integration
==================
Each particle accumulates generalized forces (forces and torques) from the contacts in which it participates. These generalized forces are then used to integrate motion equations for each particle separately; therefore, we omit $i$ indices denoting the $i$-th particle in this section.

The customary leapfrog scheme (also known as the Verlet scheme) is used, with some adjustments for rotation of non-spherical particles, as explained below. The "leapfrog" name comes from the fact that even derivatives of position/orientation are known at on-step points, whereas odd derivatives are known at mid-step points. Let us recall that we use $\prev{a}$, $\curr{a}$, $\next{a}$ for on-step values of $a$ at $t-\Dt$, $t$ and $t+\Dt$ respectively; and $\pprev{a}$, $\nnext{a}$ for mid-step values of $a$ at $t-\Dt/2$, $t+\Dt/2$.
		
Described integration algorithms are implemented in the :yref:`NewtonIntegrator` class in Yade.

Position
----------
Integrating motion consists in using current acceleration $\curraccel$ on a particle to update its position from the current value $\currpos$ to its value at the next timestep $\nextpos$. Computation of acceleration, knowing current forces $\vec{F}$ acting on the particle in question and its mass $m$, is simply

.. math:: \curraccel&=\vec{F}/m.
			
Using the 2nd order finite difference with step $\Dt$, we obtain

.. math:: \curraccel&\cong\frac{\prevpos-2\currpos+\nextpos}{\Dt^2}

from which we express

.. math::

	\nextpos&=2\currpos-\prevpos+\curraccel\Dt^2 =\\
		&=\currpos+\Dt\underbrace{\left(\frac{\currpos-\prevpos}{\Dt}+\curraccel\Dt\right)}_{(\dagger)}.

Typically, $\prevpos$ is already not known (only $\currpos$ is); we notice, however, that

.. math:: \pprevvel&\simeq\frac{\currpos-\prevpos}{\Dt},

i.e. the mean velocity during the previous step, which is known. Plugging this approximate into the $(\dagger)$ term, we also notice that mean velocity during the current step can be approximated as

.. math:: \nnextvel&\simeq\pprevvel+\curraccel\Dt,

which is $(\dagger)$; we arrive finally at

.. math:: \nextpos&=\currpos+\Dt\left(\pprevvel+\curraccel\Dt\right).

The algorithm can then be written down by first computing current mean velocity $\nnextvel$ which we need to store for the next step (just as we use its old value $\pprevvel$ now), then computing the position for the next time step $\nextpos$:

.. math::
	:label: eq-leapfrog-nextvel
	:nowrap:

	\begin{align*}
		\nnextvel&=\pprevvel+\curraccel\Dt \\
		\nextpos&=\currpos+\nnextvel\Dt.
	\end{align*}

Positions are known at times $i\Delta t$ (if $\Delta t$ is constant) while velocities are known at $i\Delta t+\frac{\Delta t}{2}$. The facet that they interleave (jump over each other) in such way gave rise to the colloquial name "leapfrog" scheme.

Orientation (spherical)
------------------------
Updating particle orientation $\curr{q}$ proceeds in an analogous way to position update. First, we compute current angular acceleration $\curraaccel$ from known current torque $\vec{T}$. For spherical particles where the inertia tensor is diagonal in any orientation (therefore also in current global orientation), satisfying $\vec{I}_{11}=\vec{I}_{22}=\vec{I}_{33}$, we can write

.. math:: \curraaccel_i&=\vec{T}_i/\vec{I}_{11},

We use the same approximation scheme, obtaining an equation analogous to :eq:`eq-leapfrog-nnextvel`

.. math:: \nnextangvel&=\pprevangvel+\Dt\curraaccel.

The quaternion $\Delta q$ representing rotation vector $\nnextangvel\Dt$ is constructed, i.e. such that

.. math::
	:nowrap:

	\begin{align*}
		(\Delta q)_{\theta}&=|\nnextangvel|, \\
		(\Delta q)_{\vec{u}}&=\normalized{\nnextangvel}
	\end{align*}

Finally, we compute the next orientation $\next{q}$ by rotation composition

.. math:: \next{q}&=\Delta q\curr{q}.

Orientation (aspherical)
------------------------
Integrating rotation of aspherical particles is considerably more complicated than their position, as their local reference frame is not inertial. Rotation of rigid body in the local frame, where inertia matrix $\mat{I}$ is diagonal, is described in the continuous form by Euler's equations ($i\in\{1,2,3\}$ and $i$, $j$, $k$ are subsequent indices):

.. math:: \vec{T}_i=\mat{I}_{ii}\dot{\vec{\omega}}_i+(\mat{I}_{kk}-\mat{I}_{jj})\vec{\omega}_j\vec{\omega}_k.

Due to the presence of the current values of both $\vec{\omega}$ and $\dot{\vec{\omega}}$, they cannot be solved using the standard leapfrog algorithm (that was the case for translational motion and also for the spherical bodies' rotation where this equation reduced to $\vec{T}=\mat{I}\dot{\vec{\omega}}$).
			
The algorithm presented here is described by [Allen1989]_ (pg. 84--89) and was designed by Fincham for molecular dynamics problems; it is based on extending the leapfrog algorithm by mid-step/on-step estimators of quantities known at on-step/mid-step points in the basic formulation. Although it has received criticism and more precise algorithms are known ([Omelyan1999]_, [Neto2006]_, [Johnson2008]_), this one is currently implemented in Yade for its relative simplicity.

.. Finchman: Leapfrog Rotational Algorithms: http://www.informaworld.com/smpp/content~content=a756872469&db=all
	Schvanberg: Leapfrog Rotational Algorithms: http://www.informaworld.com/smpp/content~content=a914299295&db=all

			
Each body has its local coordinate system based on the principal axes of inertia for that body. We use $\locframe{\bullet}$ to denote vectors in local coordinates. The orientation of the local system is given by the current particle's orientation $\curr{q}$ as a quaternion; this quaternion can be expressed as the (current) rotation matrix $\mat{A}$. Therefore, every vector $\vec{a}$ is transformed as $\locframe{\vec{a}}=q\vec{a}q^{*}=\mat{A}\vec{a}$. Since $\mat{A}$ is a rotation (orthogonal) matrix, the inverse rotation $\mat{A}^{-1}=\mat{A}^{T}$.

For given particle in question, we know

* $\loccurr{\mat{I}}$ (constant) inertia matrix; diagonal, since in local, principal coordinates,
* $\curr{\vec{T}}$ external torque,
* $\curr{q}$ current orientation (and its equivalent rotation matrix $\mat{A}$),
* $\pprev{\vec{\omega}}$  mid-step angular velocity,
* $\pprev{\vec{L}}$ mid-step angular momentum; this is an auxiliary variable that must be tracked in addition for use in this algorithm. It will be zero in the initial step.

Our goal is to compute new values of the latter three, that is $\nnext{\vec{L}}$, $\next{q}$, $\nnext{\vec{\omega}}$. We first estimate current angular momentum and compute current local angular velocity:

.. math::
	:nowrap:

	\begin{align*}
		\curr{\vec{L}}&=\pprev{\vec{L}}+\curr{\vec{T}}\frac{\Dt}{2}, &\loccurr{\vec{L}}&=\mat{A}\curr{\vec{L}}, \\
		\nnext{\vec{L}}&=\pprev{\vec{L}}+\curr{\vec{T}}\Dt, &\nnext{\locframe{\vec{L}}}&=\mat{A}\nnext{\vec{L}}, \\
		\loccurr{\vec{\omega}}&=\curr{\locframe{\mat{I}}}{}^{-1}\loccurr{\vec{L}}, \\
		\nnext{\locframe{\vec{\omega}}}&=\curr{\locframe{\mat{I}}}{}^{-1}\nnext{\locframe{\vec{L}}}. \\
	\end{align*}

Then we compute $\curr{\dot{q}}$, using $\curr{q}$ and $\loccurr{\vec{\omega}}$:

.. math::
	:label: eq-quaternion-derivative
	:nowrap:

		\begin{align*}
			\begin{pmatrix}\curr{\dot{q}}_w \\ \curr{\dot{q}}_x \\ \curr{\dot{q}}_y \\ \curr{\dot{q}}_z\end{pmatrix}&=
				\def\cq{\curr{q}}
				\frac{1}{2}\begin{pmatrix}
					\cq_w & -\cq_x & -\cq_y & -\cq_z \\
					\cq_x & \cq_w & -\cq_z & \cq_y \\
					\cq_y & \cq_z & \cq_w & -\cq_x \\
					\cq_z & -\cq_y & \cq_x & \cq_w
				\end{pmatrix}
				\begin{pmatrix} 0 \\ \loccurr{\vec{\omega}}_x \\ \loccurr{\vec{\omega}}_y \\ \loccurr{\vec{\omega}}_z	\end{pmatrix},  \\
				\nnext{q}&=\curr{q}+\curr{\dot{q}}\frac{\Dt}{2}.\\
		\end{align*}

We evaluate $\nnext{\dot{q}}$ from $\nnext{q}$ and $\nnext{\locframe{\vec{\omega}}}$ in the same way as in :eq:`eq-quaternion-derivative` but shifted by $\Dt/2$ ahead. Then we can finally compute the desired values

.. math::
	:nowrap:

	\begin{align*}
		\next{q}&=\curr{q}+\nnext{\dot{q}}\Dt, \\
		\nnext{\vec{\omega}}&=\mat{A}^{-1}\nnext{\locframe{\vec{\omega}}}
	\end{align*}

Clumps (rigid aggregates)
-------------------------
DEM simulations frequently make use of rigid aggregates of particles to model complex shapes [Price2007]_ called *clumps*, typically composed of many spheres. Dynamic properties of clumps are computed from the properties of its members: the clump's mass $m_c$ is summed over members, the inertia tensor $\mathbf{I}_c$ with respect to the clump's centroid is computed using the parallel axes theorem; local axes are oriented such that they are principal and inertia tensor is diagonal and clump's orientation is changed to compensate rotation of the local system, as to not change the clump members' positions in global space. Initial positions and orientations of all clump members in local coordinate system are stored.

In Yade (class :yref:`Clump`), clump members behave as stand-alone particles during simulation for purposes of collision detection and contact resolution, except that they have no contacts created among themselves within one clump. It is at the stage of motion integration that they are treated specially. Instead of integrating each of them separately, forces/torques on those particles $\vec{F}_i$, $\vec{T}_i$ are converted to forces/torques on the clump itself. Let us denote $r_i$ relative position of each particle with regards to clump's centroid, in global orientation. Then summary force and torque on the clump are

.. math::
	:nowrap:

	\begin{align*}
		\vec{F}_c&=\sum F_i/m_c, \\
		\vec{T}_c&=\sum r_i\times F_i + T_i.
	\end{align*}

Motion of the clump is then integrated, using aspherical rotation integration. Afterwards, clump members are displaced in global space, to keep their initial positions and orientations in the clump's local coordinate system. In such a way, relative positions of clump members are always the same, resulting in the behavior of a rigid aggregate.

Numerical damping
-----------------
In simulations of quasi-static phenomena, it it desirable to dissipate kinetic energy of particles. Since most constitutive laws (including :yref:`Law_ScGeom_FrictPhys_Basic` shown above, :ref:`sect-formulation-stress-cundall`) do not include velocity-based damping (such as one in [Addetta2001]_), it is possible to use artificial numerical damping. The formulation is described in [Pfc3dManual30]_, although our version is slightly adapted. The basic idea is to decrease forces which increase the particle velocities and vice versa by $(\Delta F)_d$, comparing the current acceleration sense and particle velocity sense. This is done by component, which makes the damping scheme clearly non-physical, as it is not invariant with respect to coordinate system rotation; on the other hand, it is very easy to compute. Cundall proposed the form (we omit particle indices $i$ since it applies to all of them separately):

.. math:: \frac{(\Delta \vec{F})_{dw}}{\vec{F}_w}=-\lambda_d\sign(\vec{F}_w\pprev{\dot{\vec{u}}}_{w}),\quad w\in\{x,y,z\}

where $\lambda_d$ is the damping coefficient. This formulation has several advantages [Hentz2003]_:

* it acts on forces (accelerations), not constraining uniform motion;
* it is independent of eigenfrequencies of particles, they will be all damped equally;
* it needs only the dimensionless parameter $\lambda_d$ which does not have to be scaled.

In Yade, we use the adapted form

.. math::
	:label: eq-damping-yade

	\frac{(\Delta\vec{F})_{dw}}{\vec{F}_w}=-\lambda_d\sign\vec{F}_w\underbrace{\left(\pprev{\dot{u}}_w+\frac{\curr{\ddot{\vec{u}}}_w\Dt}{2}\right)}_{\simeq\curr{\dot{u}}_w},

where we replaced the previous mid-step velocity $\pprev{\dot{u}}$ by its on-step estimate in parentheses. This is to avoid locked-in forces that appear if the velocity changes its sign due to force application at each step, i.e. when the particle in question oscillates around the position of equilibrium with $2\Dt$ period.

In Yade, damping :eq:`eq-damping-yade` is implemented in the :yref:`NewtonIntegrator` engine; the damping coefficient $\lambda_d$ is :yref:`NewtonIntegrator.damping`.

.. _sect-formulation-dt:
Stability considerations
------------------------------------------------------
..
	\def\Dtcr{\Dt_{\rm cr}}
	% http://en.wikipedia.org/wiki/Harmonic_oscillator
	% eigenfrequency is (1/2π)√(k/m) for simple harmonic oscillator
	% Bruno gives √(m/k) 
	% http://imechanica.org/node/7670#comment-13672: Δt_crit=2/ω_max
	The leapfrog integration scheme is conditionally stable, i.e. not magnifying errors, provided $\Dt<\Dtcr$ where $\Dtcr$ is the *critical timestep*, above which the integration is unstable. Usually, $\Dt$ is taken as a fraction of $\Dtcr$; this fraction is called the *timestep safety factor*, with meaningful values $\in\langle 0,1)$.

Critical timestep
^^^^^^^^^^^^^^^^^
In order to ensure stability for the explicit integration sceheme, an upper limit is imposed on $\Dt$:

.. math::
	:label: eq-dt-angular

	\Dtcr=\frac{2}{\omega_{\rm max}}


where $\omega_{\rm max}$ is the highest eigenfrequency within the system.

Single mass-spring system
""""""""""""""""""""""""""

Single 1D mass-spring system with mass $m$ and stiffness $K$ is governed by the equation

.. math:: m\ddot{x}=-Kx

where $x$ is displacement from the mean (equilibrium) position. The solution of harmonic oscillation is $x(t)=A\cos(\omega t+\phi)$ where phase $\phi$ and amplitude $A$ are determined by initial conditions. The angular frequency

.. math::
	:label: eq-dt-omega
	
	\omega^{(1)}=\sqrt{\frac{K}{m}}

does not depend on initial conditions. Since there is one single mass, $\omega_{\rm max}^{(1)}=\omega^{(1)}$. Plugging :eq:`eq-dt-omega` into :eq:`eq-dt-angular`, we obtain

.. math:: \Dtcr^{(1)}=2/\omega_{\rm max}^{(1)}=2\sqrt{m/K}

for a single oscillator.
	
General mass-spring system
""""""""""""""""""""""""""
				
In a general mass-spring system, the highest frequency occurs if two connected masses $m_i$, $m_j$ are in opposite motion; let us suppose they have equal velocities (which is conservative) and they are connected by a spring with stiffness $K_{i}$: displacement $\Delta x_i$ of $m_i$ will be accompained by $\Delta x_j=-\Delta x_i$ of $m_j$, giving $\Delta F_i=-K_{i}(\Delta x_i-(-\Delta x_i))=-2K_{i}\Delta x_i$. That results in apparent stiffness $K_{i}^{(2)}=2K_{i}$, giving maximum angular frequency of the whole system

.. math:: \omega_{\rm max}=\max_i\sqrt{K_i^{(2)}/m_i}.
			
The overall critical timestep is then

.. math::
	:label: eq-dtcr-global

	\Dtcr=\frac{2}{\omega_{\rm max}}=\min_i\, 2\sqrt{\frac{m_i}{K_i^{(2)}}}=\min_i\, 2\sqrt{\frac{m_i}{2K_i}}=\min_i \sqrt{2}\sqrt{\frac{m_i}{K_i}}.

This equation can be used for all 6 degrees of freedom (DOF) in translation and rotation, by considering generalized mass and stiffness matrices $M$ and $K$, and replacing fractions $\frac{m_i}{K_i}$ by eigen values of $K.M^{-1}$. The critical timestep is then associated to the eigen mode with highest frequency : 

.. math::
	:label: eq-dtcr-axes

	\Dtcr=\min {\Dtcr}_k,\quad k\in\{1,...,6\}.


				
DEM simulations
""""""""""""""""
In DEM simulations, per-particle stiffness $\vec{K}_{ij}$ is determined from the stiffnesses of contacts in which it participates [Chareyre2005]_. Suppose each contact has normal stiffness $K_{Nk}$, shear stiffness $K_{Tk}=\xi K_{Nk}$ and is oriented by normal $\vec{n}_{k}$. A translational stiffness matrix $\vec{K}_{ij}$ can be defined as the sum of contributions of all contacts in which it participates (indices $k$), as

.. math::
	:label: eq-dtcr-particle-stiffness
	
	\vec{K}_{ij}=\sum_k (K_{Nk}-K_{Tk})\vec{n}_{i}\vec{n}_{j}+K_{Tk}=\sum_j K_{Nk}\left((1-\xi)\vec{n}_{i}\vec{n}_{j}+\xi\right)

with $i$ and $j\in\{x,y,z\}$. Equations :eq:`eq-dtcr-axes` and :eq:`eq-dtcr-particle-stiffness` determine $\Dtcr$ in a simulation. A similar approach generalized to all 6 DOFs is implemented by the :yref:`GlobalStiffnessTimeStepper` engine in Yade. The derivation of generalized stiffness including rotational terms is very similar but not developed here, for simplicity. For full reference, see "PFC3D - Theoretical Background".
					
Note that for computation efficiency reasons, eigenvalues of the stiffness matrices are not computed. They are only approximated assuming than DOF's are uncoupled, and using diagonal terms of $K.M^{-1}$. They give good approximates in typical mechanical systems.

There is one important condition that $\omega_{\rm max}>0$: if there are no contacts between particles and $\omega_{\rm max}=0$, we would obtain value $\Dtcr=\infty$. While formally correct, this value is numerically erroneous: we were silently supposing that stiffness remains constant during each timestep, which is not true if contacts are created as particles collide. In case of no contact, therefore, stiffness must be pre-estimated based on future interactions, as shown in the next section.
				

.. _sect-dt-pwave:
Estimation of $\Dtcr$ by wave propagation speed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		
Estimating timestep in absence of interactions is based on the connection between interaction stiffnesses and the particle's properties. Note that in this section, symbols $E$ and $\rho$ refer exceptionally to Young's modulus and density of *particles*, not of macroscopic arrangement.
			 
In Yade, particles have associated :yref:`Material` which defines density $\rho$ (:yref:`Material.density`), and also may define (in :yref:`ElastMat` and derived classes) particle's "Young's modulus" $E$ (:yref:`ElastMat.young`). $\rho$ is used when particle's mass $m$ is initially computed from its $\rho$, while $E$ is taken in account when creating new interaction between particles, affecting stiffness $K_N$. Knowing $m$ and $K_N$, we can estimate :eq:`eq-dtcr-particle-stiffness` for each particle; we obviously neglect 

* number of interactions per particle $N_i$; for a "reasonable" radius distribution, however, there is a geometrically imposed upper limit (6 for a packing of spheres with equal radii, for instance);
* the exact relationship the between particles' rigidities $E_i$, $E_j$, supposing only that $K_N$ is somehow proportional to them.

By defining $E$ and $\rho$, particles have continuum-like quantities. Explicit integration schemes for continuum equations impose a critical timestep based on sonic speed $\sqrt{E/\rho}$; the elastic wave must not propagate farther than the minimum distance of integration points $l_{\rm min}$ during one step. Since $E$, $\rho$ are parameters of the elastic continuum and $l_{\rm min}$ is fixed beforehand, we obtain

.. math:: \Dtcr^{(c)}=l_{\rm min}\sqrt{\frac{\rho}{E}}.

For our purposes, we define $E$ and $\rho$ for each particle separately; $l_{\rm min}$ can be replaced by the sphere's radius $R_i$; technically, $l_{\rm min}=2R_i$ could be used, but because of possible interactions of spheres and facets (which have zero thickness), we consider $l_{\rm min}=R_i$ instead. Then

.. FIXME Why not by sphere's diameter instead? We could say it is because :yref:`Sphere`-:yref:`Facet` interaction which has half length?

.. math:: \Dtcr^{(p)}=\min_i R_i \sqrt{\frac{\rho_i}{E_i}}.

This algorithm is implemented in the :yref:`yade.utils.PWaveTimeStep` function.
	
Let us compare this result to :eq:`eq-dtcr-global`; this necessitates making several simplifying hypotheses:

* all particles are spherical and have the same radius $R$;
* the sphere's material has the same $E$ and $\rho$
* the average number of contacts per sphere is $N$;
* the contacts have sufficiently uniform spatial distribution around each particle;
* the $\xi=K_N/K_T$ ratio is constant for all interactions;
* contact stiffness $K_N$ is computed from $E$ using a formula of the form

  .. math:: 
     :label: eq-dt-kn

     K_N=E\pi'R',

  where $\pi'$ is some constant depending on the algorithm in use\footnote{For example, $\pi'=\pi/2$ in the concrete particle model (:yref:`Ip2_CpmMat_CpmMat_CpmPhys`), while $\pi'=2$ in the classical DEM model (:yref:`Ip2_FrictMat_FrictMat_FrictPhys`) as implemented in Yade.} and $R'$ is half-distance between spheres in contact, equal to $R$ for the case of interaction radius $R_I=1$. If $R_I=1$ (and $R'\equiv R$ by consequence), all interactions will have the same stiffness $K_N$. In other cases, we will consider $K_N$ as the average stiffness computed from average $R'$ (see below).

As all particles have the same parameters, we drop the $i$ index in the following formulas.
	
We try to express the average per-particle stiffness from :eq:`eq-dtcr-particle-stiffness`. It is a sum over all interactions where $K_{N}$ and $\xi$ are scalars that will not rotate with interaction, while $\vec{n}_w$ is $w$-th component of unit interaction normal $\vec{n}$. Since we supposed uniform spatial distribution, we can replace $\vec{n}_w^2$ by its average value $\overline{\vec{n}}_w^2$. Recognizing components of $\vec{n}$ as direction cosines, the average values of $\vec{n}_w^2$ is $1/3$. %we find the average value by integrating over all possible orientations, which are uniformly distributed in space:
			
Moreover, since all directions are equal, we can write the per-body stiffness as $K=\vec{K}_w$ for all $w\in\{x,y,z\}$. We obtain 

.. math:: K=\sum K_N\left((1-\xi)\frac{1}{3}+\xi\right)=\sum K_N\frac{1-2\xi}{3}

and can put constant terms (everything) in front of the summation. $\sum 1$ equals the number of contacts per sphere, i.e. $N$. Arriving at

.. math:: K=N K_N \frac{1-2\xi}{3},

we substitute $K$ into :eq:`eq-dtcr-global` using :eq:`eq-dt-kn`:

.. math:: \Dtcr=\sqrt{2}\sqrt{\frac{m}{K}}=\sqrt{2}\sqrt{\frac{\frac{4}{3}\pi R^3\rho}{N E\pi'R\frac{1-2\xi}{3}}}=\underbrace{R\sqrt{\frac{\rho}{E}}}_{\Dtcr^{(p)}}2\sqrt{\frac{\pi/\pi'}{N(1-2\xi)}}.

The ratio of timestep $\Dtcr^{(p)}$ predicted by the p-wave velocity and numerically stable timestep $\Dtcr$ is the inverse value of the last (dimensionless) term:

.. math:: \frac{\Dtcr^{(p)}}{\Dtcr}=2\sqrt{\frac{N(1+\xi)}{\pi/\pi'}}.

Actual values of this ratio depend on characteristics of packing $N$, $K_N/K_T=\xi$ ratio and the way of computing contact stiffness from particle rigidity. Let us show it for two models in Yade:

Concrete particle model
   computes contact stiffness from the equivalent area $A_{\rm eq}$ first :eq:`eq-strain-crosssection`,

   .. math::
      :nowrap:
   
      \begin{align*}
         A_{\rm eq}&=\pi R^2 
         K_N&=\frac{A_{\rm eq}E}{d_0}.
      \end{align*}

   $d_0$ is the initial contact length, which will be, for interaction radius :eq:`eq-strain-interaction-radius` $R_I>1$, in average larger than $2R$. For $R_I=1.5$ (sect.~\ref{sect-calibration-elastic-properties}), we can roughly estimate $\overline{d}_0=1.25\cdot2R=\frac{5}{2}R$, getting

   .. math:: K_N=E\left(\frac{2}{5}\pi\right)R

   where $\frac{2}{5}\pi=\pi'$ by comparison with :eq:`eq-dt-kn`.

   Interaction radius $R_I=1.5$ leads to average $N\approx12$ interactions per sphere for dense packing of spheres with the same radius $R$. $\xi=0.2$ is calibrated (sect.~\ref{sect-calibration-elastic-properties}) to match the desired macroscopic Poisson's ratio $\nu=0.2$.

   Finally, we obtain the ratio

   .. math:: \frac{\Dtcr^{(p)}}{\Dtcr}=2\sqrt{\frac{12(1-2\cdot0.2)}{\frac{\pi}{(2/5)\pi}}}=3.39,

   showing significant overestimation by the p-wave algorithm.
                  
Non-cohesive dry friction model
   is the basic model proposed by Cundall explained in \ref{sect-formulation-stress-cundall}. Supposing almost-constant sphere radius $R$ and rather dense packing, each sphere will have $N=6$ interactions on average (that corresponds to maximally dense packing of spheres with a constant radius). If we use the :yref:`Ip2_FrictMat_FrictMat_FrictPhys` class, we have $\pi'=2$, as $K_N=E2R$; we again use $\xi=0.2$ (for lack of a more significant value). In this case, we obtain the result

   .. math:: \frac{\Dtcr^{(p)}}{\Dtcr}=2\sqrt{\frac{6(1-2\cdot0.2)}{\pi/2}}=3.02

   which again overestimates the numerical critical timestep.

To conclude, p-wave timestep gives estimate proportional to the real $\Dtcr$, but in the cases shown, the value of about $\Dt=0.3\Dtcr^{(p)}$ should be used to guarantee stable simulation.

Non-elastic $\Dt$ constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let us note at this place that not only $\Dtcr$ assuring numerical stability of motion integration is a constraint. In systems where particles move at relatively high velocities, position change during one timestep can lead to non-elastic irreversible effects such as damage. The $\Dt$ needed for reasonable result can be lower $\Dtcr$. We have no rigorously derived rules for such cases.

.. 	http://books.google.cz/books?id=_KTsjCZtt_EC&lpg=PP1&ots=nK0B4wui2F&dq=stability%20analysis&pg=PR5#v=onepage&q&f=false

..			The leap-frog integration scheme assumes constant $\Dt$. Although this is typically the case, and even dynamic time-stepping techniques such as :yref:`GlobalStiffnessTimeStepper` change $\Dt$ infrequently, big changes in timestep could destabilize the integration [Skeel1993]_.

Periodic boundary conditions
============================
While most DEM simulations happen in $R^3$ space, it is frequently useful to avoid boundary effects by using periodic space instead. In order to satisfy periodicity conditions, periodic space is created by repetition of parallelepiped-shaped cell. In Yade, periodic space is implemented in the :yref:`Cell` class. The cell is determined by

* \item the reference size $\vec{s}$ (:yref:`Cell.refSize`), giving reference cell configuration (which is always perpendicular): axis-aligned cuboid with corners $(0,0,0)$ and $\vec{s}$;
* \item the transformation matrix $\mat{T}$ (:yref:`Cell.trsf`).

The transformation matrix $\mat{T}$ can hold arbitrary linear transformation composed of scaling, rotation and shear. Volume change of the cell is given by $\det\mat{T}$. The cell can be manipulated by directly changing its transformation matrix $\mat{T}$ and its reference size $\vec{s}$.

Additionally, we define teransformation gradient $\nabla \vec{v}$ (:yref:`Cell.velGrad`) which can be automatically integrated at every step using the Euler scheme

.. math:: \next{\mat{T}}=\curr{\mat{T}}+\nabla \vec{v} \Dt.

Along with the automatic integration of cell transformation, there is an option to homothetically displace all particles so that $\nabla \vec{v}$ is swept linearly over the whole simulation (enabled via :yref:`NewtonIntegrator.homotheticCellResize`). This avoids all boundary effects coming from change of the transformation.

Collision detection in periodic cell
------------------------------------
In usual implementations, particle positions are forced to be inside the cell by wrapping their positions if they get over the boundary (so that they appear on the other side). As we wanted to avoid abrupt changes of position (it would make particle's velocity inconsistent with step displacement change), a different method was chosen.


.. _sect-cell-approx-collision:

Approximate collision detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pass 1 collision detection (based on sweep and prune algorithm, sect.~\ref{sect-sweep-and-prune}) operates on axis-aligned bounding boxes (:yref:`Aabb`) of particles. During the collision detection phase, bounds of all :yref:`Aabb's<Aabb>` are wrapped inside the cell in the first step. At subsequent runs, every bound remembers by how many cells it was initially shifted from coordinate given by the :yref:`Aabb` and uses this offset repeatedly as it is being updated from :yref:`Aabb` during particle's motion. Bounds are sorted using the periodic insertion sort algorithm (sect.~\ref{sect-periodic-insertion-sort}), which tracks periodic cell boundary $||$.

Upon inversion of two :yref:`Aabb`'s, their collision along all three axes is checked, wrapping real coordinates inside the cell for that purpose.

This algorithm detects collisions as if all particles were inside the cell but without the need of constructing "ghost particles" (to represent periodic image of a particle which enters the cell from the other side) or changing the particle's positions.

It is required by the implementation (and partly by the algorithm itself) that particles do not span more than half of the current cell size along any axis; the reason is that otherwise two (or more) contacts between both particles could appear, on each side. Since Yade identifies contacts by :yref:`Body.id` of both bodies, they would not be distinguishable.

In presence of shear, the sweep-and-prune collider could not sort bounds independently along three axes: collision along $x$ axis depends on the mutual position of particles on the $y$ axis. Therefore, bounding boxes *are expressed in transformed coordinates* which are perpendicular in the sense of collision detection. This requires some extra computation: :yref:`Aabb` of sphere in transformed coordinates will no longer be cube, but cuboid, as the sphere itself will appear as ellipsoid after transformation. Inversely, the sphere in simulation space will have a parallelepiped bounding "box", which is cuboid around the ellipsoid in transformed axes (the :yref:`Aabb` has axes aligned with transformed cell basis). This is shown in fig. `fig-cell-shear-aabb`_.

The restriction of a single particle not spanning more than half of the transformed axis becomes stringent as :yref:`Aabb` is enlarged due to shear. Considering :yref:`Aabb` of a sphere with radius $r$ in the cell where $x'\equiv x$, $z'\equiv z$, but $\angle(y,y')=\phi$, the $x$-span of the :yref:`Aabb` will be multiplied by $1/\cos\phi$. For the infinite shear $\phi\to\pi/2$, which can be desirable to simulate, we have $1/\cos\phi \to \infty$. Fortunately, this limitation can be easily circumvented by realizing the quasi-identity of all periodic cells which, if repeated in space, create the same grid with their corners: the periodic cell can be flipped, keeping all particle interactions intact, as shown in fig. `fig-cell-flip`_. It only necessitates adjusting the :yref:`Interaction.cellDist` of interactions and re-initialization of the collider (``Collider::invalidatePersistentData``). Cell flipping is implemented in the :yref:`yade.utils.flipCell` function.

.. _fig-cell-flip:
.. figure:: fig/cell-flip.*
	
	Flipping cell (:yref:`yade.utils.flipCell`) to avoid infinite stretch of the bounding boxes' spans with growing $\phi$. Cell flip does not affect interactions from the point of view of the simulation. The periodic arrangement on the left is the same as the one on the right, only the cell is situated differently between identical grid points of repetition; at the same time $|\phi_2|<|\phi_1|$ and sphere bounding box's $x$-span stretched by $1/\cos\phi$ becomes smaller. Flipping can be repeated, making effective infinite shear possible.

This algorithm is implemented in :yref:`InsertionSortCollider` and is used whenever simulation is periodic (:yref:`Omega.isPeriodic`); individual :yref:`BoundFunctor's<BoundFunctor>` are responsible for computing sheared :yref:`Aabb's<Aabb>`; currently it is implemented for spheres and facets (in :yref:`Bo1_Sphere_Aabb` and :yref:`Bo1_Facet_Aabb` respectively).

.. _fig-cell-shear-aabb:
.. figure:: fig/cell-shear-aabb.pdf
   
   Constructing axis-aligned bounding box (:yref:`Aabb`) of a sphere in simulation space coordinates (without periodic cell -- left) and transformed cell coordinates (right), where collision detection axes $x'$, $y'$ are not identical with simulation space axes $x$, $y$. Bounds' projection to axes is shown by orange lines.

Exact collision detection
^^^^^^^^^^^^^^^^^^^^^^^^^
When the collider detects approximate contact (on the :yref:`Aabb` level) and the contact does not yet exist, it creates *potential* contact, which is subsequently checked by exact collision algorithms (depending on the combination of :yref:`Shapes<Shape>`). Since particles can interact over many periodic cells (recall we never change their positions in simulation space), the collider embeds the relative cell coordinate of particles in the interaction itself (:yref:`Interaction.cellDist`) as an *integer* vector $c$. Multiplying current cell size $\mat{T}\vec{s}$ by $c$ component-wise, we obtain particle offset $\Delta \vec{x}$ in aperiodic $R^3$; this value is passed (from :yref:`InteractionLoop`) to the functor computing exact collision (:yref:`IGeomFunctor`), which adds it to the position of the particle :yref:`Interaction.id2`.

By storing the integral offset $c$, $\Delta\vec{x}$ automatically updates as cell parameters change.

.. _sect-periodic-insertion-sort:
Periodic insertion sort algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The extension of sweep and prune algorithm (described in :ref:`sect-sweep-and-prune`) to periodic boundary conditions is non-trivial. Its cornerstone is a periodic variant of the insertion sort algorithm, which involves keeping track of the "period" of each boundary; e.g. taking period $\langle 0,10)$, then $8_1\equiv-2_2<2_2$ (subscript indicating period). Doing so efficiently (without shuffling data in memory around as bound wraps from one period to another) requires moving period boundary rather than bounds themselves and making the comparison work transparently at the edge of the container.

This algorithm was also extended to handle non-orthogonal periodic :yref:`Cell` boundaries by working in transformed rather than Cartesian coordinates; this modifies computation of :yref:`Aabb` from Cartesian coordinates in which bodies are positioned (treated in detail in :ref:`sect-cell-approx-collision`).

The sort algorithm is tracking :yref:`Aabb` extrema along all axes. At the collider's initialization, each value is assigned an integral period, i.e. its distance from the cell's interior expressed in the cell's dimension along its respective axis, and is wrapped to a value inside the cell. We put the period number in subscript.
				
Let us give an example of coordinate sequence along $x$ axis (in a real case, the number of elements would be even, as there is maximum and minimum value couple for each particle; this demonstration only shows the sorting algorithm, however.)

.. math:: \sortlines{& 4_1 &\sortInv& 12_2 &\sortSep& -1_2 &\sortInv& -2_4 &\sortInv& 5_0 &}

with cell $x$-size $s_x=10$. The $4_1$ value then means that the real coordinate $x_i$ of this extremum is $x_i+1\cdot10=4$, i.e. $x_i=-4$. The $||$ symbol denotes the periodic cell boundary.

Sorting starts from the first element in the cell, i.e. right of $||$, and inverts elements as in the aperiodic variant. The rules are, however, more complicated due to the presence of the boundary $||$:

========================== ==============================================================
($\leq$)                   stop inverting if neighbors are ordered;
($||\bullet$)              current element left of $||$ is below 0 (lower period boundary); in this case, decrement element's period, decrease its coordinate by $s_x$ and move $||$ right;
($\bullet||$)              current element right of $||$ is above $s_x$ (upper period boundary); increment element's period, increase its coordinate by $s_x$ and move $||$ left;
($\crossBound$)            inversion across $||$ must subtract $s_x$ from the left coordinate during comparison. If the elements are not in order, they are swapped, but they must have their periods changed as they traverse $||$. Apply ($||\circ$) if necessary;
($||\circ$)                if after ($\crossBound$) the element that is now right of $||$ has $x_i<s_x$, decrease its coordinate by $s_x$ and decrement its period. Do not move $||$.
========================== ==============================================================

In the first step, ($||\bullet$) is applied, and inversion with $12_2$ happens; then we stop because of ($\leq$):

.. math::

				\sortlines{
					& 4_1 &\sortInv& 12_2 &\sortSep& \currelem{-1_2} &\sortInv& -2_4 &\sortInv& 5_0, &  \\
					& 4_1 &\sortInv& 12_2 &\sortInv& \isnleq{9_1} &\sortSep& -2_4 &\sortInv& 5_0, & \\
					& 4_1 &\sortInv& \isleq{9_1} &\sortInv& 12_2 &\sortSep& -2_4 &\sortInv& 5_0. &
				}

We move to next element $\currelem{-2_4}$; first, we apply ($||\bullet$), then invert until ($\leq$):

.. math::

				\sortlines{
					& 4_1 &\sortInv& 9_1 &\sortInv& 12_2 &\sortSep& \currelem{-2_4} &\sortInv& 5_0, & \\
					& 4_1 &\sortInv& 9_1 &\sortInv& 12_2 &\sortInv& \isnleq{8_3} &\sortSep& 5_0, & \\
					& 4_1 &\sortInv& 9_1 &\sortInv& \isnleq{8_3} &\sortInv& 12_2 &\sortSep& 5_0, & \\
					& 4_1 &\sortInv& \isleq{8_3} &\sortInv& 9_1 &\sortInv& 12_2 &\sortSep& 5_0. &
				}

The next element is $\currelem{5_0}$; we satisfy ($\crossBound$), therefore instead of comparing $12_2>5_0$, we must do $(12_2-s_x)=2_3\leq5$; we adjust periods when swapping over $||$ and apply ($||\circ$), turning $12_2$ into $2_3$; then we keep inverting, until ($\leq$):

.. math::

				\sortlines{
					& 4_1 &\sortInv& 8_3 &\sortInv& 9_1 &\sortInv& 12_2 &\sortSep& \isnleq{5_0}, & \\
					& 4_1 &\sortInv& 8_3 &\sortInv& 9_1 &\sortInv& \isnleq{5_{-1}} &\sortSep& 2_3, & \\
					& 4_1 &\sortInv& 8_3 &\sortInv& \isnleq{5_{-1}} &\sortInv& 9_1 &\sortSep& 2_3, & \\
					& 4_1 &\sortInv& \isleq{5_{-1}} &\sortInv& 8_3 &\sortInv& 9_1 &\sortSep& 2_3. & \\
				}

We move (wrapping around) to $\currelem{4_1}$, which is ordered:

.. math::
				\sortlines{
					& \currelem{4_1}\ar@(dr,dl)[rrrrrrrr]|{\geq} &\sortInv& 5_{-1} &\sortInv& 8_3 &\sortInv& 9_1 &\sortSep& 2_3 &
				}

and so is the last element

.. math:: 

				\sortlines{
					& 4_1 &\sortInv& \isleq{5_{-1}} &\sortInv& 8_3 &\sortInv& 9_1 &\sortSep& 2_3. &
				}

Computational aspects
=====================
Cost
----
The DEM computation using an explicit integration scheme demands a relatively high number of steps during simulation, compared to implicit scehemes. The total computation time $Z$ of simulation spanning $T$ seconds (of simulated time), containing $N$ particles in volume $V$ depends on:

* linearly, the number of steps $i=T/(s_t \Dtcr)$, where $s_t$ is timestep safety factor; $\Dtcr$ can be estimated by p-wave velocity using $E$ and $\rho$ (sect.~\ref{sect-dt-pwave}) as $\Dtcr^{(p)}=r\sqrt{\frac{\rho}{E}}$. Therefore

  .. math:: i=\frac{T}{s_t r}\sqrt{\frac{E}{\rho}}.

* the number of particles $N$; for fixed value of simulated domain volume $V$ and particle radius $r$

  .. math:: N=p\frac{V}{\frac{4}{3}\pi r^3},
  
  where $p$ is packing porosity, roughly $\frac{1}{2}$ for dense irregular packings of spheres of similar radius.

  The dependency is not strictly linear (which would be the best case), as some algorithms do not scale linearly; a case in point is the sweep and prune collision detection algorithm introduced in :ref:`sect-sweep-and-prune`_, with scaling roughly $\bigO{N \log N}$.
	
  The number of interactions scales with $N$, as long as packing characteristics are the same.
* the number of computational cores $\numCPU$; in the ideal case, the dependency would be inverse-linear were all algorithms parallelized (in Yade, collision detection is not).

Let us suppose linear scaling. Additionally, let us suppose that the material to be simulated ($E$, $\rho$) and the simulation setup ($V$, $T$) are given in advance. Finally, dimensionless constants $s_t$, $p$ and $\numCPU$ will have a fixed value. This leaves us with one last degree of freedom, $r$. We may write

.. math:: Z\propto iN\frac{1}{\numCPU}=\frac{T}{s_t r}\sqrt{\frac{E}{\rho}} p\frac{V}{\frac{4}{3}\pi r^3} \frac{1}{\numCPU}\propto \frac{1}{r}\frac{1}{r^3}=\frac{1}{r^4}.

This (rather trivial) result is essential to realize DEM scaling; if we want to have finer results, refining the "mesh" by halving $r$, the computation time will grow $2^4=16$ times.

For very crude estimates, one can use a known simulation to obtain a machine "constant"

.. math:: \mu=\frac{Z}{Ni}

with the meaning of time per particle and per timestep (in the order of $10^{-6}\,{\rm s}$ for current machines). $\mu$ will be only useful if simulation characteristics are similar and non-linearities in scaling do not have major influence, i.e. $N$ should be in the same order of magnitude as in the reference case.
				
Result indeterminism
--------------------
It is naturally expected that running the same simulation several times will give exactly the same results: although the computation is done with finite precision, round-off errors would be deterministically the same at every run. While this is true for *single-threaded* computation where exact order of all operations is given by the simulation itself, it is not true anymore in *multi-threaded* computation which is described in detail in later sections.
				
The straight-forward manner of parallel processing in explicit DEM is given by the possibility of treating interactions in arbitrary order. Strain and stress is evaluated for each interaction independently, but forces from interactions have to be summed up. If summation order is also arbitrary (in Yade, forces are accumulated for each thread in the order interactions are processed, then summed together), then the results can be slightly different. For instance ::

   (1/10.)+(1/13.)+(1/17.)=0.23574660633484162
   (1/17.)+(1/13.)+(1/10.)=0.23574660633484165

As forces generated by interactions are assigned to bodies in quasi-random order, summary force $F_i$ on the body can be different between single-threaded and multi-threaded computations, but also between different runs of multi-threaded computation with exactly the same parameters. Exact thread scheduling by the kernel is not predictable since it depends on asynchronous events (hardware interrupts) and other unrelated tasks running on the system; and it is thread scheduling that ultimately determines summation order of force contributions from interactions.
				
Numerical damping influence
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The effect of summation order can be significantly amplified by the usage of a *discontinuous* damping function in :yref:`NewtonIntegrator` given in :eq:`eq-damping-yade` as

.. math:: \frac{(\Delta\vec{F})_{dw}}{\vec{F}_w}=-\lambda_d\sign\vec{F}_w\left(\pprev{\dot{u}}_w+\frac{\curr{\ddot{\vec{u}}}_w\Dt}{2}\right).

If the $\sign$ argument is close to zero then the least significant finite precision artifact can determine whether the equation (relative increment of $\vec{F}_w$) is $+\lambda_d$ or $-\lambda_d$. Given commonly used values of $\lambda_d=0.4$, it means that such artifact propagates from least significant place to the most significant one at once.

Getting Started
===============

Installation
------------
*TBA*


Basic Usage
-----------
helanal is used to calculate helix properties for helices with at least
9 residues.

Here, we demonstrate the basic usage and options of helanal. First, import 
MDAnalysis and helanal::

    import MDAnalysis as mda
    import helanal

For this example, we will use datafiles from the MDAnalysis tests::

    from MDAnalysis.tests.datafiles import PSF, DCD
    u = mda.Universe(PSF, DCD)

To analyse a single helix, pass a selection with one atom per residue 
(normally this will be the Cα atoms)::

    hel = helanal.HELANAL(u, select='name CA and resnum 161-187')
    hel.run()

All computed properties for each simulation frame will be available as
arrays in ``.results``, e.g.::

    helix_tilt = hel.results.global_tilts

For more details about the properties calculated and method used, 
see :doc:`properties_computed`.


Further Options
---------------

- The **simulation frames** over which analysis is performed can be specified 
  using ``start``, ``stop``, and/or ``step``, or by providing a list 
  ``frames`` with which to slice the trajectory:: 

    hel.run(start=5, step=10)

- The **reference axis** to which helix tilts and screw angles are calculated 
  can be specified using ``ref_axis``::

    hel_xaxis = helanal.HELANAL(u, select='name CA and resnum 161-187',
                                ref_axis=[1,0,0])

- To analyse **multiple helices** at once, pass in a list of selection
  strings::

    hel_multi = helanal.HELANAL(u, select=('name CA and resnum 100-160',
                                           'name CA and resnum 200-230'))

  Each property in ``.results`` will now be a list of arrays, one for 
  each helix.  

- By default, the results of a single-helix analysis are flattened, such
  that each result is a single array, rather than a list-of-arrays of 
  length 1. This behaviour can be turned off (i.e. to be consistent 
  with the behaviour for multi-helix analysis) using 
  ``flatten_single_helix``::

    hel_noflat = helanal.HELANAL(u, select='name CA and resnum 161-187',
                                 flatten_single_helix=False)

- Analysis can be carried out directly on a single set of positions using
  the :func:`helix_analysis` function::

    hel_xyz = helanal.helix_analysis(u.atoms.positions, ref_axis=[0, 0, 1])

  Results are returned as a dictionary.



Citations
---------

helanal is based on the `HELANAL algorithm`_ from  [Bansal2000]_, which itself
uses the method of [Sugeta1967]_ to characterise each local axis. Please cite 
them when using this module in published work.

.. [Sugeta1967] Sugeta, H. and Miyazawa, T. 1967. General method for
   calculating helical parameters of polymer chains from bond lengths, bond
   angles and internal rotation angles. *Biopolymers* 5 673 - 679

.. [Bansal2000] Bansal M, Kumar S, Velavan R. 2000.
   HELANAL - A program to characterise helix geometry in proteins.
   *J Biomol Struct Dyn.*  17(5):811-819.

.. _`HELANAL algorithm`:
   https://web.archive.org/web/20090226192455/http://www.ccrnp.ncifcrf.gov/users/kumarsan/HELANAL/helanal.f

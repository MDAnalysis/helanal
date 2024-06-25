Getting Started
===============

Installation
------------
*TBA*


Example use
-----------
Import MDAnalysis and helanal::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD
    import helanal

You can pass in a single selection::

    u = mda.Universe(PSF, DCD)
    hel = helanl.HELANAL(u, select='name CA and resnum 161-187')
    hel.run()

All computed properties are available in ``.results``::

    print(hel.results.summary)

Alternatively, you can analyse several helices at once by passing
in multiple selection strings::

    hel2 = helanal.HELANAL(u, select=('name CA and resnum 100-160',
                                      'name CA and resnum 200-230'))

The :func:`helix_analysis` function will carry out helix analysis on
atom positions, treating each row of coordinates as an alpha-carbon
equivalent::

    hel_xyz = helanal.helix_analysis(u.atoms.positions, ref_axis=[0, 0, 1])


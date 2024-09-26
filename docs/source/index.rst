.. helanal documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to helanal's documentation!
=========================================================

This module contains code to analyse protein helices using the
HELANAL_ algorithm ([Bansal2000]_, [Sugeta1967]_).

HELANAL_ quantifies the geometry of helices in proteins on the basis of their
Cα atoms. It can determine local structural features such as the local
helical twist and rise, local helix origins and bending angles along the
helix.

.. toctree::
   :maxdepth: 1
   :caption: Contents

   getting_started
   properties_computed
   api


References
^^^^^^^^^^

.. _HELANAL: https://pubmed.ncbi.nlm.nih.gov/10798526/

.. [Sugeta1967] Sugeta, H. and Miyazawa, T. 1967. General method for
   calculating helical parameters of polymer chains from bond lengths, bond
   angles and internal rotation angles. *Biopolymers* 5 673 - 679

.. [Bansal2000] Bansal M, Kumar S, Velavan R. 2000.
   HELANAL - A program to characterise helix geometry in proteins.
   *J Biomol Struct Dyn.*  17(5):811-819.


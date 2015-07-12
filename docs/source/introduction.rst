Introduction
============

The Python Library of Earth Spherical Geometry (PyESG) is aiming to perform spherical geometry on the earth and interpolate data from an unstructured mesh grid to another one, which is the so-called "regridding" or "remapping".
Thus two main usages of PyESG are:

* Spherical geometry
* Regridding

and the second usage (regridding) is the main goal of this library.

There are some existing implementations of regridding, such as NCL_ and EMSF_.
NCL_ mainly performs interpolation between regular mesh grids, and recently (after the version of `6.1.0 <http://www.ncl.ucar.edu/Document/Functions/ESMF/ESMF_regrid.shtml>`_) it concludes the implementation of EMSF_, which performs interpolation between any mesh grids.
However, the APIs of EMSF_ are ported to specific global models that it is not so convinient to be used in a lightweighted way.
By a lightweighted way I mean calling a function like::

    func(data, mesh_old, mesh_new)

Therefore, PyESG is developed to make regridding between any mesh grids easier.

.. _NCL: http://www.ncl.ucar.edu
.. _EMSF: http://www.earthsystemmodeling.org

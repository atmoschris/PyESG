Regridding
==========

This section will introduce the algorithm for regridding in the PyESG.

Basic Idea
----------
The basic idea of regridding between two mesh grids, including the case of unstructured mesh grids, is interpolation.
Spherical interpolation is very complicated, and it is assumed in PyESG that for a small regional area, planar interpolation is appropriage with limit error.
Since the mesh grid is unstructured, the interpolation method used is **barycentric interpolation**.

To make the API easy to use, the Python function of regridding is like below::

    pyesg.Interp.regrid(mesh_old, mesh_new)

This function dose not regrid a specific data, but calculate the regridding information when one performs regridding from ``mesh_old`` to ``mesh_new``.
If the mesh grids are two-dimensional, say ``(nlat, nlon)``, then the shape of the return is ``(3, nlat, nlon)``, which stores the location information and weights of the three points that used to calculate the interpolated result.

So how to realize this function?

* Step 1: `Search`_

To calculate the interpolated value of a specific point ``(lat, lon)`` in the ``mesh_new``, we need to find out the three points in the ``mesh_old`` used for interpolation.

* Step 2: `Interpolation`_

After those three points are found, the barycentric interpolation can be performed.

Search
------


Interpolation
-------------

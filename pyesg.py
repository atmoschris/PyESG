#!/usr/bin/env python3

# ==============================================================================
# Author: Feng Zhu
# Date: 2015-07-12 12:16:04 Version 0.0.8: the bugs in spehrical geometry fixed
# Date: 2014-12-30 18:19:41 Version 0.0.3: a quick usable version
__version__ = '0.0.8'
# ==============================================================================
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import progressbar
import click

from IPython import embed

# ==============================================================================
# Constants
# ------------------------------------------------------------------------------


class Earth:

    radius = 6371  # km
    diameter = 2 * np.pi * radius  # km


class System:

    undef = -99999.

# ==============================================================================

# ==============================================================================
# Mathematics
# ------------------------------------------------------------------------------


class Math:

    def haversine(theta):

        hsin = (1 - np.cos(theta))/2.

        return hsin


# ==============================================================================

# ==============================================================================
# Public variables
# ------------------------------------------------------------------------------


class Public:

    aa, bb = -1, -1
    grids = []

    first_check = True
    first_search = True

# ==============================================================================

# ==============================================================================
# Spherical geometry
# ------------------------------------------------------------------------------


class Point(object):
    '''
    A point on the earth defined by the latitude and longitude (unit: degree).
    '''
    def __init__(self, lat, lon):
        if lat > 90 or lat < -90:
            raise ValueError('The range of latitude should be [-90, 90]!')

        if lon > 180 or lon < -180:
            raise ValueError('The range of longitude should be [-180, 180]!')

        self.lat = np.radians(lat)
        self.lon = np.radians(lon)

    def __str__(self):
        lat_deg = np.degrees(self.lat)
        lon_deg = np.degrees(self.lon)

        return '(' + str(lat_deg) + ', ' + str(lon_deg) + ')'

    def lat_deg(self):
        ''' radians -> degree
        '''
        lat_deg = np.degrees(self.lat)

        return lat_deg

    def lon_deg(self):
        ''' radians -> degree
        '''
        lon_deg = np.degrees(self.lon)

        return lon_deg

    def spherical_coord(self):
        ''' p (lat_rad, lon_rad) -> x, y, z

        Return the (x, y, z) in a UNIT spherical coordinate system.

        Reference: `<http://en.wikipedia.org/wiki/Spherical_coordinate_system>`_
        '''

        x = np.cos(self.lat) * np.cos(self.lon)
        y = np.cos(self.lat) * np.sin(self.lon)
        z = np.sin(self.lat)

        return x, y, z

    def vector(self):
        '''
        Return the vector from the center of the Earth to the point.
        '''
        vector = np.array(self.spherical_coord())

        return vector


class Arc(object):
    '''
    An arc on the earth defined by two points p1 and p2.
    It also can be seen as the relationship between two points.
    '''
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

    def __str__(self):
        return str(self.p1) + ' <-> ' + str(self.p2)

    def delta_lat(self):
        delta_lat = self.p2.lat - self.p1.lat

        return delta_lat

    def delta_lon(self):
        delta_lon = self.p2.lon - self.p1.lon

        return delta_lon

    def delta_sigma(self):
        '''
        Used a more solid equation for the calculation of delta sigma.

        Reference: `<http://en.wikipedia.org/wiki/https://en.wikipedia.org/wiki/Great-circle_distance>`_
        '''

        # delta_sigma = 2 * np.arcsin(
            # np.sqrt(
                # np.sin(self.delta_lat()/2.)**2 +
                # np.cos(self.p1.lat) *
                # np.cos(self.p2.lat) *
                # np.sin(self.delta_lon()/2)**2
            # )
        # )

        delta_sigma = np.arctan2(
            np.sqrt(
                (np.cos(self.p2.lat)*np.sin(self.delta_lon()))**2 +
                (
                    np.cos(self.p1.lat)*np.sin(self.p2.lat) -
                    np.sin(self.p1.lat)*np.cos(self.p2.lat) *
                    np.cos(self.delta_lon())
                )**2
            ),
            np.sin(self.p1.lat)*np.sin(self.p2.lat) +
            np.cos(self.p1.lat)*np.cos(self.p2.lat)*np.cos(self.delta_lon())
        )

        return delta_sigma

    def distance(self):
        '''
        The difficulty is the calculation of delta_sigma().
        Unit: km.

        Reference: `<http://en.wikipedia.org/wiki/https://en.wikipedia.org/wiki/Great-circle_distance>`_
        '''

        distance = Earth.radius * self.delta_sigma()

        return distance

    def rad(self):
        '''
        Convert great-circle distance on the earth to radians.
        '''

        rad = 2 * np.pi * self.distance() / Earth.diameter

        return rad

    def waypoint(self, k):
        ''' Calculate the location of a selected point (lat, lon) according to:
        * the location of point 1 (lat1, lon1);
        * the location of point 2 (lat2, lon2);
        * the coefficient k decides the position between point 1 and point 2,

        e.g.:
        * when k = 0.0, (lat, lon) is point 1;
        * when k = 0.5, (lat, lon) is the mid-point;
        * when k = 1.0, (lat, lon) is point 2.

        Reference: `<http://en.wikipedia.org/wiki/Great-circle_navigation>`_

        '''
        alpha1 = np.arctan2(
            np.sin(self.delta_lon()),
            np.cos(self.p1.lat)*np.tan(self.p2.lat) -
            np.sin(self.p1.lat)*np.cos(self.delta_lon())
        )

        alpha0 = np.arcsin(np.sin(alpha1) * np.cos(self.p1.lat))

        sigma01 = np.arctan2(
            np.tan(self.p1.lat),
            np.cos(alpha1)
        )  # Napier's rules

        sigma = k*self.delta_sigma() + sigma01

        lat = np.arctan2(
            np.cos(alpha0) * np.sin(sigma),
            np.sqrt(np.cos(sigma)**2 + np.sin(alpha0)**2*np.sin(sigma)**2)
        )

        lon01 = np.arctan2(
            np.sin(alpha0) * np.sin(sigma01),
            np.cos(sigma01)
        )

        lon0 = self.p1.lon - lon01

        lon_tmp = np.arctan2(
            np.sin(alpha0) * np.sin(sigma),
            np.cos(sigma)
        )

        lon = lon_tmp + lon0

        wp = Point(np.degrees(lat), np.degrees(lon))

        return wp


class Triangle(object):
    '''
    An triangle on the earth defined by three points p1, p2, and p3.
    It also can be seen as the relationship between three points.
    '''
    def __init__(self, p1, p2, p3):

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

        arc1 = Arc(self.p2, self.p3)
        arc2 = Arc(self.p3, self.p1)
        arc3 = Arc(self.p1, self.p2)

        if (
            Check.is_same_great_circle(arc1, arc2) or
            Check.is_same_great_circle(arc1, arc3) or
            Check.is_same_great_circle(arc2, arc3)
        ):
            print(self.p1, self.p2, self.p3)
            raise ValueError('The three given points are on the same arc!')

    def __str__(self):
        return str(self.p1) + ' <-> ' + str(self.p2) + ' <-> ' + str(self.p3)

    def angles(self):
        '''
        Calculate the included angle between two sides on the earth.
        If we set a is the side p2-p3, b the side p3-p1, and c the side p1-p2.
        Then the return value A is the included angle betwen sides b and c.

        References:
            * `<https://en.wikipedia.org/wiki/Spherical_law_of_cosines>`_
            * `<https://en.wikipedia.org/wiki/Haversine_formula>`_
        '''
        # arc1 = Arc(self.p2, self.p3)
        # arc2 = Arc(self.p3, self.p1)
        # arc3 = Arc(self.p1, self.p2)

        # a = arc1.rad()
        # b = arc2.rad()
        # c = arc3.rad()

        # A = np.arccos(
        #     (np.cos(a) - np.cos(b)*np.cos(c)) /
        #     (np.sin(b)*np.sin(c))
        # )

        # B = np.arccos(
        #     (np.cos(b) - np.cos(a)*np.cos(c)) /
        #     (np.sin(a)*np.sin(c))
        # )

        # C = np.arccos(
        #     (np.cos(c) - np.cos(a)*np.cos(b)) /
        #     (np.sin(a)*np.sin(b))
        # )

        w = self.p1.vector()  # A
        v = self.p2.vector()  # B
        u = self.p3.vector()  # C

        ta = u - np.dot(w, np.dot(w, u))
        tb = v - np.dot(w, np.dot(w, v))
        A = np.arccos(
            np.dot(ta/np.linalg.norm(ta), tb/np.linalg.norm(tb))
        )

        ta = u - np.dot(v, np.dot(v, u))
        tb = w - np.dot(v, np.dot(v, w))
        B = np.arccos(
            np.dot(ta/np.linalg.norm(ta), tb/np.linalg.norm(tb))
        )

        ta = v - np.dot(u, np.dot(u, v))
        tb = w - np.dot(u, np.dot(u, w))
        C = np.arccos(
            np.dot(ta/np.linalg.norm(ta), tb/np.linalg.norm(tb))
        )

        if np.isnan(A) or np.isnan(B) or np.isnan(C):
            print('WRONG:', A, B, C)
            print('%.20f %.20f' % (self.p1.lat_deg(), self.p1.lon_deg()))
            print('%.20f %.20f' % (self.p2.lat_deg(), self.p2.lon_deg()))
            print('%.20f %.20f' % (self.p3.lat_deg(), self.p3.lon_deg()))
            embed()
            raise ValueError('WRONG arccos()!!!')

        return A, B, C

    def area(self):
        ''' Calculate the area of the triangle bounded by the sides made by the
        three points p1 (lat1, lon1), p2 (lat2, lon2), and p3 (lat3, lon3)
        according to the Girard's Theorem:
        $area = R^2 * E$,
        where R is the radius of the sphere, and E the angle excess:
        $E = A + B + C - \pi$.
        Cosine rules are used to calculate the angles A, B, and C.

        References:
            * `<http://mathworld.wolfram.com/SphericalTriangle.html>`_
            * `<http://www.princeton.edu/~rvdb/WebGL/GirardThmProof.html>`_
            * `<http://en.wikipedia.org/wiki/Spherical_trigonometry>`_
            * `<http://mathforum.org/library/drmath/view/65316.html>`_
        ]

        '''
        A, B, C = self.angles()

        E = A + B + C - np.pi

        area = Earth.radius**2 * E

        return area


class Quadrangle(object):
    '''
    An quadrangle on the earth defined by three points p1, p2, p3, and p4.
    It also can be seen as the relationship between four points.

    Note: p1 -> p2 -> p3 -> p4 should be rotative.
    '''
    def __init__(self, p1, p2, p3, p4):
        if not Check.is_rotative(p1, p2, p3, p4):
            raise ValueError('The given four points are in wrong sequence!')

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

    def __str__(self):
        return (
            str(self.p1) + ' <-> ' + str(self.p2) +
            ' <-> ' + str(self.p3) + ' <-> ' + str(self.p4)
        )

    def angles(self):
        '''
        Treated as two triangles.
        '''
        tri1 = Triangle(self.p1, self.p2, self.p4)
        tri2 = Triangle(self.p2, self.p3, self.p4)

        A, B1, D1 = tri1.angles()
        C, D2, B2 = tri2.angles()

        B = B1 + B2
        D = D1 + D2

        return A, B, C, D

    def area(self):
        '''

        Reference: `<http://mathworld.wolfram.com/SphericalPolygon.html>`_
        '''

        A, B, C, D = self.angles()

        E = A + B + C + D - 2*np.pi

        area = Earth.radius**2 * E

        return area


class Mesh(object):
    '''
    Unstructed mesh grids, which is defined by
    2 dimenional arrays lat2d and lon2d.
    '''
    def __init__(self, lat2d, lon2d):
        self.lat2d = lat2d
        self.lon2d = lon2d

    def __str__(self):
        return str(self.lat2d.shape) + ' x ' + str(self.lon2d.shape)


class Check:

    def is_close_enough(v1, v2, e=1e-8):
        '''
        Check if two values are close enough.
        '''
        if abs(v1 - v2) <= e:
            return True

        else:
            return False

    def is_equal(p1, p2):
        '''
        Check if the given two points are equal.
        '''
        if p1 is None or p2 is None:
            return False

        elif (
            Check.is_close_enough(p1.lat, p2.lat) and
            Check.is_close_enough(p1.lon, p2.lon)
        ):
            return True

        else:
            return False

    def is_one_of_grids(point, mesh, first_time=True):
        '''
        Check if a given point is one of the mesh grids.
        '''
        nlat = mesh.lat2d.shape[0]
        nlon = mesh.lat2d.shape[1]

        if first_time:
            for i in np.arange(nlat):
                for j in np.arange(nlon):
                    pm = Point(mesh.lat2d[i, j], mesh.lon2d[i, j])
                    Public.grids.append((pm.lat, pm.lon))

        Public.first_check = False

        if (point.lat, point.lon) in Public.grids:
            return True
        else:
            return False

    def is_rotative(p1, p2, p3, p4):
        '''
        Check if the given four points are rotative (positive or negative).
        '''
        vec12 = np.array([p2.lat-p1.lat, p2.lon-p1.lon])
        vec23 = np.array([p3.lat-p2.lat, p3.lon-p2.lon])
        vec34 = np.array([p4.lat-p3.lat, p4.lon-p3.lon])
        vec41 = np.array([p1.lat-p4.lat, p1.lon-p4.lon])

        n123 = np.cross(vec12, vec23)
        n234 = np.cross(vec23, vec34)
        n341 = np.cross(vec34, vec41)

        if np.dot(n123, n234) > 0 and np.dot(n234, n341) > 0:
            return True
        else:
            return False

    def is_same_great_circle(arc1, arc2):
        '''
        Check if two given arcs are on the same great_circle.

        It is tricky to determine e.
        Since 1 degree = 1e5 m,
        when e=1e-8, the precision is about 1e-7 degree = 1 cm

        Generally, this method is quite robust.
        '''
        vec11 = arc1.p1.vector()
        vec12 = arc1.p2.vector()
        vec21 = arc2.p1.vector()
        vec22 = arc2.p2.vector()

        n1 = np.cross(vec11, vec12)
        n2 = np.cross(vec21, vec22)

        t = np.cross(n1, n2)
        # mag_t = np.sqrt(np.dot(t, t))
        mag_t = np.linalg.norm(t)

        if Check.is_close_enough(mag_t, 0):
            return True
        else:
            return False

    def is_on_great_circle(point, arc):
        '''
        Check if a given point is on the great_circle defined by the given arc.
        See is_same_great_circle.
        '''
        arc1 = Arc(point, arc.p1)
        arc2 = Arc(point, arc.p2)

        if Check.is_same_great_circle(arc1, arc2):
            return True
        else:
            return False

    def is_waypoint(point, arc, method='inner'):
        '''
        Check if a given point is on the given arc.

        This method is similar to is_same_great_circle and is_intersected.
        It is tricky to determin e.

        Reference: `<http://en.wikipedia.org/wiki/Great-circle_navigation>`_
        '''
        if not Check.is_on_great_circle(point, arc):
            # print('The given point is not on the given great-circle.')
            return False

        else:
            A = arc.p1.vector()
            B = arc.p2.vector()
            P = point.vector()

            N1 = np.cross(A, P)
            N2 = np.cross(B, P)
            T = np.dot(N1, N2)
            # print('[Debug] T:', T)

            if method == 'inner':
                if T <= 0:
                    return True
                else:
                    return False
            elif method == 'outer':
                if T > 0:
                    return True
                else:
                    return False
            else:
                raise ValueError('Wrong is_waypoint method!')

    def is_intersected(arc1, arc2):
        '''
        Check if two great-circle arcs are intersected.

        Reference:
            * `<http://www.mathworks.com/matlabcentral/newsreader/view_thread/276271>`_
            * `<http://stackoverflow.com/questions/2954337/great-circle-rhumb-line-intersection>`_
        '''

        vec11 = arc1.p1.vector()
        vec12 = arc1.p2.vector()
        vec21 = arc2.p1.vector()
        vec22 = arc2.p2.vector()

        n1 = np.cross(vec11, vec12)
        n2 = np.cross(vec21, vec22)

        t = np.cross(n1, n2)
        mag_t = np.sqrt(np.dot(t, t))

        if Check.is_close_enough(mag_t, 0):
            raise ValueError('The given two arcs lie on the same great-circle!')
        else:
            p1 = t / mag_t

        p2 = -p1

        lat1_deg = np.degrees(np.arcsin(p1[2]))
        lon1_deg = np.degrees(np.arctan2(p1[1], p1[0]))
        intersected_p1 = Point(lat1_deg, lon1_deg)

        lat2_deg = np.degrees(np.arcsin(p2[2]))
        lon2_deg = np.degrees(np.arctan2(p2[1], p2[0]))
        intersected_p2 = Point(lat2_deg, lon2_deg)

        if (
            Check.is_waypoint(intersected_p1, arc1) and
            Check.is_waypoint(intersected_p1, arc2)
        ):
            return True, intersected_p1

        elif (
            Check.is_waypoint(intersected_p2, arc1) and
            Check.is_waypoint(intersected_p2, arc2)
        ):
            return True, intersected_p2

        else:
            return False, None

    def is_on_triangle(point, triangle):
        '''
        Check if a given point (lat, lon) is on one of the sides
        of the given triangle.
        '''
        arc1 = Arc(triangle.p2, triangle.p3)
        arc2 = Arc(triangle.p3, triangle.p1)
        arc3 = Arc(triangle.p1, triangle.p2)

        if Check.is_waypoint(point, arc1, method='inner'):
            return True
        elif Check.is_waypoint(point, arc2, method='inner'):
            return True
        elif Check.is_waypoint(point, arc3, method='inner'):
            return True
        else:
            return False

    def is_inside_triangle(point, triangle):
        '''
        Check if a given point (lat, lon) is inside
        (including on) the given triangle.

        Reference: `<http://math.stackexchange.com/questions/1151428/point-within-a-spherical-triangle-given-areas>`_
        '''
        arc1 = Arc(triangle.p2, triangle.p3)
        arc2 = Arc(triangle.p3, triangle.p1)
        arc3 = Arc(triangle.p1, triangle.p2)

        if Check.is_on_triangle(point, triangle):
            return True

        if Check.is_waypoint(point, arc1, method='outer'):
            return False
        elif Check.is_waypoint(point, arc2, method='outer'):
            return False
        elif Check.is_waypoint(point, arc3, method='outer'):
            return False

        north_pole = Point(90, point.lon_deg())
        arc = Arc(point, north_pole)

        num_intersected_points = 0

        check1, inter_p1 = Check.is_intersected(arc, arc1)
        check2, inter_p2 = Check.is_intersected(arc, arc2)
        check3, inter_p3 = Check.is_intersected(arc, arc3)

        if check1 is True:
            num_intersected_points += 1

        if check2 is True:
            num_intersected_points += 1

        if check3 is True:
            num_intersected_points += 1

        if Check.is_equal(inter_p1, inter_p2):
            num_intersected_points -= 1

        if Check.is_equal(inter_p1, inter_p3):
            num_intersected_points -= 1

        if Check.is_equal(inter_p2, inter_p3):
            num_intersected_points -= 1

        if num_intersected_points % 2 == 0:
            return False  # outside
        else:
            return True  # inside

    def is_inside_quadrangle(point, quadrangle):
        '''
        Check if a given point (lat, lon) is inside the given quadrangle.
        Treated as two triangles.
        '''
        tri1 = Triangle(quadrangle.p1, quadrangle.p2, quadrangle.p4)
        tri2 = Triangle(quadrangle.p2, quadrangle.p3, quadrangle.p4)

        if (
            Check.is_inside_triangle(point, tri1) or
            Check.is_inside_triangle(point, tri2)
        ):
            return True
        else:
            return False


class Search:

    def quadrangle_bm(point, mesh):
        '''
        [bench mark algorithm]
        Search the quadrangle which the point located in with the hardest way.
        '''
        aa = None
        bb = None

        if Check.is_one_of_grids(point, mesh):
            raise ValueError('The given point is one of the mesh grids!')

        nlat = mesh.lat2d.shape[0]
        nlon = mesh.lat2d.shape[1]

        for i in np.arange(nlat-1):
            for j in np.arange(nlon-1):
                pm11 = Point(mesh.lat2d[i, j], mesh.lon2d[i, j])
                pm12 = Point(mesh.lat2d[i, j+1], mesh.lon2d[i, j+1])
                pm21 = Point(mesh.lat2d[i+1, j], mesh.lon2d[i+1, j])
                pm22 = Point(mesh.lat2d[i+1, j+1], mesh.lon2d[i+1, j+1])
                quadr = Quadrangle(pm11, pm12, pm22, pm21)

                if Check.is_inside_quadrangle(point, quadr):
                    aa, bb = i, j
                    return quadr, (aa, bb)

        if aa or bb is None:
            print('%.9f %.9f' % (point.lat_deg(), point.lon_deg()))
            raise ValueError('Quadrangle searching failed!')

    def quadrangle_o1(point, mesh, debug_mode=True):
        '''
        NOTE: DEBUG!!!

        [optimal version 1 based on the bench mark algorithm]
        Search the quadrangle which the point located in.
        '''
        aa = None
        bb = None

        if Check.is_one_of_grids(point, mesh):
            raise ValueError('The given point is one of the mesh grids!')

        nlat = mesh.lat2d.shape[0]
        nlon = mesh.lat2d.shape[1]

        i_min = 1e8
        i_max = -1e8
        for i in np.arange(nlat-1):
            if (
                point.lat_deg() >= np.min(mesh.lat2d[i, :]) and
                point.lat_deg() <= np.max(mesh.lat2d[i+1, :])
            ):
                if i < i_min:
                    i_min = i

                if i > i_max:
                    i_max = i

        j_min = 1e8
        j_max = -1e8
        for j in np.arange(nlon-1):
            if (
                point.lon_deg() >= np.min(mesh.lon2d[:, j]) and
                point.lon_deg() <= np.max(mesh.lon2d[:, j+1])
            ):
                if j < j_min:
                    j_min = j

                if j > j_max:
                    j_max = j

        smesh = Mesh(
            mesh.lat2d[i_min:i_max+1, j_min:j_max+1],
            mesh.lon2d[i_min:i_max+1, j_min:j_max+1]
        )

        nlat_small = smesh.lat2d.shape[0]
        nlon_small = smesh.lat2d.shape[1]

        for i in np.arange(nlat_small-1):
            for j in np.arange(nlon_small-1):
                pm11 = Point(smesh.lat2d[i, j], smesh.lon2d[i, j])
                pm12 = Point(smesh.lat2d[i, j+1], smesh.lon2d[i, j+1])
                pm21 = Point(smesh.lat2d[i+1, j], smesh.lon2d[i+1, j])
                pm22 = Point(smesh.lat2d[i+1, j+1], smesh.lon2d[i+1, j+1])
                quadr = Quadrangle(pm11, pm12, pm22, pm21)

                if Check.is_inside_quadrangle(point, quadr):
                    aa, bb = i+i_min, j+j_min
                    return quadr, (aa, bb)

        if aa or bb is None:
            print('%.9f %.9f' % (point.lat_deg(), point.lon_deg()))
            if debug_mode:
                np.savetxt('mesh_lat2d.dat', mesh.lat2d, '%.9f')
                np.savetxt('mesh_lon2d.dat', mesh.lon2d, '%.9f')
            raise ValueError('Quadrangle searching failed!')

    def quadrangle_o2(
        point, mesh, debug_mode=True, first_guess=(0, 0)
    ):
        '''
        NOTE: DEBUG!!!

        [optimal version 2]
        Search the quadrangle which the point located in.
        In version 2, first_guess is given to test first, if not found,
        then searching will be continued around the first guess quadrangle.

        The first_guess (aa, bb) is the lower left cornor of
        '''
        quadr = Form.quadrangle(mesh, first_guess)

        if Check.is_inside_quadrangle(point, quadr):
            # print(first_guess)
            return quadr, first_guess

        else:
            # if the first_guess is wrong, then
            # 0. N = 2
            # 1. test lower_left, if wrong, then
            # 2. move right x N, move up x N, move left x N, move down x N,
            #   if still not found, then
            # 3. N = N + 2, goto step 1

            N = 0
            last_guess = first_guess

            while N <= 1e8:

                N += 2
                # print('N:', N)
                # print('last_guess:', last_guess)

                # first check lower left
                guess = Move.lower_left(last_guess)
                # print(guess)
                quadr = Form.quadrangle(mesh, guess)

                if quadr is not None:

                    if Check.is_inside_quadrangle(point, quadr):

                        return quadr, guess

                # start cycling

                # move right x N
                for i in np.arange(N):
                    guess = Move.right(guess)
                    # print(guess)
                    quadr = Form.quadrangle(mesh, guess)
                    if quadr is None:
                        continue
                    elif Check.is_inside_quadrangle(point, quadr):
                        # print(guess)
                        return quadr, guess

                # move up x N
                for i in np.arange(N):
                    guess = Move.up(guess)
                    # print(guess)
                    quadr = Form.quadrangle(mesh, guess)
                    if quadr is None:
                        continue
                    elif Check.is_inside_quadrangle(point, quadr):
                        return quadr, guess

                # move left x N
                for i in np.arange(N):
                    guess = Move.left(guess)
                    # print(guess)
                    quadr = Form.quadrangle(mesh, guess)
                    if quadr is None:
                        continue
                    elif Check.is_inside_quadrangle(point, quadr):
                        return quadr, guess

                # move down x N
                for i in np.arange(N):
                    guess = Move.down(guess)
                    # print(guess)
                    quadr = Form.quadrangle(mesh, guess)
                    if quadr is None:
                        continue
                    elif Check.is_inside_quadrangle(point, quadr):
                        return quadr, guess

                last_guess = guess

        if quadr is None:
            print('N:', N)
            print('%.9f %.9f' % (point.lat_deg(), point.lon_deg()))
            if debug_mode:
                np.savetxt('mesh_lat2d.dat', mesh.lat2d, '%.9f')
                np.savetxt('mesh_lon2d.dat', mesh.lon2d, '%.9f')
            raise ValueError('Quadrangle searching failed!')

    def triangle(point, mesh, first_time=True):
        '''
        Search the triangle which the point located in.
        '''

        if first_time:
            quadr, (aa, bb) = Search.quadrangle_bm(point, mesh)
            #  quadr, (aa, bb) = Search.quadrangle_o1(point, mesh)
            Public.aa, Public.bb = aa, bb
        else:
            quadr, (aa, bb) = Search.quadrangle_o2(
                point, mesh,
                first_guess=(Public.aa, Public.bb)
            )
            Public.aa, Public.bb = aa, bb

        # print(point, 'in', Form.quadrangle(mesh, (Public.aa, Public.bb)))

        tri1 = Triangle(quadr.p1, quadr.p2, quadr.p4)
        tri2 = Triangle(quadr.p2, quadr.p3, quadr.p4)

        cc = aa + 1
        dd = bb + 1

        Public.first_search = False

        if Check.is_inside_triangle(point, tri1):
            return tri1, (aa, bb), (aa, dd), (cc, bb)
        elif Check.is_inside_triangle(point, tri2):
            return tri2, (aa, dd), (cc, dd), (cc, bb)
        else:
            raise ValueError('Something wrong!')


class Interp:
    '''
    Interpolation algorithms.
    '''
    def barycentric(point, triangle):
        ''' (point, triangle) -> weight1, weight2, weight3

        Barycentric Interpolation: interpolation in a triangle.

        Reference: `https://classes.soe.ucsc.edu/cmps160/Fall10/resources/barycentricInterpolation.pdf>`_
        '''
        if not Check.is_inside_triangle(point, triangle):
            raise ValueError('The given point is not in the given triangle!')

        # print('Interpolating', point)
        # print(triangle)

        arc1 = Arc(triangle.p2, triangle.p3)
        arc2 = Arc(triangle.p1, triangle.p3)
        arc3 = Arc(triangle.p1, triangle.p2)

        if Check.is_waypoint(point, arc1):

            tri2 = Triangle(point, triangle.p1, triangle.p3)
            tri3 = Triangle(point, triangle.p1, triangle.p2)

            A1 = 0
            A2 = tri2.area()
            A3 = tri3.area()

        elif Check.is_waypoint(point, arc2):

            tri1 = Triangle(point, triangle.p2, triangle.p3)
            tri3 = Triangle(point, triangle.p1, triangle.p2)

            A1 = tri1.area()
            A2 = 0
            A3 = tri3.area()

        elif Check.is_waypoint(point, arc3):

            tri1 = Triangle(point, triangle.p2, triangle.p3)
            tri2 = Triangle(point, triangle.p1, triangle.p3)

            A1 = tri1.area()
            A2 = tri2.area()
            A3 = 0

        else:

            tri1 = Triangle(point, triangle.p2, triangle.p3)
            tri2 = Triangle(point, triangle.p1, triangle.p3)
            tri3 = Triangle(point, triangle.p1, triangle.p2)

            A1 = tri1.area()
            A2 = tri2.area()
            A3 = tri3.area()

        A = A1 + A2 + A3
        # A_tri = triangle.area()

        # if not Check.is_close_enough(A, A_tri):
        #     print(A_tri - A)
        #     print(Check.is_inside_triangle(point, triangle))
        #     raise ValueError(
        #               'The triangle is too big for barycentric interpolation!'
        #           )

        weight1 = A1 / A
        weight2 = A2 / A
        weight3 = A3 / A

        return weight1, weight2, weight3

    def regrid(mesh_old, mesh_new, method='standard'):
        ''' (mesh_old, mesh_new) -> matrix of weights and points indices.

        Calculate the remapping coefficients (weights)
        from an old mesh to a new mesh.
        * When method == standard,
        the search algorithm can resolve any situation but slow;
        * When method == quick,
        the situation is that mesh_old and mesh_new are very similar
        to each other with some points nudged.
        '''
        # nlat_old = mesh_old.lat2d.shape[0]
        # nlon_old = mesh_old.lat2d.shape[1]

        nlat_new = mesh_new.lat2d.shape[0]
        nlon_new = mesh_new.lat2d.shape[1]

        weights = np.ndarray(shape=(3, nlat_new, nlon_new))
        weights[:, :, :] = System.undef

        lat_index = np.ndarray(shape=(3, nlat_new, nlon_new))
        lat_index[:, :, :] = System.undef

        lon_index = np.ndarray(shape=(3, nlat_new, nlon_new))
        lon_index[:, :, :] = System.undef

        click.secho('Runing regrid...', fg='green')

        if method == 'quick':

            pbar = progressbar.ProgressBar()
            for i in pbar(np.arange(nlat_new)):
                for j in np.arange(nlon_new):

                    # serpentine scanning
                    if i % 2 == 1:
                        jm = j
                    else:
                        jm = nlon_new - 1 - j

                    p_new = Point(mesh_new.lat2d[i, jm], mesh_new.lon2d[i, jm])
                    p_old = Point(mesh_old.lat2d[i, jm], mesh_old.lon2d[i, jm])

                    if Check.is_equal(p_new, p_old):
                        continue

                    else:
                        # print('Interpolating for', p_new)
                        (
                            tri,
                            (p1_lat_index, p1_lon_index),
                            (p2_lat_index, p2_lon_index),
                            (p3_lat_index, p3_lon_index)
                        ) = Search.triangle(
                            p_new, mesh_old, first_time=Public.first_search
                        )

                        weights[:, i, jm] = Interp.barycentric(p_new, tri)
                        lat_index[:, i, jm] = (
                            p1_lat_index, p2_lat_index, p3_lat_index
                        )
                        lon_index[:, i, jm] = (
                            p1_lon_index, p2_lon_index, p3_lon_index
                        )

        elif method == 'standard':

            pbar = progressbar.ProgressBar()
            for i in pbar(np.arange(nlat_new)):
                for j in np.arange(nlon_new):

                    # serpentine scanning
                    if i % 2 == 1:
                        jm = j
                    else:
                        jm = nlon_new - 1 - j

                    p_new = Point(mesh_new.lat2d[i, jm], mesh_new.lon2d[i, jm])
                    p_old = Point(mesh_old.lat2d[i, jm], mesh_old.lon2d[i, jm])

                    if Check.is_one_of_grids(
                        p_new, mesh_old, first_time=Public.first_check
                    ):
                        continue

                    else:
                        # print('Interpolating for', p_new)
                        (
                            tri,
                            (p1_lat_index, p1_lon_index),
                            (p2_lat_index, p2_lon_index),
                            (p3_lat_index, p3_lon_index)
                        ) = Search.triangle(
                            p_new, mesh_old, first_time=Public.first_search
                        )

                        weights[:, i, jm] = Interp.barycentric(p_new, tri)
                        lat_index[:, i, jm] = (
                            p1_lat_index, p2_lat_index, p3_lat_index
                        )
                        lon_index[:, i, jm] = (
                            p1_lon_index, p2_lon_index, p3_lon_index
                        )

        else:
            raise ValueError('Wrong regrid method!')

        return weights, lat_index, lon_index


class Move:
    '''
    The class to move the quadrangle under control by
    the index of the lower left corner of the quadrangle.
    '''

    def lower_left(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa-1, bb-1)

    def lower_right(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa-1, bb+1)

    def upper_left(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa+1, bb-1)

    def upper_right(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa+1, bb+1)

    def up(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa+1, bb)

    def down(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa-1, bb)

    def left(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa, bb-1)

    def right(ll_corner):
        aa = ll_corner[0]
        bb = ll_corner[1]
        return (aa, bb+1)


class Form:

    def quadrangle(mesh, ll_corner):
        '''
        Form a quadrangle in the given mesh with the given lower left corner.
        '''
        aa, bb = ll_corner

        if aa < 0 or bb < 0:

            click.secho('The ll_corner does not exists!', fg='red')
            quadr = None

        else:

            cc = aa + 1
            dd = bb + 1

            p11 = Point(mesh.lat2d[aa, bb], mesh.lon2d[aa, bb])
            p12 = Point(mesh.lat2d[aa, dd], mesh.lon2d[aa, dd])
            p21 = Point(mesh.lat2d[cc, bb], mesh.lon2d[cc, bb])
            p22 = Point(mesh.lat2d[cc, dd], mesh.lon2d[cc, dd])

            quadr = Quadrangle(p11, p12, p22, p21)

        return quadr


class Plot:

    def points(*points, delta_deg=0.5, ms_point=300):
        '''
        Plot points.
        '''
        p0 = points[0]
        p_lon = p0.lon_deg()
        p_lat = p0.lat_deg()

        m = Basemap(
            llcrnrlon=p_lon-delta_deg, llcrnrlat=p_lat-delta_deg,
            urcrnrlon=p_lon+delta_deg, urcrnrlat=p_lat+delta_deg
        )

        for pt in points:
            m.scatter(pt.lon_deg(), pt.lat_deg(), ms_point)

    def point_mesh(point, mesh, delta_deg=0.5, ms_point=300, ms_mesh=30):
        '''
        Plot the location of a point along with a mesh.
        '''

        p_lat = point.lat_deg()
        p_lon = point.lon_deg()

        m = Basemap(
            llcrnrlon=p_lon-delta_deg, llcrnrlat=p_lat-delta_deg,
            urcrnrlon=p_lon+delta_deg, urcrnrlat=p_lat+delta_deg
        )

        m.drawcoastlines()
        m.drawmapboundary()

        m.scatter(mesh.lon2d, mesh.lat2d, ms_mesh, marker='o', color='k')
        m.scatter(p_lon, p_lat, ms_point, marker='*', color='r')

        plt.show()

    def point_quadrangle(point, quadr, delta_deg=.5, ms_point=300, ms_quadr=30):
        '''
        Plot the location of a point along with a mesh.
        '''

        p_lat, p_lon = point.lat_deg(), point.lon_deg()

        p1_lat, p1_lon = quadr.p1.lat_deg(), quadr.p1.lon_deg()
        p2_lat, p2_lon = quadr.p2.lat_deg(), quadr.p2.lon_deg()
        p3_lat, p3_lon = quadr.p3.lat_deg(), quadr.p3.lon_deg()
        p4_lat, p4_lon = quadr.p4.lat_deg(), quadr.p4.lon_deg()
        quadr_lat = np.array([p1_lat, p2_lat, p3_lat, p4_lat])
        quadr_lon = np.array([p1_lon, p2_lon, p3_lon, p4_lon])

        m = Basemap(
            llcrnrlon=p_lon-delta_deg, llcrnrlat=p_lat-delta_deg,
            urcrnrlon=p_lon+delta_deg, urcrnrlat=p_lat+delta_deg
        )

        m.drawcoastlines()
        m.drawmapboundary()

        m.scatter(quadr_lon, quadr_lat, ms_quadr, marker='o', color='k')
        m.scatter(p_lon, p_lat, ms_point, marker='*', color='r')

        plt.show()

    def point_triangle(point, tri, delta_deg=.5, ms_point=300, ms_quadr=30):
        '''
        Plot the location of a point along with a mesh.
        '''

        p_lat, p_lon = point.lat_deg(), point.lon_deg()

        p1_lat, p1_lon = tri.p1.lat_deg(), tri.p1.lon_deg()
        p2_lat, p2_lon = tri.p2.lat_deg(), tri.p2.lon_deg()
        p3_lat, p3_lon = tri.p3.lat_deg(), tri.p3.lon_deg()
        tri_lat = np.array([p1_lat, p2_lat, p3_lat])
        tri_lon = np.array([p1_lon, p2_lon, p3_lon])

        m = Basemap(
            llcrnrlon=p_lon-delta_deg, llcrnrlat=p_lat-delta_deg,
            urcrnrlon=p_lon+delta_deg, urcrnrlat=p_lat+delta_deg
        )

        m.drawcoastlines()
        m.drawmapboundary()

        m.scatter(tri_lon, tri_lat, ms_quadr, marker='o', color='k')
        m.scatter(p_lon, p_lat, ms_point, marker='*', color='r')

        plt.show()

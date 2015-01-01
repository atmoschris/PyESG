#!/usr/bin/env python3

#========================================================================================
# Author: Feng Zhu
# Date: 2014-12-30 18:19:41
#========================================================================================
import math
import numpy as np

#========================================================================================
# Constants
#----------------------------------------------------------------------------------------
class Earth():
    radius = 6371 # km
    diameter = 2 * math.pi * radius # km

#========================================================================================

#========================================================================================
# Spherical geometry
#----------------------------------------------------------------------------------------
class Point(object):
    '''
    A point on the earth defined by the latitude and longitude (unit: degree).
    '''
    def __init__(self, lat, lon):
        if lat > 90 or lat < -90:
            raise ValueError('The range of latitude should be [-90, 90]')

        if lon > 180 or lon < -180:
            raise ValueError('The range of longitude should be [-180, 180]')

        self.lat = math.radians(lat)
        self.lon = math.radians(lon)

    #def deg(self):
        #''' radians -> degree
        #'''
        #lat_deg = math.degrees(self.lat)
        #lon_deg = math.degrees(self.lon)

        #return lat_deg, lon_deg

    def __str__(self):
        ''' radians -> degree
        '''
        lat_deg = math.degrees(self.lat)
        lon_deg = math.degrees(self.lon)

        return '(' + str(lat_deg) + ', ' + str(lon_deg) + ')'

    def spherical_coord(self):
        ''' p (lat, lon) -> x, y, z

        Return the (x, y, z) in a UNIT spherical coordinate system.

        [reference: http://en.wikipedia.org/wiki/Spherical_coordinate_system]
        '''

        x = math.cos(self.lat) * math.cos(self.lon)
        y = math.cos(self.lat) * math.sin(self.lon)
        z = math.sin(self.lat)

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
        delta_sigma = 2 * math.asin(
            math.sqrt(
                math.sin(self.delta_lat()/2.)**2 + math.cos(self.p1.lat)*math.cos(self.p2.lat)*math.sin(self.delta_lon()/2)**2
                )
            )

        return delta_sigma

    def distance(self):
        '''
        Calculate the great-circle distance of the arc.

        [reference: http://en.wikipedia.org/wiki/Great-circle_distance]
        '''

        distance = Earth.radius * self.delta_sigma()

        return distance

    def rad(self):
        '''
        Convert great-circle distance on the earth to radians.
        '''

        rad =  2 * math.pi * self.distance() / Earth.diameter

        return rad

    def waypoint(self, k):
        ''' Calculate the location of a selected point (lat, lon) according to:
        a) the location of point 1 (lat1, lon1);
        b) the location of point 2 (lat2, lon2);
        c) the coefficient k decides the position between point 1 and point 2,
        e.g., when k = 0.0, (lat, lon) is point 1;
              when k = 0.5, (lat, lon) is the mid-point;
              when k = 1.0, (lat, lon) is point 2.

        [reference: http://en.wikipedia.org/wiki/Great-circle_navigation]
        '''
        alpha1 = math.atan2(
                math.sin(self.delta_lon()),
                (math.cos(self.p1.lat))*math.tan(self.p2.lat) - math.sin(self.p1.lat)*math.cos(self.delta_lon())
                )

        alpha0 = math.asin(
                math.sin(alpha1) * math.cos(self.p1.lat)
                )

        sigma01 = math.atan2(
                math.tan(self.p1.lat),
                math.cos(alpha1)
                )

        sigma02 = sigma01 + self.delta_sigma()

        sigma = k * (sigma02-sigma01) + sigma01

        lat = math.atan2(
                math.cos(alpha0) * math.sin(sigma),
                math.sqrt(
                    math.cos(sigma)**2 + math.sin(alpha0)**2*math.sin(sigma)**2
                    )
                )

        lon01 = math.atan2(
                math.sin(alpha0) * math.sin(sigma01),
                math.cos(sigma01)
                )

        lon0 = self.p1.lon - lon01

        lon_tmp = math.atan2(
                math.sin(alpha0) * math.sin(sigma),
                math.cos(sigma)
                )
        lon = lon_tmp + lon0

        wp = Point(math.degrees(lat), math.degrees(lon))

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

    def __str__(self):
        return str(self.p1) + ' <-> ' + str(self.p2) + ' <-> ' + str(self.p3)

    def angles(self):
        '''
        Calculate the included angle between two sides on the earth.
        If we set a is the side p2-p3, b the side p3-p1, and c the side p1-p2.
        Then the return value A is the included angle betwen sides b and c.
        '''
        arc1 = Arc(self.p2, self.p3)
        arc2 = Arc(self.p3, self.p1)
        arc3 = Arc(self.p1, self.p2)

        a = arc1.rad()
        b = arc2.rad()
        c = arc3.rad()

        A = math.acos(
                (math.cos(a) - math.cos(b)*math.cos(c)) /
                (math.sin(b)*math.sin(c))
                )
        B = math.acos(
                (math.cos(b) - math.cos(a)*math.cos(c)) /
                (math.sin(a)*math.sin(c))
                )
        C = math.acos(
                (math.cos(c) - math.cos(a)*math.cos(b)) /
                (math.sin(a)*math.sin(b))
                )

        return A, B, C

    def area(self):
        ''' Calculate the area of the triangle bounded by the sides made by the
        three points p1 (lat1, lon1), p2 (lat2, lon2), and p3 (lat3, lon3) according to
        the Girard's Theorem:
            area = R^2 * E,
        where R is the radius of the sphere, and E the angle excess:
            E = A + B + C - pi.
        Cosine rules are used to calculate the angles A, B, and C.

        [references:
        http://www.princeton.edu/~rvdb/WebGL/GirardThmProof.html
        http://en.wikipedia.org/wiki/Spherical_trigonometry
        http://mathforum.org/library/drmath/view/65316.html
        ]
        '''
        A, B, C = self.angles()

        E = A + B + C - math.pi

        area = Earth.radius**2 * E

        return area

class Quadrangle(object):
    '''
    An quadrangle on the earth defined by three points p1, p2, p3, and p4.
    It also can be seen as the relationship between four points.

    Note: p1 -> p2 -> p3 -> p4 should be rotative.
    '''
    def __init__(self, p1, p2, p3, p4):
        if Check.is_rotative(p1, p2, p3, p4) == False:
            raise ValueError('The given four points are in wrong sequence!')

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

    def __str__(self):
        return str(self.p1) + ' <-> ' + str(self.p2) + ' <-> ' + str(self.p3) + ' <-> ' + str(self.p4)

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

        [reference: http://mathworld.wolfram.com/SphericalPolygon.html]
        '''

        A, B, C, D = self.angles()

        E = A + B + C + D - 2*math.pi

        area = Earth.radius**2 * E

        return area

class Mesh(object):
    '''
    Unstructed mesh grids, which is defined by 2 dimenional arrays lat2d and lon2d.
    '''
    def __init__(self, lat2d, lon2d):
        self.lat2d = lat2d
        self.lon2d = lon2d

    def __str__(self):
        return str(self.lat2d.shape) + ' x ' + str(self.lon2d.shape)

class Check():

    def is_close_enough(v1, v2, e=0.0001):
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
        if p1 == None or p2 == None:
            return False

        elif p1.lat == p2.lat and p1.lon == p2.lon:
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

    def is_on_great_circle(point, arc):
        '''
        Check if a given point is on the great_circle defined by the given arc.
        '''
        arc1 = Arc(point, arc.p1)
        arc2 = Arc(point, arc.p2)

        vec11 = arc1.p1.vector()
        vec12 = arc1.p2.vector()
        vec21 = arc2.p1.vector()
        vec22 = arc2.p2.vector()

        n1 = np.cross(vec11, vec12)
        n2 = np.cross(vec21, vec22)

        t = np.cross(n1, n2)
        mag_t = math.sqrt(np.dot(t, t))

        if mag_t == 0:
            return True

        else:
            return False

    def is_waypoint(point, arc, method='inner'):
        '''
        Check if a given point is on the given arc.
        '''
        if Check.is_on_great_circle == False:
            return False

        arc1 = Arc(point, arc.p1)
        arc2 = Arc(point, arc.p2)

        d = arc.distance()
        d1 = arc1.distance()
        d2 = arc2.distance()

        if method == 'inner':
            if Check.is_close_enough(d1 + d2, d):
                return True
            else:
                return False

        elif method == 'outer':
            if Check.is_close_enough(max(d1,d2), d + min(d1,d2)):
                return True
            else:
                return False

    def is_intersected(arc1, arc2):
        '''
        Check if two great-circle arcs are intersected.

        [reference:
        http://www.mathworks.com/matlabcentral/newsreader/view_thread/276271
        http://stackoverflow.com/questions/2954337/great-circle-rhumb-line-intersection
        ]
        '''

        vec11 = arc1.p1.vector()
        vec12 = arc1.p2.vector()
        vec21 = arc2.p1.vector()
        vec22 = arc2.p2.vector()

        n1 = np.cross(vec11, vec12)
        n2 = np.cross(vec21, vec22)

        t = np.cross(n1, n2)
        mag_t = math.sqrt(np.dot(t, t))

        if mag_t != 0:
            p1 = t / mag_t
        else:
            raise ValueError('The given two arcs lie on the same great-circle!')

        p2 = -p1

        lat1_deg = math.degrees(math.asin(p1[2]))
        lon1_deg = math.degrees(math.atan2(p1[1], p1[0]))
        intersected_p1 = Point(lat1_deg, lon1_deg)

        lat2_deg = math.degrees(math.asin(p2[2]))
        lon2_deg = math.degrees(math.atan2(p2[1], p2[0]))
        intersected_p2 = Point(lat2_deg, lon2_deg)

        arc_p1_11 = Arc(intersected_p1, arc1.p1)
        arc_p2_11 = Arc(intersected_p2, arc1.p1)

        if arc_p1_11.distance() <= arc_p2_11.distance():
            intersected_p = intersected_p1
        else:
            intersected_p = intersected_p2

        arc_p_11 = Arc(intersected_p, arc1.p1)
        arc_p_12 = Arc(intersected_p, arc1.p2)
        arc_p_21 = Arc(intersected_p, arc2.p1)
        arc_p_22 = Arc(intersected_p, arc2.p2)

        d_p_11 = arc_p_11.distance()
        d_p_12 = arc_p_12.distance()
        d_p_21 = arc_p_21.distance()
        d_p_22 = arc_p_22.distance()

        d1 = arc1.distance()
        d2 = arc2.distance()

        if Check.is_close_enough(d_p_11+d_p_12, d1) and Check.is_close_enough(d_p_21+d_p_22, d2):
            is_intersected = True
            return is_intersected, intersected_p
        else:
            is_intersected = False
            return is_intersected, None

    def is_inside_triangle(point, triangle):
        '''
        Check if a given point (lat, lon) is inside the given triangle.
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

        if Check.is_waypoint(point, arc1, method='outer'):
            return False
        elif Check.is_waypoint(point, arc2, method='outer'):
            return False
        elif Check.is_waypoint(point, arc3, method='outer'):
            return False

        north_pole = Point(90, 0)
        arc = Arc(point, north_pole)

        num_intersected_points = 0

        check1, inter_p1 = Check.is_intersected(arc, arc1)
        check2, inter_p2 = Check.is_intersected(arc, arc2)
        check3, inter_p3 = Check.is_intersected(arc, arc3)

        if check1 == True:
            num_intersected_points += 1

        if check2 == True:
            num_intersected_points += 1

        if check3 == True:
            num_intersected_points += 1

        if Check.is_equal(inter_p1, inter_p2):
            num_intersected_points -= 1

        if Check.is_equal(inter_p1, inter_p3):
            num_intersected_points -= 1

        if Check.is_equal(inter_p2, inter_p3):
            num_intersected_points -= 1

        if num_intersected_points % 2 == 0:
            return False # outside
        else:
            return True # inside

    def is_inside_quadrangle(point, quadrangle):
        '''
        Check if a given point (lat, lon) is inside the given quadrangle.
        Treated as two triangles.
        '''
        tri1 = Triangle(quadrangle.p1, quadrangle.p2, quadrangle.p4)
        tri2 = Triangle(quadrangle.p2, quadrangle.p3, quadrangle.p4)

        if Check.is_inside_triangle(point, tri1) or Check.is_inside_triangle(point, tri2):
            return True
        else:
            return False

class Search():
    def quadrangle(point, mesh):
        '''
        Search the quadrangle which the point located in.
        '''
        nlat = mesh.lat2d.shape[0]
        nlon = mesh.lat2d.shape[1]

        p11_max = Point(mesh.lat2d[0,0], mesh.lon2d[0,0])
        p12_max = Point(mesh.lat2d[0,nlon-1], mesh.lon2d[0,nlon-1])
        p21_max = Point(mesh.lat2d[nlat-1,0], mesh.lon2d[nlat-1,0])
        p22_max = Point(mesh.lat2d[nlat-1,nlon-1], mesh.lon2d[nlat-1,nlon-1])

        quadr_max = Quadrangle(p11_max, p12_max, p22_max, p21_max)

        if Check.is_inside_quadrangle(point, quadr_max) == False:
            raise ValueError('The given point is not inside the mesh!')

        aa, bb = 0, 0 # left lower corner
        cc, dd = nlat-1, nlon-1 # right upper corner
        ac = cc - aa
        bd = dd - bb

        while ac > 1 or bd > 1:

            print(aa, bb, cc, dd)

            if ac != 1 and bd != 1:

                ee, ff = (aa+cc) // 2, (bb+dd) // 2

                p111 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                p112 = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                p121 = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                p122 = Point(mesh.lat2d[ee,ff], mesh.lon2d[ee,ff])

                quadr1 = Quadrangle(p111, p112, p122, p121)

                p211 = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                p212 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                p221 = Point(mesh.lat2d[ee,ff], mesh.lon2d[ee,ff])
                p222 = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])

                quadr2 = Quadrangle(p211, p212, p222, p221)

                p311 = Point(mesh.lat2d[ee,ff], mesh.lon2d[ee,ff])
                p312 = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])
                p321 = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])
                p322 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

                quadr3 = Quadrangle(p311, p312, p322, p321)

                p411 = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                p412 = Point(mesh.lat2d[ee,ff], mesh.lon2d[ee,ff])
                p421 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                p422 = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])

                quadr4 = Quadrangle(p411, p412, p422, p421)

                if Check.is_inside_quadrangle(point, quadr1):
                    cc = ee
                    dd = ff

                elif Check.is_inside_quadrangle(point, quadr2):
                    cc = ee
                    bb = ff

                elif Check.is_inside_quadrangle(point, quadr3):
                    aa = ee
                    bb = ff

                elif Check.is_inside_quadrangle(point, quadr4):
                    aa = ee
                    dd = ff

                ac = cc - aa
                bd = dd - bb

            elif ac == 1:

                ff = (bb+dd) // 2

                p111 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                p112 = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                p121 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                p122 = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])

                quadr1 = Quadrangle(p111, p112, p122, p121)

                p211 = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                p212 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                p221 = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])
                p222 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

                quadr2 = Quadrangle(p211, p212, p222, p221)

                if Check.is_inside_quadrangle(point, quadr1):
                    dd = ff

                elif Check.is_inside_quadrangle(point, quadr2):
                    bb = ff

                ac = cc - aa
                bd = dd - bb

            elif bd == 1:

                ee = (aa+cc) // 2

                p111 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                p112 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                p121 = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                p122 = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])

                quadr1 = Quadrangle(p111, p112, p122, p121)

                p211 = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                p212 = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])
                p221 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                p222 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

                quadr2 = Quadrangle(p211, p212, p222, p221)

                if Check.is_inside_quadrangle(point, quadr1):
                    cc = ee

                elif Check.is_inside_quadrangle(point, quadr2):
                    aa = ee

                ac = cc - aa
                bd = dd - bb

        print(aa, bb, cc, dd)
        p11 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
        p12 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
        p21 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
        p22 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

        return Quadrangle(p11, p12, p22, p21)

    def triangle(point, mesh):
        '''
        Search the triangle which the point located in.
        '''
        quadr = Search.quadrangle(point, mesh)

        tri1 = Triangle(quadr.p1, quadr.p2, quadr.p4)
        tri2 = Triangle(quadr.p2, quadr.p3, quadr.p4)

        if Check.is_inside_triangle(point, tri1):
            return tri1
        elif Check.is_inside_triangle(point, tri1):
            return tri2

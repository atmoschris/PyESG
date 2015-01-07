#!/usr/bin/env python3

#========================================================================================
# Author: Feng Zhu
# Date: 2014-12-30 18:19:41
__version__ = '0.0.13'
#========================================================================================
import math

import numpy as np

import progressbar
import colors

#========================================================================================
# Constants
#----------------------------------------------------------------------------------------
class Earth:
    radius = 6371 # km
    diameter = 2 * math.pi * radius # km

class System:
    undef = -99999.

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
            raise ValueError('The range of latitude should be [-90, 90]!')

        if lon > 180 or lon < -180:
            raise ValueError('The range of longitude should be [-180, 180]!')

        self.lat = math.radians(lat)
        self.lon = math.radians(lon)

    def __str__(self):
        lat_deg = math.degrees(self.lat)
        lon_deg = math.degrees(self.lon)

        return '(' + str(lat_deg) + ', ' + str(lon_deg) + ')'

    def lat_deg(self):
        ''' radians -> degree
        '''
        lat_deg = math.degrees(self.lat)

        return lat_deg

    def lon_deg(self):
        ''' radians -> degree
        '''
        lon_deg = math.degrees(self.lon)

        return lon_deg

    def spherical_coord(self):
        ''' p (lat_rad, lon_rad) -> x, y, z

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
        + the location of point 1 (lat1, lon1);
        + the location of point 2 (lat2, lon2);
        + the coefficient k decides the position between point 1 and point 2,

        e.g.:
        + when k = 0.0, (lat, lon) is point 1;
        + when k = 0.5, (lat, lon) is the mid-point;
        + when k = 1.0, (lat, lon) is point 2.
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

        arc1 = Arc(self.p2, self.p3)
        arc2 = Arc(self.p3, self.p1)
        arc3 = Arc(self.p1, self.p2)

        if Check.is_same_great_circle(arc1, arc2) or Check.is_same_great_circle(arc1, arc3) or Check.is_same_great_circle(arc2, arc3):
            print(self.p1, self.p2, self.p3)
            #import pdb; pdb.set_trace()
            raise ValueError('The three given points are on the same arc!')

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
        $area = R^2 * E$,
        where R is the radius of the sphere, and E the angle excess:
        $E = A + B + C - \pi$.
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
        if not Check.is_rotative(p1, p2, p3, p4):
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
        if p1 == None or p2 == None:
            return False

        elif Check.is_close_enough(p1.lat, p2.lat, e=1e-8) and Check.is_close_enough(p1.lon, p2.lon, e=1e-8):
            return True

        else:
            return False

    def is_one_of_grids(point, mesh):
        '''
        Check if a given point is one of the mesh grids.
        '''
        nlat = mesh.lat2d.shape[0]
        nlon = mesh.lat2d.shape[1]

        grids = []

        for i in np.arange(nlat):
            for j in np.arange(nlon):
                pm = Point(mesh.lat2d[i,j], mesh.lon2d[i,j])
                grids.append((pm.lat, pm.lon))

        #if (point.lat_deg(), point.lon_deg()) in grids:
        if (point.lat, point.lon) in grids:
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
        '''
        vec11 = arc1.p1.vector()
        vec12 = arc1.p2.vector()
        vec21 = arc2.p1.vector()
        vec22 = arc2.p2.vector()

        n1 = np.cross(vec11, vec12)
        n2 = np.cross(vec21, vec22)
        #print(n1)
        #print(n2)

        t = np.cross(n1, n2)
        mag_t = math.sqrt(np.dot(t, t))
        #print(mag_t)

        if Check.is_close_enough(mag_t, 0, e=1e-8):
            return True
        else:
            return False

    def is_on_great_circle(point, arc):
        '''
        Check if a given point is on the great_circle defined by the given arc.
        '''
        arc1 = Arc(point, arc.p1)
        arc2 = Arc(point, arc.p2)

        if Check.is_same_great_circle(arc1, arc2):
            return True
        else:
            return False

    def is_waypoint(point, arc, method='inner'):
        '''
        NOTE: When the length of the given arc is too long, the result may be wrong.

        Check if a given point is on the given arc.

        [reference: http://en.wikipedia.org/wiki/Great-circle_navigation]
        '''
        if not Check.is_on_great_circle(point, arc):
            #print('--------')
            return False

        arc1 = Arc(point, arc.p1)
        arc2 = Arc(point, arc.p2)
        #print(point)
        #print(arc.p1)
        #print(arc.p2)

        d = arc.distance()
        d1 = arc1.distance()
        d2 = arc2.distance()
        #print(Check.is_same_great_circle(arc1, arc2))
        #print(Check.is_same_great_circle(arc1, arc))
        #print('d, d1, d2', d, d1, d2)
        #print('d1 + d2 - d', d1 + d2 - d)

        if method == 'inner':
            if Check.is_close_enough(d1+d2, d, e=1e-2):
                return True
            else:
                return False

        elif method == 'outer':
            if Check.is_close_enough(max(d1,d2), d+min(d1,d2), e=1e-2):
                #print('close enough')
                return True
            else:
                #print('not close enough')
                return False

        else:
            raise ValueError('Wrong is_waypoint method!')

    def is_intersected(arc1, arc2):
        '''
        TODO: debug

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
        #print(vec11, vec12, vec21, vec22)

        n1 = np.cross(vec11, vec12)
        n2 = np.cross(vec21, vec22)

        t = np.cross(n1, n2)
        mag_t = math.sqrt(np.dot(t, t))
        #print(t)
        #print(mag_t)

        if Check.is_close_enough(mag_t, 0, e=1e-8):
            raise ValueError('The given two arcs lie on the same great-circle!')
        else:
            p1 = t / mag_t

        p2 = -p1
        #print(p1)
        #print(p2, p2[0], p2[1], p2[2])

        lat1_deg = math.degrees(math.asin(p1[2]))
        lon1_deg = math.degrees(math.atan2(p1[1], p1[0]))
        intersected_p1 = Point(lat1_deg, lon1_deg)
        #print(intersected_p1)
        #print(Check.is_on_great_circle(intersected_p1, arc1))
        #print(Check.is_on_great_circle(intersected_p1, arc2))
        #print(Check.is_waypoint(intersected_p1, arc1))
        #print(Check.is_waypoint(intersected_p1, arc2))

        lat2_deg = math.degrees(math.asin(p2[2]))
        lon2_deg = math.degrees(math.atan2(p2[1], p2[0]))
        intersected_p2 = Point(lat2_deg, lon2_deg)
        #print(intersected_p2)
        #print(arc1)
        #print(arc2)
        #print(Check.is_on_great_circle(intersected_p2, arc1))
        #print(Check.is_on_great_circle(intersected_p2, arc2))
        #print(Check.is_waypoint(intersected_p2, arc1))
        #print(Check.is_waypoint(intersected_p2, arc2))

        arc_p1_11 = Arc(intersected_p1, arc1.p1)
        arc_p2_11 = Arc(intersected_p2, arc1.p1)

        if arc_p1_11.distance() <= arc_p2_11.distance():
            intersected_p = intersected_p1
            #print('p = p1')
        else:
            intersected_p = intersected_p2
            #print('p = p2')

        #print(intersected_p)

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

        if Check.is_close_enough(d_p_11+d_p_12, d1, e=1e-2) and Check.is_close_enough(d_p_21+d_p_22, d2, e=1e-2):
            return True, intersected_p
        else:
            return False, None

    def is_on_triangle(point, triangle):
        '''
        Check if a given point (lat, lon) is on one of the sides of the given triangle.
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
        Check if a given point (lat, lon) is inside (including on) the given triangle.
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

        #north_pole = Point(90, 0)
        north_pole = Point(90, point.lon_deg())
        arc = Arc(point, north_pole)

        num_intersected_points = 0

        check1, inter_p1 = Check.is_intersected(arc, arc1)
        check2, inter_p2 = Check.is_intersected(arc, arc2)
        check3, inter_p3 = Check.is_intersected(arc, arc3)

        if check1 == True:
            #print('Intersected with arc1:', arc1, 'at', inter_p1)
            num_intersected_points += 1

        if check2 == True:
            #print('Intersected with arc2:', arc2, 'at', inter_p2)
            num_intersected_points += 1

        if check3 == True:
            #print('Intersected with arc3:', arc3, 'at', inter_p3)
            num_intersected_points += 1

        if Check.is_equal(inter_p1, inter_p2):
            #print('Same intersected point with arc1 and arc2.')
            num_intersected_points -= 1

        if Check.is_equal(inter_p1, inter_p3):
            #print('Same intersected point with arc1 and arc3.')
            num_intersected_points -= 1

        if Check.is_equal(inter_p2, inter_p3):
            #print('Same intersected point with arc2 and arc3.')
            num_intersected_points -= 1

        #print(num_intersected_points)
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
        #tri3 = Triangle(quadrangle.p1, quadrangle.p3, quadrangle.p4)
        #tri4 = Triangle(quadrangle.p1, quadrangle.p2, quadrangle.p3)

        if Check.is_inside_triangle(point, tri1) or Check.is_inside_triangle(point, tri2):
        #if Check.is_inside_triangle(point, tri1) or Check.is_inside_triangle(point, tri2) or Check.is_inside_triangle(point, tri3) or Check.is_inside_triangle(point, tri4):
            return True
        else:
            return False

class Search:
    def quadrangle(point, mesh):
        '''
        TODO: debug

        Search the quadrangle which the point located in.
        '''
        #print(point)
        #print(mesh)
        nlat = mesh.lat2d.shape[0]
        nlon = mesh.lat2d.shape[1]

        if Check.is_one_of_grids(point, mesh):
            raise ValueError('The given point is one of the mesh grids!')

        aa, bb = 0, 0 # lower left corner
        cc, dd = nlat-1, nlon-1 # upper right corner

        p11_max = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
        p12_max = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
        p21_max = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
        p22_max = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

        quadr_max = Quadrangle(p11_max, p12_max, p22_max, p21_max)

        if not Check.is_inside_quadrangle(point, quadr_max):
            raise ValueError('The given point is not inside the mesh!')

        ac = cc - aa
        bd = dd - bb

        while ac > 1 or bd > 1:

            last_ac = ac
            last_bd = bd

            if ac != 1 and bd != 1:

                ee, ff = (aa+cc) // 2, (bb+dd) // 2

                p111 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                p112 = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                p121 = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                p122 = Point(mesh.lat2d[ee,ff], mesh.lon2d[ee,ff])

                quadr1 = Quadrangle(p111, p112, p122, p121)
                #print('quadr1:',quadr1)
                #print(Check.is_inside_quadrangle(point, quadr1))

                p211 = p112
                p212 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                p221 = p122
                p222 = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])

                quadr2 = Quadrangle(p211, p212, p222, p221)
                #print('quadr2:',quadr2)
                #print(Check.is_inside_quadrangle(point, quadr2))

                p311 = p122
                p312 = p222
                p321 = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])
                p322 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

                quadr3 = Quadrangle(p311, p312, p322, p321)
                #print('quadr3:',quadr3)
                #print(Check.is_inside_quadrangle(point, quadr3))

                p411 = p121
                p412 = p122
                p421 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                p422 = p321

                quadr4 = Quadrangle(p411, p412, p422, p421)
                #print('quadr4:', quadr4)
                #print(Check.is_inside_quadrangle(point, quadr4))

                if Check.is_inside_quadrangle(point, quadr1):
                    #print('Point is inside quadr1:', quadr1)
                    cc = ee
                    dd = ff

                elif Check.is_inside_quadrangle(point, quadr2):
                    #print('Point is inside quadr2:', quadr2)
                    cc = ee
                    bb = ff

                elif Check.is_inside_quadrangle(point, quadr3):
                    #print('Point is inside quadr3:', quadr3)
                    aa = ee
                    bb = ff

                elif Check.is_inside_quadrangle(point, quadr4):
                    #print('Point is inside quadr4:', quadr4)
                    aa = ee
                    dd = ff

                else:
                    print(colors.red('ac != 1 and bd != 1, need to search around...'))

                    if not Check.is_waypoint(p112, Arc(p111,p212)):
                        tri1 = Triangle(p111, p112, p212)
                        if Check.is_inside_triangle(point, tri1):
                            print('tri1')
                            case = 1

                    if not Check.is_waypoint(p222, Arc(p212,p322)):
                        tri2 = Triangle(p212, p222, p322)
                        if Check.is_inside_triangle(point, tri2):
                            print('tri2')
                            case = 2

                    if not Check.is_waypoint(p321, Arc(p421,p322)):
                        tri3 = Triangle(p421, p321, p322)
                        if Check.is_inside_triangle(point, tri3):
                            print('tri3')
                            case = 3

                    if not Check.is_waypoint(p121, Arc(p421,p111)):
                        tri4 = Triangle(p421, p121, p111)
                        if Check.is_inside_triangle(point, tri4):
                            print('tri4')
                            case = 4

            elif ac == 1:

                ff = (bb+dd) // 2

                p111 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                p112 = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                p121 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                p122 = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])

                quadr1 = Quadrangle(p111, p112, p122, p121)

                p211 = p112
                p212 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                p221 = p122
                p222 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

                quadr2 = Quadrangle(p211, p212, p222, p221)

                if Check.is_inside_quadrangle(point, quadr1):
                    dd = ff

                elif Check.is_inside_quadrangle(point, quadr2):
                    bb = ff

                else:
                    print(colors.red('ac == 1, need to search around...'))

                    if not Check.is_waypoint(p112, Arc(p111,p212)):
                        tri1 = Triangle(p111, p112, p212)
                        if Check.is_inside_triangle(point, tri1):
                            print('tri1')
                            case = 1

                    if not Check.is_waypoint(p122, Arc(p121,p222)):
                        tri3 = Triangle(p122, p121, p222)
                        if Check.is_inside_triangle(point, tri3):
                            print('tri3')
                            case = 3

            elif bd == 1:

                ee = (aa+cc) // 2

                p111 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                p112 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                p121 = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                p122 = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])

                quadr1 = Quadrangle(p111, p112, p122, p121)

                p211 = p121
                p212 = p122
                p221 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                p222 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

                quadr2 = Quadrangle(p211, p212, p222, p221)

                if Check.is_inside_quadrangle(point, quadr1):
                    cc = ee

                elif Check.is_inside_quadrangle(point, quadr2):
                    aa = ee

                else:
                    print(colors.red('bd == 1, need to search around...'))

                    if not Check.is_waypoint(p122, Arc(p112,p222)):
                        tri2 = Triangle(p122, p112, p222)
                        if Check.is_inside_triangle(point, tri2):
                            print('tri2')
                            case = 2

                    if not Check.is_waypoint(p121, Arc(p111,p221)):
                        tri4 = Triangle(p121, p111, p221)
                        if Check.is_inside_triangle(point, tri4):
                            print('tri4')
                            case = 4

            #print((aa, bb), (cc, dd), (ee, ff))
            ac = cc - aa
            bd = dd - bb

            if ac == last_ac and bd == last_bd:

                print(colors.red('Need to search around...'))

                if case == 1:

                    cc = aa
                    last_aa = aa
                    for k in np.arange(1, last_aa):
                        print(k)
                        aa_tmp = aa - k

                        p111_tmp = Point(mesh.lat2d[aa_tmp,bb], mesh.lon2d[aa_tmp,bb])
                        p112_tmp = Point(mesh.lat2d[aa_tmp,ff], mesh.lon2d[aa_tmp,ff])
                        p121_tmp = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                        p122_tmp = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])
                        quadr1_tmp = Quadrangle(p111_tmp, p112_tmp, p122_tmp, p121_tmp)

                        p211_tmp = Point(mesh.lat2d[aa_tmp,ff], mesh.lon2d[aa_tmp,ff])
                        p212_tmp = Point(mesh.lat2d[aa_tmp,dd], mesh.lon2d[aa_tmp,dd])
                        p221_tmp = Point(mesh.lat2d[cc,ff], mesh.lon2d[cc,ff])
                        p222_tmp = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])
                        quadr2_tmp = Quadrangle(p211_tmp, p212_tmp, p222_tmp, p221_tmp)

                        if Check.is_inside_quadrangle(point, quadr1_tmp):
                            dd = ff
                            aa = aa_tmp
                            break

                        elif Check.is_inside_quadrangle(point, quadr2_tmp):
                            bb = ff
                            aa = aa_tmp
                            break

                elif case == 2:

                    bb = dd
                    last_dd = dd
                    for k in np.arange(1, nlon-last_dd):
                        print(k)
                        dd_tmp = dd + k

                        p111_tmp = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                        p112_tmp = Point(mesh.lat2d[aa,dd_tmp], mesh.lon2d[aa,dd_tmp])
                        p121_tmp = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                        p122_tmp = Point(mesh.lat2d[ee,dd_tmp], mesh.lon2d[ee,dd_tmp])
                        quadr1_tmp = Quadrangle(p111_tmp, p112_tmp, p122_tmp, p121_tmp)

                        p211_tmp = Point(mesh.lat2d[ee,bb], mesh.lon2d[ee,bb])
                        p212_tmp = Point(mesh.lat2d[ee,dd_tmp], mesh.lon2d[ee,dd_tmp])
                        p221_tmp = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                        p222_tmp = Point(mesh.lat2d[cc,dd_tmp], mesh.lon2d[cc,dd_tmp])
                        quadr2_tmp = Quadrangle(p211_tmp, p212_tmp, p222_tmp, p221_tmp)

                        if Check.is_inside_quadrangle(point, quadr1_tmp):
                            cc = ee
                            dd = dd_tmp
                            break

                        elif Check.is_inside_quadrangle(point, quadr2_tmp):
                            aa = ee
                            dd = dd_tmp
                            break

                elif case == 3:

                    aa = cc
                    last_cc = cc
                    for k in np.arange(1, nlat-last_cc):
                        print(k)
                        cc_tmp = cc + k

                        p111_tmp = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                        p112_tmp = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                        p121_tmp = Point(mesh.lat2d[cc_tmp,bb], mesh.lon2d[cc_tmp,bb])
                        p122_tmp = Point(mesh.lat2d[cc_tmp,ff], mesh.lon2d[cc_tmp,ff])
                        quadr1_tmp = Quadrangle(p111_tmp, p112_tmp, p122_tmp, p121_tmp)

                        p211_tmp = Point(mesh.lat2d[aa,ff], mesh.lon2d[aa,ff])
                        p212_tmp = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                        p221_tmp = Point(mesh.lat2d[cc_tmp,ff], mesh.lon2d[cc_tmp,ff])
                        p222_tmp = Point(mesh.lat2d[cc_tmp,dd], mesh.lon2d[cc_tmp,dd])
                        quadr2_tmp = Quadrangle(p211_tmp, p212_tmp, p222_tmp, p221_tmp)

                        if Check.is_inside_quadrangle(point, quadr1_tmp):
                            dd = ff
                            cc = cc_tmp
                            break

                        elif Check.is_inside_quadrangle(point, quadr2_tmp):
                            bb = ff
                            cc = cc_tmp
                            break

                elif case == 4:

                    dd = bb
                    last_bb = bb
                    for k in np.arange(1, last_bb):
                        print(k)
                        bb_tmp = bb - k

                        p111_tmp = Point(mesh.lat2d[aa,bb_tmp], mesh.lon2d[aa,bb_tmp])
                        p112_tmp = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                        p121_tmp = Point(mesh.lat2d[ee,bb_tmp], mesh.lon2d[ee,bb_tmp])
                        p122_tmp = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])
                        quadr1_tmp = Quadrangle(p111_tmp, p112_tmp, p122_tmp, p121_tmp)

                        p211_tmp = Point(mesh.lat2d[ee,bb_tmp], mesh.lon2d[ee,bb_tmp])
                        p212_tmp = Point(mesh.lat2d[ee,dd], mesh.lon2d[ee,dd])
                        p221_tmp = Point(mesh.lat2d[cc,bb_tmp], mesh.lon2d[cc,bb_tmp])
                        p222_tmp = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])
                        quadr2_tmp = Quadrangle(p211_tmp, p212_tmp, p222_tmp, p221_tmp)

                        if Check.is_inside_quadrangle(point, quadr1_tmp):
                            cc = ee
                            bb = bb_tmp
                            break

                        elif Check.is_inside_quadrangle(point, quadr2_tmp):
                            aa = ee
                            bb = bb_tmp
                            break
                        else:
                            raise ValueError('Wrong case!')

                ac = cc - aa
                bd = dd - bb

        #print(aa, bb, cc, dd)
        p11 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
        p12 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
        p21 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
        p22 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

        quadr = Quadrangle(p11, p12, p22, p21)
        #print(quadr)

        return quadr, aa, bb, cc, dd
        #return aa, bb, cc, dd, quadr

    def quadrangle_accurate(point, mesh):
        '''
        Search the quadrangle which the point located in. This method is accurate but slow.
        '''
        nlat = mesh.lat2d.shape[0]
        nlon = mesh.lat2d.shape[1]

        if Check.is_one_of_grids(point, mesh):
            raise ValueError('The given point is one of the mesh grids!')

        for i in np.arange(nlat-1):
            for j in np.arange(nlon-1):
                aa = i
                bb = j
                cc = i+1
                dd = j+1

                p11 = Point(mesh.lat2d[aa,bb], mesh.lon2d[aa,bb])
                p12 = Point(mesh.lat2d[aa,dd], mesh.lon2d[aa,dd])
                p21 = Point(mesh.lat2d[cc,bb], mesh.lon2d[cc,bb])
                p22 = Point(mesh.lat2d[cc,dd], mesh.lon2d[cc,dd])

                quadr = Quadrangle(p11, p12, p22, p21)

                if Check.is_inside_quadrangle(point, quadr):
                    return quadr, aa, bb, cc, dd

    def triangle(point, mesh):
        '''
        Search the triangle which the point located in.
        '''
        quadr, aa, bb, cc, dd = Search.quadrangle(point, mesh)
        #quadr, aa, bb, cc, dd = Search.quadrangle_accurate(point, mesh)
        #print(aa, bb, cc, dd)

        tri1 = Triangle(quadr.p1, quadr.p2, quadr.p4)
        tri2 = Triangle(quadr.p2, quadr.p3, quadr.p4)
        #print(point)
        #print(tri1)
        #print(tri2)

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

        [reference: https://classes.soe.ucsc.edu/cmps160/Fall10/resources/barycentricInterpolation.pdf]
        '''
        if not Check.is_inside_triangle(point, triangle):
            raise ValueError('The given point is not in the given triangle!')

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
        A_tri = triangle.area()

        #if not Check.is_close_enough(A, A_tri):
            #print(A_tri - A)
            #print(Check.is_inside_triangle(point, triangle))
            #raise ValueError('The triangle is too big for barycentric interpolation!')

        weight1 = A1 / A
        weight2 = A2 / A
        weight3 = A3 / A

        return weight1, weight2, weight3

    def regrid(mesh_old, mesh_new, method='standard'):
        ''' (mesh_old, mesh_new) -> matrix of weights and points indices.

        Calculate the remapping coefficients (weights) from an old mesh to a new mesh.
        + When method == standard, the search algorithm can resolve any situation but slow;
        + When method == quick, the situation is that mesh_old and mesh_new are very similar
        to each other with some points nudged.
        '''
        nlat_old = mesh_old.lat2d.shape[0]
        nlon_old = mesh_old.lat2d.shape[1]


        nlat_new = mesh_new.lat2d.shape[0]
        nlon_new = mesh_new.lat2d.shape[1]

        weights = np.ndarray(shape=(3,nlat_new,nlon_new))
        weights[:,:,:] = System.undef

        lat_index = np.ndarray(shape=(3,nlat_new,nlon_new))
        lat_index[:,:,:] = System.undef

        lon_index = np.ndarray(shape=(3,nlat_new,nlon_new))
        lon_index[:,:,:] = System.undef

        print(colors.green('Runing regrid...'))

        if method == 'quick':

            pbar = progressbar.ProgressBar()
            for i in pbar(np.arange(nlat_new)):
                for j in np.arange(nlon_new):
                    p_new = Point(mesh_new.lat2d[i,j], mesh_new.lon2d[i,j])
                    p_old = Point(mesh_old.lat2d[i,j], mesh_old.lon2d[i,j])

                    if Check.is_equal(p_new, p_old):
                        continue

                    else:
                        #print('Interpolating for', p_new)
                        tri, (p1_lat_index, p1_lon_index), (p2_lat_index, p2_lon_index), (p3_lat_index, p3_lon_index) = Search.triangle(p_new, mesh_old)

                        weights[:,i,j] = Interp.barycentric(p_new, tri)
                        lat_index[:,i,j] = p1_lat_index, p2_lat_index, p3_lat_index
                        lon_index[:,i,j] = p1_lon_index, p2_lon_index, p3_lon_index

        elif method == 'standard':

            for i in np.arange(nlat_new):
                print(colors.green('Processing the (' + str(i+1) + ' of ' + str(nlat_new) + ') row... ' + '{:3.0f}'.format((i+1)/nlat_new*100) + '%'))
                pbar = progressbar.ProgressBar()
                for j in pbar(np.arange(nlon_new)):
                    p_new = Point(mesh_new.lat2d[i,j], mesh_new.lon2d[i,j])
                    p_old = Point(mesh_old.lat2d[i,j], mesh_old.lon2d[i,j])

                    if Check.is_one_of_grids(p_new, mesh_old):
                        continue

                    else:
                        #print('Interpolating for', p_new)
                        tri, (p1_lat_index, p1_lon_index), (p2_lat_index, p2_lon_index), (p3_lat_index, p3_lon_index) = Search.triangle(p_new, mesh_old)

                        weights[:,i,j] = Interp.barycentric(p_new, tri)
                        lat_index[:,i,j] = p1_lat_index, p2_lat_index, p3_lat_index
                        lon_index[:,i,j] = p1_lon_index, p2_lon_index, p3_lon_index

        else:
            raise ValueError('Wrong regrid method!')

        return weights, lat_index, lon_index

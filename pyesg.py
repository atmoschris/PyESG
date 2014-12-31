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

    def is_inside(point, triangle):
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
            is_inside = False # outside
        else:
            is_inside = True # inside

        return is_inside

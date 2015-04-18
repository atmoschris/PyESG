#!/usr/bin/env python3
from pyesg import *
import numpy as np

p1 = Point(26.179439544677734, -79.83404541015625)
p2 = Point(26.19446258730258, -79.685885669798)
p3 = Point(26.17958576623122, -79.83482906697526)

arc1 = Arc(p1, p2)
arc2 = Arc(p2, p3)

print(Check.is_same_great_circle(arc1, arc2))

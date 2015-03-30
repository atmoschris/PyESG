#!/usr/bin/env python3
from pyesg import *
import numpy as np

p1 = Point(27.0, -77.2)
p2 = Point(27.963257, -76.791260)

arc = Arc(p1, p2)
gcd = arc.distance()

print(gcd)

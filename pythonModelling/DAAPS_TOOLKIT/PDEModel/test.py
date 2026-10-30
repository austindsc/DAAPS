import os
import numpy as np
from scipy.interpolate import interp1d

from .OpticsDataLoader import *
from . import MathEquations as ME

"""
call
import PDEhelpers.test as test
or whatever else you want

I was using this for testing pathes and getting the package working together
"""

## get the relative path to the data so we can run from anywhere
script_dir = os.path.dirname(os.path.realpath(__file__))
print("OpticsDataLoader running from", script_dir)  

dhTest = OpticsData()
print("got an optics data thing")
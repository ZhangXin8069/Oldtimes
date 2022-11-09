from __future__ import division
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import nest_asyncio
import matplotlib.pyplot as plt
from scipy import*
from sympy import *
from sympy.plotting import *
from zx_h import *
from zx0 import *
import sklearn as skl
x, y, z, t = symbols('x y z t')
k, m, n = symbols('k m n', integer=True)
f, g, h = symbols('f g h', cls=Function)
nest_asyncio.apply()
zx = 'nb'
cy = 'cy'

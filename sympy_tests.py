
# 1. pip install mpmath
# 2. pip install sympy

from sympy import *

# test working
x = Symbol('x')
limit(sin(x)/x,x,0) # 1
integrate(1/x,x) # log(x)


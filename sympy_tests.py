
# 1. pip install mpmath
# 2. pip install sympy

from sympy import *
import pyperclip

# test working
x = Symbol('x')
limit(sin(x)/x,x,0) # 1
integrate(1/x,x) # log(x)


x, y, z, theta = symbols('x y x theta')
diff(cos(x), x) # >>> -sin(x)

sigma_x, sigma_y, tau_xy = symbols('sigma_x sigma_y tau_xy')
sigma_theta = 1/2 * (sigma_x + sigma_y) + 1/2 * (sigma_x - sigma_y)* cos(2*theta) + tau_xy * sin (2* theta)
print(latex(sigma_theta))

diff(sigma_theta, theta)
print(latex(diff(sigma_theta, theta)))

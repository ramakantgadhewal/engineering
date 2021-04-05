
# pip install fluids thermo pint
# https://fluids.readthedocs.io/Examples

import fluids

from fluids.units import  *

from thermo.units import Stream
from math import *

# Example 7-1 Smooth Pipe (Plastic)
# Given: Water at 30°C is flowing through 20 metres of 50 mm standard wall 
# plastic pipe (smooth wall) at a rate of 200 litres per minute.
# Find: The Reynolds number and friction factor.

Q = 200*u.l/u.min
T = 30*u.degC
P = 2*u.bar # assumed
water = Stream('water', T=T, P=P, Q=Q)

NPS, Di, Do, t = nearest_pipe(Di=50*u.mm)
v = Q/(pi/4*Di**2)
Re = Reynolds(D=Di, rho=water.rho, mu=water.mu, V=v)
print('Reynolds number = %s' %Re)
fd = friction_factor(Re=Re, eD=_roughness['Glass']/Di)
print('Darcy friction factor = %s' %fd)

# Example 7.2 L, L over D and K from Kv for Conventional Type Valves
# A 150 mm class 125 Y-pattern globe valve with a Kv=500 flow coefficient is given.
# Calculate the resistance coefficient K, L/D equivalent, and the length for complete turbulence in the flow.
# Use schedule 40, 150 mm pipe as a reference.

D = .154*u.m
Kv = 500*u.m**3/u.hour

K = Kv_to_K(Kv, D)
L_D = L_equiv_from_K(K, fd=.015)
L = D*L_D

print('Loss coefficient = %s' %K)
print('Equivalent length = %s' % L_D)
print('Length for complete turbulence = %s' %L)

# Example 7.3 L, L over D, K, and Kv for Conventional Type Valves
# A 100 mm class 600 steel angle valve, has a full area seat.
# Calculate its resistance coefficient K, flow coefficient Kv, the equivalent length of it L/D, and the length for complete turbulent L.

NPS, Di, Do, t = nearest_pipe(Do=0.103*u.m, schedule='80')
fd = 0.0165 # provided - note equivalent length is proportional to this value
d = 0.0972*u.m # diameter of seat
K = K_angle_valve_Crane(D1=d, D2=Di, fd=fd, style=1)
Kv = K_to_Kv(K, d)
L_D = L_equiv_from_K(K, fd)
L = L_D*d

print('Loss coefficient = %s' %K)
print('Valve flow coefficient = %s' %Kv)
print('Equivalent length = %s' % L_D)
print('Length for complete turbulence = %s' %L)

# Example 7.18 Gas

# A natural gas pipeline (schedule 20 14””) is 100 miles long. 
# Inlet pressure is 1300 psia, and outlet pressure is 300 psia.
# The average temperature is 40 deg F. The gas composition is 75% methane, 
# 21% ethane, and 4 % propane.
# Find the flow rate in MMscfd.

from fluids.units import *
from math import pi
P1 = 1300*u.psi
P2 = 300*u.psi
T = 40*u.degF
L = 100*u.miles

mu = 1.1e-5*u.Pa #
# mu = 1.915E-5*u.Pa # A more correct value, but hinders matching the problem
NPS, Di, Do, t = nearest_pipe(NPS=14, schedule='20')
A = 0.25*pi*Di**2

from thermo import *
from thermo import PRMIX
Tcs = [190.564, 305.32, 369.83]
Pcs = [4599000.0, 4872000.0, 4248000.0]
omegas = [0.008, 0.098, 0.152]
MWs = [16.04246, 30.06904, 44.09562]
zs = [0.75, 0.21, .04]
MW = sum(zs[i]*MWs[i] for i in range(3))*u.g/u.mol

eos_1 = PRMIX(T=T.to(u.K).magnitude, P=P1.to(u.Pa).magnitude, zs=zs, Tcs=Tcs,
         Pcs=Pcs, omegas=omegas)
eos_2 = PRMIX(T=T.to(u.K).magnitude, P=P2.to(u.Pa).magnitude, zs=zs, Tcs=Tcs, 
        Pcs=Pcs, omegas=omegas)
eos_std = PRMIX(T=288.15, P=101325.0, zs=zs, Tcs=Tcs, Pcs=Pcs, omegas=omegas)

Vm1 = eos_1.V_g*u.m**3/u.mol
rho1 = (Vm1)**-1*MW

Vm2 = eos_2.V_g*u.m**3/u.mol
rho2 = (Vm2)**-1*MW

Vm_std = eos_std.V_g*u.m**3/u.mol
rho_std = (Vm_std)**-1*MW

roughness = 0.0018*u.inch

rho = 0.5*(rho1 + rho2)

v = 10.0 # Initial guess for velocity
Re = rho*v*Di/mu
fd = friction_factor(Re=Re, eD=roughness/Di)

for i in range(8):
    # Solve for velocity with sequential substitution
    m = isothermal_gas(rho, fd, P1=P1, P2=P2, L=L, D=Di)
    v = m/(A*rho)
    Re = rho*v*Di/mu
    fd = friction_factor(Re=Re, eD=roughness/Di)
Q = v*A
correction = rho_std/rho

# pint does not support mmscfd
mmscfd = Q.to(u.ft**3/u.day).magnitude/1e6/correction
print('According to the full isothermal equation, the flowrate is %g MMscfd' %(mmscfd))

### For reference only
# Note Z_avg should be used in the Panhandle and Weymouth equations
# However Crane omits them; they are not used here to match Crane.
Z_avg = 0.5*(eos_1.Z_g + eos_2.Z_g)
MW_air = 28.966*u.g/u.mol
SG = MW/MW_air

# Crane does not use the efficiency term on Weymouth
Q_Weymouth = Weymouth(SG, Tavg=T, L=L, D=Di, P1=P1, P2=P2, Zavg=1, E=1)

mmscfd = Q_Weymouth.to(u.ft**3/u.day).magnitude/1e6
print('According to the Weymouth equation, the flowrate is %g MMscfd' %(mmscfd))

Q_panhandle = Panhandle_A(SG, Tavg=T, L=L, D=Di, P1=P1, P2=P2, Zavg=1, E=0.92)

mmscfd = Q_panhandle.to(u.ft**3/u.day).magnitude/1e6
print('According to the Panhandle A equation, the flowrate is %g MMscfd' 
  %(mmscfd))

# Besides the simplifications Crane makes, the methods are fair approximations. 
# The values in Crane are 107.8, 105.1, and 128.2 mmscfd respectively. 
# The isothermal calculation employed by Crane is a simplified one, 
# explaining the difference - it does not account for expansion rigorously.
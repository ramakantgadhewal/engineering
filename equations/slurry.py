import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import plotly as py
import plotly.graph_objs as go

import plotly.express as px
import doctest


dataset = 2 # 1 or Example 3-8, 2 or Weir Tables 5.2 / 5.3

# Slurry Systems Handbook Example 3-8
# Figure 3-10 is showing Tau - Tau 0 but labelled the axis Tau??

if dataset == 1:
    rate_of_shear = [100, 150, 200, 300, 400, 500, 600, 700, 800]
    shear_stress = [10.93, 12.27, 13.49, 15.68, 17.66, 19.49, 21.2, 22.84, 24.43]
    tau_tau0 = [4.11, 5.45, 6.67, 8.87, 10.85, 12.67, 14.39, 16.03, 17.61]
elif dataset == 2:
    shear_stress = [20.8, 21.6, 23.24, 33.20]
    rate_of_shear = [35.5, 56.4, 99.2, 126.8]
else:
    pass


ax = sns.lineplot(x=rate_of_shear, y=shear_stress)
ax.set(ylim=(0,30))


## Interplolate and Extrapolate

points = np.array([(1, 1), (2, 4), (3, 1), (9, 3)])

# get x and y vectors
x = rate_of_shear # points[:,0]
y = shear_stress # points[:,1]
x_lam = rate_of_shear[:-1]
y_lam = shear_stress[:-1]

# calculate polynomial
z = np.polyfit(x_lam, y_lam, 2)
f = np.poly1d(z)

# other x values
other_x = np.array([1.2, 1.34, 1.57, 1.7, 3.6, 3.8, 3.9, 4.0, 5.4, 6.6, 
            7.2, 7.3, 7.7, 8, 8.9, 9.1, 9.3])
other_y = f(other_x)

# calculate new x's and y's
x_new = np.linspace(0, 30, 5)
y_new = f(x_new)

# Creating the dataset, and generating the plot
trace1 = go.Scatter(
    x=x,
    y=y,
    mode='markers',
    name='Data',
    marker=dict(
        size=12
    )
)

trace2 = go.Scatter(
    x=x_new, #other_x,
    y=y_new, #other_y,
    name='Interpolated/Extrapolated Data',
    mode='markers',
    marker=dict(
        symbol='square-open',
        size=12
    )
)

layout = go.Layout(
    title='Interpolation and Extrapolation of Y From X',
)

data2 = [trace1, trace2]
fig2 = go.Figure(data=data2, layout=layout)
fig2.show()

### Slurry Systems Handbook 
# Example 4-1
# A sample of slurry is sieved for particle size. The data is collected in the 
# laboratory (see Table 4-1). Plot the data on a logarithmic graph and 
# determine the d50.
# Solution
# The data is plotted in Figure 4-6; the d50 is determined to be 145 micros.

size_micron = [425, 300, 212, 150, 106, 75, 53, 45, 38, -38]
cumul_pass_prt = [97.2, 87.1, 68.3, 51.3, 35.9, 20.5, 14.5, 11.8, 10.8, 0]

px.line(x=size_micron, y=cumul_pass_prt, log_x=True, log_y=True)



def durand_velocity(F_L=None, Di=None, solid_SG=None, liquid_SG=1., 
                    Cv=None, d50=None):
    """Calculates the transitional velocity based Durand and Condolios (1952) 
       equation for uniformly sized sand and gravel. If no Durand factor is
       provided then it will be calculated per Schiller and Herbich (1991) 
       using volume concentration and d50 particle sizes.

    Args:
        F_L (float, optional): Durands factor. Defaults to None.
        Di (float, optional): Pipe inside diameter (m). Defaults to None.
        solid_SG (float, optional): Solids specific gravity. Defaults to None.
        liquid_SG (float, optional): Liquids specific gravity. Defaults to 1..
        Cv (float, optional): Volumetric concentration (0-1). Defaults to None.
        d50 (float, optional): Particle d50 size in micron. Defaults to None.
    
    Notes
    -----
    Example 4-2
    A slurry mixture has a d50 of 300 micron. The slurry is pumped in a 30 in pipe 
    with an ID of 28.28􏲑. The volumetric concentration is 0.27. Using Equations 
    4-4 and 4-2, deter- mine the speed of deposition for a sand–water mixture 
    if the specific gravity of sand is 2.65.

    Examples
    --------
    Example 4-2 from Abulnaga[1]
    >>> durand_velocity(d50=300, Di=28.25*0.0254, Cv=27, solid_SG=2.65)
    (0.9644571340435376, 4.647574091490368)\n

    Example A1-5 from Weir[2]
    >>> durand_velocity(d50=50, Di=0.20, Cv=.05, solid_SG=2.65, F_L=1.34)
    (1.34, 3.4090793021576955)\n

    References
    ----------
    [1]: Abulnaga, Baha E. Slurry Systems Handbook, 2002, McGraw Hill
    [2]: Weir Slurry Pumping Manual, 2002.

    """
    if Cv > 1:
        Cv /= 100
    if Di > 10:
        Di /= 1000 
    
    if F_L:
       pass
    else:
      # No F_L provided so calculate. 
      F_L = (1.3 * Cv**0.125) * (1 - np.e**(-6.9 * d50/1000)) 
    # Transitional velocity calculated
    V_D = F_L * np.sqrt(2 * scipy.constants.g * Di * 
                    (solid_SG-liquid_SG)/liquid_SG)
    return F_L, V_D      


F_L, V_D = durand_velocity(d50=500, Di=0.20, Cv=30, solid_SG=2.65)
doctest.testmod()

def slurry_concentrations(rho_s, rho_l, rho_m=None, Cw=None, Cv=None, m_s=None):
    """Calculates the slurry concentration parameters. 
    

    Args:
        Cw (float): Concentration by weight percent (0 - 1)
        Cv (float): Concentration by volume percent (0 - 1)
        rho_s (float): Solids specific gravity.
        rho_l (float): Liquids specific gravity.
        m_s (float, optional): Solids mass per unit of time. Defaults to None.

    Returns:
        tuple: Of parameters.

    Examples:
    >>> slurry_concentrations(Cw=0.3, rho_s=2.65, rho_l=1)
    (1.2296983758700697, 0, 0, 0.1392111368909513, 0.3)\n

    Example 1-A[1]
    >>> slurry_concentrations(Cw=0.3, rho_s=2.65, rho_l=1, m_s=140)
    (1.2296983758700697, 466.6666666666667, 379.49685534591197, 0.1392111368909513, 0.3)\n
    
    Example C11.3[3]
    >>> slurry_concentrations(rho_m=1.167, rho_s=1.4, rho_l=1, m_s=None)
    (1.167, 0, 0, 0.4175000000000002, 0.5008568980291347)\n

    Notes
    -----
    Example 1-A and Warman Slurry Handbook Section 3.12.
    
    References
    ----------
    [1]: Abulnaga, Baha E. Slurry Systems Handbook, 2002, McGraw Hill
    [2]: Weir Slurry Pumping Manual, 2000, 2002.
    [3]: Nayyar, Mohinder L. Piping Handbook, 7th Edition. (McGraw-Hill: 2000)
    """
    if Cw:
        if Cw > 1:
            Cw /= 100

    if rho_m:
        Cw = (rho_s * (rho_m - rho_l)) / (rho_m * (rho_s - rho_l))
        Cv = Cw * rho_m / rho_s


    if Cv:
        if Cv > 1:
            Cv /= 100
        _rho_m = ( (Cv * rho_s) + (1-Cv) )
        Cw = Cv * rho_s / _rho_m

    if m_s:    
        # Mass of equivalent volume of liquid equivalent to the solid content 
        m_eq_l_eq_s = m_s * rho_l / rho_s
        
        # Mass of liquid in the slurry mixture at a concentration by weight
        m_l = m_s * (1 - Cw)/Cw
        
        # Total mass of equivalent volume of liquid
        m_total_l_eq = m_l + m_eq_l_eq_s
        
        # Total mass of slurry mixture transported per unit of time 
        m_mix = m_s + m_l
        
        # Volumetric flow per unit of time
        q = m_total_l_eq / rho_l
        
        # Density of the slurry mixture
        rho_m = m_mix / q

        # Concentration of Solids by Volume
        Cv = Cw * rho_m / rho_s
        
        return rho_m, m_mix, q, Cv, Cw 
    
    # Density of the slurry mixture
    rho_m = 1 / ( (Cw / rho_s) + (1 - Cw)/rho_l)

    # Concentration of Solids by Volume
    Cv = Cw * rho_m / rho_s


    return rho_m, 0, 0, Cv, Cw


rho_m, m_mix, q, Cv, Cw = slurry_concentrations(Cw=50, rho_m=None, rho_s=1.4, 
                                            rho_l=1, m_s=100)
doctest.testmod()
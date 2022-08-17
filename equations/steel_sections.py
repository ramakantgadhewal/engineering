
import numpy as np
import pint
ureg = pint.UnitRegistry()


import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section



# Variables: https://sectionproperties.readthedocs.io/en/latest/rst/api.html#sectionproperties-class



# Channel Properties
# https://www.engineersedge.com/standard_material/Steel_channel_properties.htm

def unit_conversions_universal_metric(size: str, weight: float, area: float, 
        d: float, b_f: float, t_f: float, t_w: float, i_xx: float, i_yy: float, 
        x_b: float, inputUnits='u'):
    
    u_u = 'lb/ft' # Uniform load per length
    u_m = 'kg/m'
    a_u = 'in**2'
    a_m = 'mm**2'
    l_u = 'in'
    l_m = 'mm'
    i_u = 'in**4'
    i_m = 'mm**4'

    if inputUnits == 'u':

        size_weight = weight * ureg(u_u).to(u_m)
        size_area = area * ureg(a_u).to(a_m)
        size_d = d * ureg(l_u).to(l_m)
        size_bf = b_f * ureg(l_u).to(l_m)
        size_tf = t_f * ureg(l_u).to(l_m)
        size_tw = t_w * ureg(l_u).to(l_m)
        size_ixx = i_xx * ureg(i_u).to(i_m)
        size_iyy = i_yy * ureg(i_u).to(i_m)
        size_xb = x_b * ureg(l_u).to(l_m)

    elif inputUnits == 'm':

        size_weight = weight * ureg(u_m).to(u_u)
        size_area = area * ureg(a_m).to(a_u)
        size_d = d * ureg(l_m).to(l_u)
        size_bf = b_f * ureg(l_m).to(l_u)
        size_tf = t_f * ureg(l_m).to(l_u)
        size_tw = t_w * ureg(l_m).to(l_u)
        size_ixx = i_xx * ureg(i_m).to(i_u)
        size_iyy = i_yy * ureg(i_m).to(i_u)
        size_xb = x_b * ureg(l_u).to(l_u)


    print(f"Size {size} conversions:")
    print(f"Mass per unit length (weight): {size_weight:.1f~P}")
    print(f"Section area (area): {size_area:.1f~P}")
    print(f"Depth of section (d): {size_d:.1f~P}")
    print(f"Flange width (b_f): {size_bf:.1f~P}")
    print(f"Flange thickness (t_f): {size_tf:.1f~P}")
    print(f"Web thickness (t_w): {size_tw:.1f~P}")
    print(f"Second moment of area about centroidal x-axis (i_xx): {size_ixx:e~P}")
    print(f"Second moment of area about centroidal y-axis (i_yy): {size_iyy:e~P}")
    print(f"x_b: {size_xb:.1f~P}")
    


def get_steel_properties(d, b, t_f, t_w, r, n_r):

    geometry = steel_sections.channel_section(d, b, t_f, t_w, r, n_r) 
    geometry.create_mesh(mesh_sizes=[2.5])  # Adds the mesh to the geometry

    section = Section(geometry)
    section.calculate_geometric_properties()
    section.calculate_plastic_properties()
    section.calculate_warping_properties()

    section.plot_centroids()
    section.display_results()

if __name__ == '__main__':

    #### Unit Conversions from American to Metric

    # https://www.engineersedge.com/standard_material/Steel_channel_properties.htm
    # ('C8x13.75', 13.75, 4.04, 8, 2.343, 0.39, 0.303, 36.1, 1.53, 0.553)
    unit_conversions_universal_metric(size='C8x13.75', weight=13.75, area=4.04, 
            d=8, b_f=2.343, t_f=0.39, t_w=0.303, i_xx=36.1, i_yy=1.53, 
            x_b=0.553, inputUnits = 'u')


    #### STEEL PROPERTIES

    # AS/NZS 3679.1:2016 Figure D4 Parallel-Flange Channels
    # (d=200, b=75, t_f=12, t_w=6, r=12, n_r=8) # PFC200
    # (d=150, b=75, t_f=9.5, t_w=6, r=10, n_r=8) # PFC150

    # get_steel_properties(d=8, b=2.343, t_f=0.39, t_w=0.303, r=.2, n_r=8) # C8x13.75 Ixx = 36.1, Iyy = 1.53
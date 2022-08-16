
import numpy as np
import pint
ureg = pint.UnitRegistry()


import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section



# Variables: https://sectionproperties.readthedocs.io/en/latest/rst/api.html#sectionproperties-class



# Channel Properties
# https://www.engineersedge.com/standard_material/Steel_channel_properties.htm

def unit_conversions_universal_metric(size: str, weight: float, area: float, 
        d: float, bf: float, tf: float, tw: float, ixx: float, iyy: float, 
        xb: float):

    size_weight = weight * ureg('lb/ft').to('kg/m')
    size_area = area * ureg('in**2').to('mm**2')
    size_d = d * ureg('in').to('mm')
    size_bf = bf * ureg('in').to('mm')
    size_tf = tf * ureg('in').to('mm')
    size_tw = tw * ureg('in').to('mm')
    size_ixx = ixx * ureg('in**4').to('mm**4')
    size_iyy = iyy * ureg('in**4').to('mm**4')
    size_xb = xb * ureg('in').to('mm')

    print(f"Size {size} conversions:")
    print(f"Mass per unit length: {size_weight:.1f~P}")
    print(f"Section area: {size_area:.1f~P}")
    print(f"Depth of section: {size_d:.1f~P}")
    print(f"Flange width: {size_bf:.1f~P}")
    print(f"Flange thickness: {size_tf:.1f~P}")
    print(f"Web thickness: {size_tw:.1f~P}")
    print(f"Second moment of area about centroidal x-axis: {size_ixx:e~P}")
    print(f"Second moment of area about centroidal y-axis: {size_iyy:e~P}")
    print(f"{size_xb:.1f~P}")
    


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
    unit_conversions_universal_metric('C8x13.75', 13.75, 4.04, 8, 2.343, 0.39, 0.303, 
                36.1, 1.53, 0.553)


    #### STEEL PROPERTIES

    # AS/NZS 3679.1:2016 Figure D4 Parallel-Flange Channels
    # (d=200, b=75, t_f=12, t_w=6, r=12, n_r=8) # PFC200
    # (d=150, b=75, t_f=9.5, t_w=6, r=10, n_r=8) # PFC150

    get_steel_properties(d=8, b=2.343, t_f=0.39, t_w=0.303, r=.2, n_r=8) # C8x13.75 Ixx = 36.1, Iyy = 1.53
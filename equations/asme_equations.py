import numpy as np
import pint
ureg = pint.UnitRegistry()
# ASME BPVC-VIII_D1


def ug_131_e_2(P, P_d, W_A=None, do=None, k_D=0.64, C=2407, 
                solve_for="W_A", A=None, M=None, T=None, w= 1000, 
                W_T=None, Z=None):
    """
    UG-127 Nonreclosing Pressure Relief Devices (a) Rupture Disk Devices. 
    Subsection UG-127(a)(2)(-a)(-1)(+c) referring to UG-131(e)(2).

    Args:
        do (float: [mm] diameter
        P (float): [kPa] Set pressure plus atmospheric
        P_d (float): [kPa] pressure at discharge from device
        k_discharge (float, optional): 0.64 per UG-127(a)(2)(-a)(-1)(+c).
        C (integer): constant for gas or vapour. 
        A (float, optional): [mm2] actual discharge area through the device. Defaults to None.
        M (float, optional): molecular weight. Defaults to None.
        T (float, optional): [C] absolute temperature. Defaults to None.
        w (float, optional): [kg/m3] specific weight of water. Defaults to 1000 (62.4 lbf/ft3).
        W_T (float, optional): [kg/hr] theoretical flow. Defaults to None.
        Z (float, optional): compressibility factor. Defaults to None.

    Returns:
        do float: [mm] Diameter based on Actual flow.
        W_A float: [kg/h] actual flow.
    """

    # convert to US Customary units for equation

    P = P * ureg.kilopascal
    P_ = P.to('psi')

    P_d = P_d * ureg.kilopascal
    P_d_ = P_d.to('psi')

    w = w * ureg("kilogram / meter ** 3")
    w_ = w.to("pounds / foot ** 3")

    # solve for actual flow
    if solve_for == "W_A":
        do = do * ureg.millimeters
        do_ = do.to('inch')

        A = np.pi / 4 * do**2
        A_ = A.to('inch ** 2')
    

        W_T_ = C * A_ * np.sqrt( (P_ - P_d_) * w_) # lb/hr
        W_T_ = W_T_.magnitude * ureg("lb / hour")

        W_T = W_T_.to("kg / hour") 
        W_A = W_T * k_D 

        return do, W_A

    # solve for diameter for actual flow provided.
    if solve_for == "do_a":
        W_T = W_A / k_D * ureg("kg / hour")
        W_A = W_T * k_D
        W_T_ = W_T.to("lbs / hour")

        A_ = W_T_ / ( C * np.sqrt( (P_ - P_d_) * w_) ) # lb/hr
        A_ = A_.magnitude * ureg("inch **2")
        do_ = np.sqrt( A_ * 4 / np.pi )
        do = do_.to('mm')

        return do, W_A

do, W_A = ug_131_e_2(1780, 101, do=None, W_A=800e3, solve_for="do_a")
# do, W_A = ug_131_e_2(1780, 101, do=80, W_A=None, solve_for="W_A")

print(f"Actual flow: {W_A}")
print(f"Diameter for actual flow: {do}")

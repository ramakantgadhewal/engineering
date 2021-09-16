import numpy as np
import pint
ureg = pint.UnitRegistry()
# ASME BPVC-VIII_D1


def asme_viii_d1_ug_131_e_2(P, P_d, W_A=None, do=None, k_D=0.64, C=2407, 
                solve_for="W_A", A=None, M=None, T=None, w= 1000, 
                W_T=None, Z=None):
    """
    NOTE: This equations has moved from ASME VIII to ASME VIII 2021.
    UG-127 Nonreclosing Pressure Relief Devices (a) Rupture Disk Devices. 
    Subsection UG-127(a)(2)(-a)(-1)(+c) referring to UG-131(e)(2).
    Refer to Mandatory Appendix 11 Capacity Conversions for Safety Valves

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

do, W_A = asme_viii_d1_ug_131_e_2(1780, 101, do=None, W_A=800e3, solve_for="do_a")
do, W_A = asme_viii_d1_ug_131_e_2(1780, 101, do=80, W_A=None, solve_for="W_A")

print(f"Actual flow: {W_A}")
print(f"Diameter for actual flow: {do}")


def asme_xiii_9_7_6_4_d_3(P, P_d, rho_l, k_D=0.64, C=5.092, W_A = None, do=None,
                            A=None, W_T=None, solve_for="W_A"):
    """
    ASME Section XIII 2021 Rules for Overpressure Protection.  9.7.6 Device
    Family Certification by the Coefficient of Discharge Method. 
    (.4 (d) For Tests With Water or Other Incompressible Fluids (3) For nozzle)



    Args:
        P (float): [MPaa] absolute relieving pressure 
        P_d (float): [MPaa] absolute discharge pressure
        rho_l (float): [kg/m3] density of fluid at device inlet condition 
        k_D (float, optional): Defaults to 0.64.
        C (float, optional): constant for water SI units. Defaults to 5.092.
        W_A (float, optional): Defaults to None.
        do (float, optional): actual diameter of discharge. Defaults to None.
        A (float, optional): [mm2] actual discharge area through the device
                            at developed lift. Defaults to None.
        W_T (float, optional): [kg/h] theoretical relieving capacity. Defaults to None.
        solve_for (str, optional): "W_A" solve for actual relieving capacity. 
                                   "do_a" solve for actual relieving diameter. Defaults to "W_A".

    Returns:
        [float]: [mm] diameter, [kg/hr] actual relieving capacity
    """

    # solve for actual flow
    if solve_for == "W_A":
        do = do * ureg("mm")

        A = np.pi / 4 * do**2

        W_T = C * A * np.sqrt( (P - P_d) * rho_l) 
        W_T = W_T.magnitude * ureg("kg / hour")

        W_A = W_T * k_D 

        return do, W_A

    # solve for diameter for actual flow provided.
    if solve_for == "do_a":
        W_T = W_A / k_D * ureg("kg / hour")
        W_A = W_T * k_D

        A = W_T / ( C * np.sqrt( (P - P_d) * rho_l) ) 
        A = A.magnitude * ureg("mm **2")
        do = np.sqrt( A * 4 / np.pi )

        return do, W_A


do, W_A = asme_xiii_9_7_6_4_d_3(1.780, 0.101, rho_l=1000, W_A=800e3, solve_for="do_a")
do, W_A = asme_xiii_9_7_6_4_d_3(1.780, 0.101, rho_l=1000, do=80, W_A=None, solve_for="W_A")


print(f"Actual flow: {W_A}")
print(f"Diameter for actual flow: {do}")
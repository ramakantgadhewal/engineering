import numpy as np

# ASME BPVC-VIII_D1


def ug_131_e_2(do, P, P_d, k_discharge=0.64, C=2407, A=None, M=None, T=None, 
                w= 62.4, W_T=None, Z=None):
    """
    UG-127 Nonreclosing Pressure Relief Devices (a) Rupture Disk Devices. 
    Subsection UG-127(a)(2)(-a)(-1)(+c) referring to UG-131(e)(2).

    Args:
        do (float: [in] diameter, [in]
        P (float): [psi] Set pressure plus atmospheric
        P_d (float): [psi] pressure at discharge from device
        k_discharge (float, optional): 0.64 per UG-127(a)(2)(-a)(-1)(+c).
        C (integer): constant for gas or vapor. 
        A (float, optional): [in2] actual discharge area through the device. Defaults to None.
        M (float, optional): molecular weight. Defaults to None.
        T (float, optional): [F] absolute temperature. Defaults to None.
        w (float, optional): [lbf/ft3] specifc weight of water. Defaults to 62.4.
        W_T (float, optional): [lb/hr] theoretical flow. Defaults to None.
        Z (float, optional): compressibility factor. Defaults to None.

    Returns:
        float: [kg/h] theroetical flow.
    """

    A = np.pi / 4 * do**2

    W_T = C * A * np.sqrt( (P - P_d) * w) # lb/hr
    W_T_si = W_T / 2.2046 # kg/h

    if k_discharge:
        W_T_si = W_T_si * k_discharge


    return W_T_si

W_T_si = ug_131_e_2(4, 275, 14.7, C=2407, k_discharge=0.6)
print(f"Pressure relief system shall not exceed: {W_T_si /1000:.1f} t/h of water")
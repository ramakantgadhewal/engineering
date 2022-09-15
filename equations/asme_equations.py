import numpy as np
import pint
ureg = pint.UnitRegistry()

from scipy import interpolate

import fluids

from fluids.units import  nearest_pipe, u
from mpmath import sec
from math import radians, cos

from misc import get_ideal_number_of_bins_per_sturgis_rule

###############################################################################
#                                                                             #
#                      ASME BPVC I Power Boilers                              #
#                                                                             #
###############################################################################



def asme_bpvc_I_pg_56(D, t_lug, t, b=1, S_a=None, S=None, e=None, l=None, 
                        W=None, W_ang=None):
    """
    ASME BPVC I Power Boilers PG-56 Loading on Structural Attachments. 
    2021 version page 49 (pdf page 111). 

    Args:
        e (float) eccentricity of W (see Figure PG-56.1.2), in. (mm)
        l (float): length of attachment of tube, in. (mm)
        W (float): eccentric load applied to lug, lb (N)
        W_ang (float): elastic load angle applied to lug, degree
        D (float): Outside diameter of the pipe, in. (mm)
        t (float): Attachment thickness, in. (mm)
        S_a (float): allowable stress value from Section II, Part D, Subpart 1, 
                        Table 1A, psi (MPa)
        S = pressure stress in tube determined by the equation in PG-27.2.1, 
                psi (MPa)

    Returns:
        K (float): tube attachment width design factor from Table PG-56.2, 
                dimensionless
        X (float): a parameter used to determine Lf, dimensionless, where
            b = unit width = 1.0 in. (25.4 mm)
        W_r (float): load component normal to tube axis, lb (N)
        L_f_c (float): a compression load factor, dimensionless
        L_f_t (float): a tension load factor, dimensionless
        S_t (float): portion of allowable stress available for attachment 
            loading from Step 3, psi (MPa)
        L_a_c (float): maximum allowable unit load in compression, lb/in (N/mm) 
            of attachment from PG-56.2. 
        L_a_t (float): maximum allowable unit load in tension, lb/in (N/mm) 
            of attachment from PG-56.2. 
        L_t (float): actual unit load calculated from PG-56.1.1, lb/in (N/mm)
        L_c (float): actual unit load calculated from PG-56.1.1, lb/in (N/mm)
        S_b_c
        S_b_t
    """


    # Step 1 Determine K from Table PG-56.2.

    # Table PG-56.2 Tube Attachment Angle Design Factor, K
    angle_of_attachment = '0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 \
                            80 85 90'
    K_values = '1.000 1.049 1.108 1.162 1.224 1.290 1.364 1.451 1.545 1.615 \
                1.730 1.836 1.949 2.076 2.221 2.341 2.513 2.653 2.876'

    angle_deg = np.arcsin(t_lug/(D/2)) * 180 / np.pi
    
    K =  mf.interpolate_1d_from_strings(angle_of_attachment, K_values, angle_deg )

    X = b * D / t**2

    # Step 2. Determine load factor, Lf,

    # (a) Compression Loading
    L_f_c = 1.618 * X **(-1.020-0.014 * np.log10(X) + 0.005 * (np.log10(X))**2)

    # (b) Tension Loading

    L_f_t = 49.937 * X**(-2.978 + 0.898 * np.log10(X) - 0.139*(np.log10(X)**2))

    # Step 3. Determine available stress, St.

    S_t = 2.0 * S_a - S


    # Step 4. Using values obtained in Steps 1 through 3, determine maximum 
    # allowable unit load, La.
    
    L_a_c = K * b * L_f_c * S_t
    L_a_t = K * b * L_f_t * S_t

    # PG-56.1.2 Calculate the actual unit load lb/in. (N/mm)
    
    W_r = W * np.sin(W_ang * np.pi / 180) # lb

    We = W * np.cos(W_ang * np.pi / 180) * e # lb in

    # W inward to pipe is negative for compression and positive for tension.
    # The second part is +/- so need to get absolute value. 
    L_t = W_r / l + abs(6 * We / l**2)
    L_c = W_r / l - abs(6 * abs(We) / l**2)

    # PG-56.1.1 Check if actual unit load is less than allowable

    if L_t < L_a_t:
        print(f'Actual unit load to allowable in tension ratio: \
                    {L_t/L_a_t:.2f}')
    else:
        print(f'Warning: Actual unit load greater than allowable in tension:\
                     {L_t/L_a_t:.2f}')
    if abs(L_c) < abs(L_a_c):
        print(f'Actual unit load less to allowable in compression ratio: \
                    {abs(L_c/L_a_c):.2f}')
    else:
        print(f'Warning: Actual unit load greater than allowable in compression:\
                     {abs(L_c/L_a_c):.2f}')

    # Calculate the bending stresses due to the lug load

    S_b_c = L_c / (K * b * L_f_c)
    S_b_t = L_t / (K * b * L_f_t)

    return K, X, L_f_c, L_f_t, S_t, L_a_c, L_a_t, L_t, L_c, S_b_c, S_b_t


###############################################################################
#                                                                             #
#      ASME BPVC VIII D1 Rules for Construction of Pressure Vessels           #
#                                                                             #
###############################################################################


def asme_bpvc_viii_d1_ug_131_e_2(P, P_d, W_A=None, do=None, k_D=0.64, C=2407, 
                solve_for="W_A", A=None, M=None, T=None, w= 1000, 
                W_T=None, Z=None):
    """
    NOTE: This equations has moved from ASME VIII to ASME XIII 2021.
    UG-127 Nonreclosing Pressure Relief Devices (a) Rupture Disk Devices. 
    Subsection UG-127(a)(2)(-a)(-1)(+c) referring to UG-131(e)(2).
    Refer to Mandatory Appendix 11 Capacity Conversions for Safety Valves

    Args:
        do (float: [mm] diameter
        P (float): [kPa] Set pressure plus atmospheric
        P_d (float): [kPa] pressure at discharge from device
        k_discharge (float, optional): 0.64 per UG-127(a)(2)(-a)(-1)(+c).
        C (integer): constant for gas or vapour. 
        A (float, optional): [mm2] actual discharge area through the device. 
                            Defaults to None.
        M (float, optional): molecular weight. Defaults to None.
        T (float, optional): [C] absolute temperature. Defaults to None.
        w (float, optional): [kg/m3] specific weight of water. 
                            Defaults to 1000 (62.4 lbf/ft3).
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


###############################################################################
#                                                                             #
#             ASME BPVC XIII Rules for Overpressure Protection                #
#                                                                             #
###############################################################################

def asme_bpvc_xiii_9_7_6_4_d_3(P, P_d, rho_l, K=0.62, C=5.092, W_A = None, 
                                 do=None, A=None, W_T=None, solve_for="W_A"):
    """
    ASME Section XIII 2021 Rules for Overpressure Protection.  9.7.6 Device
    Family Certification by the Coefficient of Discharge Method. 
    (.4 (d) For Tests With Water or Other Incompressible Fluids (3) For nozzle)



    Args:
        P (float): [MPaa] absolute relieving pressure 
        P_d (float): [MPaa] absolute discharge pressure
        rho_l (float): [kg/m3] density of fluid at device inlet condition 
        K (float, optional): coefficient of discharge per 4.1.3.2(a)(4) 
                            Defaults to 0.62.
        C (float, optional): constant for water SI units. Defaults to 5.092.
        W_A (float, optional): Defaults to None.
        do (float, optional): actual diameter of discharge. Defaults to None.
        A (float, optional): [mm2] actual discharge area through the device
                            at developed lift. Defaults to None.
        W_T (float, optional): [kg/h] theoretical relieving capacity. 
                                Defaults to None.
        solve_for (str, optional): "W_A" solve for actual relieving capacity. 
                                   "do_a" solve for actual relieving diameter. 
                                    Defaults to "W_A".

    Returns:
        [float]: [mm] diameter, [kg/hr] actual relieving capacity
    """

    # solve for actual flow
    if solve_for == "W_A":
        do = do * ureg("mm")

        A = np.pi / 4 * do**2

        W_T = C * A * np.sqrt( (P - P_d) * rho_l) 
        W_T = W_T.magnitude * ureg("kg / hour")

        W_A = W_T * K 

        return do, W_A

    # solve for diameter for actual flow provided.
    if solve_for == "do_a":
        W_T = W_A / K * ureg("kg / hour")
        W_A = W_T * K

        A = W_T / ( C * np.sqrt( (P - P_d) * rho_l) ) 
        A = A.magnitude * ureg("mm **2")
        do = np.sqrt( A * 4 / np.pi )

        return do, W_A


###############################################################################
#                                                                             #
#                       ASME B31.3 Process Piping                             #
#                                                                             #
###############################################################################

def asme_b31_3_304_1_2(P=None, D=None, E=0.85, W=0.80, Y=0.4, S= None, t=None,
                        c=None, d=None, solve_for='thk'):
    """
    ASME B31.3 304.1.2 Straight Pipe Under Internal Pressure. 
    Equation 3a for t < D/6.

    Args:
        P (float): internal design gage pressure (MPa).
        D (float): outside diameter of pipe as listed in tables of standards 
                    or specifications or as measured. (mm)
        S (float): stress value for material from Table A-1 or Table A-1M
                    A106 Gr.B = 138 MPa < 200°C.
                    A312 TP304L = 115 MPa < 150°C.
        E (float): quality factor from Table A-1A or Table A-1B
                    Note: Ranges from 0.85 ERW to 1.00 seamless.
        W (float): weld joint strength reduction factor in accordance with 
                    para. 302.3.5(e)
                    Note: Ranges from 0.60 to 1.00
        Y (float): coefficient from Table 304.1.1, valid for t < D/6 and for 
                    materials shown. The value of Y may be interpolated for 
                    intermediate temperatures. For t ≥ D/6, 
                    Y = (d + 2c) / (D + d + 2c)
                    Note: 482°C and below Y = 0.4 for all material except 
                    Gray iron is 0.0. 
        d (float): inside diameter of pipe. For pressure design calculation, 
                    the inside diameter of the pipe is the maximum value 
                    allowable under the purchase specification. (mm)
        c (float): sum of the mechanical allowances (thread or groove depth) 
                    plus corrosion and erosion allowances. (mm)
    
    Returns:
        t (float): pressure design thickness, as calculated in accordance with 
                    para. 304.1.2 for internal pressure or as determined in 
                    accordance with para. 304.1.3 for external pressure
                    

    Examples:
        asme_b31_3_304_1_2(2000,4,1,1,0.4,t=0.3) # 12533 psi Peng Peng pg 197
        asme_b31_3_304_1_2(2000/145.04,4*25.4,1,1,0.4,S=12533/145.04, t=1, solve_for='t') # 7.62 mm
        asme_b31_3_304_1_2(2000/145.04,4*25.4,1,1,0.4,S=12533/145.04, t=7.62, solve_for='P') # 13.79 MPa
        asme_b31_3_304_1_2(2000/145.04,4*25.4,1,1,0.4,S=12533/145.04, t=7.62, solve_for='S') # 86.41 MPa

    Setup:
        from sympy import var, Eq, solve, N, nsolve, latex
        var('t P D S E W Y')
        eq = Eq(t, P * D / (2 * (S * E * W + P * Y) ) )
        print(latex(eq))
        solve(eq, P) # Copy and paste below
        eq_P = Eq(P, 2*E*S*W*t/(D - 2*Y*t) )
        print(latex(eq_P))
        solve(eq, S) # Copy and paste below
        eq_S = Eq(S, P*(D - 2*Y*t)/(2*E*W*t))
        print(latex(eq_S))

    """

    if solve_for.lower() in  ['t', 'tk', 'thk', 'thickness']:
        print('yes')
        t = P * D / (2 * (S * E * W + P * Y) )
        if t < D/6:
            return t
        else:
            Y = (d + 2 * c) / (D + d + 2 *c)
            print(f'Thickness to diameter not satisfied. Y={Y}')
    
    if solve_for.lower() in ['s', 'stress' ]:
        S = P*(D - 2*Y*t)/(2*E*W*t)
        return S

    if solve_for.lower() in ['p', 'pressure']:
        P = 2*E*S*W*t/(D - 2*Y*t)
        return P





###############################################################################
#                                                                             #
#     ASME B31.4 Pipeline Transportation Systems for Liquids and Slurries     #
#                                                                             #
###############################################################################



def asme_b31_4_402_3(D, P_i, t, S_H=None):
    """
    ASME B31.4 402.3 Stress from Internal Pressure. For both restrained and 
    unrestrained pipelines, the circumferential (hoop) shress due to internal 
    pressure is calculated as: 
    S_H = P_i D / (2t)

    Args:
        D (float): outside diameter of pipe, in. (mm) 
        P_i (float): internal design gage pressure, psi (bar) 
        S_H (float): circumferential (hoop) stress due to internal pressure, 
        psi (bar) 
        t (float): wall thickness of pipe, in. (mm) 

    The above equation may not be applicable for pipe D/t less than 20.

    """

    if D/t < 20:
        print(f"The above equation is only applicable for D/t > 20 and your\
            ratio is {D/t}")
    
    else:
        # Changed equation from 20 to 2 and return bar and not MPa
        # This is done so pint can be used correctly.
        S_H = (P_i * D) / (2 * t)
        
        return S_H

def asme_b31_4_402_5_1(E, alpha, T_1, T_2, S_E = None):
    """
    ASME B31.4 402.5 Stress From Thermal Expansion. 
    402.5.1 Restrained Pipe. Thermal expansion stress in restrained pipe is
    calculated as: 

    S_E = E * alpha * (T_1 - T_2)

    Args:
        E (float): moduli of elasticity. 402.2.2 Moduli of Elasticity. 
          Flexibility calculations shall be based on the modulus of elasticity
          at ambient temperature. 
        alpha (float): coefficient of thermal expansion, in./in./°F (mm/mm/°C)
          402.2.1 Coefficient of Thermal Expansion. The linear coefficient of 
          thermal expansion for carbon and low alloy high tensile steel 
          may be taken as 6.5 x 10-6 in./in.°F for temperatures up to 250°F 
          (11.7 X 10-6 mm/mm/°C for temperatures up to 120°C.
        T_1 (float): temperature of the pipe at installation or completion of 
        final tie-in, °F (°C) 
        T_2 (float): operating temperature, °F (°C) 
        S_E (float, optional): thermal expansion stress, psi (MPa) . 
                                Defaults to None.
    
      In the above equation, compressive stress is a negative value.
    """


    S_E = E * alpha * (T_1 - T_2)

    return S_E


def asme_b31_4_402_6_1( S_E, S_H, Fa=0, M=0, Z=1, A=1, v=0.3):
    """
    ASME B31.4 402.6.1 Restrained Pipe. Longitudinal stress in restrained pipe
    is calculated as: 

    S_L = S_E + v S_H + M/Z + Fa/A

    Args:
        A (_type_): metal area of nominal pipe cross section, in^2 (cm^2)
        Fa (_type_): axial force, such as weight on a riser, lb (N)
        M (_type_): bending moment, in.-lb (N m)
        Z (_type_): section modulus of the pipe, in^3 (cm^3)
        S_E (_type_): thermal expansion stress, psi (MPa)
        S_H (_type_): circumferential (hoop) stress due to internal pressure, 
                        psi (MPa) 
        v (_type_): Poisson's ratio
          402.2.3 Poisson's ratio, v. Poisson's ratio shall be taken as 0.3 for
          steel
    """
    
    
    S_L = S_E + v * S_H + M/Z + Fa/A
    
    return S_L


def asme_b31_4_402_7(S_L, S_H, S_t= None):
    """
    402 Calculation of Stresses. 402.7 Combining of Stresses.
    In restrained pipe, the Longitudinal and circumferential (hoop) stress
    are combined in accordance with the maximum shear stress theory as follows:
    Args:
        S_L (float): _description_
        S_H (_type_): _description_
        S_t (_type_, optional): _description_. Defaults to None.
    """

    if S_t:
        S_eq = 2. * ( (S_L - S_H)**2 + S_t**2)**0.5
    else:
        S_eq = abs(S_L - S_H)
        print(S_eq)
    

    S_eq_det = (((S_H)**2 - S_H*S_L + S_L**2))**0.5

    return S_eq


def asme_b31_4_403_2_1(P_i, D, F=0.72, E=1, A=0, S=None, S_y=None, t_n= None):
    """
    403.2 Criteria for Pipe Wall Thickness and Allowances. Returns
    403.2.1 Criteria. The nominal wall thickness of straight sections of steel 
    pipe shall be equal to or greater than t11 determined in accordance with 
    the following equation: 

        t_n >= t + A


    Args:
        A (float): sum of allowances for threading, grooving, corrosion, and 
        erosion as required in paras. 403.2.2 through 403.2.4, and increase in 
        wall thickness if used as protective measure in para. 403.1 
        tn (float): nominal wall thickness satisfying requirements for pressure 
        and allowances pressure design wall thickness as calculated in inches 
        (millimeters) in accordance with the following equations: 

        
        t = P_i * D / (2S)

        t (float): 
        P_i (_type_): internal design gage pressure, psi (bar)
        D (_type_): outside diameter of pipe, in. (mm) 
        S (_type_): applicable allowable stress value, psi (MPa), as determined 
        by the following equation: 

        S = F X E x Sy

        F (_type_): design factor based on nominal wall thickness. Generally
                    < 0.72 or 0.80 for slurry.
        E (float): weld joint factor as specified in Table 403.2.1-1 
        S_y (_type_): specified minimum yield strength of the pipe, psi (MPa) 

    """

    # if S_y provided. Calculate the nominal wall thickness required.
    if S_y:
        S = F * E * S_y

        t = P_i * D / (2 * S)

        t_n = t + A

    
    # To calculate the stress based on a thickness t_n, set S_y to None.
    if t_n:
        S = P_i * D / (2 * (t_n - A))
    
    return t_n, S


###############################################################################
#                                                                             #
#               PIPE STRESS ENGINEERING, Peng and Peng                        #
#                                                                             #
###############################################################################

def peng_peng_eq_10_15(D, t, E, alpha, T_1, T_2, v, S_H):
    """
    10.2 Behavior of Long Pipeline. 10.2.2 Anchor Force. The first step in 
    finding the potential movement of the pipeline is to determine the 
    force required to stop the movement of the pipeline. From the total
    expansion given in Eq. (10.14) the force and stress can be derived per 
    Eq. (10.15). (See page 335-336, and example 10.5 pg 346). 

    F = AS = AEe = A{E alpha(T-2-T_1) + (0.5 - v) S_H}

    Args:
        A (_type_): cross sectional area of the pipe material
        E (_type_): Youngs modulus.
        alpha (_type_): Thermal expansion coefficient. 
        T_1 (_type_): Installation or Construction temperature.
        T_2 (_type_): Operating temperature. 
        v (_type_): Poisson's ratio.
        S_H (_type_): Pressure hoop stress based on thickness of pipe used.
    
    Returns:
        F (float): The anchor force. Total axial driving force. AKA the 
                    potential pipeline expansion force.
    """

    A = np.pi * D * t
    F = A * (E * alpha * (T_2 - T_1) + (0.5 - v) * S_H) 
    return F




def api_spec_6a_d_3(D, P, sigma, E, f, K, S):
    """
    Calculate the bolt loads and torque per API Specification 6A / ISO
    10423 Specifcation for Wellhead and Christmas Tree Equipment. 2005 edition.
    Annex D Recommended flange bolt torque. 2021 edition has moved to Annex H
    Recommended Assembly of Closure Bolting.

    Args:
        D (_type_): is the thread major diameter (mm)
        P (_type_): is the thread pith, equal to 1 / number of threads per unit 
                    length (mm)
        sigma (_type_): is the stress in the bolt (MPa)
        E (_type_): is the pitch diameter of thread (mm)
        f (_type_): is the friction coefficient, 0.07 and 0.13
        K (_type_): is the nut internal chamfer, equal to 3.175 mm
        S (_type_): is the secant of 30 degrees. 1.1547

    Call: 
    D = [12.7, 16.7, 20.84, 24.84] * ureg('mm') # Nominal diameter
    P = [1.75, 2, 2.5, 3] * ureg('mm') # Pitch of thread
    sigma = 120 * ureg('MPa') # Bolt stress 
    E = [10.8, 14.7, 18.3, 22] * ureg('mm') # Pitch diameter
    f = 0.13 # friction coefficientH = 
    K = 3.175 * ureg('mm') # Nut internal chamfer


    for k, v in enumerate(D):
        api_spec_6a_d_3(D=D[k], P=P[k], sigma=sigma, E = E[k], f=f, K=K, S=S)


    """

    A_s = np.pi / 4 * ((D - 0.974 * P))**2 # stress area (mm2)

    S = 1/cos(radians(30)) # sec(radians(30)) == 1/cos(radians(30))

    A_s_ = A_s.magnitude
    D_ = D.magnitude
    P_ = P.magnitude
    sigma_ = sigma.magnitude
    E_ = E.magnitude
    K_ = K.magnitude

    F = sigma * A_s # Force per stud
    F_ = F.magnitude

    H_ = 1.5 * D_ + 3.175 # hex size (nut)

    T = F_ * E_ * (P_ + np.pi * f * E_ * S) / (2 * (np.pi * E_ - P_ * f * S)) \
        + (F_ * f * (H_ + D_ + K_)/4 ) # Torque
    T = T * ureg('N mm')
    T.to('N m')
    print(f"{D} - Bolt load in tension: {F.to('kN'):.0f} \
        and Torque {T.to('N m'):.0f}")






###############################################################################
#                                                                             #
#                            VALIDATION CHECKS                                #
#                                                                             #
###############################################################################

def validation_checks():


    #--- Peng & Peng Example 10.5 Example Calculations of Basic Pipeline Behaviors

    # (a) Minimum wall thickness required. API-5L Grade X52 ERW. 
    S_y = 52000 # psi
    t, S = asme_b31_4_403_2_1(P_i= 1200, D= 20, F= 0.72, E = 1.0, A = 0, 
                    S=None, S_y= S_y, t_n= None) # 0.321 in, 37,440 psi

    # (b) Pressure hoop stress - Based on selected nominal thickness

    t, S_H = asme_b31_4_403_2_1(P_i= 1200, D= 20, F= 0.72, E = 1.0, A = 0, 
        S=None, S_y= None, t_n=0.344) # 0.344 in, 34,833 psi

    # (c) Longitudinal stress at the fully restrained portion.

    S_E = asme_b31_4_402_5_1(E= 29.5E6, alpha= 6.5E-6, T_1 = 50, 
                                            T_2= 170) #-23,010 psi
    # Check S_E to Table 403.3.1-1. Restrained pipe S_E <= 0.9 S_y
    S_E <= 0.9 * S_y

    S_L = asme_b31_4_402_6_1(S_E= S_E, S_H= S_H) # -12,545 psi
    # Check S_E to Table 403.3.1-1. Restrained pipe S_L <= 0.9 S_y
    S_L <= 0.9 * S_y

    # (d) Combined equivalent stress at the fully restrained portion

    S_eq = asme_b31_4_402_7(S_L, S_H) # 47,428 psi
    # Check S_eq to Table 403.3.1-1. Restrained pipe S_eq <= 0.9 S_y
    S_eq <= 0.9 * S_y

    # (e) Design wall wall thickness

    t, S_H = asme_b31_4_403_2_1(P_i= 1200, D= 20, F= 0.72, E = 1.0, A = 0, 
            S=None, S_y= None, t_n=0.375) # 0.375 in, 32,000 psi
    S_E = asme_b31_4_402_5_1(E= 29.5E6, alpha= 6.5E-6, T_1 = 50, T_2= 170)
    S_L = asme_b31_4_402_6_1(S_E= S_E, S_H= S_H) 
    S_eq = asme_b31_4_402_7(S_L, S_H) # 45,410 psi
    S_eq  <= 0.9 * S_y # B31.4 satisfied

    # (f) Anchor load. 

    F = peng_peng_eq_10_15(20, 0.375, 29.5E6, 6.5E-6, 50, 170, 0.3, S_H) # 692,900 lb


    # Metric check

    # (a) Minimum wall thickness required. API-5L Grade X52 ERW. 
    S_y_si = (52000 * u.psi).to('bar') # psi
    P_i_si = (1200 * u.psi).to('bar')
    D_si = (20 * u.inch).to('mm')
    t_si = (0.321 * u.inch).to('mm') # 8.15 mm
    S_si = (37440 * u.psi).to('bar') # 2,581 bar
    E_si = (29.5E6 * u.psi).to('bar')
    T_1_si = (50 * u.degF).to('degC')
    T_2_si = (170 * u.degF).to('degC')


    t, S = asme_b31_4_403_2_1(P_i= P_i_si, D= D_si, F= 0.72, E = 1.0, A = 0, 
                    S=None, S_y= S_y, t_n= None) # 8.15 mm, 2,581 bar

    # (b) Pressure hoop stress - Based on selected nominal thickness

    t, S_H = asme_b31_4_403_2_1(P_i= P_i_si, D= D_si, F= 0.72, E = 1.0, A = 0, 
        S=None, S_y= None, t_n=8.74 * u.mm) # 8.74 mm, 2,402 bar


    # (c) Longitudinal stress at the fully restrained portion.

    S_E = asme_b31_4_402_5_1(E= E_si, alpha= 11.7E-6, T_1 = T_1_si, 
                                            T_2= T_2_si) #-1,586 bar
    S_E = S_E.magnitude * u("bar")


    S_L = asme_b31_4_402_6_1(S_E= S_E, S_H= S_H) # -865 bar

    # (d) Combined equivalent stress at the fully restrained portion

    S_eq = asme_b31_4_402_7(S_L, S_H) # 3,270 bar 
    # Check S_eq to Table 403.3.1-1. Restrained pipe S_eq <= 0.9 S_y
    S_eq <= 0.9 * S_y_si


    # (e) Design wall wall thickness

    t, S_H = asme_b31_4_403_2_1(P_i= P_i_si, D= D_si, F= 0.72, E = 1.0, A = 0, 
            S=None, S_y= None, t_n=9.53 * u.mm) # 9.53 mm, 2,206 bar
    S_L = asme_b31_4_402_6_1(S_E= S_E, S_H= S_H)
    S_eq = asme_b31_4_402_7(S_L, S_H) # 3,131 bar
    S_eq  <= 0.9 * S_y_si # B31.4 satisfied

    # (f) Anchor load. 

    F = peng_peng_eq_10_15(D_si.magnitude, t.magnitude, E_si.magnitude, 11.7E-6, 
                T_1_si.magnitude, T_2_si.magnitude, 0.3, S_H.magnitude) # 3,082 kN
    F = (F / 10 / 1000) * u.kilonewton # work out the units with bar goes to 1/10 N
    F.to('lbfs')




    ###############################################################################

    ### Validate Lug calculation per Peng and Peng Example page 197. 


    P = 2000        # Pipe design pressure, psi
    D = 4           # pipe diameter, in
    t = 0.3         # pipe thickness, in

    t_lug = 0.5     # lug thickness, in
    e = 2           # lug height, in
    l = 4.5         # lug length, in
    b = 1           # 1 in or 25.4 mm
    W = -1000        # load applied to the lug, lb. (-) compression. (+) tension
    W_ang = 50      # load angle applied to the lug, degree

    S_a = 15000     # allowable stress, psi
    S = 12533       # pressure hoop stress due to design pressure, psi

    S = asme_b31_3_304_1_2(P=P,D = D, E= 1, W= 1, Y = 0.4, t= t)

    K, X, L_f_c, L_f_t, S_t, L_a_c, L_a_t, L_t, L_c, S_b_c, S_b_t =\
    asme_bpvc_I_pg_56(D, t_lug, t, b, S_a, S, e, l, W, W_ang) # -14645 / 4499 psi

    st_bend_comp = -14656 * ureg('psi')
    st_bend_ten = 4999. * ureg('psi')

    P = 2000*6.895E-3       # Pipe design pressure, MPa (2000 psi)
    D = 4 *25.4             # pipe diameter, mm (4 in)
    t = 0.3 *25.4           # pipe thickness, mm (0.3 in)

    t_lug = 0.5*24.4        # lug thickness, mm (0.5 in)
    e = 2*25.4              # lug height, mm (2 in)
    l = 4.5*25.4            # lug length, mm, (4.5 in)
    b = 1*25.4              # 1 in or 25.4 mm
    W = -1000*4.448         # load applied to the lug, N (-1000 lb). 
                                # (-) compression. (+) tension
    W_ang = 50              # load angle applied to the lug, degree

    S_a = 15000 * 6.896E-3  # allowable stress, MPa (15000 psi)
    S = 12533*6.896E-3      # pressure hoop stress due to design pressure, MPa 
                                # (12533 psi)

    S = asme_b31_3_304_1_2(P=P,D = D, E= 1, W= 1, Y = 0.4, t= t)

    K, X, L_f_c, L_f_t, S_t, L_a_c, L_a_t, L_t, L_c, S_b_c, S_b_t =\
    asme_bpvc_I_pg_56(D, t_lug, t, b, S_a, S, e, l, W, W_ang) # -101 / 31 MPa

    st_bend_comp.to('MPa')
    st_bend_ten.to('MPa')


    ###############################################################################




def functions():
    """Dummy functions for use in the outline in VS Code"""

###############################################################################
#                                                                             #
#                            FUNCTION CALLS                                   #
#                                                                             #
###############################################################################



###############################################################################
#                   API Bolt Load and Torque Calculation                      #
###############################################################################

def api_bolt_load_calc():

    D = [12.7, 16.7, 20.84, 24.84] * ureg('mm') # Nominal diameter
    P = [1.75, 2, 2.5, 3] * ureg('mm') # Pitch of thread
    sigma = 120 * ureg('MPa') # Bolt stress 
    E = [10.8, 14.7, 18.3, 22] * ureg('mm') # Pitch diameter
    f = 0.07 # friction coefficientH = 
    K = 3.175 * ureg('mm') # Nut internal chamfer

    for k, v in enumerate(D):
        api_spec_6a_d_3(D=D[k], P=P[k], sigma=sigma, E = E[k], f=f, K=K, S=S)




###############################################################################
#                   Lug Calculation                                           #
###############################################################################

def lug_calcs():

    P = 4.800               # Pipe design pressure, MPa 
    D = 114.3               # pipe diameter, mm 
    t = 8.564               # pipe thickness, mm 

    # H300061 PS1-SP00 Line Stop
    t_lug = 10              # lug thickness, mm 
    e = 100 / 8             # lug height, mm. This is the eccentricity.
    l = 100                 # lug length, mm
    b = 1*25.4              # 1 in or 25.4 mm
    W = 32000               # load applied to the lug, N 
                            # (-) compression. (+) tension
    W_ang = 0               # load angle applied to the lug, degree

    S_a = 138               # allowable stress, MPa A106 Gr. B
                                
    # pressure hoop stress due to design pressure, MPa 
    S = asme_b31_3_304_1_2(P = P, D = D, E= 1, W= 1, Y = 0.4, t= t)  


    K, X, L_f_c, L_f_t, S_t, L_a_c, L_a_t, L_t, L_c, S_b_c, S_b_t =\
    asme_bpvc_I_pg_56(D, t_lug, t, b, S_a, S, h, l, W, W_ang)



# ----------------------------------------------------------------

def misc_calcs():

    NPS, Di, Do, t = nearest_pipe(NPS=4, schedule='80')
    P = 4800 * u.kPa

    S_H_bar = asme_b31_4_402_3(D=Do.to('mm'), P_i = P.to('bar')  , t=t.to('mm'))
    S_H = S_H_bar.to('MPa')




    do, W_A = asme_bpvc_viii_d1_ug_131_e_2(1780, 101, do=None, W_A=800e3, solve_for="do_a")
    do, W_A = asme_bpvc_viii_d1_ug_131_e_2(1780, 101, do=80, W_A=None, solve_for="W_A")

    print(f"Actual flow: {W_A}")
    print(f"Diameter for actual flow: {do}")



    do, W_A = asme_bpvc_xiii_9_7_6_4_d_3(1.780, 0.101, rho_l=1000, W_A=800e3, solve_for="do_a")
    do, W_A = asme_bpvc_xiii_9_7_6_4_d_3(1.780, 0.101, rho_l=1000, do=80, W_A=None, solve_for="W_A")


    print(f"Actual flow: {W_A}")
    print(f"Diameter for actual flow: {do}")


# 


if __name__ == "__main__":

    def main():

        # OD32 Stainless Steel Pipe. Check maximum pressure.
        d=32
        NPS, Di, Do, t = nearest_pipe(Do=d*u('mm'), schedule='10')
        asme_b31_3_304_1_2(2, Do,S=115,t=t, solve_for='P') # 13.89 MPa
        asme_b31_3_304_1_2(13.89, Do, S=115, solve_for='t') # 2.77 mm 
        asme_b31_3_304_1_2(13.89, Do, t=t, solve_for='S') # 115 MPa


    main()


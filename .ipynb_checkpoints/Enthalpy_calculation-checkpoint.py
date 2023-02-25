import pandas as pd

import pandas as pd

def H_miedema (A: str, B: str, c_A: float, c_B: float, 
       type_of_compound: int, miedema_coifficients_path=None)->float:
    '''
    Takes the parameters of a chemical compound and returns H, i.e. the enthalpy of the miedema
    
    A, B: name of elements
    c_A, c_B: atomic fractions
    type_of_compound: 1-disordered solid solution, 2-ordered compound, 3-amorphous phase
    miedema_coefficients_path: path to file with corresponding coefficients
    
    Example:
    >> H('Ca', 'Co', 1, 3, 2, 'miedema_coefficients.xlsx') 
    >> 11.999574005596852
    '''
    miedema_coeff_df = pd.read_excel(miedema_coifficients_path, index_col=0)
    
    names_of_params = ['electronegativity', 'discontinuity', 'volume', 'transitivity', 'radius']
    
    QP = 9.4
    
    RP = 0

    phi_A, n_A, V_A, T_A, r_A = miedema_coeff_df.loc[A, names_of_params]
    phi_B, n_B, V_B, T_B, r_B = miedema_coeff_df.loc[B, names_of_params]

    delta_phi = phi_A - phi_B
    delta_n = n_A - n_B
    
    # surface concentration
    c_A_s = round(c_A * V_A / (c_A * V_A + c_B * V_B), 3)
    c_B_s = round(c_B * V_B / (c_A * V_A + c_B * V_B), 3)
    
    match type_of_compound:
        case 1: f = c_A_s * c_B_s
        case 2: f = c_A_s * c_B_s * (1 + 5 * (c_A_s * c_B_s)**2)
        case _: f = c_A_s * c_B_s * (1 + 8 * (c_A_s * c_B_s)**2)
    
    match T_A, T_B:
        case (1, 1): P = 14.1
        case (1, 0) | (0, 1): P = 12.3
        case _: P = 10.6
    
    H = 2 * f * (c_A * V_A + c_B * V_B) / (1 / n_A + 1 / n_B) * (-P * delta_phi**2 + QP * P * delta_n**2)
    Sc = 1 - c_B_s * (V_A - V_B) / (c_A_s * V_A + c_B_s * V_B)
    H *= Sc
    return H, r_A, r_B

def L_J_potential(r, U0, b):
    return 4*U0*((b/r)**12 - (b/r)**6)

def H(A, B, c_A, c_B, 
      type_of_compound, miedema_coifficients_path, 
      multiplicity_A, multiplicity_B, 
      occupancy_A, occupancy_B,
      distance_AB, cnt=1, mode='L-J'
     ):
    N_avog = 6.02 * 10**23 
    match mode:
        case 'L-J': potential=L_J_potential
        case _: raise ValueError
    
    H_m, r_A, r_B = H_miedema(A, B, c_A, c_B, type_of_compound, miedema_coifficients_path)
    
    H_m /= ( 96.5 * cnt * occupancy_A * occupancy_B)
    
    
    enthalpy = -96.5 / 2 * N_avog * (cnt  #26 -> число атомов в ячейке
                           * occupancy_A * occupancy_B 
                           *  H_m *  potential(distance_AB / (r_A + r_B), 1, 1))
    return enthalpy
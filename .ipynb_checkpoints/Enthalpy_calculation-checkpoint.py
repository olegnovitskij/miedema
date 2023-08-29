import numpy as np
import multiprocessing as mp
import pandas as pd
from functools import partial
import itertools

def unfold(l):
    return [item for sublist in l for item in sublist]

def H_miedema (A: str, B: str, 
       type_of_compound: int, num: int=51, 
               mode: str='minmax',
               radius_type: str='atom',
               miedema_coifficients_path='miedema_coefficients.xlsx', **kwargs)->float:
    '''
    Takes the parameters of a chemical compound and returns H, i.e. the enthalpy of the miedema
    
    A, B: name of elements
    type_of_compound: 1-disordered solid solution, 2-ordered compound, 3-amorphous phase
    miedema_coefficients_path: path to file with corresponding coefficients
    
    Example:
    >> H('Ca', 'Co', 1, 3, 2, 'miedema_coefficients.xlsx') 
    >> 11.999574005596852
    '''
    
    miedema_coeff_df = pd.read_excel(miedema_coifficients_path, index_col=0)
    match radius_type:
        case 'ion': name_radius_column = 'radius_ion'
        case 'atom' : name_radius_column = 'radius_atom'
        
    names_of_params = ['electronegativity', 'discontinuity', 'volume', 'transitivity', name_radius_column, 'R']
    
    QP = 9.4
    


    phi_A, n_A, V_A, T_A, r_A, R_A = miedema_coeff_df.loc[A, names_of_params]
    phi_B, n_B, V_B, T_B, r_B, R_B = miedema_coeff_df.loc[B, names_of_params]

    RP = 0
    match T_A, T_B:
        case (1, 0) | (0, 1): RP = R_A * R_B
    
    
    delta_phi = phi_A - phi_B
    delta_n = n_A - n_B
    
    match T_A, T_B:
        case (1, 1): P = 14.1
        case (1, 0) | (0, 1): P = 12.3
        case _: P = 10.6
    # surface concentration
    H_ar = []
    for c_A in np.linspace(0, 1, num=num):
        c_B = 1 - c_A    
        c_A_s = round(c_A * V_A / (c_A * V_A + c_B * V_B), 3)
        c_B_s = round(c_B * V_B / (c_A * V_A + c_B * V_B), 3)

        match type_of_compound:
            case 1: f = c_A_s * c_B_s
            case 2: f = c_A_s * c_B_s * (1 + 8 * (c_A_s * c_B_s)**2)
            case _: f = c_A_s * c_B_s * (1 + 5 * (c_A_s * c_B_s)**2)

        H = 2 * P * f * (c_A * V_A + c_B * V_B) / (1 / n_A + 1 / n_B) * (-delta_phi**2 + QP * delta_n**2 - RP)
        # Sc = 1 - c_B_s * (V_A - V_B) / (c_A_s * V_A + c_B_s * V_B)
        # H *= Sc
        H_ar.append(H)
        
    if abs(min(H_ar)) > abs(max(H_ar)):
        global_optim = min(H_ar)
    else:
        global_optim = max(H_ar)
    match mode:
        case 'minmax': res = global_optim
        case 'mean': res = np.mean(H_ar)
    return res, r_A, r_B

def Morse_potential(x, n=12, m=6):
    '''
    U0 ~ H_miedema
    
    '''
    gamma = (n * m / 2)**0.5
    return (np.exp(-2 * gamma * (x-1)) - 2 * np.exp(-gamma *(x-1)))

def L_J_potential(x, n=12, m=6):
    '''
    U0 ~ H_miedema
    
    '''
    return ((m / (n - m)) * (x**(-n) - (n / m) * x**(-m)))


def make_calculation_for_one_row(dist_arr_row, all_params_combinations=None, type_compound=2, is_abs_potential=False):
    A_ar, site_A, B_ar, site_B, dist, occ_A_ar, occ_B_ar, cnt = dist_arr_row
    all_calculations = []
    for i, A in enumerate(A_ar):
        for j, B in enumerate(B_ar):
            cur_calculation = [A, site_A, B, site_B, occ_A_ar[i], occ_B_ar[j], dist, cnt]

            for mode, potential_type, radius_type in all_params_combinations:

                U0, r_A, r_B = H_miedema(A, B, type_compound, mode=mode, radius_type=radius_type)

                match potential_type:
                    case 'L-J': potential_func = L_J_potential
                    case 'Morse': potential_func = Morse_potential

                if is_abs_potential:               
                    potential = abs(potential_func(dist / (r_A + r_B)))
                else:
                    potential = potential_func(dist / (r_A + r_B))

                E = cnt * occ_A_ar[i] * occ_B_ar[j] * abs(U0) * potential

                if site_A == site_B:
                    E /= 2

                cur_calculation.extend([U0, r_A, r_B, potential, E])
            all_calculations.append(cur_calculation)

    return all_calculations
        
def calculate_enthalpy(dist_arr: np.array, 
                       type_compound: int=2,
                       modes=['mean', 'minmax'],
                       potential_types=['L-J', 'Morse'],
                       radius_types=['atom', 'ion'],
                       n_jobs=1,
                       is_abs_potential=False,
                       output_file='calculation.xlsx'):
    
    all_params_combinations = np.array(np.meshgrid(modes, potential_types, radius_types)).T.reshape(-1,3)
    names_of_all_params_combinations = list(map(lambda x: '_'.join(x), all_params_combinations))
    name_of_values = ['Hm', 'r_A', 'r_B', 'potential', 'E']
    name_of_additional_columns = [value + '_' + suffix for suffix in names_of_all_params_combinations for value in name_of_values]
    
    if n_jobs == -1:
        num_processes = mp.cpu_count() - 1
    else:
        num_processes = min(mp.cpu_count() - 1, n_jobs)
        
    all_calculations = []
    kwargs_for_mp = {'all_params_combinations': all_params_combinations,
                    'type_compound': type_compound,
                    'is_abs_potential': is_abs_potential}
    
    with mp.Pool(processes=num_processes) as pool:
        all_calculations_folded = pool.map(partial(make_calculation_for_one_row, **kwargs_for_mp), dist_arr)
    
    
    all_calculations = unfold(all_calculations_folded)
                    
    calculation_result = pd.DataFrame(all_calculations, columns=['A', 'site_A', 'B', 'site_B', 'occ_A', 'occ_B', 'dist', 'cnt'] + name_of_additional_columns)
    
    calculation_result.to_excel(output_file)

def convert_to_float(frac_str):
    frac_str.replace(',', '.')
    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac
    

import numpy as np
from scipy.stats import t
from DataExtractor import extract_energies_and_n_values
from ExplicativeVariableMatrixBuilder import build_explicative_variable_matrix

file_path = 'HFCO_exp.dat'
 
def run_adjustment(is_constant_term, is_cubic_term):

    energies, n_values = extract_energies_and_n_values(file_path)
    energies_input = energies.reshape(-1,1)

    X, parameters_name =  build_explicative_variable_matrix(n_values, is_cubic_term, is_constant_term)
    
    params_array = np.array([str(key) + " " + str(value) for key, value in parameters_name.items()])

    option_name = "_is_constant_term_" + str(is_constant_term) + "_is_cubic_term_" + str(is_cubic_term)
    filename = f"parameters_name{option_name}.dat"
    np.savetxt(filename , params_array, fmt="%s")

    X_T = X.T
    XtX = X_T @ X
    XtX_inv = np.linalg.inv(XtX)
    
    U, S, Vt = np.linalg.svd(XtX, full_matrices=False)
        
    SVD_XtX_inv = Vt.T @ np.diag(1/S) @ U.T

    print("Ecart maximum entre XtX_inv obtenue avec np.linalg.inv et np.linalg.svd :", np.max(SVD_XtX_inv - XtX_inv))
    np.savetxt("singular_values" + option_name  + ".dat", S, fmt='%.15g')

    parameter_estimated = XtX_inv @ X_T @ energies_input 
    energies_estimated = X @ parameter_estimated
    energies_error = energies_input - energies_estimated

    parameter_estimated_with_SVD = SVD_XtX_inv @ X_T @ energies_input 
    energies_estimated_with_SVD = X @ parameter_estimated_with_SVD

   
    np.savetxt("parameter_estimated" + option_name  + ".dat", parameter_estimated, fmt='%.15g')
    np.savetxt("energies_error" + option_name  + ".dat", energies_error, fmt='%.15g')
    np.savetxt("parameter_SVD_effect" + option_name  + ".dat", abs(parameter_estimated - parameter_estimated_with_SVD), fmt='%.15g')
    np.savetxt("energies_SVD_effect" + option_name  + ".dat", abs(energies_estimated  - energies_estimated_with_SVD), fmt='%.15g')
    
    print("SVD effet parametre max difference:", np.max(abs(parameter_estimated - parameter_estimated_with_SVD)))
    print("SVD effet energie max difference:", np.max(abs(energies_estimated  - energies_estimated_with_SVD)))
    

    rmsd = np.std(energies_estimated - energies_input)
    print("rmsd wihthout SVD :", rmsd )    
    rmsd_with_SVD = np.std(energies_estimated_with_SVD - energies_input)
    print("rmsd with SVD :", rmsd_with_SVD)

    N = len(energies_estimated)
    M = XtX_inv.shape[0] 
    sigma_estimator = np.sqrt(np.sum((energies_estimated - energies_input)**2) / (N-M))

    print(sigma_estimator)

    sqrt_XtX_inv_diag = np.sqrt(XtX_inv.diagonal())
    covariance_matrix_norm = np.zeros((M, M))  
    for k in range(M):
        for l in range(M):
            covariance_matrix_norm[k, l] = XtX_inv[k, l] / (sqrt_XtX_inv_diag[k] * sqrt_XtX_inv_diag[l])


    np.savetxt("matrice_covariance" + option_name + ".dat", covariance_matrix_norm, fmt="%.15f", delimiter=" ")
    
    threshold = 0.8
    indices = np.argwhere((np.abs(covariance_matrix_norm) > threshold) & (np.abs(covariance_matrix_norm) < 1.0 - 1e-15) )
    for i, j in indices:
        if i < j :
            print(f"Indices: ({i+1}, {j+1}), Valeur: {covariance_matrix_norm[i, j]}")

    sqrt_XtX_inv_diag_column = sqrt_XtX_inv_diag.reshape(-1,1)

    alpha = 0.01
    t_student = t.ppf(1 - alpha / 2, df= N - M)
    ic_parameter_estimated = t_student * sigma_estimator *  sqrt_XtX_inv_diag_column
    
    np.savetxt("ic_parameter_estimated" + option_name  +".dat", ic_parameter_estimated, fmt='%.15g')

    parameter_bump_up = parameter_estimated

    sigma_with_bump_up = np.zeros(M)
    sigma_with_bump_down = np.zeros(M)

    for k in range(M):
        parameter_bump_up = parameter_estimated.copy()  
        parameter_bump_up[k] *= 1.1
        energies_estimated_with_bump_up = X @ parameter_bump_up.reshape(-1,1)
        sigma_with_bump_up[k] = np.sqrt(np.sum((energies_estimated_with_bump_up - energies_input )**2) / N)

    np.savetxt("sigma_with_bump_up" + option_name  +".dat", sigma_with_bump_up.reshape(-1,1), fmt='%.15g')

    for k in range(M):
        parameter_bump_down = parameter_estimated.copy()  
        parameter_bump_down[k] *= 0.9
        energies_estimated_with_bump_down = X @ parameter_bump_down.reshape(-1,1)
        sigma_with_bump_down[k] = np.sqrt(np.sum((energies_estimated_with_bump_down - energies_input )**2) / N)

    np.savetxt("sigma_with_bump_down" + option_name  +".dat", sigma_with_bump_down.reshape(-1,1), fmt='%.15g')

    mean = energies_error.mean()
    print("Erreur moyenne :" +  str(energies_error.mean()))

    lower_bound = mean - rmsd 
    upper_bound = mean + rmsd 

    count_in_range = np.sum((energies_error >= lower_bound) & (energies_error <= upper_bound))
    percentage_in_range = (count_in_range / len(energies_error)) * 100

    print("Pourcentage à +/- rmsd :" + str(percentage_in_range))

    lower_bound = mean - sigma_estimator
    upper_bound = mean + sigma_estimator

    count_in_range = np.sum((energies_error >= lower_bound) & (energies_error <= upper_bound))
    percentage_in_range = (count_in_range / len(energies_error)) * 100

    print("Pourcentage à +/- sigma_estimator :" + str(percentage_in_range))

    # Paramètres
    index_start = 7
    threshold = 0.4

    # Extraire la sous-matrice
    sub_matrix = covariance_matrix_norm[index_start:, index_start:]


    diagonal = np.eye(sub_matrix.shape[0], dtype=bool)
    below_threshold = (abs(sub_matrix) < threshold) & ~diagonal
    percentage = np.sum(below_threshold) / (sub_matrix.size - sub_matrix.shape[0]) * 100

    print(f"Pourcentage des termes inférieurs à {threshold}  : {percentage:.2f}%")

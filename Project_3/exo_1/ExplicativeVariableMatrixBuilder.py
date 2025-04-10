import numpy as np

def build_explicative_variable_matrix(n_values, is_cubic_term, is_constant_term):

    zero_singular_value_threshold = 1e-13
    n_samples, n_modes = n_values.shape
    
    linear_terms = n_values + 0.5
    parameters_name = {}
    parameter_counter = 1
    if is_constant_term == True:
        X = np.ones((n_samples, 1))
        X = np.hstack((X, linear_terms))
        parameters_name[parameter_counter] = "T_0" 
        parameter_counter += 1 
    
    else :
        X = linear_terms
    
    for alpha in range(n_modes):
        parameters_name[parameter_counter] = "omega_" + str(alpha + 1)
        parameter_counter += 1 


    # Ajout des termes x_alpha,beta 
    for alpha in range(n_modes):
        for beta in range(alpha, n_modes):

            is_quadratic_term = check_is_quadratic_term(n_values, alpha, beta)
            
            if (is_quadratic_term) :
                quadratic_term = (linear_terms[:, alpha] * linear_terms[:, beta]).reshape(-1, 1)
                X_new = np.hstack((X, quadratic_term))
                _, S, _ = np.linalg.svd(X_new, full_matrices=False)

                if np.min(S) > zero_singular_value_threshold:
                    X = np.hstack((X, quadratic_term))
                    parameters_name[parameter_counter] = "x_" + str(alpha + 1) + str(beta + 1)
                    parameter_counter += 1

    if is_cubic_term == False:
        return X, parameters_name

    # Ajout des termes y_alpha, beta, gamma 
    for alpha in range(n_modes):
        for beta in range(alpha, n_modes):
            for gamma in range(beta, n_modes):
                
                is_cubic_term = check_is_cubic_term(n_values, alpha, beta, gamma)
                
                if (is_cubic_term):
                    cubic_term = (linear_terms[:, alpha] * linear_terms[:, beta] * linear_terms[:, gamma]).reshape(-1, 1)
                    X_new = np.hstack((X, cubic_term))
                    _, S, _ = np.linalg.svd(X_new, full_matrices=False)
                  
                    if np.min(S) > zero_singular_value_threshold:
                        X = np.hstack((X, cubic_term))
                        parameters_name[parameter_counter] = "y_" + str(alpha + 1) + str(beta + 1) + str(gamma + 1)
                        parameter_counter += 1 
                   
    return X, parameters_name


def check_is_quadratic_term(n_values, alpha, beta):
    for n_value in n_values:
        if (n_value[alpha] != 0) and (n_value[beta] !=0):
            return True
    return False

def check_is_cubic_term(n_values, alpha, beta, gamma):
    for n_value in n_values:
        if (n_value[alpha] != 0) and (n_value[beta] !=0) and  (n_value[gamma] !=0):
            return True
    return False
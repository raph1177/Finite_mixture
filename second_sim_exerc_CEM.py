# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#### Second simulation exercise ####

### Bias-reduced estimation of finite mixtures ###

### By Raphaël Langevin ###
#%%

import numpy as np
import scipy
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.rcParams['figure.figsize'] = [14, 7]
import os
os.chdir('C:\\Users\\raph1\\.spyder-py3\\Work') ### Directory to save the results ###

#%%
##### Generate simulated assignment data for K components with non-constant time average #####
def generate_normal_data(seed, K, alpha_K, range_beta, range_cov, N_ind, N_period, N_covariates, random_var, AR_process,sigma2_X):
    #### Creation of the state probabilities and categorical sequence #####
    np.random.seed(seed)
    alpha = [alpha_K/K] * K
    Initial_state_prob = np.random.RandomState(seed).dirichlet(alpha,size=N_ind)
    Transition_matrix = np.random.dirichlet(alpha,size=K)     
    State_prob = np.zeros(shape=(K,N_period,N_ind),dtype=float, order='C')
    assignment_mat = np.zeros(shape=(N_ind,N_period,K), dtype=int, order='C')
    for t in range(N_period):
        if AR_process == 1: 
            State_prob[0:,t] = np.transpose(Initial_state_prob @ np.linalg.matrix_power(Transition_matrix,t))
        else: ## No AR(1) process
            State_prob[0:,t] = np.transpose(np.random.RandomState(seed).dirichlet(alpha,size=N_ind))
        for i in range(N_ind):
            assignment_mat[i,t] = np.random.RandomState((t*N_ind+i)+seed).multinomial(n=1, pvals= State_prob[0:,t,i])   
    
    ### Creation of the matrix of covariates ####
    sigma_mat1 = np.identity(n=N_covariates) 
    mean_covariates = np.zeros(shape=(K,N_covariates), dtype=float, order='C')
    mean_cov = np.zeros(shape=(K,N_covariates), dtype=float, order='C')
    if K==1: mean_cov=[0]
    else:
        for j in range(K):
            if (range_cov!=0): mean_cov[j] = N_covariates*[np.arange(-range_cov,range_cov+(2*range_cov)/(K-1), (2*range_cov)/(K-1))[j]]
            mean_covariates[j,0:] = np.transpose(scipy.stats.multivariate_normal.rvs(mean=mean_cov[j], cov=sigma_mat1, size=1, random_state=j*(seed+1)))
    mat_covariates = np.zeros(shape =(K,N_ind,N_period,N_covariates), dtype=float, order='C')
    sigma_mat_K =  np.zeros(shape=(K,N_covariates,N_covariates), dtype=float, order='C')  
    for j in range(K):
        cov1 = np.diag(N_covariates*[sigma2_X])
        for h in range(N_covariates):
            for l in range(h):
                cov1[h,l] = np.random.normal(loc=0,scale=1,size=1)[0]       
        sigma_mat_K[j] = cov1@np.transpose(cov1)
        for i in range(N_ind):
            mat_covariates[j,i] = np.random.RandomState((j*N_ind+i)+seed).multivariate_normal(mean= mean_covariates[j], cov=sigma_mat_K[j],size=N_period)
    
    ### Selection of covariates and creation of time/cross-sectional averages ###
    X_obs = np.zeros(shape =(N_ind,N_period,N_covariates), dtype=float, order='C')
    X_unit_mean_assign = np.zeros(shape =(N_period,K), dtype=float, order='C')
    for t in range(N_period):
        for i in range(N_ind):
            X_obs[i,t] =  mat_covariates[0:,i,t].T@assignment_mat[i,t]      
        for g in range(K):
            X_unit_mean_assign[t,g] = np.average(mat_covariates[g,0:,t,0], weights= assignment_mat[0:,t,g], axis=0)           
    X_time_mean = np.average(X_obs, axis=1) 

    #### Compute real pi_k of each component ####
    pi_star = np.sum(assignment_mat, axis=(0,1))/(N_period*N_ind)
    ### Pi_k for each period ###
    pi_tstar = np.sum(assignment_mat, axis=0)/N_ind

    #### Generation of the vector of parameters ######
    if range_beta==0: mean_beta=[0]*K
    else: mean_beta = np.arange(-range_beta,range_beta+(2*range_beta)/(K-1), (2*range_beta)/(K-1))
    Beta_star = np.zeros(shape=(2,K), dtype=float, order='C')
    for i in range(K):
        Beta_star[:2,i] = np.random.RandomState(seed+i).normal(size=2, loc=mean_beta[i],scale=1)
    Varcov_random_eff = np.diagflat(np.arange(start=1,stop=K+1,step=1,dtype=float))
    ID_random_eff = scipy.stats.multivariate_normal.rvs(mean=None , cov=Varcov_random_eff, size=N_ind, random_state=seed)
    if K==1:
        ID_random_eff = np.reshape(ID_random_eff, shape=(N_ind,K))
    Random_eff = np.random.RandomState(seed).normal(size=N_period*N_ind, loc=0, scale=np.sqrt(random_var))
    Random_eff = np.reshape(Random_eff, (N_ind,N_period))
    
    ### Generation of the outcome ####
    Y_obs = np.zeros(shape=(N_ind,N_period), dtype=float, order='C')
    Y_obs_assign = np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C') 
    Time_fixed_eff = np.zeros(shape=(N_period,K), dtype=float, order='C')
    X_obs_time_mean = np.zeros(shape=(N_ind,N_period,N_covariates), dtype=float, order='C')
    
    for t in range(N_period):
        Time_fixed_eff[t] = np.random.RandomState(seed+t).normal(size=K, loc=X_unit_mean_assign[t] ,scale=1) 
        X_obs_time_mean[0:,t] = X_time_mean[0:]
        for i in range(N_ind):
            Y_obs[i,t] = X_obs[i,t,0] * (Beta_star[0]@assignment_mat[i,t]) + X_obs_time_mean[i,t,0]*(Beta_star[1] @ assignment_mat[i,t])+ ID_random_eff[i] @ assignment_mat[i,t] + Time_fixed_eff[t] @ assignment_mat[i,t] + Random_eff[i,t]
            
    return Y_obs, X_obs, Beta_star, assignment_mat, pi_star, ID_random_eff, Time_fixed_eff, Varcov_random_eff, Transition_matrix, random_var, N_ind, N_period, N_covariates, mean_covariates, Y_obs_assign, X_obs_time_mean, pi_tstar, sigma_mat_K

#%%
#### Classification error rate at the true values ####
def optimal_err_class(K_real,mahalanobis,dens_min = 1e-320):
    K=K_real
    mean_covariates = generated_data[13]
    sigma_mat_K = generated_data[17]
    Beta_star = generated_data[2]
    assignment_mat = generated_data[3]
    Time_fixed_eff = generated_data[6]
    X_obs_time_mean = generated_data[15]
    Varcov_random_eff = generated_data[7]
    Y_obs = generated_data[0]
    pi_star = generated_data[4]
    X_obs = generated_data[1]
    Time_fixed_eff_mat = N_ind*[np.diag(np.ones(shape = N_period, dtype=float))]
    Time_fixed_eff_mat = np.reshape(Time_fixed_eff_mat, shape=(N_ind,N_period, N_period))
    X_obs_total = np.c_[X_obs, X_obs_time_mean, Time_fixed_eff_mat] 
    X_model =  np.c_[X_obs_total[0:,0:,0:2], X_obs_total[0:,0:,np.shape(X_obs_total)[2]-N_period:np.shape(X_obs_total)[2]]]
    Beta_star_total = np.zeros(shape=(K,2+N_period), dtype=float, order='C')
    Beta_star_total = np.transpose(np.vstack((Beta_star[0],Beta_star[1],Time_fixed_eff)))
 
    ### Distance computation ###
    joint_density, density_y = np.zeros(shape=(N_ind,N_period,K), dtype=float, order='C'), np.zeros(shape=(N_ind,N_period,K), dtype=float, order='C')
    res_iter, res_X_iter = np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C'), np.zeros(shape=(K,N_ind,N_period,N_covariates), dtype=float, order='C'),    
    for j in range(K):
        for i in range(N_ind):
            res_iter[j,i] = X_model[i] @ Beta_star_total[j] - Y_obs[i]
            res_X_iter[j,i] = X_obs[i] - mean_covariates[j]
            density_y[i,0:,j] = scipy.stats.norm.pdf(x=res_iter[j,i], loc=0, scale=np.sqrt(random_var+Varcov_random_eff[j,j]))
            density_y[i,np.where(density_y[i,0:,j]<dens_min),j]=dens_min
            if (mahalanobis==1):
                joint_density[i,0:,j] = -np.diag(res_X_iter[j,i]@np.linalg.inv(sigma_mat_K[j])@res_X_iter[j,i].T)
            else:
                joint_density[i,0:,j] = scipy.stats.multivariate_normal.pdf(x=res_X_iter[j,i], mean=None, cov=sigma_mat_K[j])*density_y[i,0:,j]    
        
    ## Classification error at the true values, with covariates ##
    resp_max = np.argmax(joint_density, axis=2)
    resp_class = np.zeros(shape=(N_ind, N_period, K), dtype=int, order='C')
    class_error = np.zeros(shape=(N_ind, N_period, K), dtype=int, order='C')
    for i in range(N_ind):
        for t in range(N_period):
            resp_class[i,t,resp_max[i,t]] = 1
    class_error = np.abs(resp_class - assignment_mat)/2
    error_class_with_cov = np.sum(class_error)/(N_ind*N_period)

    ### Log-likelihood at the true parameter values ###
    true_log_lik = 0
    true_MC_log_lik = 0

    for i in range(N_ind):
        true_log_lik += np.sum(np.log(np.sum(pi_star*density_y[i], axis=1)))
        if (mahalanobis==0):
            true_MC_log_lik += np.sum(np.log(np.sum(assignment_mat[i]*joint_density[i], axis=1)))
        else:
            true_MC_log_lik += np.sum(np.log(np.sum(assignment_mat[i]*density_y[i], axis=1)))
    return error_class_with_cov, true_log_lik, true_MC_log_lik

#%% 
##### Generate initial values ##### 
def Initial_values(seed, K, var_init, dens_min, alpha_K):
    np.random.seed(seed)
    N_ind  = generated_data[10]
    N_period = generated_data[11]
    N_covariates = generated_data[12]
    X_obs = generated_data[1]
    X_obs_time_mean = generated_data[15]
    
    alpha = [alpha_K/K] * K 
    State_prob = np.zeros(shape=(K,N_period,N_ind),dtype=float, order='C')
    assignment_mat_init = np.zeros(shape=(K,N_ind,N_period), dtype=int, order='C')
    for i in range(N_ind):
        for t in range(N_period):
            State_prob[0:,t,i] = np.random.RandomState(seed).dirichlet(alpha,size=1)
            assignment_mat_init[0:,i,t] = np.random.RandomState((t*N_ind+i)+seed).multinomial(n=1, pvals= State_prob[:K,t,i])     
    
    Initial_beta = np.zeros(shape=(K, 2+N_period), dtype=float, order='C')
    Initial_mu = np.zeros(shape=(K, N_covariates), dtype=float, order='C')
    Time_fixed_eff_mat = N_ind*[np.diag(np.ones(shape = N_period, dtype=float))]
    Time_fixed_eff_mat = np.reshape(Time_fixed_eff_mat, shape=(N_ind,N_period, N_period))
    X_model = np.c_[X_obs[0:,0:,0:1], X_obs_time_mean[0:,0:,0:1], Time_fixed_eff_mat]

    resp_init = assignment_mat_init 
    Initial_precision = K*[np.linalg.inv(np.diagflat([var_init]*(N_period)))]
    return X_model, resp_init, Initial_beta, Initial_mu, Initial_precision, Time_fixed_eff_mat
#%%
##### EM with covariates density ######
def EM_without_cov_dens(K, Max_iter, tolerance, true_init=0, var_max=1000, var_min=1e-8, dens_min=1e-300):
    N_ind  = generated_data[10]
    N_period = generated_data[11]
    log_lik_total = list()
    list.append(log_lik_total,0)
    var_cov_y = np.linalg.inv(Init_values[4])
    precision_matrix_y = Init_values[4]
    Y_obs = generated_data[0]
    if (true_init==1): resp_iter1 = np.reshape(generated_data[3], shape=(K,N_ind, N_period))
    else: resp_iter1 = Init_values[1]
    X_model = Init_values[0]
    pi_g = K*[1/K]

    ### Iteration for the whole estimation procedure with the EM algorithm ###
    for l in range(Max_iter):
        ### GLS computing (M-STEP) ###
        S_xx = np.zeros(shape =(K, 2+N_period, 2+N_period), dtype=float, order='C')
        S_xy = np.zeros(shape =(K, 2+N_period), dtype=float, order='C')
        X_obs_tilde_g, Y_obs_tilde_g = np.zeros(shape =(N_ind, N_period, 2+N_period), dtype=float, order='C'),np.zeros(shape =(N_ind, N_period), dtype=float, order='C')
        beta_iter = np.zeros(shape =(K, 2+N_period), dtype=float, order='C')
        PX_g, PY_g = np.zeros(shape =(N_ind, N_period, 2+N_period), dtype=float, order='C'),  np.zeros(shape =(N_ind, N_period), dtype=float, order='C')
        for g in range(K):
            Y_obs_tilde_g = Y_obs*resp_iter1[g]
            X_obs_tilde_g = np.transpose(np.transpose(resp_iter1[g])*np.transpose(X_model))            
            PX_g = precision_matrix_y[g] @ X_obs_tilde_g
            PY_g = (precision_matrix_y[g] @ Y_obs_tilde_g.T).T           
            for k in range(2+N_period):
                S_xy[g,k] = np.sum(np.diag((X_obs_tilde_g.T @ PY_g)[k]))
            for i in range(N_ind):
                S_xx[g] += X_obs_tilde_g[i].T@PX_g[i]
            beta_iter[g] = np.linalg.inv(S_xx[g]) @ S_xy[g]
        
        #### E-STEP #####
        res_iter =  np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C')
        resp_iter1 = np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C')
        density_y = np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C')
        for j in range(K):
            res_iter[j] = X_model @ beta_iter[j] - Y_obs
            density_y[j] = pi_g[j]*scipy.stats.norm.pdf(res_iter[j], loc=0, scale=np.sqrt(np.diag(var_cov_y[j])))
            arg = np.where(density_y[j]<dens_min)
            density_y[j,arg[0],arg[1]] = dens_min
        resp_iter1 = density_y/np.sum(density_y, axis=0)
        sum_NT = np.sum(resp_iter1, axis=(1,2))
        
        ### Covariance matrix updates ###
        var_cov_y = np.zeros(shape=(K, N_period, N_period), dtype=float, order='C')
        var_total, var_unit, var_random = np.zeros(shape=(K), dtype=float, order='C'), np.zeros(shape=(K), dtype=float, order='C'), np.zeros(shape=(K), dtype=float, order='C')
        precision_matrix_y = np.zeros(shape=(K, N_period, N_period), dtype=float, order='C')
        for j in range(K):
            var_total[j] = np.sum(resp_iter1[j]*res_iter[j]**2)/(sum_NT[j]-2-N_period)
            mean_unit_error = np.average(res_iter[j], weights=resp_iter1[j], axis=1)
            time_weights = np.sum(resp_iter1[j], axis=1)      
            var_unit[j] = np.average((mean_unit_error - np.average(mean_unit_error, weights = time_weights))**2, weights = time_weights)
            var_random[j]  =  var_total[j] - var_unit[j]
            var_cov_y[j] = np.diag(N_period*[var_random[j]]) + var_unit[j]
            eigen_cov_y = np.linalg.eigh(var_cov_y[j])
            eigen_cov_y[0][np.where(np.abs(eigen_cov_y[0])< var_min)] = var_min
            var_cov_y[j] = eigen_cov_y[1] @ np.diag(eigen_cov_y[0]) @ np.transpose(eigen_cov_y[1]) 
            precision_matrix_y[j] = np.linalg.inv(var_cov_y[j])

        ### Log-likelihood computation ###
        pi_g = np.zeros(K, dtype=float, order='C')
        penalty = np.zeros(K, dtype=float, order='C')
        for j in range(K):
            pi_g[j] = np.sum(resp_iter1[j])/(N_ind*N_period)
            penalty[j] = -(var_random[j]+var_unit[j])**(-1)-np.log(var_random[j]+var_unit[j])
        log_lik = np.sum(np.log(np.sum(density_y, axis=0)))+(1/np.sqrt(N_ind*N_period))*np.sum(penalty)
        list.append(log_lik_total, log_lik)
        
        ### Break mechanism ###
        if (l==0): log_lik_total = log_lik_total[1:]
        elif ((1-(log_lik_total[l]/log_lik_total[l-1]))< tolerance): break
    
        ## Save last results ##
        resp_iter_last = resp_iter1
        beta_iter_last = beta_iter
        var_total_last = var_total
        pi_g_last = pi_g
        var_unit_last = var_unit
        var_cov_y_last = var_cov_y
    
    return log_lik_total, l, resp_iter_last, beta_iter_last, var_total_last, pi_g_last, var_unit_last, var_cov_y_last

#%%
##### CEM with covariates density ######
def CEM_with_cov_dens(K, Max_iter, mahalanobis, var_max=1000, var_min=1e-8, var_min_X=1e-8, dens_min=1e-300, var_cov_x_init=0):
    N_ind  = generated_data[10]
    N_period = generated_data[11]
    N_covariates = generated_data[12]
    log_lik_total = list()
    sum_changes = list()
    list.append(log_lik_total,0)
    var_cov_X = np.array(K*[var_init*np.identity(N_covariates)])
    var_cov_y = np.linalg.inv(Init_values[4])
    precision_matrix_y = Init_values[4]
    Y_obs = generated_data[0]
    X_obs = generated_data[1]
    resp_iter = Init_values[1]
    Mu_tilde = Init_values[3]
    X_model = Init_values[0]
    resp_class_last = 0
    changes = 0
    beta_iter_last = np.zeros(shape =(K, 2+N_period), dtype=float, order='C')
    ### First set of covariate mean ###
    for j in range(K):
        for p in range(N_covariates):
            Mu_tilde[j,p] = np.sum(np.diag(np.transpose(X_obs[0:,0:,p]) @ resp_iter[j]))/np.sum(resp_iter[j])

    for l in range(Max_iter):
        ### GLS computing (M-STEP) ###
        S_xx = np.zeros(shape =(K, 2+N_period, 2+N_period), dtype=float, order='C')
        S_xy = np.zeros(shape =(K, 2+N_period), dtype=float, order='C')
        X_obs_tilde_g, Y_obs_tilde_g = np.zeros(shape =(N_ind, N_period, 2+N_period), dtype=float, order='C'),np.zeros(shape =(N_ind, N_period), dtype=float, order='C')
        beta_iter = np.zeros(shape =(K, 2+N_period), dtype=float, order='C')
        PX_g, PY_g = np.zeros(shape =(N_ind, N_period, 2+N_period), dtype=float, order='C'),  np.zeros(shape =(N_ind, N_period), dtype=float, order='C')
        for g in range(K):
            Y_obs_tilde_g = Y_obs*resp_iter[g]
            X_obs_tilde_g = np.transpose(np.transpose(resp_iter[g])*np.transpose(X_model))            
            PX_g = precision_matrix_y[g] @ X_obs_tilde_g
            PY_g = (precision_matrix_y[g] @ Y_obs_tilde_g.T).T           
            for k in range(2+N_period):
                S_xy[g,k] = np.sum(np.diag((X_obs_tilde_g.T @ PY_g)[k]))
            for i in range(N_ind):
                S_xx[g] += X_obs_tilde_g[i].T@PX_g[i]
            if (np.linalg.det(S_xx[g])==0): beta_iter[g] = beta_iter_last[g]
            else: beta_iter[g] = np.linalg.inv(S_xx[g]) @ S_xy[g]

        #### E-STEP #####
        density_y = np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C')
        res_X_iter, res_iter = np.zeros(shape=(K,N_ind,N_period,N_covariates), dtype=float, order='C'), np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C')
        resp_iter = np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C')
        Mu_tilde1 = np.zeros(shape=(K,N_covariates), dtype=float, order='C')
        joint_density = np.zeros(shape=(K,N_ind,N_period), dtype=float, order='C')
        for j in range(K):
            res_iter[j] = X_model @ beta_iter[j] - Y_obs
            res_X_iter[j] = X_obs - Mu_tilde[j]
            density_y[j] = scipy.stats.norm.pdf(x=res_iter[j], loc=0, scale=np.sqrt(np.diag(var_cov_y[j])))
            if (mahalanobis==1): 
                for i in range(N_ind):
                    joint_density[j,i] = -np.diag(res_X_iter[j,i]@np.linalg.inv(var_cov_X[j])@res_X_iter[j,i].T)
            else:
                joint_density[j] = scipy.stats.multivariate_normal.pdf(x=res_X_iter[j], mean=None, cov=var_cov_X[j])*density_y[j]
        resp_class  = np.argmax(joint_density, axis=0)
        for j in range(K):
            arg_resp = np.where(resp_class==j)
            resp_iter[j,arg_resp[0],arg_resp[1]] = 1    
            for p in range(N_covariates):
                Mu_tilde1[j,p] = np.sum(np.diag(np.transpose(X_obs[0:,0:,p]) @ resp_iter[j]))
        sum_NT = np.sum(resp_iter, axis=(1,2))
        
        ### Covariance matrix updates ###
        Mu_tilde = np.zeros(shape=(K,N_covariates), dtype=float, order='C')
        res_X_iter1 = np.zeros(shape=(K,N_ind,N_covariates,N_period), dtype=float, order='C')
        var_cov_y, var_cov_X = np.zeros(shape=(K, N_period, N_period), dtype=float, order='C'), np.zeros(shape=(K, N_covariates, N_covariates), dtype=float, order='C')
        precision_matrix_y = np.zeros(shape=(K, N_period, N_period), dtype=float, order='C')
        var_total, var_unit, var_random = np.zeros(shape=(K), dtype=float, order='C'), np.zeros(shape=(K), dtype=float, order='C'), np.zeros(shape=(K), dtype=float, order='C')        
        for j in range(K):
            Mu_tilde[j] = Mu_tilde1[j]/np.sum(resp_iter[j])   
            var_total[j] = np.sum(resp_iter[j]*(res_iter[j]**2))/(sum_NT[j]-2-N_period)
            res_iter1 = res_iter[j,np.sum(resp_iter[j], axis=1)>1] ###
            resp_iter1 = resp_iter[j,np.sum(resp_iter[j], axis=1)>1]
            mean_unit_error = np.average(res_iter1, weights=resp_iter1, axis=1)
            time_weights = np.sum(resp_iter1, axis=1)
            var_unit[j] = np.average((mean_unit_error - np.average(mean_unit_error, weights = time_weights))**2, weights = time_weights)
            var_random[j]  =  var_total[j] - var_unit[j]
            var_cov_y[j] = np.diag(N_period*[var_random[j]]) + var_unit[j]
            for t in range(N_period):
                res_X_iter1[j,0:,0:,t] = X_obs[0:,t,0:] - Mu_tilde[j]
                var_cov_X[j] += ((resp_iter[j,0:,t]*np.transpose(res_X_iter1[j,0:,0:,t]))@(res_X_iter1[j,0:,0:,t]))
            var_cov_X[j] = var_cov_X[j]/(np.sum(resp_iter[j])-N_covariates)
            eigen_cov_x = np.linalg.eigh(var_cov_X[j])
            eigen_cov_x[0][np.where(np.abs(eigen_cov_x[0])< var_min_X)] = var_min_X
            var_cov_X[j] = eigen_cov_x[1] @ np.diag(eigen_cov_x[0]) @ np.transpose(eigen_cov_x[1])            
            eigen_cov_y = np.linalg.eigh(var_cov_y[j])
            eigen_cov_y[0][np.where(np.abs(eigen_cov_y[0])< var_min)] = var_min
            var_cov_y[j] = eigen_cov_y[1] @ np.diag(eigen_cov_y[0]) @ np.transpose(eigen_cov_y[1]) 
            precision_matrix_y[j] = np.linalg.inv(var_cov_y[j])
        
        ### Log-likelihood computation ###
        pi_g = np.zeros(K, dtype=float, order='C')
        for j in range(K):
            pi_g[j] = np.sum(resp_iter[j])/(N_ind*N_period)   
        if (mahalanobis==0):
            log_lik = np.sum(np.log(np.sum(joint_density*resp_iter, axis=0)))
        else:
            log_lik = np.sum(np.log(np.sum(density_y*resp_iter, axis=0)))
        list.append(log_lik_total, log_lik)
        if (l!=0): changes = np.sum(resp_class_last!=resp_class)
        list.append(sum_changes, changes)

        ### Break mechanism ###
        if (l==0): 
            log_lik_total = log_lik_total[1:]
            sum_changes = sum_changes[1:]
        elif (sum_changes[l-2]==0):break
        
        ## Save last results ##
        resp_class_last = resp_class
        resp_iter_last = resp_iter
        beta_iter_last = beta_iter
        var_total_last = var_total
        pi_g_last = pi_g
        var_unit_last = var_unit
        var_cov_y_last = var_cov_y
        Mu_tilde_last = Mu_tilde   

    return log_lik_total, l, resp_iter_last, beta_iter_last, var_total_last, pi_g_last, var_unit_last, var_cov_y_last, Mu_tilde_last, sum_changes
    
#%%
def Relabelling_EM(K, Max_iter, tolerance, true_init):
    result_EM = EM_without_cov_dens(K, Max_iter, tolerance=tolerance, true_init=0, var_max=1000, var_min=1e-8, dens_min=1e-300)
    resp_iter1 = result_EM[2]
    beta_iter = result_EM[3]
    var_total = result_EM[4]
    pi_g = result_EM[5]
    var_unit = result_EM[6]
    X_model = Init_values[0]
    Y_obs = generated_data[0]
    Beta_star = generated_data[2]
    Varcov_random_eff = generated_data[7]
    assignment_mat = generated_data[3]
    pi_star = generated_data[4]
    Time_fixed_eff = generated_data[6]
    random_var = generated_data[9]      
            
    #### Relabelling of the component ####
    Beta_star_total = np.c_[np.transpose(Beta_star[0]),np.transpose(Beta_star[1]),np.transpose(Time_fixed_eff)]
    Var_star_total = np.zeros(shape=(np.shape(np.diag(Varcov_random_eff))[0]), dtype=float, order='C')
    for g in range(np.shape(np.diag(Varcov_random_eff))[0]):
        Var_star_total[g] = random_var + np.diag(Varcov_random_eff)[g]
    Var_star_assign =  np.zeros(shape=(K), dtype=float, order='C')
    Beta_star_assign =  np.zeros(shape=(K, 2+N_period), dtype=float, order='C')
    pi_star_assign = np.zeros(shape=(K), dtype=float, order='C')
    resp_max = np.argmax(resp_iter1, axis=0)
    resp_class = np.zeros(shape=(N_ind, N_period, K), dtype=int, order='C')
    for i in range(N_ind):
        for t in range(N_period):
            resp_class[i,t,resp_max[i,t]] = 1   
    assignment_mat1 = np.zeros(shape=(N_ind,N_period,K), dtype=float, order='C')
    assign_comp = np.zeros(shape=(K_real, K), dtype=float, order='C')
    for i in range(K):
        for j in range(K_real):
            assign_comp[j,i] = np.sum(np.abs(Beta_star_total[j]-beta_iter[i]))
    for i in range(K_real):
        min_position = np.where(assign_comp == np.min(assign_comp))
        Var_star_assign[min_position[1]] = Var_star_total[min_position[0]]
        Beta_star_assign[min_position[1]] = Beta_star_total[min_position[0]]
        assignment_mat1[0:,0:,min_position[1]] = assignment_mat[0:,0:,min_position[0]]
        pi_star_assign[min_position[1]] = pi_star[min_position[0]]
        assign_comp[min_position[0],0:] = (N_ind*N_period)*9999
        assign_comp[0:,min_position[1]] = (N_ind*N_period)*9999
    
    ## Bias and classification error (for the last iteration)
    bias_cov = np.zeros(shape=(K, 2), dtype=float, order='C') 
    bias_beta = np.zeros(shape=(K, 2+N_period), dtype=float, order='C') 
    bias_cov_w = np.zeros(shape=(K, 2), dtype=float, order='C') 
    bias_beta_w = np.zeros(shape=(K, 2+N_period), dtype=float, order='C')
    bias_pi = np.zeros(shape=(K), dtype=float, order='C')
    for g in range(K):
        bias_beta_w[g] = pi_star_assign[g]*(beta_iter[g] - Beta_star_assign[g])
        bias_cov_w[g] = [pi_star_assign[g]*((var_total[g]-var_unit[g]) - random_var), pi_star_assign[g]*(var_unit[g] - (Var_star_assign[g]-random_var))]
        bias_beta[g] = beta_iter[g] - Beta_star_assign[g]
        bias_cov[g] = [(var_total[g]-var_unit[g]) - random_var, var_unit[g] - (Var_star_assign[g]-random_var)]  
    bias_pi = pi_g - pi_star_assign      
    bias_beta = np.reshape(bias_beta, shape=(K*(2+N_period)))
    bias_cov = np.reshape(bias_cov, shape=(K*2))
    bias_beta_w = np.reshape(bias_beta_w, shape=(K*(2+N_period)))
    bias_cov_w = np.reshape(bias_cov_w, shape=(K*2))
    class_error = np.zeros(shape=(N_ind, N_period, K), dtype=int, order='C')
    pred_Y  = np.zeros(shape=(K,N_ind, N_period), dtype=float, order='C')
    pred_error = np.zeros(shape=(N_ind, N_period), dtype=float, order='C')
    for g in range(K):
        pred_Y[g] = X_model @ beta_iter[g]    
    pred_error = np.average(Y_obs - pred_Y, weights= resp_iter1, axis=0)
    class_error = np.abs(resp_class - assignment_mat1)/2
    class_error = np.sum(class_error)/(N_ind*N_period)

    ## Prediction error ##
    R_2 = 1-(np.sum(pred_error**2)/np.sum((Y_obs-np.average(Y_obs))**2))
    R_2a = 1-((1-R_2)*(N_ind*N_period-1)/(N_ind*N_period-K*(N_covariates*2+N_period)))
    RMSE = np.sqrt(np.sum(bias_beta**2)/bias_beta.shape[0])
    RMSE_var = np.sqrt(np.sum(bias_cov**2)/bias_cov.shape[0])
    RMSE_w = np.sqrt(np.sum(bias_beta_w**2)/(bias_beta_w.shape[0]/K))
    RMSE_var_w = np.sqrt(np.sum(bias_cov_w**2)/(bias_cov_w.shape[0]/K))
    return pi_star_assign, class_error, bias_pi, bias_beta,bias_cov,bias_beta_w,bias_cov_w, resp_class, R_2a, R_2, RMSE, RMSE_w, RMSE_var, RMSE_var_w, Var_star_assign, Beta_star_assign, beta_iter

#%%
def Relabelling_CEM(K, Max_iter):
    result_CEM = CEM_with_cov_dens(K=K, Max_iter=Max_iter, mahalanobis=mahalanobis, var_max=1000, var_min=1e-8, dens_min=1e-300)
    resp_iter1 = result_CEM[2]
    beta_iter = result_CEM[3]
    var_total = result_CEM[4]
    pi_g = result_CEM[5]
    Mu_tilde = result_CEM[8]
    var_unit = result_CEM[6]
    X_model = Init_values[0]
    Y_obs = generated_data[0]
    Beta_star = generated_data[2]
    Varcov_random_eff = generated_data[7]
    assignment_mat = generated_data[3]
    pi_star = generated_data[4]
    Time_fixed_eff = generated_data[6]
    random_var = generated_data[9]
    mean_covariates = generated_data[13]

    #### Relabelling of the component ####
    Beta_star_total = np.c_[np.transpose(Beta_star[0]),np.transpose(Beta_star[1]),np.transpose(Time_fixed_eff)]
    Var_star_total = np.zeros(shape=(np.shape(np.diag(Varcov_random_eff))[0]), dtype=float, order='C')
    for g in range(np.shape(np.diag(Varcov_random_eff))[0]):
        Var_star_total[g] = random_var + np.diag(Varcov_random_eff)[g]
    Var_star_assign =  np.zeros(shape=(K), dtype=float, order='C')
    Beta_star_assign =  np.zeros(shape=(K, 2+N_period), dtype=float, order='C')
    pi_star_assign = np.zeros(shape=(K), dtype=float, order='C')
    resp_max = np.argmax(resp_iter1, axis=0)
    resp_class = np.zeros(shape=(N_ind, N_period, K), dtype=int, order='C')
    for i in range(N_ind):
        for t in range(N_period):
            resp_class[i,t,resp_max[i,t]] = 1   
    assignment_mat1 = np.zeros(shape=(N_ind,N_period,K), dtype=float, order='C')
    assign_comp = np.zeros(shape=(K_real, K), dtype=float, order='C')
    covariates_star_assign = np.zeros(shape=(K,N_covariates), dtype=float, order='C')
    for i in range(K):
        for j in range(K_real):
            assign_comp[j,i] = np.sum(np.abs(Beta_star_total[j]-beta_iter[i]))
    for i in range(K_real):
        min_position = np.where(assign_comp == np.min(assign_comp))
        Var_star_assign[min_position[1]] = Var_star_total[min_position[0]]
        Beta_star_assign[min_position[1]] = Beta_star_total[min_position[0]]
        assignment_mat1[0:,0:,min_position[1]] = assignment_mat[0:,0:,min_position[0]]
        pi_star_assign[min_position[1]] = pi_star[min_position[0]]
        covariates_star_assign[min_position[1]] = mean_covariates[min_position[0]]
        assign_comp[min_position[0],0:] = (N_ind*N_period)*9999
        assign_comp[0:,min_position[1]] = (N_ind*N_period)*9999

    ## Bias and classification error (for the last iteration)
    bias_cov = np.zeros(shape=(K, 2), dtype=float, order='C') 
    bias_beta = np.zeros(shape=(K, 2+N_period), dtype=float, order='C') 
    bias_cov_w = np.zeros(shape=(K, 2), dtype=float, order='C') 
    bias_beta_w = np.zeros(shape=(K, 2+N_period), dtype=float, order='C')
    bias_pi = np.zeros(shape=(K), dtype=float, order='C')
    bias_covariates_w = np.zeros(shape=(K, N_covariates), dtype=float, order='C') 
    bias_covariates = np.zeros(shape=(K, N_covariates), dtype=float, order='C') 
    for g in range(K):
        bias_beta_w[g] = pi_star_assign[g]*(beta_iter[g] - Beta_star_assign[g])
        bias_cov_w[g] = [pi_star_assign[g]*((var_total[g]-var_unit[g]) - random_var), pi_star_assign[g]*(var_unit[g] - (Var_star_assign[g]-random_var))]
        bias_beta[g] = beta_iter[g] - Beta_star_assign[g]
        bias_cov[g] = [(var_total[g]-var_unit[g]) - random_var, var_unit[g] - (Var_star_assign[g]-random_var)]
        bias_covariates[g] =  covariates_star_assign[g] - Mu_tilde[g]
        bias_covariates_w[g] = pi_star_assign[g]*(covariates_star_assign[g] - Mu_tilde[g])       
    bias_pi = pi_g - pi_star_assign      
    bias_beta = np.reshape(bias_beta, shape=(K*(2+N_period)))
    bias_cov = np.reshape(bias_cov, shape=(K*2))
    bias_beta_w = np.reshape(bias_beta_w, shape=(K*(2+N_period)))
    bias_cov_w = np.reshape(bias_cov_w, shape=(K*2))
    class_error = np.zeros(shape=(N_ind, N_period, K), dtype=int, order='C')
    pred_Y  = np.zeros(shape=(K,N_ind, N_period), dtype=float, order='C')
    pred_error = np.zeros(shape=(N_ind, N_period), dtype=float, order='C')
    for g in range(K):
        pred_Y[g] = X_model @ beta_iter[g]    
    pred_error = np.average(Y_obs - pred_Y, weights= resp_iter1, axis=0)
    class_error = np.abs(resp_class - assignment_mat1)/2
    class_error = np.sum(class_error)/(N_ind*N_period)

    ## Prediction error ##
    R_2 = 1-(np.sum(pred_error**2)/np.sum((Y_obs-np.average(Y_obs))**2))
    R_2a = 1-((1-R_2)*(N_ind*N_period-1)/(N_ind*N_period-K*(N_covariates*2+N_period)))
    RMSE = np.sqrt(np.sum(bias_beta**2)/bias_beta.shape[0])
    RMSE_var = np.sqrt(np.sum(bias_cov**2)/bias_cov.shape[0])
    RMSE_w = np.sqrt(np.sum(bias_beta_w**2)/(bias_beta_w.shape[0]/K))
    RMSE_var_w = np.sqrt(np.sum(bias_cov_w**2)/(bias_cov_w.shape[0]/K))
    return pi_star_assign, class_error, bias_pi, bias_beta,bias_cov,bias_beta_w,bias_cov_w, resp_class, R_2a, R_2, RMSE, RMSE_w, RMSE_var, RMSE_var_w, Var_star_assign, Beta_star_assign, beta_iter, bias_covariates, bias_covariates_w
#%%
### Check the Data-Generating Process #####
### Continuous outcome ###
seed=1 ## Random generation seed
K=3 ## Number of estimated components
K_real=3 ## True number of components  
N_ind=500 ## Number of individuals in the panel 
N_period=5 ## Number of periods
N_covariates=10  ## Number of covariates (all continuous)
range_beta=0 ## Increase for better EM performance ##
range_cov=0 ## Increase for better CEM performance ##
random_var=1 ## Variance of the idiosyncratic error term, default = 1
AR_process = 1  ### =1 group membership is an AR(1) process, otherwise random
alpha_K_gen = 10  ## Parameter the assignment process, default = 10
sigma2_X = 1 ### Diagonal term of P_g, default = 1
mahalanobis = 0 ### =1 use the Mahalanobis distance for classification in the CEM, otherwise joint density classifier with a normal multivariate distribution


generated_data = generate_normal_data(seed=seed, K=K_real, alpha_K=alpha_K_gen, range_beta=range_beta, sigma2_X = sigma2_X, range_cov=range_cov, N_ind=N_ind, N_period=N_period, N_covariates=N_covariates, random_var=random_var, AR_process=AR_process)
optim = optimal_err_class(K_real=K_real, mahalanobis=mahalanobis)
print("pi_star=",generated_data[4],"class_error_with_cov=",optim[0],"loglik_true_value=",optim[1], "class_loglik_true_value=",optim[2]) 
print("pi_star_t=",generated_data[-2])


#%%
## Monte carlo setup - Continuous outcome ##
n_rep = 10 ## Number of replications
n_b= 10 ## Number of initial values per replication (random initialization)
N_ind = 500 ## Number of individuals in the panel 
N_period= 5 ## Number of periods
N_covariates=5 ## Number of covariates (all continuous)
K=2 ## Number of estimated components
K_real = 2 ## True number of components  
range_beta=0 ### Range of values for the beta^0_g, default = 0, increase for better EM performance 
range_cov=0 ### Range of values for the mean mu^0_g, default = 0, Increase for better CEM performance
random_var = 1 ## Variance of the idiosyncratic error term, default = 1
sigma2_X=1 ### Diagonal term of P_g, default = 1
var_init = 100 ## Initial outcome variance for both algorithms, default =100 ##
min_eigen = 1e-8 ### Minimum eigenvalue for covariance matrices
tolerance = 0.0001 ## Tolerance for convergence of EM algorithm ##
Max_iter = 100 ## Maximum nummber of iterations ##
AR_process = 1 ### =1 group membership is an AR(1) process, otherwise random
alpha_K_gen = 10 ## Parameter the assignment process, default = 10
alpha_K = 5 ## Parameter for the assignement process of the initial values, default = 5
mahalanobis = 0 ### =1 use the Mahalanobis distance for classification in the CEM, otherwise joint density classifier with a normal multivariate distribution
true_init = 0 ### =1 initialize with the true parameter values for the EM, otherwise random

### Store results ###
class_error_CEM = np.zeros(shape=(n_rep),dtype=float, order='C')
class_error_EM = np.zeros(shape=(n_rep),dtype=float, order='C')
bias_beta_CEM = np.zeros(shape=(n_rep, (2+N_period)*K),dtype=float, order='C')
bias_beta_EM = np.zeros(shape=(n_rep, (2+N_period)*K),dtype=float, order='C')
bias_var_CEM = np.zeros(shape=(n_rep, K*2),dtype=float, order='C')
bias_var_EM = np.zeros(shape=(n_rep, K*2),dtype=float, order='C')
bias_pi_CEM = np.zeros(shape=(n_rep, K),dtype=float, order='C')
bias_pi_EM = np.zeros(shape=(n_rep, K),dtype=float, order='C')
pi_star_rep = np.zeros(shape=(n_rep, K),dtype=float, order='C') 

for k in range(n_rep):
    generated_data = generate_normal_data(seed=k, K=K_real, alpha_K=alpha_K_gen, range_beta=range_beta, sigma2_X = sigma2_X, range_cov=range_cov, N_ind=N_ind, N_period=N_period, N_covariates=N_covariates, random_var=random_var, AR_process=AR_process)
    loglik_total = np.zeros(shape=(n_b, 3),dtype=float, order='C')
    pi_star_rep[k] = generated_data[4]
    for m in range(n_b):
        Init_values =  Initial_values(seed=m, K=K, alpha_K = alpha_K, var_init=var_init, dens_min=1e-320)
        try:
            result_CEM = CEM_with_cov_dens(K=K, Max_iter=Max_iter, mahalanobis=mahalanobis, var_max=1000, var_min=min_eigen, dens_min=1e-320)
            result_EM = EM_without_cov_dens(K=K, Max_iter=Max_iter, tolerance=tolerance, true_init=true_init, var_max=1000, var_min=min_eigen, dens_min=1e-320)
        except:
            pass
        loglik_total[m,0] = m
        loglik_total[m,1] = result_CEM[0][-1]
        loglik_total[m,2] = result_EM[0][-1] 
        plt.plot(result_CEM[0])
        plt.plot(result_EM[0])
        plt.show()

    ## Select the starting values that maximizes the log-likelihood and check the biases ##
    Init_values = Initial_values(seed=np.argmax(loglik_total[0:,1]), K=K, alpha_K = alpha_K, var_init=var_init, dens_min=1e-320)
    bias_CEM = Relabelling_CEM(K=K, Max_iter=Max_iter)
    class_error_CEM[k] = bias_CEM[1]
    bias_pi_CEM[k] = bias_CEM[2]
    bias_beta_CEM[k] = bias_CEM[3]
    bias_var_CEM[k] = bias_CEM[4]
    print(bias_CEM[11])
    print(class_error_CEM[k])
    Init_values = Initial_values(seed=np.argmax(loglik_total[0:,2]), K=K, alpha_K = alpha_K, var_init=var_init, dens_min=1e-320)
    bias_EM = Relabelling_EM(K=K, Max_iter=Max_iter, tolerance=tolerance, true_init=true_init)        
    class_error_EM[k] = bias_EM[1]
    bias_pi_EM[k] = bias_EM[2]
    bias_beta_EM[k] = bias_EM[3]
    bias_var_EM[k] = bias_EM[4]
    print(bias_EM[11])
    print(class_error_EM[k])

## Bias computation ##
print(np.average(bias_beta_CEM, axis=0))
print(np.average(bias_beta_EM, axis=0))
print(np.average(bias_var_CEM, axis=0))
print(np.average(bias_var_EM, axis=0))
print(np.average(bias_pi_CEM, axis=0))
print(np.average(bias_pi_EM, axis=0))
plt.hist(class_error_CEM, bins=5)
plt.hist(class_error_EM, bins=5)
plt.show()

### Results of function (reminder) ###
#EM: pi_star_assign, class_error, bias_pi, bias_beta,bias_cov,bias_beta_w,bias_cov_w, resp_class, R_2a, R_2, RMSE, RMSE_w, RMSE_var, RMSE_var_w, Var_star_assign, Beta_star_assign, beta_iter
#CEM: pi_star_assign, class_error, bias_pi, bias_beta,bias_cov,bias_beta_w,bias_cov_w, resp_class, R_2a, R_2, RMSE, RMSE_w, RMSE_var, RMSE_var_w, Var_star_assign, Beta_star_assign, beta_iter, bias_covariates, bias_covariates_w

#%%
np.savetxt('class_error_EM_cov5_K2_500.txt', class_error_EM, delimiter=",", fmt='%.12g') 
np.savetxt('class_error_CEM_cov5_K2_500.txt', class_error_CEM, delimiter=",", fmt='%.12g') 
np.savetxt('bias_beta_CEM_cov5_K2_500.txt', bias_beta_CEM, delimiter=",", fmt='%.12g') 
np.savetxt('bias_beta_EM_cov5_K2_500.txt', bias_beta_EM, delimiter=",", fmt='%.12g') 
np.savetxt('bias_var_CEM_cov5_K2_500.txt', bias_var_CEM, delimiter=",", fmt='%.12g') 
np.savetxt('bias_var_EM_cov5_K2_500.txt', bias_var_EM, delimiter=",", fmt='%.12g') 
np.savetxt('bias_pi_CEM_cov5_K2_500.txt', bias_pi_CEM, delimiter=",", fmt='%.12g') 
np.savetxt('bias_pi_EM_cov5_K2_500.txt', bias_pi_EM, delimiter=",", fmt='%.12g') 
np.savetxt('pi_star_rep_cov5_K2_500.txt', pi_star_rep, delimiter=",", fmt='%.12g')


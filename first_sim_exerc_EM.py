# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 11:38:58 2025

@author: raph1177
"""
### First simulation exercise ###

### Bias-reduced estimation of finite mixtures ###

### By Raphaël Langevin ###
#%%
import pandas as pd
import numpy as np
import math
import scipy
import matplotlib.pyplot as plt
import random
import os
import time
np.version.version
os.chdir('C:\\Users\\raph1\.spyder-py3\Work')
#%%
### Normal distributions with G=2 - EM algorithm ###
def mix_normal_EM(nb_simul, d, max_iter, penalty, mu1_0, mu2_0, sigma1_0, sigma2_0, pi1_0, dens_min = 1e-320, min_sd=1e-10, tolerance = 1e-10):
    t = time.process_time()
    quant = (0.9999,0.9999,0.9995,0.9995,0.999,0.999,0.995,0.995,0.99,0.99,0.98,0.98,0.97,0.97,0.95,0.95)
    log_lik_out_init = np.zeros(shape=(d,max_iter,nb_simul,len(quant)), order="C")
    log_lik_true_init = np.zeros(shape=(d,max_iter,nb_simul), order="C")
    mu1_high = np.zeros(shape=(d-1,nb_simul), order="C")
    mu2_high = np.zeros(shape=(d-1,nb_simul), order="C")
    sd1_high = np.zeros(shape=(d-1,nb_simul), order="C")
    sd2_high = np.zeros(shape=(d-1,nb_simul), order="C")
    pi1_high = np.zeros(shape=(d-1,nb_simul), order="C")
    pi2_high = np.zeros(shape=(d-1,nb_simul), order="C")
    
    ## Start simulations ##
    for k in range(nb_simul):
        for n in range(1,d):
            N_1 = int(pi1_0*10*(10**n))
            N_2 = int((1-pi1_0)*10*(10**n))
            sum_N = N_1+N_2
            pi2_0 = 1-pi1_0
            np.random.seed(k)
            data1 = np.random.normal(loc=mu1_0, scale=sigma1_0, size=N_1)
            data2 = np.random.normal(loc=mu2_0, scale=sigma2_0, size=N_2)
            data_total = np.hstack([data1,data2])
            del data1,data2
            
            ## Outliers initiation ##
            j_list = list()
            mu1_list = list()
            mu2_list = list()
            sd1_list = list()
            sd2_list = list()
            pi1_list = list()
            pi2_list = list()
            for b in range(len(quant)):
                if (b%2==0): 
                    quant_1 = np.quantile(data_total, q=quant[b])
                    data_outlier  = data_total[np.where(data_total>=quant_1)]
                    data_outlier1 = data_total[np.where(data_total<quant_1)]
                if (b%2==1): 
                    quant_1 = np.quantile(data_total, q=(1-quant[b]))
                    data_outlier  = data_total[np.where(data_total<=quant_1)]
                    data_outlier1  = data_total[np.where(data_total>quant_1)]
                mu1 = np.array(np.mean(data_outlier))
                sd1 = np.array(np.std(data_outlier))
                mu2 = np.array(np.mean(data_outlier1))
                sd2 = np.array(np.std(data_outlier1))
                if (math.isnan(mu1)): mu1 = np.max(data_total)
                if (math.isnan(mu2)): mu2 = np.min(data_total)
                sd1 = np.nanmax([sd1,min_sd])
                sd2 = np.nanmax([sd2,min_sd])
                dens1 = scipy.stats.norm.pdf(data_total, loc=mu1, scale=sd1)
                dens2 = scipy.stats.norm.pdf(data_total, loc=mu2, scale=sd2) 
                dens1[np.where(dens1<dens_min)]=dens_min
                dens2[np.where(dens2<dens_min)]=dens_min
                post_prob1 = (dens1)/(dens1+dens2)
                post_prob2 = 1-post_prob1
                pi1 = [np.shape(data_outlier)[0]/sum_N]
                pi2 = [np.shape(data_outlier1)[0]/sum_N]
                N1 = np.sum(post_prob1)
                N2 = np.sum(post_prob2)
                penalty1 = -sd1**(-2)-np.log(sd1**2)
                penalty2 = -sd2**(-2)-np.log(sd2**2)
                if penalty==1: log_lik_out_init[n,0,k,b] = np.sum(np.log(pi1*dens1+pi2*dens2))+1/np.sqrt(sum_N)*(penalty1+penalty2)
                else: log_lik_out_init[n,0,k,b] = np.sum(np.log(pi1*dens1+pi2*dens2))
                
                ## Outliers initiation - Iterative algorithm ##
                for j in range(1,max_iter):
                    mu1 = np.hstack([mu1,np.average(data_total, weights=post_prob1)])
                    mu2 = np.hstack([mu2,np.average(data_total, weights=post_prob2)])
                    N1 = np.hstack([N1, np.sum(post_prob1)])
                    N2 = np.hstack([N2, np.sum(post_prob2)])
                    sd1 = np.hstack([sd1,np.max([np.sqrt(1/(N1[j])*np.sum(post_prob1*(data_total-mu1[j])**2)),min_sd])])
                    sd2 = np.hstack([sd2,np.max([np.sqrt(1/(N2[j])*np.sum(post_prob2*(data_total-mu2[j])**2)),min_sd])])
                    dens1 = scipy.stats.norm.pdf(data_total, loc=mu1[j], scale=sd1[j])
                    dens2 = scipy.stats.norm.pdf(data_total, loc=mu2[j], scale=sd2[j])
                    dens1[np.where(dens1<dens_min)]=dens_min
                    dens2[np.where(dens2<dens_min)]=dens_min
                    post_prob1 = (pi1[j-1]*dens1)/(pi1[j-1]*dens1+pi2[j-1]*dens2)
                    post_prob2 = 1-post_prob1
                    if (np.mean(post_prob2)==0.0):
                        j = j-1
                        break
                    if (np.mean(post_prob1)==0.0):
                        j = j-1
                        break
                    pi1.append(np.mean(post_prob1))
                    pi2.append(np.mean(post_prob2))
                    total_dens = (pi1[j])*dens1+(pi2[j])*dens2
                    penalty1 = -sd1[j]**(-2)-np.log(sd1[j]**2)
                    penalty2 = -sd2[j]**(-2)-np.log(sd2[j]**2)
                    if penalty==1: log_lik_out_init[n,j,k,b] = np.sum(np.log(total_dens))+1/np.sqrt(sum_N)*(penalty1+penalty2) ##1/np.sqrt(sum_N)*
                    else: log_lik_out_init[n,j,k,b] = np.sum(np.log(total_dens))
                    del total_dens, dens1, dens2
                    if j>1: #### Neglect the potential initial drop in the likelihood value (this is a caveat of the EM algorithm) ###
                        if ((log_lik_out_init[n,j,k,b]-log_lik_out_init[n,j-1,k,b])<tolerance): break        

                
                ### Select the parameters with the highest likelihood value ###
                j_list.append(j-1)
                mu1_list.append(mu1[-2])
                mu2_list.append(mu2[-2])
                sd1_list.append(sd1[-2])
                sd2_list.append(sd2[-2])
                pi1_list.append(pi1[-2])
                pi2_list.append(pi2[-2])
                if b>0:
                    if (log_lik_out_init[n,j_list[b-1],k,b-1] > log_lik_out_init[n,j_list[b],k,b]): 
                        log_lik_out_init[n,0:,k,b] = log_lik_out_init[n,0:,k,b-1]
                        j_list[b] = j_list[b-1]
                        mu1_list[b] = mu1_list[b-1]
                        mu2_list[b] = mu2_list[b-1]
                        sd1_list[b] = sd1_list[b-1]
                        sd2_list[b] = sd2_list[b-1]
                        pi1_list[b] = pi1_list[b-1]
                        pi2_list[b] = pi2_list[b-1]
            
            ## Initiation with the true parameter values ###
            mu1_t, mu2_t, sd1_t, sd2_t = mu1_0, mu2_0, sigma1_0, sigma2_0
            dens1_t = scipy.stats.norm.pdf(data_total, loc=mu1_t, scale=sd1_t)
            dens2_t = scipy.stats.norm.pdf(data_total, loc=mu2_t, scale=sd2_t) 
            dens1_t[np.where(dens1_t<dens_min)]=dens_min
            dens2_t[np.where(dens2_t<dens_min)]=dens_min
            post_prob1_t = (pi1_0*dens1_t)/(pi1_0*dens1_t+pi2_0*dens2_t)
            post_prob2_t = 1-post_prob1_t
            pi1_t = [pi1_0]
            pi2_t = [pi2_0]
            N1_t = np.sum(post_prob1_t)
            N2_t = np.sum(post_prob2_t)
            penalty1_t = -sd1_t**(-2)-np.log(sd1_t**2)
            penalty2_t = -sd2_t**(-2)-np.log(sd2_t**2)
            if penalty==1: log_lik_true_init[n,0,k] = np.sum(np.log(pi1_t*dens1_t+pi2_t*dens2_t))+1/np.sqrt(sum_N)*(penalty1_t+penalty2_t) ##1/np.sqrt(sum_N)*
            else: log_lik_true_init[n,0,k] = np.sum(np.log(pi1_t*dens1_t+pi2_t*dens2_t))
            
            ## True values initiation - Iterative algorithm ###
            for m in range(1,max_iter):
                mu1_t = np.hstack([mu1_t,np.average(data_total, weights=post_prob1_t)])
                mu2_t = np.hstack([mu2_t,np.average(data_total, weights=post_prob2_t)])
                N1_t = np.hstack([N1_t, np.sum(post_prob1_t)])
                N2_t = np.hstack([N2_t, np.sum(post_prob2_t)])
                sd1_t = np.hstack([sd1_t,np.max([np.sqrt(1/(N1_t[m])*np.sum(post_prob1_t*(data_total-mu1_t[m])**2)),min_sd])])
                sd2_t = np.hstack([sd2_t,np.max([np.sqrt(1/(N2_t[m])*np.sum(post_prob2_t*(data_total-mu2_t[m])**2)),min_sd])])
                dens1_t = scipy.stats.norm.pdf(data_total, loc=mu1_t[m], scale=sd1_t[m])
                dens2_t = scipy.stats.norm.pdf(data_total, loc=mu2_t[m], scale=sd2_t[m])
                dens1_t[np.where(dens1_t<dens_min)]=dens_min
                dens2_t[np.where(dens2_t<dens_min)]=dens_min
                post_prob1_t = (pi1_t[m-1]*dens1_t)/(pi1_t[m-1]*dens1_t+pi2_t[m-1]*dens2_t)
                post_prob2_t = 1-post_prob1_t
                if (np.mean(post_prob2_t)==0.0): 
                    m = m-1
                    break
                if (np.mean(post_prob1_t)==0.0): 
                    m = m-1
                    break
                pi1_t.append(np.mean(post_prob1_t))
                pi2_t.append(np.mean(post_prob2_t))
                total_dens_t = pi1_t[m]*dens1_t+pi2_t[m]*dens2_t            
                penalty1_t = -sd1_t[m]**(-2)-np.log(sd1_t[m]**2)
                penalty2_t = -sd2_t[m]**(-2)-np.log(sd2_t[m]**2)
                if penalty==1: log_lik_true_init[n,m,k] = np.sum(np.log(total_dens_t))+1/np.sqrt(sum_N)*(penalty1_t+penalty2_t) ##1/np.sqrt(sum_N)*
                else: log_lik_true_init[n,m,k] = np.sum(np.log(total_dens_t))
                del total_dens_t, dens1_t, dens2_t
                if m>1: #### Neglect the potential initial drop in the likelihood value (this is a caveat of the EM algorithm) ###
                    if ((log_lik_true_init[n,m,k]-log_lik_true_init[n,m-1,k])<tolerance): break
            
            #### Select the final parameters with the highest likelihood value ####
            #### I deal with label switching by ordering the estimated mean values ###
            ### Except when mu1_0 = mu2_0, then switching is avoided by ordering standard deviations ###
            if (log_lik_true_init[n,m-1,k] >= log_lik_out_init[n,j_list[b],k,b]):
                if mu1_t[m-1] > mu2_t[m-1]: ### Change by sd1_t[m-1] < sd2_t[m-1] if mu1_0 = mu2_0  ###
                    mu1_high[n-1,k] = mu1_t[m-1]
                    mu2_high[n-1,k] = mu2_t[m-1]
                    sd1_high[n-1,k] = sd1_t[m-1]
                    sd2_high[n-1,k] = sd2_t[m-1]
                    pi1_high[n-1,k] = pi1_t[m-1]
                    pi2_high[n-1,k] = pi2_t[m-1]
                else: 
                    mu1_high[n-1,k] = mu2_t[m-1]
                    mu2_high[n-1,k] = mu1_t[m-1]
                    sd1_high[n-1,k] = sd2_t[m-1]
                    sd2_high[n-1,k] = sd1_t[m-1]
                    pi1_high[n-1,k] = pi2_t[m-1]
                    pi2_high[n-1,k] = pi1_t[m-1]
            if (log_lik_true_init[n,m-1,k] < log_lik_out_init[n,j_list[b],k,b]):
                if mu1_list[b] > mu2_list[b]: ### Change by sd1_list[b] < sd2_list[b] if mu1_0 = mu2_0  ###
                    mu1_high[n-1,k] = mu1_list[b]
                    mu2_high[n-1,k] = mu2_list[b]
                    sd1_high[n-1,k] = sd1_list[b]
                    sd2_high[n-1,k] = sd2_list[b]
                    pi1_high[n-1,k] = pi1_list[b]
                    pi2_high[n-1,k] = pi2_list[b]
                else: 
                    mu1_high[n-1,k] = mu2_list[b]
                    mu2_high[n-1,k] = mu1_list[b]
                    sd1_high[n-1,k] = sd2_list[b]
                    sd2_high[n-1,k] = sd1_list[b]
                    pi1_high[n-1,k] = pi2_list[b]
                    pi2_high[n-1,k] = pi1_list[b]        
    elapsed_time = time.process_time() - t
    return mu1_high, mu2_high, sd1_high, sd2_high, pi1_high, pi2_high, elapsed_time
#%%
## Run simulations ##
d = 4 ## Maximum sample size, where N = 10**n, with n = 2,...,d
nb_simul = 1000 ## Number of replications
max_iter = 100 ## Maxmimum number of iterations
penalty = 1 ## =1 for adding the penalty term
mu1_0 = 0  ## mu1_0 > mu2_0
mu2_0 = 0
sigma1_0 = 0.75 # sd1 < sd2
sigma2_0 = 1.25
pi1_0= 0.7
results_normal = mix_normal_EM(nb_simul=nb_simul, d=d, max_iter=max_iter, penalty=penalty, mu1_0=mu1_0, mu2_0=mu2_0, sigma1_0=sigma1_0, sigma2_0=sigma2_0, pi1_0=pi1_0, dens_min = 1e-320, min_sd=1e-10)
mu1_high = results_normal[0]
mu2_high = results_normal[1]
sd1_high = results_normal[2]
sd2_high = results_normal[3]
pi1_high = results_normal[4]
pi2_high = results_normal[5]
print(results_normal[6]/60)


#%%
### Plot and save the results ####
## Mu_1 ##
plt.plot(np.mean(mu1_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(mu1_high, axis=1), yerr=1.96*np.std(mu1_high,ddof=1, axis=1))
plt.plot((d-1)*[mu1_0])
plt.title("Mu_1 average")
plt.show()
plt.hist(mu1_high[d-2,0:], bins=40)
plt.title("Mu_1 distribution")
plt.show()
np.savetxt('norm_mu1_high_0_sd_1_07.txt', mu1_high, delimiter=",", fmt='%.12g')

## Mu_2 ##
plt.plot(np.mean(mu2_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(mu2_high, axis=1), yerr=1.96*np.std(mu2_high,ddof=1, axis=1))
plt.plot((d-1)*[mu2_0]) 
plt.title("Mu_2 average")
plt.show()
plt.hist(mu2_high[d-2,0:], bins=40)
plt.title("Mu_2 distribution")
plt.show()
np.savetxt('norm_mu2_high_0_sd_1_07.txt', mu2_high, delimiter=",", fmt='%.12g') 

## Sd_1 ##
plt.plot(np.mean(sd1_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(sd1_high, axis=1), yerr=(1.96*np.std(sd1_high,ddof=1, axis=1)))
plt.plot((d-1)*[sigma1_0])
plt.title("Sd_1 average")
plt.show()
plt.hist(sd1_high[d-2,0:], bins=40)
plt.title("Sd_1 distribution")
plt.show()
np.savetxt('norm_sd1_high_0_sd_1_07.txt', sd1_high, delimiter=",", fmt='%.12g') 

## Sd_2 ##
plt.plot(np.mean(sd2_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(sd2_high, axis=1), yerr=1.96*np.std(sd2_high,ddof=1, axis=1))
plt.plot((d-1)*[sigma2_0])
plt.title("Sd_2 average")
plt.show()
plt.hist(sd2_high[d-2,0:], bins=40)
plt.title("Sd_2 distribution")
plt.show()
np.savetxt('norm_sd2_high_0_sd_1_07.txt', sd2_high, delimiter=",", fmt='%.12g') 

## Pi_1 ##
plt.plot(np.mean(pi1_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(pi1_high, axis=1), yerr=1.96*np.std(pi1_high,ddof=1, axis=1))
plt.plot((d-1)*[pi1_0])
plt.title("Pi_1 average") 
plt.show()
plt.hist(pi1_high[d-2,0:], bins=40)
plt.title("Pi_1 distribution")
plt.show()
np.savetxt('norm_pi1_high_0_sd_1_07.txt', pi1_high, delimiter=",", fmt='%.12g') 

## Pi_2 ##
plt.plot(np.mean(pi2_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(pi2_high, axis=1), yerr=1.96*np.std(pi2_high,ddof=1, axis=1))
plt.plot((d-1)*[1-pi1_0])
plt.title("Pi_2 average") 
plt.show()
plt.hist(pi2_high[d-2,0:], bins=40)
plt.title("Pi_2 distribution") 
plt.show()
np.savetxt('norm_pi2_high_0_sd_1_07.txt', pi2_high, delimiter=",", fmt='%.12g') 

## Scatter plots ##
plt.scatter(mu1_high[d-2,0:],mu2_high[d-2,0:])
plt.xlabel("Estiamted Mu_1 value, N=10,000")
plt.ylabel("Estiamted Mu_2 value, N=10,000")
plt.show()
plt.scatter(sd1_high[d-2,0:],sd2_high[d-2,0:])
plt.xlabel("Estiamted Sd_1 value, N=10,000")
plt.ylabel("Estiamted Sd_2 value, N=10,000")
plt.show()
plt.scatter(pi1_high[d-2,0:],pi2_high[d-2,0:])
plt.xlabel("Estiamted Pi_1 value, N=10,000")
plt.ylabel("Estiamted Pi_2 value, N=10,000")
plt.show()

#%%
### Poisson distributions with G=2 ##
def mix_poisson_EM(nb_simul, d, max_iter, mu1_0, mu2_0, pi1_0, dens_min = 1e-320, min_sd=1e-10, tolerance = 1e-10):
    t = time.process_time()
    quant = (0.9999,0.9995,0.999,0.995,0.99,0.98,0.97,0.95)
    log_lik_true_init = np.zeros(shape=(d,max_iter,nb_simul), order="C")
    log_lik_out_init = np.zeros(shape=(d,max_iter,nb_simul,len(quant)), order="C")
    mu1_high = np.zeros(shape=(d-1,nb_simul), order="C")
    mu2_high = np.zeros(shape=(d-1,nb_simul), order="C")
    pi1_high = np.zeros(shape=(d-1,nb_simul), order="C")
    pi2_high = np.zeros(shape=(d-1,nb_simul), order="C")
    for k in range(nb_simul):
        for n in range(1,d):
            N_1 = int(pi1_0*10*(10**n))
            N_2 = int((1-pi1_0)*10*(10**n))
            sum_N = N_1+N_2
            pi2_0 = 1-pi1_0
            np.random.seed(k)
            data1 = np.random.poisson(lam=mu1_0, size=N_1)
            data2 = np.random.poisson(lam=mu2_0, size=N_2)
            data_total = np.hstack([data1,data2])
            
            ## Outliers initiation ##
            j_list = list()
            mu1_list = list()
            mu2_list = list()
            pi1_list = list()
            pi2_list = list()
            for b in range(len(quant)):
                quant_1 = np.quantile(data_total, q=quant[b])
                data_outlier  = data_total[np.where(data_total>=quant_1)]
                data_outlier1 = data_total[np.where(data_total<quant_1)]
                mu1 = np.array(np.mean(data_outlier1))
                mu2 = np.array(np.mean(data_outlier))
                if (math.isnan(mu2)): mu2 = np.min(data_total)
                dens1 = scipy.stats.poisson.pmf(data_total, mu=mu1, loc=0)
                dens2 = scipy.stats.poisson.pmf(data_total, mu=mu2, loc=0) 
                dens1[np.where(dens1<dens_min)]=dens_min
                dens2[np.where(dens2<dens_min)]=dens_min
                post_prob1 = (dens1)/(dens1+dens2)
                post_prob2 = 1-post_prob1
                pi1 = [np.shape(data_outlier)[0]/sum_N]
                pi2 = [np.shape(data_outlier1)[0]/sum_N]
                N1 = np.sum(post_prob1)
                N2 = np.sum(post_prob2)
                log_lik_out_init[n,0,k,b] = np.sum(np.log(pi1*dens1+pi2*dens2))
                
                ## Outliers initiation - Iterative algorithm ##
                for j in range(1,max_iter):
                    mu1 = np.hstack([mu1,np.average(data_total, weights=post_prob1)])
                    mu2 = np.hstack([mu2,np.average(data_total, weights=post_prob2)])
                    N1 = np.hstack([N1, np.sum(post_prob1)])
                    N2 = np.hstack([N2, np.sum(post_prob2)])
                    dens1 = scipy.stats.poisson.pmf(data_total, mu=mu1[j], loc=0)
                    dens2 = scipy.stats.poisson.pmf(data_total, mu=mu2[j], loc=0)
                    dens1[np.where(dens1<dens_min)]=dens_min
                    dens2[np.where(dens2<dens_min)]=dens_min
                    post_prob1 = (pi1[j-1]*dens1)/(pi1[j-1]*dens1+pi2[j-1]*dens2)
                    post_prob2 = 1-post_prob1
                    if (np.mean(post_prob2)==0.0):
                        j = j-1
                        break
                    if (np.mean(post_prob1)==0.0):
                        j = j-1
                        break
                    pi1.append(np.mean(post_prob1))
                    pi2.append(np.mean(post_prob2))
                    total_dens = pi1[j]*dens1+pi2[j]*dens2
                    log_lik_out_init[n,j,k,b] = np.sum(np.log(total_dens))
                    del total_dens, dens1, dens2
                    if j>1: #### Neglect the potential initial drop in the likelihood value (this is a caveat of the EM algorithm) ###
                        if ((log_lik_out_init[n,j,k,b]-log_lik_out_init[n,j-1,k,b])<tolerance): break        
                j_list.append(j-1)
                mu1_list.append(mu1[-2])
                mu2_list.append(mu2[-2])
                pi1_list.append(pi1[-2])
                pi2_list.append(pi2[-2])
                if b>0:
                    if (log_lik_out_init[n,j_list[b-1],k,b-1] > log_lik_out_init[n,j_list[b],k,b]): 
                        log_lik_out_init[n,0:,k,b] = log_lik_out_init[n,0:,k,b-1]
                        j_list[b] = j_list[b-1]
                        mu1_list[b] = mu1_list[b-1]
                        mu2_list[b] = mu2_list[b-1]
                        pi1_list[b] = pi1_list[b-1]
                        pi2_list[b] = pi2_list[b-1]
    
            ## True values initiation ###
            mu1_t, mu2_t = mu1_0, mu2_0
            dens1_t = scipy.stats.poisson.pmf(data_total, mu=mu1_t, loc=0)
            dens2_t = scipy.stats.poisson.pmf(data_total, mu=mu2_t, loc=0) 
            dens1_t[np.where(dens1_t<dens_min)]=dens_min
            dens2_t[np.where(dens2_t<dens_min)]=dens_min
            post_prob1_t = (pi1_0*dens1_t)/(pi1_0*dens1_t+pi2_0*dens2_t)
            post_prob2_t = 1-post_prob1_t
            pi1_t = [pi1_0]
            pi2_t = [pi2_0]
            N1_t = np.sum(post_prob1_t)
            N2_t = np.sum(post_prob2_t)
            log_lik_true_init[n,0,k] = np.sum(np.log(pi1_t*dens1_t+pi2_t*dens2_t))
            
            ## True values initiation - Iterative algorithm ###
            for m in range(1,max_iter):
                mu1_t = np.hstack([mu1_t,np.average(data_total, weights=post_prob1_t)])
                mu2_t = np.hstack([mu2_t,np.average(data_total, weights=post_prob2_t)])
                dens1_t = scipy.stats.poisson.pmf(data_total, mu=mu1_t[m], loc=0)
                dens2_t = scipy.stats.poisson.pmf(data_total, mu=mu2_t[m], loc=0)
                dens1_t[np.where(dens1_t<dens_min)]=dens_min
                dens2_t[np.where(dens2_t<dens_min)]=dens_min
                post_prob1_t = (pi1_t[m-1]*dens1_t)/(pi1_t[m-1]*dens1_t+pi2_t[m-1]*dens2_t)
                post_prob2_t = 1-post_prob1_t
                N1_t = np.hstack([N1_t, np.sum(post_prob1_t)])
                N2_t = np.hstack([N2_t, np.sum(post_prob2_t)])
                if (np.mean(post_prob2_t)==0.0): 
                    m = m-1
                    break
                if (np.mean(post_prob1_t)==0.0): 
                    m = m-1
                    break
                pi1_t.append(np.mean(post_prob1_t))
                pi2_t.append(np.mean(post_prob2_t))
                total_dens_t = pi1_t[m]*dens1_t+pi2_t[m]*dens2_t            
                log_lik_true_init[n,m,k] = np.sum(np.log(total_dens_t))
                del total_dens_t, dens1_t, dens2_t
                if m>1: #### Neglect the potential initial drop in the likelihood value (this is a caveat of the EM algorithm) ###
                    if ((log_lik_true_init[n,m,k]-log_lik_true_init[n,m-1,k])<tolerance): break
            
            if (log_lik_true_init[n,m-1,k] >= log_lik_out_init[n,j_list[b],k,b]):
                if mu1_t[m-1] > mu2_t[m-1]:
                    mu1_high[n-1,k] = mu1_t[m-1]
                    mu2_high[n-1,k] = mu2_t[m-1]
                    pi1_high[n-1,k] = pi1_t[m-1]
                    pi2_high[n-1,k] = pi2_t[m-1]
                else: 
                    mu1_high[n-1,k] = mu2_t[m-1]
                    mu2_high[n-1,k] = mu1_t[m-1]
                    pi1_high[n-1,k] = pi2_t[m-1]
                    pi2_high[n-1,k] = pi1_t[m-1]
            if (log_lik_true_init[n,m-1,k] < log_lik_out_init[n,j_list[b],k,b]):
                if mu1_list[b] > mu2_list[b]:
                    mu1_high[n-1,k] = mu1_list[b]
                    mu2_high[n-1,k] = mu2_list[b]
                    pi1_high[n-1,k] = pi1_list[b]
                    pi2_high[n-1,k] = pi2_list[b]
                else: 
                    mu1_high[n-1,k] = mu2_list[b]
                    mu2_high[n-1,k] = mu1_list[b]
                    pi1_high[n-1,k] = pi2_list[b]
                    pi2_high[n-1,k] = pi1_list[b]
    elapsed_time = time.process_time() - t
    return mu1_high, mu2_high, pi1_high, pi2_high, elapsed_time

#%%
d = 4 ## Maximum sample size, where N = 10**n, with n = 2,...,d
nb_simul = 1000 ## Number of replications
max_iter = 100 ## Maxmimum number of iterations
mu1_0 = 5  ## mu1_0 > mu2_0
mu2_0 = 3
pi1_0=0.5
results_poisson = mix_poisson_EM(nb_simul=nb_simul, d=d, max_iter=max_iter, mu1_0=mu1_0, mu2_0=mu2_0, pi1_0=pi1_0, dens_min = 1e-320, min_sd=1e-10, tolerance = 1e-10)
mu1_high = results_poisson[0]
mu2_high = results_poisson[1]
pi1_high = results_poisson[2]
pi2_high = results_poisson[3]
print(results_poisson[4]/60)

#%%
### Plot and save the results ####
## Mu_1 ##
plt.plot(np.mean(mu1_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(mu1_high, axis=1), yerr=1.96*np.std(mu1_high,ddof=1, axis=1))
plt.plot((d-1)*[mu1_0])
plt.title("Mu_1 average")
plt.show()
plt.hist(mu1_high[d-2,0:], bins=40)
plt.title("Mu_1 distribution")
plt.show()
np.savetxt('pois_mu1_high_5_3_1.txt', mu1_high, delimiter=",", fmt='%.12g')

## Mu_2 ##
plt.plot(np.mean(mu2_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(mu2_high, axis=1), yerr=1.96*np.std(mu2_high,ddof=1, axis=1))
plt.plot((d-1)*[mu2_0]) 
plt.title("Mu_2 average")
plt.show()
plt.hist(mu2_high[d-2,0:], bins=40)
plt.title("Mu_2 distribution")
plt.show()
np.savetxt('pois_mu2_high_5_3_1.txt', mu2_high, delimiter=",", fmt='%.12g') 

## Pi_1 ##
plt.plot(np.mean(pi1_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(pi1_high, axis=1), yerr=1.96*np.std(pi1_high,ddof=1, axis=1))
plt.plot((d-1)*[pi1_0])
plt.title("Pi_1 average") 
plt.show()
plt.hist(pi1_high[d-2,0:], bins=40)
plt.title("Pi_1 distribution")
plt.show()
np.savetxt('pois_pi1_high_5_3_1.txt', pi1_high, delimiter=",", fmt='%.12g') 

## Pi_2 ##
plt.plot(np.mean(pi2_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(pi2_high, axis=1), yerr=1.96*np.std(pi2_high,ddof=1, axis=1))
plt.plot((d-1)*[1-pi1_0])
plt.title("Pi_2 average")
plt.show()
plt.hist(pi2_high[d-2,0:], bins=40)
plt.title("Pi_2 distribution") 
plt.show()
np.savetxt('pois_pi2_high_5_3_1.txt', pi2_high, delimiter=",", fmt='%.12g') 

## Scatter plots ##
plt.scatter(mu1_high[d-2,0:],mu2_high[d-2,0:])
plt.xlabel("Estiamted Mu_1 value, N=10,000")
plt.ylabel("Estiamted Mu_2 value, N=10,000")
plt.show()
plt.scatter(pi1_high[d-2,0:],pi2_high[d-2,0:])
plt.xlabel("Estiamted Pi_1 value, N=10,000")
plt.ylabel("Estiamted Pi_2 value, N=10,000")
plt.show()

#%%
### Exponential distributions with G=2 ##
def mix_expo_EM(nb_simul, d, max_iter, mu1_0, mu2_0, pi1_0, dens_min = 1e-320, min_sd=1e-10, tolerance = 1e-10):
    t = time.process_time()
    quant = (0.9999,0.9995,0.999,0.995,0.99,0.98,0.97,0.95)
    log_lik_out_init = np.zeros(shape=(d,max_iter,nb_simul,len(quant)), order="C")
    log_lik_true_init = np.zeros(shape=(d,max_iter,nb_simul), order="C")
    mu1_high = np.zeros(shape=(d-1,nb_simul), order="C")
    mu2_high = np.zeros(shape=(d-1,nb_simul), order="C")
    pi1_high = np.zeros(shape=(d-1,nb_simul), order="C")
    pi2_high = np.zeros(shape=(d-1,nb_simul), order="C")
    
    for k in range(nb_simul):
        for n in range(1,d):
            N_1 = int(pi1_0*10*(10**n))
            N_2 = int((1-pi1_0)*10*(10**n))
            sum_N = N_1+N_2
            pi1_0 = N_1/sum_N 
            pi2_0 = N_2/sum_N
            np.random.seed(k)
            data1 = np.random.exponential(scale=mu1_0, size=N_1)
            data2 = np.random.exponential(scale=mu2_0, size=N_2)
            data_total = np.hstack([data1,data2])
            
            ## Outliers initiation ##
            j_list = list()
            mu1_list = list()
            mu2_list = list()
            pi1_list = list()
            pi2_list = list()
            for b in range(len(quant)):
                quant_1 = np.quantile(data_total, q=quant[b])
                data_outlier  = data_total[np.where(data_total>=quant_1)]
                data_outlier1 = data_total[np.where(data_total<quant_1)]
                mu1 = np.array(np.mean(data_outlier1))
                mu2 = np.array(np.mean(data_outlier))
                if (math.isnan(mu2)): mu2 = np.min(data_total)
                dens1 = scipy.stats.expon.pdf(data_total, scale=mu1, loc=0)
                dens2 = scipy.stats.expon.pdf(data_total, scale=mu2, loc=0) 
                dens1[np.where(dens1<dens_min)]=dens_min
                dens2[np.where(dens2<dens_min)]=dens_min
                post_prob1 = (dens1)/(dens1+dens2)
                post_prob2 = 1-post_prob1
                pi1 = [np.shape(data_outlier)[0]/sum_N]
                pi2 = [np.shape(data_outlier1)[0]/sum_N]
                N1 = np.sum(post_prob1)
                N2 = np.sum(post_prob2)
                log_lik_out_init[n,0,k,b] = np.sum(np.log(pi1*dens1+pi2*dens2))
                
                ## Outliers initiation - Iterative algorithm ##
                for j in range(1,max_iter):
                    mu1 = np.hstack([mu1,np.average(data_total, weights=post_prob1)])
                    mu2 = np.hstack([mu2,np.average(data_total, weights=post_prob2)])
                    dens1 = scipy.stats.expon.pdf(data_total, scale=mu1[j], loc=0)
                    dens2 = scipy.stats.expon.pdf(data_total, scale=mu2[j], loc=0)
                    dens1[np.where(dens1<dens_min)]=dens_min
                    dens2[np.where(dens2<dens_min)]=dens_min
                    post_prob1 = (pi1[j-1]*dens1)/(pi1[j-1]*dens1+pi2[j-1]*dens2)
                    post_prob2 = 1-post_prob1
                    if (np.mean(post_prob2)==0.0):
                        j = j-1
                        break
                    if (np.mean(post_prob1)==0.0):
                        j = j-1
                        break
                    N1 = np.hstack([N1, np.sum(post_prob1)])
                    N2 = np.hstack([N2, np.sum(post_prob2)])
                    pi1.append(np.mean(post_prob1))
                    pi2.append(np.mean(post_prob2))
                    total_dens = pi1[j]*dens1+pi2[j]*dens2
                    log_lik_out_init[n,j,k,b] = np.sum(np.log(total_dens))
                    del total_dens, dens1, dens2
                    if j>1: #### Neglect the potential initial drop in the likelihood value (this is a caveat of the EM algorithm) ###
                        if ((log_lik_out_init[n,j,k,b]-log_lik_out_init[n,j-1,k,b])<tolerance): break
    
                j_list.append(j-1)
                mu1_list.append(mu1[-2])
                mu2_list.append(mu2[-2])
                pi1_list.append(pi1[-2])
                pi2_list.append(pi2[-2])
                if (b>0):
                    if (log_lik_out_init[n,j_list[b-1],k,b-1] > log_lik_out_init[n,j_list[b],k,b]): 
                        log_lik_out_init[n,0:,k,b] = log_lik_out_init[n,0:,k,b-1]
                        j_list[b] = j_list[b-1]
                        mu1_list[b] = mu1_list[b-1]
                        mu2_list[b] = mu2_list[b-1]
                        pi1_list[b] = pi1_list[b-1]
                        pi2_list[b] = pi2_list[b-1]
    
            ## True values initiation ###
            mu1_t, mu2_t = mu1_0, mu2_0
            dens1_t = scipy.stats.expon.pdf(data_total, scale=mu1_t, loc=0)
            dens2_t = scipy.stats.expon.pdf(data_total, scale=mu2_t, loc=0) 
            dens1_t[np.where(dens1_t<dens_min)]=dens_min
            dens2_t[np.where(dens2_t<dens_min)]=dens_min
            post_prob1_t = (pi1_0*dens1_t)/(pi1_0*dens1_t+pi2_0*dens2_t)
            post_prob2_t = 1-post_prob1_t
            pi1_t = [pi1_0]
            pi2_t = [pi2_0]
            N1_t = np.sum(post_prob1_t)
            N2_t = np.sum(post_prob2_t)
            log_lik_true_init[n,0,k] = np.sum(np.log(pi1_t*dens1_t+pi2_t*dens2_t))
            
            ## True values initiation - Iterative algorithm ###
            for m in range(1,max_iter):
                mu1_t = np.hstack([mu1_t,np.average(data_total, weights=post_prob1_t)])
                mu2_t = np.hstack([mu2_t,np.average(data_total, weights=post_prob2_t)])
                dens1_t = scipy.stats.expon.pdf(data_total, scale=mu1_t[m], loc=0)
                dens2_t = scipy.stats.expon.pdf(data_total, scale=mu2_t[m], loc=0)
                dens1_t[np.where(dens1_t<dens_min)]=dens_min
                dens2_t[np.where(dens2_t<dens_min)]=dens_min
                post_prob1_t = (pi1_t[m-1]*dens1_t)/(pi1_t[m-1]*dens1_t+pi2_t[m-1]*dens2_t)
                post_prob2_t = 1-post_prob1_t
                if (np.mean(post_prob2_t)==0.0): 
                    m = m-1
                    break
                if (np.mean(post_prob1_t)==0.0): 
                    m = m-1
                    break
                N1_t = np.hstack([N1_t, np.sum(post_prob1_t)])
                N2_t = np.hstack([N2_t, np.sum(post_prob2_t)])
                pi1_t.append(np.mean(post_prob1_t))
                pi2_t.append(np.mean(post_prob2_t))
                total_dens_t = pi1_t[m]*dens1_t+pi2_t[m]*dens2_t            
                log_lik_true_init[n,m,k] = np.sum(np.log(total_dens_t))
                del total_dens_t, dens1_t, dens2_t
                if m>1: #### Neglect the potential initial drop in the likelihood value (this is a caveat of the EM algorithm) ###
                    if ((log_lik_true_init[n,m,k]-log_lik_true_init[n,m-1,k])<tolerance): break
            
            if (log_lik_true_init[n,m-1,k] >= log_lik_out_init[n,j_list[b],k,b]):
                if (mu1_t[m-1] > mu2_t[m-1]):
                    mu1_high[n-1,k] = mu1_t[m-1]
                    mu2_high[n-1,k] = mu2_t[m-1]
                    pi1_high[n-1,k] = pi1_t[m-1]
                    pi2_high[n-1,k] = pi2_t[m-1]
                else:
                    mu1_high[n-1,k] = mu2_t[m-1]
                    mu2_high[n-1,k] = mu1_t[m-1]
                    pi1_high[n-1,k] = pi2_t[m-1]
                    pi2_high[n-1,k] = pi1_t[m-1]
            if (log_lik_true_init[n,m-1,k] < log_lik_out_init[n,j_list[b],k,b]):
                if (mu1_list[b] > mu2_list[b]):
                    mu1_high[n-1,k] = mu1_list[b]
                    mu2_high[n-1,k] = mu2_list[b]
                    pi1_high[n-1,k] = pi1_list[b]
                    pi2_high[n-1,k] = pi2_list[b]
                else: 
                    mu1_high[n-1,k] = mu2_list[b]
                    mu2_high[n-1,k] = mu1_list[b]
                    pi1_high[n-1,k] = pi2_list[b]
                    pi2_high[n-1,k] = pi1_list[b]
    elapsed_time = time.process_time() - t
    return mu1_high, mu2_high, pi1_high, pi2_high, elapsed_time

#%%
d = 4 ## Maximum sample size, where N = 10**n, with n = 2,...,d
nb_simul = 1000 ## Number of replications
max_iter = 100 ## Maxmimum number of iterations
mu1_0 = 1  ## mu1_0 > mu2_0
mu2_0 = 0.6
pi1_0=0.5
results_expo = mix_expo_EM(nb_simul=nb_simul, d=d, max_iter=max_iter, mu1_0=mu1_0, mu2_0=mu2_0, pi1_0=pi1_0, dens_min = 1e-320, min_sd=1e-10, tolerance = 1e-10)
mu1_high = results_expo[0]
mu2_high = results_expo[1]
pi1_high = results_expo[2]
pi2_high = results_expo[3]
print(results_expo[4]/60)

#%%
### Plot and save the results ####
## Mu_1 ##
plt.plot(np.mean(mu1_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(mu1_high, axis=1), yerr=1.96*np.std(mu1_high,ddof=1, axis=1))
plt.plot((d-1)*[mu1_0])
plt.title("Mu_1 average")
plt.show()
plt.hist(mu1_high[d-2,0:], bins=40)
plt.title("Mu_1 distribution")
plt.show()
np.savetxt('expo_mu1_high_1_06_1.txt', mu1_high, delimiter=",", fmt='%.12g')

## Mu_2 ##
plt.plot(np.mean(mu2_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(mu2_high, axis=1), yerr=1.96*np.std(mu2_high,ddof=1, axis=1))
plt.plot((d-1)*[mu2_0]) 
plt.title("Mu_2 average")
plt.show()
plt.hist(mu2_high[d-2,0:], bins=40)
plt.title("Mu_2 distribution")
plt.show()
np.savetxt('expo_mu2_high_1_06_1.txt', mu2_high, delimiter=",", fmt='%.12g') 

## Pi_1 ##
plt.plot(np.mean(pi1_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(pi1_high, axis=1), yerr=1.96*np.std(pi1_high,ddof=1, axis=1))
plt.plot((d-1)*[pi1_0])
plt.title("Pi_1 average") 
plt.show()
plt.hist(pi1_high[d-2,0:], bins=40)
plt.title("Pi_1 distribution")
plt.show()
np.savetxt('expo_pi1_high_1_06_1.txt', pi1_high, delimiter=",", fmt='%.12g') 

## Pi_2 ##
plt.plot(np.mean(pi2_high, axis=1))
plt.errorbar(x=range(0,d-1), y=np.mean(pi2_high, axis=1), yerr=1.96*np.std(pi2_high,ddof=1, axis=1))
plt.plot((d-1)*[1-pi1_0])
plt.title("Pi_2 average") 
plt.show()
plt.hist(pi2_high[d-2,0:], bins=40)
plt.title("Pi_2 distribution") 
plt.show()
np.savetxt('expo_pi2_high_1_06_1.txt', pi2_high, delimiter=",", fmt='%.12g') 

## Scatter plots ##
plt.scatter(mu1_high[d-2,0:],mu2_high[d-2,0:])
plt.xlabel("Estiamted Mu_1 value, N=10,000")
plt.ylabel("Estiamted Mu_2 value, N=10,000")
plt.show()
plt.scatter(pi1_high[d-2,0:],pi2_high[d-2,0:])
plt.xlabel("Estiamted Pi_1 value, N=10,000")
plt.ylabel("Estiamted Pi_2 value, N=10,000")

plt.show()

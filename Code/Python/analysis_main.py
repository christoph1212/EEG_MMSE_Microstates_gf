
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 16:30:55 2022

@author: Jonas A. Thiele
"""

#%% Imports

import pandas as pd
import numpy as np
import pingouin as pg
from matplotlib import pyplot as plt
from sklearn import linear_model
import scipy.stats
from sklearn.model_selection import KFold
from scipy.ndimage import measurements
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
import pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import minmax_scale
from factor_analyzer import FactorAnalyzer

#%% Load data

beh_main = pd.read_csv('beh_data_main.csv', delimiter=';', encoding= 'unicode_escape')
beh_main = beh_main[beh_main['SQ_1'].notna()]
beh_main = beh_main[beh_main['APM_Anz_correct_JT'].notna()]

intell_main = beh_main['APM_Anz_correct_JT'].to_numpy() # RAPM scores
age_main = beh_main['SQ_1'].to_numpy() # age
subs_beh = beh_main['Teilnehmer'].to_numpy() # subject IDs corresponding to behavioral data

df_complexity = pd.read_pickle('./df_complexity.pkl') # load complexity measures  
subs_neuro = np.array(df_complexity.IDs_subjects_neuro) # subject IDs corresponding to complexity data
epochs_rejected = np.array(df_complexity.epochs_rejected) # epochs removed during preprocessing

# channel names
with open("channel_names", "rb") as fp:
    channel_names = pickle.load(fp)

#%% Sort behavioral data according to order of subjects in neuro data

# get indexes where subjects from behavioral data are found in subjects of neuro data to sort behavioral data according to order of neuro data  
idx_neuro = []
no_data_idx_neuro = []
for s in list(subs_neuro):

    idx = np.where(subs_beh == s)[0]
    if idx.shape[0] == 0:
        no_data_idx_neuro.append(np.where(subs_neuro==s))
    else:
        idx_neuro.append(np.where(subs_beh == s)[0])
    
# sort age and intelligence (RAPM) scores according to order of neuro data    
age_sorted = age_main[np.concatenate(idx_neuro)]
intel_sorted = intell_main[np.concatenate(idx_neuro)]


# visualize distribution of RAPM scores
plt.figure()
ax=sns.histplot(intel_sorted,bins=20, color = 'red', alpha = 0.5)
plt.ylabel('frequency', fontsize = '16') 
plt.xlabel('RAPM', fontsize = '16') 
plt.yticks(fontsize = '16')
plt.xticks(fontsize = '16')
plt.xlim(11,39)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('RAPM_main.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# visualize distribution of ages
plt.figure()
ax=sns.histplot(age_sorted,bins=20, color = 'red', alpha = 0.5, )
plt.ylabel('frequency', fontsize = '16') 
plt.xlabel('age', fontsize = '16') 
plt.yticks(fontsize = '16')
plt.xticks(fontsize = '16')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('age_main.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# visualize distribution of rejected epochs
plt.figure()
ax=sns.histplot(epochs_rejected, color = 'red', alpha = 0.5)
plt.ylabel('frequency', fontsize = '16') 
plt.xlabel('rejected epochs', fontsize = '16') 
plt.yticks(fontsize = '16')
plt.xticks(fontsize = '16')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('epochs_rejected_main.jpg', format='jpg', dpi = 1200, bbox_inches='tight')


#%% Factor analysis

def regress_confounds(feature, confounds):
    X = confounds
    y = feature
    regr = linear_model.LinearRegression()
    regr.fit(X, y)
    prediction = regr.predict(X)
    residuals = (y - prediction)
    residuals = scipy.stats.zscore(residuals)
    
    return residuals


df_vars_fa = df_complexity.copy()

# remove variables not needed for factor analysis
idx_ID = np.where(df_complexity.columns == 'IDs_subjects_neuro')[0] 
idx_epochs = np.where(df_complexity.columns == 'epochs_rejected')[0] 

n_cluster = 8
n_scale = 4
idx_mse_cluster = []
for ch in range(n_cluster):
    for sc in range(n_scale):
            
       idx_mse_cluster.append(np.where(df_complexity.columns == 'MSE_cluster_channel_'+ str(ch) + '_scales_' + str(sc))[0])
idx_mse_cluster = np.array(idx_mse_cluster).ravel()

idx_posHoc = []
idx_posHoc.append(np.where(df_complexity.columns == 'similarity_subject_micosates')[0])
idx_posHoc.append(np.where(df_complexity.columns == 'gev_subject')[0])
idx_posHoc.append(np.where(df_complexity.columns == 'gev_group')[0])
idx_posHoc = np.array(idx_posHoc).ravel()

idx_to_remove = np.sort(np.concatenate([idx_ID, idx_epochs, idx_mse_cluster, idx_posHoc]))[::-1]

df_vars_fa = df_vars_fa.drop(df_vars_fa.columns[idx_to_remove], axis = 1) 


feature_names = list(df_vars_fa.columns) # names of variables for factor analysis
mat_var = np.array(df_vars_fa) # all variables for factor analysis
confounds = np.vstack((age_sorted, epochs_rejected)).T


# control for confounds via linear regression
intelligence_reg = regress_confounds(intel_sorted, confounds)

vars_reg = []
for el in range(mat_var.shape[1]):
    
    res = regress_confounds(mat_var[:,el], confounds)
    vars_reg.append(res)

# save table with variables
mat_var_table = pd.DataFrame(np.array(vars_reg).T, columns = df_vars_fa.columns)
mat_var_table.to_csv ('complexity_vars_controlled.csv', index = False, header=True)

# perform factor analysis with 17 factors
X = np.array(vars_reg).T
fa = FactorAnalyzer(n_factors = 17) # number of factors estimated by parallel analysis
fa.fit(X)

# visualize result of factor analysis
vari = fa.get_factor_variance()
fig,a = plt.subplots()
plt.bar(range(17),vari[1],color = 'gray', edgecolor = 'black', alpha = 0.5)
plt.plot(vari[2],color = 'black')
plt.xticks(np.arange(17), np.arange(17)+1, fontsize = '11')
plt.xlabel('factor', fontsize='14')
plt.ylabel('explained variance', fontsize='14')
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.savefig('fa_main.jpg', format='jpg', dpi = 1200, bbox_inches='tight')
plt.close()

# get top 10 variables of the 17 factors
loadings = fa.loadings_
names_all = []
idx_sort_loading_all=[]
for n in range(17):
    
    idx_sort_loading = np.argsort(loadings[:,n])
    idx_sort_loading = idx_sort_loading[::-1]
    
    
    idx_sort_loading_all.append(idx_sort_loading)
    names_all.append(np.array(feature_names)[idx_sort_loading[0:10]])

# create table with the top 10 variables of each factor sorted by their loadings
columns = []
for i in range(17):
    columns.append('Factor ' + str(i))
    
df_fa = pd.DataFrame(np.array(names_all).T, columns=columns)   
df_fa.to_csv('factor_analysis.csv', sep=',')


#%% Compute associations between complexity values and RAPM scores 

def get_association(X, y, C):
    
    C = np.array(C)
    X = minmax_scale(X, feature_range=(0, 1)) # scaling between 0 and 1 for better comparision between mean and SD 
    data = np.vstack((X, y, C)).T
    
    if C.shape[0] == 2:
        df_data = pd.DataFrame(data, columns=['X','y','conf1','conf2'])
        res = pg.partial_corr(data=df_data, x='X', y='y', covar=['conf1','conf2']) # Pearson partial correlation
    
    elif C.shape[0] == 3:
        df_data = pd.DataFrame(data, columns=['X','y','conf1','conf2', 'conf3'])
        res = pg.partial_corr(data=df_data, x='X', y='y', covar=['conf1','conf2','conf3'])
    else: 
        raise ValueError('get_associations: not implemented for your number of confounds')
        
    return res, X
        
C = [age_sorted, epochs_rejected] # confounds   
y = intel_sorted # RAPM scores

# lists to save results in
results_p_raw = []
results_p_adj = []
results_corr = []
results_mean =[]
results_std = []
results_min = []
results_max = []
results_names = []

# association between number of GFP peaks and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == 'n_gfp_peaks')[0]
X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()

res, X = get_association(X, y, C)

p = res['p-val'].to_numpy()
corr = res['r'].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p) 
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append('n_gfp_peaks')

# associations between microstate measures and intelligence
n_MS = 5
n_measures = 4

#coverage
corr_ms = np.zeros((n_MS*n_measures))
p_ms = np.zeros((n_MS*n_measures))
cnt = 0
for k in range(n_MS):
           
    idx_cols_vars = np.where(df_complexity.columns == 'coverage_' + str(k))[0]
    X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)

    p_ms[cnt] = res['p-val'].to_numpy()
    corr_ms[cnt] = res['r'].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append('ms_' + 'coverage_' + str(k))
    cnt += 1
    
# lifespan
for k in range(n_MS):
           
    idx_cols_vars = np.where(df_complexity.columns == 'lifespan_' + str(k))[0]
    X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)

    p_ms[cnt] = res['p-val'].to_numpy()
    corr_ms[cnt] = res['r'].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append('ms_' + 'lifespan_' + str(k))
    cnt += 1
    
# lifespan GFP peaks
for k in range(n_MS):
           
    idx_cols_vars = np.where(df_complexity.columns == 'lifespan_peaks_' + str(k))[0]
    X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)

    p_ms[cnt] = res['p-val'].to_numpy()
    corr_ms[cnt] = res['r'].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append('ms_' + 'lifespan_peaks_' + str(k))
    cnt += 1
    
# frequency    
for k in range(n_MS):
           
    idx_cols_vars = np.where(df_complexity.columns == 'frequence_' + str(k))[0]
    X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)

    p_ms[cnt] = res['p-val'].to_numpy()
    corr_ms[cnt] = res['r'].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append('ms_' + 'frequence_' + str(k))
    cnt += 1
  
reject, p_val_adj = fdrcorrection(p_ms.flatten().T, alpha=0.05)    
results_corr.append(corr_ms.flatten().T)
results_p_raw.append(p_ms.flatten())
results_p_adj.append(p_val_adj)  


# associations between transition probabilities of microstates and RAPM scores
corr_ms1_ms2 = np.zeros((n_MS,n_MS))
p_ms1_ms2 = np.zeros((n_MS,n_MS))
for ms1 in range(n_MS):
    for ms2 in range(n_MS):
        
        idx_cols_vars = np.where(df_complexity.columns == 'transition_probability_'+ str(ms1) + '_' + str(ms2))[0]
        X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
        res, X = get_association(X, y, C)

        p_ms1_ms2[ms1,ms2] = res['p-val'].to_numpy()
        corr_ms1_ms2[ms1,ms2] = res['r'].to_numpy()
        results_mean.append(np.mean(X))
        results_std.append(np.std(X))
        results_min.append(np.min(X))
        results_max.append(np.max(X))
        results_names.append('ms_transition_' + str(ms1) + '_to_' + str(ms2))
        
reject, p_val_adj = fdrcorrection(p_ms1_ms2.flatten(), alpha=0.05)    
results_corr.append(corr_ms1_ms2.flatten())
results_p_raw.append(p_ms1_ms2.flatten())
results_p_adj.append(p_val_adj)             


# associations between transition probabilities of microstates GFP points only and RAPM scores
corr_ms1_ms2 = np.zeros((n_MS,n_MS))
p_ms1_ms2 = np.zeros((n_MS,n_MS))
for ms1 in range(n_MS):
    for ms2 in range(n_MS):
           
        idx_cols_vars = np.where(df_complexity.columns == 'transition_probability_peaks_'+ str(ms1) + '_' + str(ms2))[0]
        X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
        res, X = get_association(X, y, C)
        
        p_ms1_ms2[ms1,ms2] = res['p-val'].to_numpy()
        corr_ms1_ms2[ms1,ms2] = res['r'].to_numpy()
        results_mean.append(np.mean(X))
        results_std.append(np.std(X))
        results_min.append(np.min(X))
        results_max.append(np.max(X))
        results_names.append('ms_transition_peak_' + str(ms1) + '_to_' + str(ms2))
        

reject, p_val_adj = fdrcorrection(p_ms1_ms2.flatten(), alpha=0.05)    
results_corr.append(corr_ms1_ms2.flatten())
results_p_raw.append(p_ms1_ms2.flatten())
results_p_adj.append(p_val_adj)   


# association between Shannon entropy of GFP signals and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == 'shannon_gfp')[0]
X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
res, X = get_association(X, y, C)
p = res['p-val'].to_numpy()
corr= res['r'].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p) 
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append('shannon_gfp')

# association betweem Fuzzy entropy of GFP signals and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == 'fuzzy_gfp')[0]
X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
res, X = get_association(X, y, C)
p = res['p-val'].to_numpy()
corr= res['r'].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p) 
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append('fuzzy_gfp')


# assosiation between Shannon entropy of each channel and RAPM scores
corr_ch = np.zeros((28,1))
p_ch = np.zeros((28,1))
for ch in range(28):
           
    idx_cols_vars = np.where(df_complexity.columns == 'shannon_ch_' + str(channel_names[ch]))[0]
    X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)
    p_ch[ch] = res['p-val'].to_numpy()
    corr_ch[ch] = res['r'].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append('shannon_' + channel_names[ch])
    
reject, p_val_adj = fdrcorrection(p_ch.flatten(), alpha=0.05)    
results_corr.append(corr_ch.flatten())
results_p_raw.append(p_ch.flatten())
results_p_adj.append(p_val_adj) 


# associations between Fuzzy entropy of each channel and RAPM scores
corr_ch = np.zeros((28,1))
p_ch = np.zeros((28,1))
for ch in range(28):
           
    idx_cols_vars = np.where(df_complexity.columns == 'fuzzy_ch_' + str(channel_names[ch]))[0]
    X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)
    p_ch[ch] = res['p-val'].to_numpy()
    corr_ch[ch] = res['r'].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append('fuzzy_' + channel_names[ch])
    
reject, p_val_adj = fdrcorrection(p_ch.flatten(), alpha=0.05)    
results_corr.append(corr_ch.flatten())
results_p_raw.append(p_ch.flatten())
results_p_adj.append(p_val_adj)     

        
# associations between MSE and RAPM scores
corr_ch_sc = np.zeros((28,20))
p_ch_sc = np.zeros((28,20))
for ch in range(28):
    for sc in range(20):
            
        idx_cols_vars = np.where(df_complexity.columns == 'MSE_'+ channel_names[ch] + '_' + str(sc))[0]
        X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
        res, X = get_association(X, y, C)
        
        p_ch_sc[ch,sc] = res['p-val'].to_numpy()
        corr_ch_sc[ch, sc] = res['r'].to_numpy()

        results_mean.append(np.mean(X))
        results_std.append(np.std(X))
        results_min.append(np.min(X))
        results_max.append(np.max(X))
        results_names.append('mse_ch_' + channel_names[ch] + '_sc_' + str(sc))
        


results_corr.append(corr_ch_sc.flatten())
results_p_raw.append(p_ch_sc.flatten())

# measure clustersize of significant correlations
signi_ch_sc = np.zeros((28,20))
signi_ch_sc[np.where(p_ch_sc < 0.05)] = 1
lw_real, num_real = measurements.label(signi_ch_sc)
area_real = measurements.sum(signi_ch_sc, lw_real, index=np.arange(1,lw_real.max() + 1)) # clustersizes 

# permutation test (null models) for computing p-value of clustersizes in real data      
areas_all = []
n_permutations = 100
for p in range(n_permutations):
    
    print(p)
    corr_ch_sc = np.zeros((28,20))
    p_ch_sc = np.zeros((28,20))
    intel_sorted_perm = np.random.permutation(intel_sorted) # Permuted RAPM scores
    
    for ch in range(28):
        for sc in range(20):
            
            idx_cols_vars = np.where(df_complexity.columns == 'MSE_'+ channel_names[ch] + '_' + str(sc))[0]
            X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
            res, X = get_association(X, intel_sorted_perm, C)
                        
            p_ch_sc[ch,sc] = res['p-val'].to_numpy()
            corr_ch_sc[ch, sc] = res['r'].to_numpy()

    signi_ch_sc = np.zeros((28,20))
    signi_ch_sc[np.where(p_ch_sc < 0.05)] = 1
    lw, num = measurements.label(signi_ch_sc)
    area = measurements.sum(signi_ch_sc, lw, index=np.arange(1,lw.max() + 1))
    areas_all.append(area)

# cluster sizes of permutated data
size_cluster_perm = np.concatenate(areas_all).astype(int)
count_cluster_sizes = np.bincount(size_cluster_perm)

CI95 = np.percentile(size_cluster_perm, 95) # 95th percentile of cluster sizes of null models

# compute p-values for each cluster size of real data
p_cluster_size = []
for n in range(1,int(area_real.max())+1):
    permute = size_cluster_perm
    real = n
    sum_bad = np.sum(permute >= real) # instances of null models have larger cluster size than real data
    p = sum_bad/permute.shape[0]
    p_cluster_size.append(p)
    
p_cluster_size = np.array(p_cluster_size)

# assign adjusted p-values to clusters of real data
p_adj = np.ones(p_ch_sc.shape)    
for c in range(num_real):
    
    c_idx = c + 1
    area_c = area_real[c]
    p_c = p_cluster_size[int(area_c)-1]
    p_adj[np.where(lw_real==c_idx)]=p_c
    
p_val_adj = np.array(p_adj).flatten()
results_p_adj.append(p_val_adj)



# associations between clustered MSE and RAPM scores
n_cluster = 8
n_scale = 4
corr_ch_sc = np.zeros((n_cluster,n_scale))
p_ch_sc = np.zeros((n_cluster,n_scale))
for ch in range(n_cluster):
    for sc in range(n_scale):
            
        idx_cols_vars = np.where(df_complexity.columns == 'MSE_cluster_channel_'+ str(ch) + '_scales_' + str(sc))[0]
        X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
        res, X = get_association(X, y, C)
        
        p_ch_sc[ch,sc] = res['p-val'].to_numpy()
        corr_ch_sc[ch, sc] = res['r'].to_numpy()
        
        
        results_mean.append(np.mean(X))
        results_std.append(np.std(X))
        results_min.append(np.min(X))
        results_max.append(np.max(X))
        results_names.append('mse_clustered_cluster_' + str(ch) + '_' + str(sc))
        

results_corr.append(corr_ch_sc.flatten())
results_p_raw.append(p_ch_sc.flatten())
reject, p_val_adj = fdrcorrection(p_ch_sc.flatten(), alpha=0.05)    
results_p_adj.append(p_val_adj)  


# associations between MSE of GFP and RAPM scores
rho_sc = np.zeros((20))
p_sc = np.zeros((20))
for sc in range(20):
            
    idx_cols_vars = np.where(df_complexity.columns == 'MSE_gfp_scale' + str(sc))[0]
    X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
    res, X = get_association(X, y, C)
    
    p_sc[sc] = res['p-val'].to_numpy()
    rho_sc[sc] = res['r'].to_numpy()
    results_mean.append(np.mean(X))
    results_std.append(np.std(X))
    results_min.append(np.min(X))
    results_max.append(np.max(X))
    results_names.append('mse_gfp_scale_' + str(sc))
    
reject, p_val_adj = fdrcorrection(p_sc.flatten(), alpha=0.05)    
results_corr.append(rho_sc.flatten())
results_p_raw.append(p_sc.flatten())
results_p_adj.append(p_val_adj)     

# Post-hoc analyses
# association between similarity of individual (subject-specific) microstates and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == 'similarity_subject_micosates')[0]
X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
res, X = get_association(X, y, C)
p = res['p-val'].to_numpy()
corr = res['r'].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p) 
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append('similarity_subject_micosates')

# association between explained variance from individual (subject-specific) microstates and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == 'gev_subject')[0]
X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
res, X = get_association(X, y, C)
p = res['p-val'].to_numpy()
corr = res['r'].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p) 
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append('gev_subject')

# association between explained variance from group microstates and RAPM scores
idx_cols_vars = np.where(df_complexity.columns == 'gev_group')[0]
X = np.array(df_complexity.iloc[:,idx_cols_vars]).ravel()
res, X = get_association(X, y, C)
p = res['p-val'].to_numpy()
corr = res['r'].to_numpy()
results_p_raw.append(p)
results_p_adj.append(p) 
results_corr.append(corr)
results_mean.append(np.mean(X))
results_std.append(np.std(X))
results_min.append(np.min(X))
results_max.append(np.max(X))
results_names.append('gev_group')


# combine all results in data frame and save
res1 = np.concatenate(results_corr)
res2 = np.concatenate(results_p_raw)
res3 = np.concatenate(results_p_adj)
res4 = np.array(results_mean)
res5 = np.array(results_std)
res6 = np.array(results_min)
res7 = np.array(results_max)

data = np.vstack((res1, res2, res3, res4, res5, res6, res7))
df_results_main = pd.DataFrame(data, columns = np.array(results_names), index = ['corr','p','p_adj','M','SD', 'min','max'])

df_results_main.to_pickle('df_results_main')


#%% Prepare variables for prediction models

df_vars_prediction = df_complexity.copy()

# remove columns not needed for predicition
idx_mse = []
for ch in range(28):
    for sc in range(20):
            
        idx_mse.append(np.where(df_complexity.columns == 'MSE_'+ channel_names[ch] + '_' + str(sc))[0])
idx_mse = np.array(idx_mse).ravel()

idx_ID = np.where(df_complexity.columns == 'IDs_subjects_neuro')[0] 
idx_epochs = np.where(df_complexity.columns == 'epochs_rejected')[0] 

idx_posHoc = []
idx_posHoc.append(np.where(df_complexity.columns == 'similarity_subject_micosates')[0])
idx_posHoc.append(np.where(df_complexity.columns == 'gev_subject')[0])
idx_posHoc.append(np.where(df_complexity.columns == 'gev_group')[0])
idx_posHoc = np.array(idx_posHoc).ravel()

idx_to_remove = np.sort(np.concatenate([idx_mse, idx_ID, idx_epochs, idx_posHoc]))[::-1]

df_vars_prediction = df_vars_prediction.drop(df_vars_prediction.columns[idx_to_remove], axis = 1) 


feature_names = list(df_vars_prediction.columns) # names of variables for the prediction model
mat_var = np.array(df_vars_prediction) # all variables for the prediction model
confounds = np.vstack((age_sorted, epochs_rejected)).T

#%% Model 1 - Explained variance by predictors

def regress_confounds(feature, confounds):
    X = confounds
    y = feature
    regr = linear_model.LinearRegression()
    regr.fit(X, y)
    prediction = regr.predict(X)
    residuals = (y - prediction)
    residuals = scipy.stats.zscore(residuals)
    
    return residuals

# controlling for confounds (via linear regression)
intelligence_reg = regress_confounds(intel_sorted, confounds)
vars_reg = []
for el in range(mat_var.shape[1]):
    
    res = regress_confounds(mat_var[:,el], confounds)
    vars_reg.append(res)


X = np.array(vars_reg).T
y = intelligence_reg   


# select features by correlation with RAPM
p_f = []
r_f = []
for f in range(X.shape[1]):
    p_f.append(scipy.stats.pearsonr(X[:,f],y)[1])
    r_f.append(scipy.stats.pearsonr(X[:,f],y)[0])
    
ind_p = np.where(np.array(p_f) < 0.05) # get feature indexes with p<0.05

ind_pos = np.where(np.array(r_f) > 0) # indexes of positively correlating features
ind_neg = np.where(np.array(r_f) < 0) # indexes of negatively correlating features

ind_p_pos = np.intersect1d(ind_p, ind_pos)
ind_p_neg = np.intersect1d(ind_p, ind_neg)

X_sel_pos = X[:,ind_p_pos] # positively correlating variables (p<0.05)
X_sel_neg = X[:,ind_p_neg] # negatively correlating variables (p<0.05)
 
y_pred_pos = np.mean(X_sel_pos, axis = 1) # mean of positively correlating variables (p<0.05)
y_pred_neg = np.mean(X_sel_neg, axis = 1) # mean of negatively correlating variables (p<0.05)

y_pred = y_pred_pos - y_pred_neg 

print('----------------------------------------------')
print('Model 1 explained y vs y_pred --- corr, p: ')
print(scipy.stats.pearsonr(y_pred,y))
print('----------------------------------------------')

plt.figure('model main')
ax = sns.regplot(x=y, y=y_pred, ci = None, color = 'black', marker = '.')
plt.xlabel('observed intelligence score', fontsize = '16')
plt.ylabel('predictied intelligence score', fontsize = '16') 
plt.yticks(fontsize = '16')
plt.xticks(fontsize = '16')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('model_main.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# save features for replication
features_pos = np.array(feature_names)[ind_p_pos]
features_neg = np.array(feature_names)[ind_p_neg]
np.save('features_pos_main.npy', features_pos)  
np.save('features_neg_main.npy', features_neg)  


#%% Model 2: K-fold with permutation for testing statistical significance of predicition

def create_folds(df, n_s=5, n_grp=None):
    df['Fold'] = -1
    
    if n_grp is None:
        skf = KFold(n_splits=n_s)
        target = df.target
    else:
        skf = StratifiedKFold(n_splits=n_s, shuffle = True)
        df['grp'] = pd.cut(df.target, n_grp, labels=False)
        target = df.grp
    
    for fold_no, (t, v) in enumerate(skf.split(target, target)):
        df.loc[v, 'Fold'] = fold_no
    return df


def regress_confounds_fold(score, confounds, train_index, test_index):
    
    regr = linear_model.LinearRegression()
    regr.fit(confounds[train_index,:], score[train_index])
    
    fittedvalues = regr.predict(confounds)
    residuals = (score - np.ravel(fittedvalues))
    residuals = scipy.stats.zscore(residuals)
    
    residuals_train = residuals[train_index]
    residuals_test = residuals[test_index]
    
    return residuals_train, residuals_test

n_folds = 10
n1_permutations = 100 # number of "true" models
n2_permutations = 1000 # number of null models
model_corr_real = [] # store correlataions of models with true RAPM score order
model_corr_rand = [] # store correlations of models with permutated RAPM scores
y_pred_all_perms = []
y_best = []
y_pred_test = []

# true model
for p in range(n1_permutations):
    
    print(p)
    
    # create stratified folds
    df = pd.DataFrame(intelligence_reg, columns = ['target'])
    df = create_folds(df, n_s=n_folds, n_grp=30)
    
    
    y_pred_all = []
    y_test_all = []
    
    # k-fold
    for k in range(n_folds):
        
        
        test_index = np.where(df.Fold == k)[0]
        train_index = np.where(df.Fold != k)[0]
        
        confounds_fold = confounds
        X_fold = mat_var
        y_fold = intel_sorted
        

        # linear regression of confounds (regression model trained on training data and applied to test data respectively)
        y_train, y_test = regress_confounds_fold(score = y_fold, confounds = confounds_fold, 
                                                 train_index = train_index, test_index = test_index)

        X_train = []
        X_test = []
        for el in range(mat_var.shape[1]):
            
            var_train, var_test = regress_confounds_fold(score = X_fold[:,el], confounds = confounds_fold, 
                                                         train_index = train_index, test_index = test_index)
            X_train.append(var_train)
            X_test.append(var_test)

        X_train = np.array(X_train).T
        X_test = np.array(X_test).T
        y_test_all.append(y_test)

        # select features by correlation with RAPM in training sample
        p_f = []
        r_f = []
        for f in range(X_train.shape[1]):
            p_f.append(scipy.stats.pearsonr(X_train[:,f],y_train)[1])
            r_f.append(scipy.stats.pearsonr(X_train[:,f],y_train)[0])
        
        ind_p = np.where(np.array(p_f) < 0.05) # get feature indexes with p<0.05
        ind_pos = np.where(np.array(r_f) > 0) # indexes of positively correlating features
        ind_neg = np.where(np.array(r_f) < 0) # indexes of negatively correlating features
        
        ind_p_pos = np.intersect1d(ind_p, ind_pos) # indexes positively correlating variables (p<0.05)
        ind_p_neg = np.intersect1d(ind_p, ind_neg) # indexes negatively correlating variables (p<0.05)
        
        # Apply selected features to test sample
        # if only positively or negatively correlated feature were found, apply only these
        if ind_p_pos.shape[0] == 0:
            
            X_sel_test_neg = X_test[:,ind_p_neg] # negatively correlating variables (p<0.05)
            y_pred_neg  = np.mean(X_sel_test_neg, axis = 1)
            y_pred = - y_pred_neg
            
        elif ind_p_neg.shape[0] == 0: 
            
            X_sel_test_pos = X_test[:,ind_p_pos] # positively correlating variables (p<0.05)
            y_pred_pos = np.mean(X_sel_test_pos, axis = 1)
            y_pred = y_pred_pos 
            
        else:
            
            X_sel_test_neg = X_test[:,ind_p_neg] # negatively correlating variables (p<0.05)
            X_sel_test_pos = X_test[:,ind_p_pos] # positively correlating variables (p<0.05)
            y_pred_neg  = np.mean(X_sel_test_neg, axis = 1)
            y_pred_pos = np.mean(X_sel_test_pos, axis = 1)
            y_pred = y_pred_pos - y_pred_neg
    
       
        y_pred_all.append(y_pred)
        
    y_pred_all_perms.append(np.concatenate(y_pred_all))
    
    model_corr = scipy.stats.pearsonr(np.concatenate(y_pred_all),np.concatenate(y_test_all)) # model accuracy of predictions across folds
    model_corr_real.append(model_corr)
    
    # save result of best model
    if model_corr[0] == np.array(model_corr_real)[:,0].max():
        
        y_best = np.concatenate(y_test_all)
        y_pred_best = np.concatenate(y_pred_all)
        
# visualize best model        
model_corr = scipy.stats.pearsonr(y_best, y_pred_best)    
plt.figure('model 10-k fold example')
ax = sns.regplot(x=y_best, y=y_pred_best, ci = None, color = 'black', marker = '.')
plt.xlabel('true intelligence score', fontsize = '16')
plt.ylabel('predictied intelligence score', fontsize = '16') 
plt.yticks(fontsize = '16')
plt.xticks(fontsize = '16')
plt.title('k-fold example - corr = ' + str(model_corr))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('model 10-k fold best.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

#%% Null models  
 
for p in range(n2_permutations):
    print(p)
    
    y_p = np.random.permutation(intel_sorted) # permutated RAPM scores
    y_p_res = regress_confounds(y_p, confounds)
    
    # create folds
    df = pd.DataFrame(y_p_res, columns = ['target'])
    df = create_folds(df, n_s=n_folds, n_grp=30)
    
    
    # k-fold
    y_pred_all = []
    y_test_all = []
    invalid = 0 # if no features are selected (no features correlating with permutated RAPM scores)
    
    for k in range(n_folds):
        
        test_index = np.where(df.Fold == k)[0]
        train_index = np.where(df.Fold != k)[0]
        
        
        confounds_fold = confounds
        X_fold = np.array(mat_var)
        y_fold = y_p
        
        # linear regression of confounds (regression model trained on training data and applied to test data respectively)
        y_train, y_test = regress_confounds_fold(score = y_fold, confounds = confounds_fold, 
                                                 train_index = train_index, test_index = test_index)

        X_train = []
        X_test = []
        for el in range(mat_var.shape[1]):
            
            var_train, var_test = regress_confounds_fold(score = X_fold[:,el], confounds = confounds_fold, 
                                                         train_index = train_index, test_index = test_index)
            X_train.append(var_train)
            X_test.append(var_test)

        X_train = np.array(X_train).T
        X_test = np.array(X_test).T
        
        y_test_all.append(y_test)
    
        # select features by correlation with RAPM in training sample
        p_f = []
        r_f = []
        for f in range(X_train.shape[1]):
            p_f.append(scipy.stats.pearsonr(X_train[:,f],y_train)[1])
            r_f.append(scipy.stats.pearsonr(X_train[:,f],y_train)[0])
        
        ind_p = np.where(np.array(p_f) < 0.05)  # get feature indexes with p<0.05
        
        ind_pos = np.where(np.array(r_f) > 0) # indexes of positively correlating features
        ind_neg = np.where(np.array(r_f) < 0) # indexes of negatively correlating features
        
        ind_p_pos = np.intersect1d(ind_p, ind_pos) # indexes positively correlating variables (p<0.05)
        ind_p_neg = np.intersect1d(ind_p, ind_neg) # indexes negatively correlating variables (p<0.05)
        
        # apply selected features in test sample
        # if only positively or negatively correlated feature were found, apply only these
        # if no features are correlated at all, set accuracy (correlation(y, y_pred)) of the model to zero
    
        if ind_p_pos.shape[0] + ind_p_neg.shape[0] == 0:
            
            invalid = 1 # if no correlation found
        
        elif ind_p_pos.shape[0] == 0:
            
            X_sel_test_neg = X_test[:,ind_p_neg] # negatively correlating variables (p<0.05)
            y_pred_neg  = np.mean(X_sel_test_neg, axis = 1)
            y_pred = - y_pred_neg
            
        elif ind_p_neg.shape[0] == 0:
           
            X_sel_test_pos = X_test[:,ind_p_pos] # positively correlating variables (p<0.05)
            y_pred_pos = np.mean(X_sel_test_pos, axis = 1)
            y_pred = y_pred_pos 
          
        else:
            
            X_sel_test_neg = X_test[:,ind_p_neg] # negatively correlating variables (p<0.05)
            X_sel_test_pos = X_test[:,ind_p_pos] # positively correlating variables (p<0.05)
            y_pred_neg  = np.mean(X_sel_test_neg, axis = 1)
            y_pred_pos = np.mean(X_sel_test_pos, axis = 1)
            y_pred = y_pred_pos - y_pred_neg
        
       
        y_pred_all.append(y_pred)
     
        
    if invalid == 1:
        model_corr = 0 # if no features are correlated at all, set accuracy (correlation(y, y_pred)) to zero
    else:
        model_corr = scipy.stats.pearsonr(np.concatenate(y_pred_all),np.concatenate(y_test_all))[0] # model accuracy of predictions across folds
    
    model_corr_rand.append(model_corr)
    

# 95th percentile of model performance values of null models    
CI95 = np.percentile(np.array(model_corr_rand), 95)

# mean value of performances of true models
real_mean = np.mean(np.array(model_corr_real)[:,0])

# p-value for k-fold model 
permute = np.array(model_corr_rand)
real = real_mean
sum_bad = np.sum(permute >= real) # sum of permutations where null model outperforms real model
p_real = sum_bad/permute.shape[0]

# visualize histogram of null model accuracies vs mean of "true" models
plt.figure()
ax=sns.histplot(np.array(model_corr_rand), color = 'red', alpha = 0.5)
plt.axvline(real_mean)
plt.ylabel('frequency', fontsize = '16') 
plt.xlabel('r', fontsize = '16') 
plt.yticks(fontsize = '16')
plt.xticks(fontsize = '16')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.locator_params(axis='x', nbins=5)
plt.savefig('hist_main_permutation_test.jpg', format='jpg', dpi = 1200, bbox_inches='tight')


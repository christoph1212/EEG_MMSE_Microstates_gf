# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 16:00:11 2022

@author: Jonas A. Thiele
"""
# Plot results main and replication sample

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import mne
from matplotlib import ticker
import pickle


df_main = pd.read_pickle('df_results_main') # associations main sample
df_repli = pd.read_pickle('df_results_repli') # associations replication sample

column_names_main = np.array(df_main.columns)
column_names_repli = np.array(df_repli.columns)
check = (column_names_main[0:-3]==column_names_repli).all() # last three columns of column_names_main are post-hoc analyses done in main sample only
print(check) # check if column names are identical between main and replication sample
column_names = column_names_repli

with open("channel_names", "rb") as fp:
    channel_names = pickle.load(fp)

raw = mne.io.read_raw_fif('info_raw.fif') # raw for topomap plotting

maps_group = np.load('microstates_group.npy').T # group microstates

#%% Associations between single complexity measures and intelligence in main and replication sample

# Shannon entropy and intelligence - main and replication sample
idx_vars = []
for ch in channel_names:
    var_name = 'shannon_' + ch 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)    
    
fig,a = plt.subplots()
plt.bar(range(idx_vars.size), np.array(df_main)[0,idx_vars], color = 'red', edgecolor = 'black', alpha=0.5)
plt.bar(range(idx_vars.size), np.array(df_repli)[0,idx_vars], color = 'gray', edgecolor = 'black', alpha=0.4)

plt.ylim([-0.35, 0.15])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '11')
plt.xticks(np.arange(idx_vars.size), np.array(channel_names), fontsize = '11', rotation = 90)
#plt.xlabel('Channel')
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('shannon_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')



# Fuzzy entropy and intelligence - main and replication sample
idx_vars = []
for ch in channel_names:
    var_name = 'fuzzy_' + ch 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)    
    
fig,a = plt.subplots()
plt.bar(range(idx_vars.size), np.array(df_main)[0,idx_vars], color = 'red', edgecolor = 'black', alpha=0.5)
plt.bar(range(idx_vars.size), np.array(df_repli)[0,idx_vars], color = 'gray', edgecolor = 'black', alpha=0.4)
plt.ylim([-0.35, 0.15])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '11')
plt.xticks(np.arange(idx_vars.size), np.array(channel_names), fontsize = '11', rotation = 90)
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('fuzzy_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')



# MSE of GFP and intelligence - main and replication sample
idx_vars = []
for sc in range(20):
    var_name = 'mse_gfp_scale_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)    
    
fig,a = plt.subplots()
plt.bar(range(idx_vars.size), np.array(df_main)[0,idx_vars], color = 'red', edgecolor = 'black', alpha = 0.5)
plt.bar(range(idx_vars.size), np.array(df_repli)[0,idx_vars], color = 'gray', edgecolor = 'black', alpha = 0.4)
plt.ylim([-0.25, 0.25])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '11')
plt.xticks([0,4,9,14,19], [1,5,10,15,20], fontsize = '11')
plt.xlabel('scale')
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('mse_gfp_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# MSE and intelligence - main sample
idx_vars = []
for ch in range(28):
    for sc in range(20):
            var_name = 'mse_ch_' + channel_names[ch] + '_sc_' + str(sc)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

mse = np.array(df_main)[0,idx_vars]
mse = np.reshape(mse,(28,20))  
    
fig,a = plt.subplots()
divnorm=colors.TwoSlopeNorm(vmin=-0.31, vcenter=0, vmax=0.36)
plt.imshow(mse, cmap="bwr", norm=divnorm)
plt.colorbar()
plt.yticks(np.arange(mse.shape[0]), np.array(channel_names), fontsize = '9')
plt.xticks([0,4,9,14,19], [1,5,10,15,20], fontsize = '9')
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.savefig('mse_main_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# plot map of significant cluster MSE main sample
p_adj = np.array(df_main)[2,idx_vars]
p_adj = np.reshape(p_adj,(28,20))  
p_sig = np.zeros((28,20))
idx, idy = np.where(p_adj<0.05)
p_sig[idx,idy]=1    
fig,a = plt.subplots()
plt.imshow(p_sig, cmap="bwr", norm=divnorm)

# MSE and intelligence - replication sample
idx_vars = []
for ch in range(28):
    for sc in range(20):
            var_name = 'mse_ch_' + channel_names[ch] + '_sc_' + str(sc)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

mse = np.array(df_repli)[0,idx_vars]
mse = np.reshape(mse,(28,20))  
    
fig,a = plt.subplots()
divnorm=colors.TwoSlopeNorm(vmin=-0.31, vcenter=0, vmax=0.36)
plt.imshow(mse, cmap="bwr", norm=divnorm)
plt.colorbar()
plt.yticks(np.arange(mse.shape[0]), np.array(channel_names), fontsize = '9')
plt.xticks([0,4,9,14,19], [1,5,10,15,20], fontsize = '9')
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.savefig('mse_repli_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# plot map of significant clusters MSE replication sample
idx_vars = []
for ch in range(28):
    for sc in range(20):
            var_name = 'mse_ch_' + channel_names[ch] + '_sc_' + str(sc)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

p_adj = np.array(df_repli)[2,idx_vars]
p_adj = np.reshape(p_adj,(28,20))  
p_sig = np.zeros((28,20))
idx, idy = np.where(p_adj<0.05)
p_sig[idx,idy]=1    
fig,a = plt.subplots()
plt.imshow(p_sig, cmap="bwr", norm=divnorm)



# Microstates
idx_order_ms = [3,0,4,2,1] # microstates to order A,B,C,D,E
ms_maps = maps_group[:,idx_order_ms] # microstates with right order

# plot microstate topomaps
fig, axes = plt.subplots(nrows=1, ncols=5)
cnt=0
for axes_row in axes:
    print(axes_row)

    mne.viz.plot_topomap(ms_maps[:,cnt], raw.info, axes=axes_row, show=False)
    axes_row.spines['right'].set_visible(False)
    axes_row.spines['top'].set_visible(False)
    #ax.set_title('%s %s' % (ch_type.upper(), extr), fontsize=14)
    cnt += 1
fig.tight_layout()
fig.savefig('microstates_topomaps.jpg',format='jpg', dpi = 1200, bbox_inches='tight')


# coverage and intelligence - main and replication sample
idx_vars = []
for sc in range(5):
    var_name = 'ms_coverage_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)   

fig,a = plt.subplots()
ms_vars = np.array(df_main)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
ms_vars = np.array(df_repli)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'gray', edgecolor = 'black', alpha = 0.4)
plt.ylim([-0.4, 0.35])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xlabel('microstates')
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('coverage_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# lifespan and intelligence - main and replication sample
idx_vars = []
for sc in range(5):
    var_name = 'ms_lifespan_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)   

fig,a = plt.subplots()
ms_vars = np.array(df_main)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
ms_vars = np.array(df_repli)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'gray', edgecolor = 'black', alpha = 0.4)
plt.ylim([-0.4, 0.35])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xlabel('microstates')
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('lifespan_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# lifespan at GFP peaks and intelligence - main and replication sample
idx_vars = []
for sc in range(5):
    var_name = 'ms_lifespan_peaks_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)   

fig,a = plt.subplots()
ms_vars = np.array(df_main)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
ms_vars = np.array(df_repli)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'gray', edgecolor = 'black', alpha = 0.4)
plt.ylim([-0.4, 0.35])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xlabel('microstates')
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('lifespan_peaks_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# frequency and intelligence - main and replication sample
idx_vars = []
for sc in range(5):
    var_name = 'ms_frequence_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)   

fig,a = plt.subplots()
ms_vars = np.array(df_main)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
ms_vars = np.array(df_repli)[0,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'gray', edgecolor = 'black', alpha = 0.4)
plt.ylim([-0.4, 0.35])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xlabel('microstates')
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('frequency_RAPM.jpg', format='jpg', dpi = 1200, bbox_inches='tight')


# transition probabilities and intelligence - main sample
idx_vars = []
for ch in range(5):
    for sc in range(5):
            ch_sort = idx_order_ms[ch]
            sc_sort = idx_order_ms[sc]
            var_name = 'ms_transition_' + str(ch_sort) + '_to_' + str(sc_sort)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

ms_trans = np.array(df_main)[0,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,ax = plt.subplots()
divnorm=colors.TwoSlopeNorm(vmin=-0.41, vcenter=0, vmax=0.41)
plt.imshow(ms_trans, cmap="bwr", norm=divnorm)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')

cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('transitions_RAPM_main.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# transition probabilities at GFP peaks and intelligence - main sample
idx_vars = []
for ch in range(5):
    for sc in range(5):
            ch_sort = idx_order_ms[ch]
            sc_sort = idx_order_ms[sc]
            var_name = 'ms_transition_peak_' + str(ch_sort) + '_to_' + str(sc_sort)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

ms_trans = np.array(df_main)[0,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,ax = plt.subplots()
divnorm=colors.TwoSlopeNorm(vmin=-0.41, vcenter=0, vmax=0.41)
plt.imshow(ms_trans, cmap="bwr", norm=divnorm)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')

cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.savefig('transitions_peak_RAPM_main.jpg', format='jpg', dpi = 1200, bbox_inches='tight')


# transition probabilities and intelligence - replication sample
idx_vars = []
for ch in range(5):
    for sc in range(5):
            ch_sort = idx_order_ms[ch]
            sc_sort = idx_order_ms[sc]
            var_name = 'ms_transition_' + str(ch_sort) + '_to_' + str(sc_sort)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

ms_trans = np.array(df_repli)[0,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,a = plt.subplots()
divnorm=colors.TwoSlopeNorm(vmin=-0.41, vcenter=0, vmax=0.41)
plt.imshow(ms_trans, cmap="bwr", norm=divnorm)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')

cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
a.spines['right'].set_visible(False)
a.spines['top'].set_visible(False)
plt.savefig('transitions_RAPM_repli.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# transition probabilities at GFP peaks and intelligence - replication sample
idx_vars = []
for ch in range(5):
    for sc in range(5):
            ch_sort = idx_order_ms[ch]
            sc_sort = idx_order_ms[sc]
            var_name = 'ms_transition_peak_' + str(ch_sort) + '_to_' + str(sc_sort)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

ms_trans = np.array(df_repli)[0,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,ax = plt.subplots()
divnorm=colors.TwoSlopeNorm(vmin=-0.41, vcenter=0, vmax=0.41)
plt.imshow(ms_trans, cmap="bwr", norm=divnorm)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '16')

cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('transitions_peak_RAPM_repli.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# Find significant associations (after correcting for multiple comparisons)
signi_main = column_names_main[np.where(np.array(df_main)[2,:]<0.05)]
signi_repli = column_names_repli[np.where(np.array(df_repli)[2,:]<0.05)]

#%% Means and standard deviations of measures

# Shannon entropy
idx_vars = []
for ch in channel_names:
    var_name = 'shannon_' + ch 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)    

fig,(ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(np.array(df_main)[3,idx_vars], raw.info, axes=ax1, show=False, vmin = 0.1, vmax = 0.9, contours=12)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt.savefig('shannon_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

fig,(ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(np.array(df_main)[4,idx_vars], raw.info, axes=ax1, show=False, vmin = 0.1, vmax = 0.25)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt.savefig('shannon_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')



# Fuzzy entropy
idx_vars = []
for ch in channel_names:
    var_name = 'fuzzy_' + ch 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)    
    
fig,(ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(np.array(df_main)[3,idx_vars], raw.info, axes=ax1, show=False, vmin = 0.1, vmax = 0.9)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt.savefig('fuzzy_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

fig,(ax1) = plt.subplots(ncols=1)
im, cm = mne.viz.plot_topomap(np.array(df_main)[4,idx_vars], raw.info, axes=ax1, show=False, vmin = 0.1, vmax = 0.25)
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=20)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt.savefig('fuzzy_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')



# MSE of GFP
idx_vars = []
for sc in range(20):
    var_name = 'mse_gfp_scale_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)    
    
fig,ax = plt.subplots()
plt.bar(range(idx_vars.size), np.array(df_main)[3,idx_vars], color = 'red', edgecolor = 'black', alpha = 0.5)
#plt.ylim([-0.25, 0.1])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks([0,4,9,14,19], [1,5,10,15,20], fontsize = '20')
plt.xlabel('scale', fontsize = '20')
plt.ylabel('mean MSE',fontsize = '20')
plt.ylim((0,0.9))
plt.locator_params(axis='y', nbins=5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('mse_gfp_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

fig,ax = plt.subplots()
plt.bar(range(idx_vars.size), np.array(df_main)[4,idx_vars], color = 'red', edgecolor = 'black', alpha = 0.5)
#plt.ylim([-0.25, 0.1])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks([0,4,9,14,19], [1,5,10,15,20], fontsize = '20')
plt.xlabel('scale', fontsize = '20')
plt.ylabel('SD MSE', fontsize = '20')
plt.ylim((0,0.25))
plt.locator_params(axis='y', nbins=4)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('mse_gfp_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')



# MSE
idx_vars = []
for ch in range(28):
    for sc in range(20):
            var_name = 'mse_ch_' + channel_names[ch] + '_sc_' + str(sc)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

var = np.array(df_main)[3,idx_vars]
var = np.reshape(var,(28,20))  
  
scale = np.array([0,4,9,14,19])
fig,(ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5)
vmin = 0.1
vmax = 0.9
im, cm = mne.viz.plot_topomap(var[:,scale[0]], raw.info, axes=ax1, show=False,vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[1]], raw.info, axes=ax2, show=False, vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[2]], raw.info, axes=ax3, show=False, vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[3]], raw.info, axes=ax4, show=False, vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[4]], raw.info, axes=ax5, show=False, vmin=vmin,vmax=vmax)
cbar_ax = fig.add_axes([0.21, 0.35, 0.6, 0.05])
cb = fig.colorbar(im, cax=cbar_ax, orientation = "horizontal")
cb.ax.tick_params(labelsize=12)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
plt.savefig('MSE_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')


var = np.array(df_main)[4,idx_vars]
var = np.reshape(var,(28,20))  
  
scale = np.array([0,4,9,14,19])
fig,(ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5)
vmin = 0.1
vmax = 0.25
im, cm = mne.viz.plot_topomap(var[:,scale[0]], raw.info, axes=ax1, show=False,vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[1]], raw.info, axes=ax2, show=False, vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[2]], raw.info, axes=ax3, show=False, vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[3]], raw.info, axes=ax4, show=False, vmin=vmin,vmax=vmax)
im, cm = mne.viz.plot_topomap(var[:,scale[4]], raw.info, axes=ax5, show=False, vmin=vmin,vmax=vmax)
cbar_ax = fig.add_axes([0.21, 0.35, 0.6, 0.05])
cb = fig.colorbar(im, cax=cbar_ax, orientation = "horizontal")
cb.ax.tick_params(labelsize=12)
tick_locator = ticker.MaxNLocator(nbins=7)
cb.locator = tick_locator
cb.update_ticks()
plt.savefig('MSE_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')



# Microstates
idx_order_ms = [3,0,4,2,1] # microstates to order A,B,C,D,E

# Coverage
idx_vars = []
for sc in range(5):
    var_name = 'ms_coverage_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)  
ms_vars = np.array(df_main)[3,idx_vars]
ms_vars = ms_vars[idx_order_ms]

fig,ax = plt.subplots()
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.9])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('coverage_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

ms_vars = np.array(df_main)[4,idx_vars]
ms_vars = ms_vars[idx_order_ms]

fig,ax = plt.subplots()
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.25])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=4)
plt.savefig('coverage_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# Lifespan
idx_vars = []
for sc in range(5):
    var_name = 'ms_lifespan_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)   
ms_vars = np.array(df_main)[3,idx_vars]
ms_vars = ms_vars[idx_order_ms]

fig,ax = plt.subplots()
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.9])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('lifespan_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

fig,ax = plt.subplots()
ms_vars = np.array(df_main)[4,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.25])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=4)
plt.savefig('lifespan_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# Lifespan GFP peaks
idx_vars = []
for sc in range(5):
    var_name = 'ms_lifespan_peaks_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)   
ms_vars = np.array(df_main)[3,idx_vars]
ms_vars = ms_vars[idx_order_ms]

fig,ax = plt.subplots()
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.9])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('lifespan_peak_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

fig,ax = plt.subplots()
ms_vars = np.array(df_main)[4,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.25])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=4)
plt.savefig('lifespan_peak_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# Frequency
idx_vars = []
for sc in range(5):
    var_name = 'ms_frequence_' + str(sc) 
    idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)   
ms_vars = np.array(df_main)[3,idx_vars]
ms_vars = ms_vars[idx_order_ms]


fig,ax = plt.subplots()
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.9])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=5)
plt.savefig('frequency_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

fig,ax = plt.subplots()
ms_vars = np.array(df_main)[4,idx_vars]
ms_vars = ms_vars[idx_order_ms]
plt.bar(range(idx_vars.size), ms_vars, color = 'red', edgecolor = 'black', alpha = 0.5)
plt.ylim([0, 0.25])
plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.yticks(fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xlabel('microstates')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.locator_params(axis='y', nbins=4)
plt.savefig('frequency_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')


# transition probabilities all time points
idx_vars = []
for ch in range(5):
    for sc in range(5):
            ch_sort = idx_order_ms[ch]
            sc_sort = idx_order_ms[sc]
            var_name = 'ms_transition_' + str(ch_sort) + '_to_' + str(sc_sort)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

ms_trans = np.array(df_main)[3,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,ax = plt.subplots()
plt.imshow(ms_trans, cmap="Reds", vmin=0.1, vmax=0.9)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('transitions_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

ms_trans = np.array(df_main)[4,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,ax = plt.subplots()
plt.imshow(ms_trans, cmap="Reds", vmin=0.099, vmax=0.25)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('transitions_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

# transition probabilities GFP peaks 
idx_vars = []
for ch in range(5):
    for sc in range(5):
            ch_sort = idx_order_ms[ch]
            sc_sort = idx_order_ms[sc]
            var_name = 'ms_transition_peak_' + str(ch_sort) + '_to_' + str(sc_sort)
            idx_vars.append(np.where(column_names==var_name)[0][0])
idx_vars = np.array(idx_vars)

ms_trans = np.array(df_main)[3,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,ax = plt.subplots()
plt.imshow(ms_trans, cmap="Reds", vmin=0.1, vmax=0.9)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('transitions_peak_mean.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

ms_trans = np.array(df_main)[4,idx_vars]
ms_trans = np.reshape(ms_trans,(5,5))  
    
fig,ax = plt.subplots()
plt.imshow(ms_trans, cmap="Reds", vmin=0.099, vmax=0.25)

plt.yticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
plt.xticks(np.arange(5), ['A','B','C','D','E'], fontsize = '20')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=16)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.grid(color='gray', linestyle='dashed', linewidth=0.5)
plt.savefig('transitions_peak_SD.jpg', format='jpg', dpi = 1200, bbox_inches='tight')

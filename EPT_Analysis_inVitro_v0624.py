# -*- coding: utf-8 -*-
"""
Created on Fri May 24 17:24:44 2024

@author: Jessica A. Martinez
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import skimage.morphology as morpho
import toml
import tkinter as tk
from tkinter import filedialog
import sys
import time
import subprocess
import warnings
from scipy.ndimage import gaussian_filter

warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN slice encountered")

sys.stdout.flush()
np.seterr(invalid='ignore')

plt.close('all')

os.chdir('Z:\Documents\Projects\MR_Safety\Subject_Specific_SAR_maps_Project\Code')

import EPTMods2 as ept
import Mods as mods

# % Write EPT settings
###################################

sg_shape = "Sphere"#"Cube"#"Cross"#
sg_degree = 2
bc_epsr = 1.0
sg_size = [3,3,3]
bc_sigma = 0.6
ad_coeff = 0.0001
volume_tomography = True

method = "Helmholtz EPT"
# method = "Convection-reaction EPT"
#######################################
#%
# 1. Load H5 file usually store in the forEPT folder 

address=ept.select_HDF5file() 
ad=os.path.dirname(address)
file_name = os.path.splitext(os.path.basename(address))[0]
mask_loop=[]
with h5py.File(address, 'r') as f:
    for name in f.keys():
        if name=='mask_0' or name=='mask_1' or name=='mask_2':
            mask_loop.append(name)
        data = f[name][()]
        globals()[name] = data
if not mask_loop:
    mask_loop.append('mask')
    mask=mods.tissue_extractor(mag, 150, solidity_thresh=0.2)
    mask[mask==1]=np.nan
    mask[mask==0]=1
#%%
# 1.1 Open H5 file
step = [2e-3, 2e-3, 2e-3]
freq = 64e6
# mask=np.ones(np.shape(sigma))

M=mask#mods.tissue_extractor(mag, 150, solidity_thresh=0.2)
trx_phase=(ph_3)

#freq=float(int(freq))
step=np.round(step,3).tolist()


sigma_values_mask=np.zeros((len(mask_loop),2))
sigma_median = np.zeros([*np.shape(trx_phase), len(sigma_values_mask)+1])
m=0
tx_sens=np.transpose(mag,(1,0,2))
for mn in mask_loop:
    print(mn)
    trx_phase=  np.transpose(ph_7/np.max(ph_7)*2*np.pi,(1,0,2))*np.rot90(locals()[mn],0)
    trx_phase=np.nan_to_num(trx_phase,0)
    trx_phase= gaussian_filter(trx_phase, sigma=0.1)
    trx_phase[trx_phase==0]=np.nan

#%% Setup additional data variables

    os.chdir(ad)
    nn = tx_sens.shape
    size = nn[::-1]

    title = "EPTlib tutorial"
    description = "Let's learn how to use EPTlib through Python"

#%% Write input h5 file
    addr=os.path.join(ad,"Settings",file_name+"-input.h5")
    ept.h5_input(file_name,addr,tx_sens,trx_phase)

    # create all the configurations with H-EPT
    if method=="Helmholtz EPT":
        title = "Helmholtz EPT tutorial"
        configurations = {}
        key = file_name+"_hept_"
        description = "Phase-based approximation -Helmholtz"
        addr_tx_sens = ""
        addr_trx_phase = addr+":/"+file_name+"-input/"+"trx-phase"
        ofname=os.path.join(ad,key+"result"+".h5")
        addr_sigma = ofname+":/"+key.replace("-","/")+"/"+"sigma"
        addr_epsr = ""
        ept.config_hept_settings(configurations, title,description,size,step,freq,method,sg_size,sg_shape,sg_degree,addr_tx_sens,addr_trx_phase,addr_sigma,addr_epsr)
        print("Configuration for Helmholtz EPT...Done")
    elif method=="Convection-reaction EPT":
        # create all the configurations with CR-EPT
        title = "Convection-reaction EPT tutorial"
        configurations = {}
        key = file_name+"_crept_"
        description = "Phase-based approximation -CR"
        addr_tx_sens = ""
        addr_trx_phase = addr+":/"+file_name+"-input/"+"trx-phase"
        ofname=os.path.join(ad,key+"result"+".h5")
        addr_sigma = ofname+":/"+key.replace("-","/")+"/"+"sigma"
        addr_epsr = ""
        ept.config_crept_settings(configurations, title,description,size,step,freq,method,sg_size,sg_shape,sg_degree,bc_sigma,bc_epsr,volume_tomography,ad_coeff,addr_tx_sens,addr_trx_phase,addr_sigma,addr_epsr)
        print("Configuration for C-R EPT...Done")
    
    # create the configuration files    
    fname = os.path.join("Settings",key+".toml")
    with open(fname,"w") as ofile:
        toml.dump(configurations,ofile)

#%% RUN EPTLIB
    
    # Construct paths
    t = os.path.join(ad, "EPTlib_0.3.3/bin")
    
    executable = os.path.join(t, "EPTlib_app.exe") 
    config = os.path.join(ad, "settings", key + ".toml")
    log = os.path.join(ad+"/logs", key + ".log")
    
    command=''.join([executable, ' ',config])
    
    os.system(command)
    
    #%%

    with h5py.File(os.path.join(ad,ofname),'r') as ept_res:
        sigma_data = ept_res[key+'/sigma'][:]
    
    sigma_values_mask[m,0]=np.nanmedian(sigma_data)
    sigma_values_mask[m,1]=np.nanpercentile(sigma_data.flatten(), 75)-np.nanpercentile(sigma_data.flatten(), 25)

    
    sigma_median[:,:,:,m]=locals()[mn]*round(np.nanmedian(sigma_data),2)
    
    
    m=m+1
sigma_median=np.nan_to_num(sigma_median,0)    
sigma_median[:,:,:,m]=np.sum(sigma_median, axis=3)
sigma_median[sigma_median==0]=np.nan

with h5py.File(os.path.join(ad, ofname), 'a') as h5file:
    h5file.create_dataset('sigma_val', data=sigma_data)
    h5file.create_dataset('values_per_tissue', data=sigma_values_mask)

#%% scroll through images
fig, ax = plt.subplots(1, 1)
tracker = ept.IndexTracker(ax, sigma_data[:,:,:] ,'jet',0.3,1)
fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
plt.title('Conductivity is :'+ str(round(np.nanmedian(sigma_data),2))+' S/m')
plt.show()

    #%% Creating Report
    # os.chdir(ad)
    # from reportlab.pdfgen import canvas
    # from io import BytesIO
    # from reportlab.lib.pagesizes import letter
    # import matplotlib.pyplot as plt
    # import numpy as np
    
    
    # # Create a PDF file
    # buffer = BytesIO()
    # c = canvas.Canvas(buffer, pagesize=letter)
    # c.drawString(100, 750, "Report {} Data ({})".format(method, file_name))
    
    # fig, axs = plt.subplots(1, 2, figsize=(10, 4))  # Adjust figure size and subplot layout
    # plt.sca(axs[0])
    # plt.imshow(phase_0[:,:,sl], vmin=-np.pi, vmax=np.pi, cmap="turbo")
    # plt.title("Phase")
    # plt.sca(axs[1])
    # plt.imshow(trx_phase[:,:,sl], vmin=-np.pi, vmax=np.pi, cmap="turbo")
    # plt.title("Denoised Phase")
    # plt.savefig("plot0.png")
    # c.drawImage("plot0.png", 100, 500, width=400, height=200)
    
    
    # # Plot Median Values per Slice
    # plt.figure(figsize=(5, 2))
    # plt.plot(np.nanmedian(sigma_data, axis=(0, 1)))
    # plt.ylim(0, 1)
    # plt.xlabel('Slice number')
    # plt.ylabel('Conductivity [S/m]')
    # plt.title('Median Values per Slice')
    # plt.savefig("plot1.png")
    # c.drawImage("plot1.png", 100, 270,width=400, height=200)
    
    # # Add Median Value at Center Slice
    # median_center_slice = np.nanmedian(sigma_data[:,:,sl])
    
    # # Display Conductivity Map for the Center Slice
    # plt.figure(figsize=(5, 2))
    # plt.imshow(sigma_data[:,:,sl], vmin=0, vmax=1, cmap="jet")
    # plt.title("Median Conductivity value center slice: {:.2f} S/m".format(median_center_slice)+ "S/m" )
    # plt.savefig("plot2.png")
    # c.drawImage("plot2.png", 100, 70,width=400, height=200)
    
    # plt.close(plt.gcf())
    # c.save()
    
    # with open(key+'_'+tissue+'.pdf', 'wb') as f:
    #     buffer.seek(0)
    #     f.write(buffer.read())
    
    # # Delete the temporary PNG files
    # os.remove("plot0.png")
    # os.remove("plot1.png")
    # os.remove("plot2.png")
    # plt.close('all')
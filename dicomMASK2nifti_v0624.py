#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:05:22 2024

@author: Jessica A. Martinez 
"""
#dicom to numpy

import numpy as np
import pydicom
import os
import matplotlib.pyplot as plt
import nibabel as nib
import cv2

from brainextractor import BrainExtractor

os.chdir('/Users/jam48/Documents/Projects/MR_Safety/Subject_Specific_SAR_maps_Project\Code')

import Mods as mods

Dicom_path=mods.select_foldr()
# items = os.listdir(Dicom_path)
folders = [item for item in os.listdir(Dicom_path) if os.path.isdir(os.path.join(Dicom_path, item))]
name=os.path.split(str(os.path.split(Dicom_path)[1]))[1]
for fol in folders:
    Dicom_spath=os.path.join(Dicom_path,fol)
    
    files=[item for item in os.listdir(Dicom_spath) if item.endswith('.dcm') and os.path.isfile(os.path.join(Dicom_spath, item))]
    files.sort()
    # files=files[:-1] 

    for dic in range (0, len(files)):
        dcm_name=Dicom_spath+"/"+files[dic]
        ds = pydicom.dcmread(dcm_name)
        if dic==0:
            globals()[fol[:-5]] = np.zeros((ds.Columns, ds.Rows, len(files)), dtype=float) 
            rescale=np.zeros((len(files),2))
        globals()[fol[:-5]][:, :, dic] = ds.pixel_array-1
    globals()[fol[:-5]][globals()[fol[:-5]]==0]=np.nan
        #rescale[dic,0]=ds.RescaleSlope
        #rescale[dic,1]=ds.RescaleIntercept
    


#%%
WM_Mask=WM_Mask*WB_Mask
CSF_Mask=CSF_Mask*WB_Mask
GM_Mask=WB_Mask-(np.nan_to_num(WM_Mask)+np.nan_to_num(CSF_Mask))
GM_Mask[GM_Mask==0]=np.nan
Voxels=np.nan_to_num(WM_Mask)*2+np.nan_to_num(CSF_Mask)*3+np.nan_to_num(GM_Mask)
Voxels[Voxels==0]=np.nan
#%%
def get_path_up_to(directory, full_path):
    # Split the full path into its components
    path_parts = full_path.split(os.sep)
    
    # Find the index of the specified directory
    if directory in path_parts:
        index = path_parts.index(directory)
        # Join the path components up to and including the specified directory
        return os.sep.join(path_parts[:index + 1])
    
NIFTIPaths = get_path_up_to('NIFTIs', Dicom_path)


nib.save(nib.Nifti1Image(WB_Mask, np.eye(4)), os.path.join(NIFTIPaths,name+'_brain.nii.gz'))
nib.save(nib.Nifti1Image(GM_Mask, np.eye(4)), os.path.join(NIFTIPaths,name+'_brain_pve_0.nii.gz'))
nib.save(nib.Nifti1Image(WM_Mask, np.eye(4)), os.path.join(NIFTIPaths,name+'_brain_pve_1.nii.gz'))
nib.save(nib.Nifti1Image(CSF_Mask, np.eye(4)), os.path.join(NIFTIPaths,name+'_brain_pve_2.nii.gz'))
nib.save(nib.Nifti1Image(Voxels, np.eye(4)), os.path.join(NIFTIPaths,name+'_brain_pve_3.nii.gz'))
#%%
output_path=os.path.join(NIFTIPaths,name+'_brain_pve_0.nii.gz')
mask = nib.load(output_path).get_fdata()

# mask[mask>0]=1
# mask[mask==0]=np.nan
#%%  
    
fig, ax = plt.subplots(1, 1)
tracker = mods.IndexTracker(ax, 
    mask, 'jet',0,3)
fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
plt.show()
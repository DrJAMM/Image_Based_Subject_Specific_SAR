#ImportDicom.py

import numpy as np
import pydicom
import os
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('Agg')
import nibabel as nib
import cv2
import h5py
import toml

#from brainextractor import BrainExtractor
os.chdir('Z:\Documents\Projects\MR_Safety\Subject_Specific_SAR_maps_Project\Code')

import Mods as mods


# Define the directory path
mpl_config_dir = "./Process0/temp/matplotlib"

# Create the directory if it doesn't exist
os.makedirs(mpl_config_dir, exist_ok=True)

# Set the environment variable
os.environ['MPLCONFIGDIR'] = mpl_config_dir

####3

Dicom_path=mods.select_foldr()

path_root=os.path.dirname(os.path.dirname(Dicom_path))

files=os.listdir(Dicom_path)
files.sort()
files=files[:-1] 
#
for dic in range (0, len(files)):
    dcm_name=Dicom_path+"/"+files[dic]
    ds = pydicom.dcmread(dcm_name)
    if dic==0:
        data = np.zeros((ds.Columns, ds.Rows, len(files)), dtype=float) 
        rescale=np.zeros((len(files),2))
        s_type=[]
    data[:, :, dic] = ds.pixel_array
    
    rescale[dic,0]=ds.RescaleSlope
    rescale[dic,1]=ds.RescaleIntercept
    s_type.append(list(ds.ImageType)[3])
    step=[float(ds.PixelSpacing[0])*1e-3,float(ds.PixelSpacing[1])*1e-3,float(ds.SpacingBetweenSlices)*1e-3]
    freq=float(ds.ImagingFrequency)*1e6
    
#%%    
# Organice Data
mag=np.zeros((data.shape[0],data.shape[1],round(data.shape[2]/4)))
real=np.zeros((data.shape[0],data.shape[1],round(data.shape[2]/4)))
imag=np.zeros((data.shape[0],data.shape[1],round(data.shape[2]/4)))
phase=np.zeros((data.shape[0],data.shape[1],round(data.shape[2]/4)))

v0 = 0
v = 0
thresholds = [round(data.shape[2] / 4) * i for i in range(1, 4)]

for st in s_type:
    if st == 'M':
        mag[:, :, v] = data[:, :, v0]
        v += 1
    elif st == 'R':
        real[:, :, v] = data[:, :, v0]
        v += 1
    elif st == 'I':
        imag[:, :, v] = data[:, :, v0]
        v += 1
    elif st == 'P':
        phase[:, :, v] = data[:, :, v0]
        v += 1
    
    v0 += 1
    
    if v0 in thresholds:
        v = 0

#%%
from skimage.restoration import unwrap_phase


phase_mf_2=np.zeros(np.shape(phase))
phase_mf_3=np.zeros(np.shape(phase))
phase_mf_4=np.zeros(np.shape(phase))
phase_mf_5=np.zeros(np.shape(phase))
phase_mf_7=np.zeros(np.shape(phase))
phase_mf_9=np.zeros(np.shape(phase))
phase_mf=np.zeros(np.shape(phase))
for sl in range (1,np.size(phase,2)):
    print(f"\rMedian Filter, slice: {sl}... Percentage Progress: {round((sl / np.size(phase,2)) * 100, 1)}%", end='', flush=True)
    phase_mf_2[:,:,sl]=mods.median_nonan(phase[:,:,sl],sl, 2)
    phase_mf_3[:,:,sl]=mods.median_nonan(phase[:,:,sl],sl, 3)
    phase_mf_4[:,:,sl]=mods.median_nonan(phase[:,:,sl],sl, 4)
    phase_mf_5[:,:,sl]=mods.median_nonan(phase[:,:,sl],sl, 5)
    phase_mf_7[:,:,sl]=mods.median_nonan(phase[:,:,sl],sl, 7)
    phase_mf_9[:,:,sl]=mods.median_nonan(phase[:,:,sl],sl, 9)
print(f"\rMedian Filter, slice: {sl}... Percentage Progress: {round( 100, 1)}%")

#%%
tail=os.path.basename(Dicom_path)
parts = tail.split('_')
new_tail_parts = [parts[-1]] + parts[:-1]
new_tail = '_'.join(new_tail_parts)
file_name = new_tail.lower() 
output_path=os.path.join(path_root,"NIFTIs")

# output_path = "//Mac/Home/Documents/Projects/MR_Safety/Subject_Specific_SAR_maps_Project/Test_Brain_Ept_Jamm/Unnamed - 758915522/T1_TSE_COR_MASKS"
#%%
if not os.path.isfile(output_path+'/'+file_name+'_brain.nii.gz'):
    print('The Data has no Masks... ',end='')
    mask=[]
    mask_0=[]
    mask_1=[]
    mask_2=[]
else:

    mask=nib.load(output_path+'/'+file_name+'_brain.nii.gz').get_fdata()
    mask_0 = nib.load(output_path+'/'+file_name+'_brain_pve_0.nii.gz').get_fdata()
    mask_1 = nib.load(output_path+'/'+file_name+'_brain_pve_1.nii.gz').get_fdata()
    mask_2 = nib.load(output_path+'/'+file_name+'_brain_pve_2.nii.gz').get_fdata()

    mask[mask>0]=1
    mask[mask==0]=np.nan
    
    mask_0[mask_0>0]=1
    mask_0[mask_0==0]=np.nan
    
    mask_1[mask_1>0]=1
    mask_1[mask_1==0]=np.nan
    
    mask_2[mask_2>0]=1
    mask_2[mask_2==0]=np.nan

#%% Denoising with the ept method
# import EPTMods2 as ept
# phase_dn = ept.denoise_mri_volume(phase)

#%%
# Load pre-made Masks using FSL
tail=os.path.basename(Dicom_path)
parts = tail.split('_')
new_tail_parts = [parts[-1]] + parts[:-1]
new_tail = '_'.join(new_tail_parts)
file_name = new_tail.lower() 
#%%
# Create HDF5 File
ept_path=path_root+'/forEPT'
if not os.path.exists(ept_path):
    os.makedirs(ept_path)
hf = h5py.File(ept_path+'/'+os.path.basename(Dicom_path)+'_phase_br.h5', 'w')
hf.create_dataset('mag', data=mag)
hf.create_dataset('ph', data=phase)
hf.create_dataset('ph_2', data=phase_mf_2)
hf.create_dataset('ph_3', data=phase_mf_3)
hf.create_dataset('ph_4', data=phase_mf_4)
hf.create_dataset('ph_5', data=phase_mf_5)
hf.create_dataset('ph_7', data=phase_mf_7)
hf.create_dataset('ph_9', data=phase_mf_9)
hf.create_dataset('mask', data=mask)
hf.create_dataset('mask_0', data=mask_0)
hf.create_dataset('mask_1', data=mask_1)
hf.create_dataset('mask_2', data=mask_2)
hf.create_dataset('step',data=step)
hf.create_dataset('freq',data=freq)

hf.close()

print('HDF5 File Created!')




# %%


# Example usage
# import matplotlib.pyplot as plt

# input_img = mag  # Example 3D image
# threshold = 150  # Example threshold

# brain_mask = mods.tissue_extractor(input_img, threshold, solidity_thresh=0.2)

# #%%
# fig, ax = plt.subplots(1, 1)
# tracker = mods.IndexTracker(ax, 
#     brain_mask, 'jet',0,1)
# fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
# plt.show()
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
from skimage.restoration import unwrap_phase

#from brainextractor import BrainExtractor
os.chdir('Z:\Documents\Projects\MR_Safety\Subject_Specific_SAR_maps_Project\Code')

import Mods as mods

#%% Load Data
# Define the directory path
mpl_config_dir = "./Process0/temp/matplotlib"

# Create the directory if it doesn't exist
os.makedirs(mpl_config_dir, exist_ok=True)

# Set the environment variable
os.environ['MPLCONFIGDIR'] = mpl_config_dir


Dicom_path_m=mods.select_foldr()

name_mag=os.path.basename(Dicom_path_m)

name_end=name_mag.split('_')[-1]
name_end=name_end.replace(name_end[-1],str(int(name_end[-1])+1))
name_ph=name_mag[:-len(name_end)]+name_end



path_root=os.path.dirname(os.path.dirname(Dicom_path_m))

Dicom_path_ph=Dicom_path_m.replace(name_mag,name_ph)


files_m=os.listdir(Dicom_path_m)
files_m.sort()
files_m=files_m[:-1] 

files_ph=os.listdir(Dicom_path_ph)
files_ph.sort()
files_ph=files_ph[:-1] 

#%%
# Getting Magnitude
for dic in range (0, len(files_m)):
    dcm_name=Dicom_path_m+"/"+files_m[dic]
    ds = pydicom.dcmread(dcm_name)
    if dic==0:
        data = np.zeros((ds.Columns, ds.Rows, len(files_m)), dtype=float) 
    data[:, :, dic] = ds.pixel_array
    step=[float(ds.PixelSpacing[0])*1e-3,float(ds.PixelSpacing[1])*1e-3,float(ds.SliceThickness)*1e-3]
    freq=float(ds.ImagingFrequency)*1e6
mag = data

# Getting Phase
for dic in range (0, len(files_ph)):
    dcm_name=Dicom_path_ph+"/"+files_ph[dic]
    ds = pydicom.dcmread(dcm_name)
    if dic==0:
        data = np.zeros((ds.Columns, ds.Rows, len(files_ph)), dtype=float) 
    data[:, :, dic] = ds.pixel_array
phase_r = data

#%% Unwrap Phase
phase_r_mf=np.zeros(np.shape(phase_r))
for sl in range (1,np.size(phase_r,2)):
    print(f"\rMedian Filter, slice: {sl}... Percentage Progress: {round((sl / np.size(phase_r,2)) * 100, 1)}%", end='', flush=True)
    phase_r_mf[:,:,sl]=mods.median_nonan(phase_r[:,:,sl],sl, 3)
print(f"\rMedian Filter, slice: {sl}... Percentage Progress: {round( 100, 1)}%", end='')
uphase=unwrap_phase(phase_r_mf,rng=1000)
uphase2=unwrap_phase(phase_r,rng=1000)

#%% Create HDF5 File
ept_path=path_root+'/forEPT'
if not os.path.exists(ept_path):
    os.makedirs(ept_path)
hf = h5py.File(ept_path+'/'+name_mag[:-len(name_end)-1]+'_'+os.path.basename(path_root)+'.h5', 'w')
hf.create_dataset('tx_sens', data=mag)
hf.create_dataset('phase', data=phase_r)

hf.create_dataset('phase_mf', data=phase_r_mf)
hf.create_dataset('phase_uw', data=uphase2)
hf.create_dataset('phase_mf_uw', data=uphase)


hf.create_dataset('step',data=step)
hf.create_dataset('freq',data=freq)

hf.close()

print('HDF5 File Created!')

#%%
# fig, ax = plt.subplots(1, 1)
# tracker = mods.IndexTracker(ax, 
#     phase_r, 'jet',0,2600)
# fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
# plt.show()
# %%

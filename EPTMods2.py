# -*- coding: utf-8 -*-
"""
Created on Fri May 24 17:28:48 2024

@author: Jessica A. Martinez
"""
import tkinter as tk
from tkinter import filedialog
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import label, find_objects
from skimage.measure import regionprops
import os
import skimage.measure  
import matplotlib.pyplot as plt  
import nibabel as nib
import subprocess
from PyQt5.QtWidgets import QApplication, QFileDialog
import sys
import cv2 as cv 

def select_HDF5file():
    # root = tk.Tk()
    # root.withdraw()  # Hide the root window
    # address= filedialog.askopenfilename  (title="Select HDF5 File")
    # # root.destroy()  # Destroy the root window
    import sys

    try:
        app = QApplication(sys.argv)  # Initialize the application with sys.argv
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog  # Uncomment this line if you face issues with native dialogs
        file_name, _ = QFileDialog.getOpenFileName(None, "Select HDF5 File", "", "HDF5 Files (*.h5 *.hdf5);;All Files (*)", options=options)
        app.quit()  # Quit the application after file selection
        return file_name
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def read_HDF5_data(file_path):
    try:
        with h5py.File(file_path, "r") as ifile:
            tx_sens = ifile["tx_sens"][()]
            trx_phase = ifile["trx_phase"][()]
            mask = ifile["mask"][()]
            mask_wm = ifile["mask_wm"][()]
            mask_gm = ifile["mask_gm"][()]
            mask_csf = ifile["mask_csf"][()]
            sigma = ifile["sigma"][()]
            step= ifile["step"][()]
            freq= ifile["freq"][()]
        return tx_sens, trx_phase, mask,mask_wm,mask_gm,mask_csf, sigma,step,freq
    except KeyError as e:
        print(f"Dataset not found: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)

def read_HDF5_ph_data(file_path):
    try:
        with h5py.File(file_path, "r") as ifile:
            
            if "tx-sens" in ifile:   
                tx_sens = ifile["tx-sens"][()]
            elif "tx_sens" in ifile:
                tx_sens = ifile["tx_sens"][()]
            else:
                tx_sens=0
                
            if "trx-phase" in ifile:
                trx_phase = ifile["trx-phase"][()]
            elif "trx_phase" in ifile:
                trx_phase = ifile["trx_phase"][()]
            else:
                trx_phase=0;   
                 
            if "mask-wm" in ifile:
                mask = ifile["mask-wm"][()]
            elif "mask" in ifile:
                mask = ifile["mask"][()]
            else:
                mask=0;    
                
            sigma = ifile["sigma"][()]
            if "step" in ifile:
                step= ifile["step"][()]
            else:
                step=0
                
            if "freq" in ifile:
                freq= ifile["freq"][()]
            else:
                freq=0
            
        return tx_sens, trx_phase, mask, sigma,step,freq
    except KeyError as e:
        print(f"Dataset not found: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
####
def create_table_conditionally(config, label):
    if (label not in config.keys() or not isinstance(config[label],dict)):
        config[label] = {}
        
def config_common_settings(config, title,description,size,step,freq):
    config["title"] = title
    config["description"] = description
    create_table_conditionally(config,"mesh")
    config["mesh"]["size"] = size
    config["mesh"]["step"] = step
    create_table_conditionally(config,"input")
    config["input"]["frequency"] = freq
    config["input"]["tx-channels"] = 1
    config["input"]["rx-channels"] = 1
    config["input"]["wrapped-phase"] = False
    
def config_io_settings(config, in_tx_sens,in_trx_phase,out_sigma,out_epsr):
    create_table_conditionally(config,"input")
    config["input"]["tx-sensitivity"] = in_tx_sens
    config["input"]["reference-image"] = in_tx_sens
    config["input"]["trx-phase"] = in_trx_phase
    create_table_conditionally(config,"output")
    config["output"]["electric-conductivity"] = out_sigma
    config["output"]["relative-permittivity"] = out_epsr


#######
    
def h5_input(file_name,addr,tx_sens,trx_phase):
    with h5py.File(addr,"w") as ofile:
        group_noiseless = ofile.create_group(file_name+"-input")
        group_noiseless.create_dataset("tx-sens",data=tx_sens)
        group_noiseless.create_dataset("trx-phase",data=trx_phase)
#%%

methods = {
    "Helmholtz EPT": 0,
    "Convection-reaction EPT": 1,
    "Gradient EPT": 2
}
# define a function configuring the method setting
def set_option_from_dictionary(config, label,value,dictionary):
    if isinstance(value,str):
        config[label] = dictionary[value]
    else:
        config[label] = value
def config_method_settings(config, method):
    set_option_from_dictionary(config, "method",method,methods)
    
shapes = {
    "Cross": 0,
    "Sphere": 1,
    "Cube": 2
}

def config_sg_settings(config, size,shape,degree):
    create_table_conditionally(config,"parameter")
    create_table_conditionally(config["parameter"],"savitzky-golay")
    config["parameter"]["savitzky-golay"]["size"] = size
    set_option_from_dictionary(config["parameter"]["savitzky-golay"], "shape",shape,shapes)
    config["parameter"]["savitzky-golay"]["degree"] = degree
def config_bc_settings(config, sigma,epsr):
    create_table_conditionally(config,"parameter")
    create_table_conditionally(config["parameter"],"dirichlet")
    config["parameter"]["dirichlet"]["electric-conductivity"] = sigma
    config["parameter"]["dirichlet"]["relative-permittivity"] = epsr
# define a function configuring the problem dimensions
def config_volume_settings(config, volume_tomography,slice=None):
    create_table_conditionally(config,"parameter")
    config["parameter"]["volume-tomography"] = volume_tomography
    if (not volume_tomography and slice):
        config["parameter"]["imaging-slice"] = slice
# define a function configuring the artificial diffusion
def config_artificial_diffusion_settings(config, coeff):
    create_table_conditionally(config,"parameter")
    if (coeff>0):
        config["parameter"]["artificial-diffusion"] = True
        config["parameter"]["artificial-diffusion-coefficient"] = coeff
def config_hept_settings(config, title,description,size,step,freq, method, sg_size,sg_shape,sg_degree, in_tx_sens,in_trx_phase,out_sigma,out_epsr):
    config_common_settings(config, title,description,size,step,freq)
    config_method_settings(config, method)
    config_sg_settings(config, sg_size,sg_shape,sg_degree)
    config_io_settings(config, in_tx_sens,in_trx_phase,out_sigma,out_epsr)   

def config_crept_settings(config, title,description,size,step,freq, method, sg_size,sg_shape,sg_degree, bc_sigma,bc_epsr, volume_tomography, ad_coeff, in_tx_sens,in_trx_phase,out_sigma,out_epsr):
    config_common_settings(config, title,description,size,step,freq)
    config_method_settings(config, method)
    config_sg_settings(config, sg_size,sg_shape,sg_degree)
    config_bc_settings(config, bc_sigma,bc_epsr)
    config_volume_settings(config, volume_tomography)
    config_artificial_diffusion_settings(config, ad_coeff)
    config_io_settings(config, in_tx_sens,in_trx_phase,out_sigma,out_epsr)
#####
# Denoise

#1. map values to match the denoise_mri_slices function
def map_phase_to_uint8(phase):
    """
    Map phase values from the range [-pi, pi] to [0, 255] for denoising.
    
    :param phase: Phase image represented as a NumPy array.
    :return: Phase image mapped to uint8 in the range [0, 255].
    """
    # Map phase values to the range [0, 255]
    phase_mapped = ((phase + np.pi) * (255 / (2 * np.pi))).astype(np.uint8)
    return phase_mapped
#2. Denoise :)
def denoise_mri_slices(volume, target_slice_index, temporal_window_size=5, h=3):
    """
    Apply Non-Local Means denoising on a sequence of MRI slices.
    
    :param volume: 3D NumPy array representing the MRI volume (dimensions: height x width x number_of_slices).
    :param target_slice_index: Index of the target slice to be denoised.
    :param temporal_window_size: Temporal window size for denoising.
    :param h: Parameter regulating filter strength.
    :return: Denoised slice as a 2D NumPy array.
    """
    # Convert 3D volume to a list of 2D slices
    slices = [volume[:, :, i] for i in range(volume.shape[2])]
    
    # Denoise the target slice
    denoised_slice = cv.fastNlMeansDenoisingMulti(
        srcImgs=slices,
        imgToDenoiseIndex=target_slice_index,
        temporalWindowSize=temporal_window_size,
        h=h
    )
    
    return denoised_slice

def denoise_mri_volume(volume, temporal_window_size=5, h=3):
    """
    Apply Non-Local Means denoising to an MRI volume.
    
    :param volume: 3D NumPy array representing the MRI volume (dimensions: height x width x number_of_slices).
    :param temporal_window_size: Temporal window size for denoising.
    :param h: Parameter regulating filter strength.
    :return: Denoised volume as a 3D NumPy array.
    """
    # Convert volume to uint8
    volume_uint8 = map_phase_to_uint8(volume)
    
    # Ensure that temporal window size is within bounds
    temporal_window_size = min(temporal_window_size, volume.shape[2])
    
    denoised_volume = np.zeros_like(volume)  # Initialize denoised volume

    # Denoise each slice separately
    for i in range(volume.shape[2]):
        # Skip denoising if the slice is all zeros
        if np.all(volume[:, :, i] == 0):
            denoised_volume[:, :, i] = volume[:, :, i]  # Keep the original slice
        else:
            # Denoise the slice
            denoised_slice_uint8 = cv.fastNlMeansDenoising(
                volume_uint8[:, :, i],
                None,  # Noisy image (set to None as the source image is already provided)
                h=h,
                templateWindowSize=7,  # Size of the window used to compute the weighted average for a pixel
                searchWindowSize=21     # Size of the window used to search for a pixel with the most similar neighborhood
            )
            # Map denoised values back to phase range
            denoised_volume[:, :, i] = map_uint8_to_phase(denoised_slice_uint8)

    return denoised_volume


#3. Back to phase
def map_uint8_to_phase(phase_uint8):
    """
    Map phase values from uint8 range [0, 255] back to the range [-pi, pi].
    
    :param phase_uint8: Denoised phase image mapped to uint8.
    :return: Phase image mapped back to the range [-pi, pi].
    """
    # Map uint8 values back to the range [-pi, pi]
    phase = (phase_uint8 * (2 * np.pi / 255)) - np.pi
    return phase
###

def m(tissue,mask_b, mask_wm, mask_gm, mask_csf):
    if tissue=='CSF':
        mask=mask_csf
    elif tissue=='WM':
        mask=mask_wm
    elif tissue=='GM':
        mask=mask_gm
    elif tissue=='Whole':
        mask=mask_b
    return mask
    




###Image Viewer


class IndexTracker:
    def __init__(self, ax, X,comap,minv, maxv):
        self.ax = ax
        # ax.set_title('slice %s' % self.ind)

        self.X = X
        rows, cols, self.slices = X.shape
        self.ind = self.slices//2

        self.im = ax.imshow(self.X[:, :, self.ind],
                            vmin=minv, vmax=maxv,cmap=comap, aspect='auto')
        plt.colorbar(self.im)
        plt.axis('off')
        self.update()

    def on_scroll(self, event):
        # print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.ax.set_title('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()
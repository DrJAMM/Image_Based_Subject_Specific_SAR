"""
Created on Fri Apr  5 15:05:22 2024

@author: Jessica A. Martinez 
"""

# Some easy functions to run the stuff smooth.

#1. Load Data path

import tkinter as tk
from tkinter import filedialog
from skimage.measure import shannon_entropy
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
from skimage.measure import label, regionprops
from skimage.segmentation import clear_border
from skimage.morphology import binary_dilation, disk
from skimage.segmentation import flood_fill
#import easygui

#%% Scroll through Images
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

#%%% Select Files
def select_foldr():
    from PyQt5.QtWidgets import QApplication, QFileDialog
    app = QApplication([])
    folder_path = QFileDialog.getExistingDirectory(None, "Select Folder")
    app.quit()
    return folder_path    

def select_folder():
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    
    folder_path = filedialog.askdirectory(title="Select Folder")
    
    root.destroy()  # Destroy the root window
    
    return folder_path

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

#%% Filters
# def median_nonan(data, filter_size):
    
#     indexer = filter_size // 2
    
#     new_image = data.copy()
    
#     nrow, ncol = data.shape
        
#     E1 = skimage.measure.shannon_entropy(data[:,:])
#     E2=0;
#     diff_entropy=0;#abs(E2-E1)
#     itera=0
    
#     while diff_entropy<=0.0001 and itera<100:
#         tmp_image = new_image.copy()
#         mask=np.ones((nrow,ncol))
#         for i in range(nrow):
#             for j in range(ncol):
#                 if mask[i,j]==1:
                
#                         itmp=0;
#                         tmp=[];
#                         # create the window
#                         for k in range(i-indexer, i+indexer+1):           
#                             for m in range(j-indexer, j+indexer+1):
                           
#                                 if (k > -1) and (k < nrow):
#                                   if (m > -1) and (m < ncol) :                      
#                                       tmp.append(tmp_image[k,m])
#                                       itmp=itmp+1
              
#                         # window created
#                         if (itmp>0):
                            
#                             contrast=np.nanmax(tmp)-np.nanmin(tmp)
#                             b1=contrast/2
#                             b2=np.nanmedian(tmp)+(contrast/2)
#                             if data[i,j]>b1 and  data[i,j]<b2:
#                                 mask[i,j]=0;
#                             else:    
#                                 new_image[i,j] = np.nanmedian(tmp)                                    
#                         else:
#                             new_image[i,j] =np.nan
                            
#             E2 = skimage.measure.shannon_entropy(new_image[:,:])
#             diff_entropy=abs(E2-E1)
#             itera=itera+1
#             # print(diff_entropy)
#     return new_image


def median_nonan(data, sl,filter_size):
    indexer = filter_size // 2
    new_image = data.copy()
    nrow, ncol = data.shape
    
    E1 = shannon_entropy(data)
    E2 = 0
    diff_entropy = 0
    itera = 0

    while diff_entropy <= 0.0001 and itera < 100:
        tmp_image = new_image.copy()
        mask = np.ones((nrow, ncol), dtype=bool)
        
        padded_image = np.pad(tmp_image, indexer, mode='edge')
        # print(f"\rMedian Filter, slice: {sl}...)%", end='')
        
        for i in range(nrow):
            # print(f"\rMedian Filter, slice: {sl}... Percentage Progress: {round((i / nrow) * 100, 1)}%", end='', flush=True)
            for j in range(ncol):
                if mask[i, j]:
                    window = padded_image[i:i + filter_size, j:j + filter_size]
                    contrast = np.nanmax(window) - np.nanmin(window)
                    b1 = contrast / 2
                    b2 = np.nanmedian(window) + (contrast / 2)
                    
                    if b1 < data[i, j] < b2:
                        mask[i, j] = False
                    else:
                        new_image[i, j] = np.nanmedian(window)
        # print(f"\rMedian Filter, slice: {sl}... Percentage Progress: {round(p/100, 1)}%", end='', flush=True)
        E2 = shannon_entropy(new_image)
        diff_entropy = abs(E2 - E1)
        itera += 1
    
    return new_image

#%% Make mask
def fill_holes(bwImg):
    """
    Fills holes in binary images using morphological operations.

    Parameters:
    - bwImg: array_like, Input binary image.

    Returns:
    - filledImg: array_like, Binary image with holes filled.
    """
    from scipy.ndimage import binary_fill_holes
    return binary_fill_holes(bwImg)

def brainExtractor(inputImg, threshold, **kwargs):
    """
    Description: Simple automated brain extraction.

    Inputs:
    - inputImg: matrix, Input image to segment (e.g. PD map)
    - threshold: int, Threshold value to use

    Outputs:
    - brainMask: matrix, Mask corresponding to the brain tissue

    Based on code from https://www.mathworks.com/matlabcentral/answers/42820-skull-striping-without-affecting-tumor-region
    Change Log:
    """
    # Default values
    solidityThresh = kwargs.get('solidityThresh', 0.04)

    # Threshold the image
    bwImg = inputImg > threshold
    bwImg = fill_holes(bwImg)

    # Pre-allocate
    brainMask = np.zeros_like(inputImg)

    # Select the brain region for each slice
    for ii in range(inputImg.shape[2]):
        # Apply labels based on thresholded image
        currLabel, num_labels = label(bwImg[:, :, ii])

        props = regionprops(currLabel, intensity_image=inputImg[:, :, ii])
        if props:
            # solidity is the percentage "filled" of an area. For the skull,the solidity will be really low.
            solidity = [prop.solidity for prop in props]
            area = [prop.area for prop in props]
            hiSolid = np.array(solidity) > solidityThresh  # get only high solidity objects
            maxArea = max(np.array(area)[hiSolid], default=0)
            brainLabel = np.argmax(np.array(area) == maxArea) + 1  # label of brain
            brain = currLabel == brainLabel  # b/w image of brain
            brain=fill_holes(brain)
            # Set the brain mask of the current slice to be all pixels corresponding to that label
            brainMask[:, :, ii] = brain.astype(np.uint8)
        # Set the brain mask of the current slice to be all pixels corresponding to that label
        

    return brainMask

def run_bet(input_image, output_image, options=''):
    """
    Run the FSL's 'bet' command from a NumPy array.

    Parameters:
        input_array (numpy.ndarray): Input 3D numpy array representing the brain image.
        output_path (str): Path to save the output brain image.
        options (str): Additional options for the bet command.

    Returns:
        None
    """
    # Convert NumPy array to NIfTI image
    # input_nifti = nib.Nifti1Image(input_array, np.eye(4))

    # Save the NIfTI image temporarily
    temp_input_path = '/Users/jam48/temp_input.nii.gz'
    # nib.save(input_nifti, temp_input_path)


    # Change current directory to the FSL bin directory

    os.chdir('/Users/jam48/fsl/bin')
    # Specify the relative path to the bet executable
    bet_path = './bet'

    # Run FSL's 'bet' command with the specified options
    command = [bet_path, input_image, output_image] + options.split()
    print(command)
    new_env = os.environ.copy()
    subprocess.run(command,env=new_env)
    
# tissue_extractor
def binary_fill_holes(binary_image):
    # Invert the binary image
    inverted_image = np.logical_not(binary_image)
    
    # Label the inverted image
    labeled_inverted_image, num_labels = label(inverted_image, return_num=True)
    
    # Select the label that corresponds to the background (assuming it is the label 0)
    filled_image = np.logical_not(labeled_inverted_image == 0)
    
    return filled_image

def tissue_extractor(input_img, threshold, **kwargs):
    """
    Simple automated brain extraction.

    Parameters:
        input_img (ndarray): Input image to segment (e.g. PD map).
        threshold (int): Threshold value to use.
        **kwargs: Optional parameters.
            solidity_thresh (float): Threshold for solidity. Default is 0.2.

    Returns:
        brain_mask (ndarray): Mask corresponding to the brain tissue.
    """

    # Default values
    solidity_thresh = kwargs.get('solidity_thresh', 0.2)

    # Threshold the image
    bw_img = input_img > threshold
    bw_img = binary_fill_holes(bw_img)

    # Pre-allocate
    brain_mask = np.zeros_like(input_img)

    # Select the brain region for each slice
    for ii in range(input_img.shape[2]):
        # Apply labels based on thresholded image
        curr_label = label(bw_img[:, :, ii])

        props = regionprops(curr_label)
        # Solidity is the percentage "filled" of an area. For the skull, the solidity will be really low.
        solidity = np.array([prop.solidity for prop in props])
        area = np.array([prop.area for prop in props])
        hi_solid = solidity > solidity_thresh  # Get only high solidity objects
        max_area = np.max(area[hi_solid])
        brain_label = [prop.label for prop in props if prop.area == max_area][0]  # Label of brain
        brain = curr_label == brain_label  # Binary image of brain

        # Set the brain mask of the current slice to be all pixels corresponding to that label
        brain_mask[:, :, ii] = brain

    return brain_mask
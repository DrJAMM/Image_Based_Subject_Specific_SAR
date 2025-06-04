## Introduction

Welcome to the Image_Based_Subject_Specific_SAR repository! Here you’ll find everything you need to recreate the, practical approach we developed in <a href="https://doi.org/10.1002/mrm.30547" target="_blank" rel="noopener noreferrer">https://doi.org/10.1002/mrm.30547</a> for generating subject-specific brain SAR maps from routine MRI data. Instead of relying on complex electromagnetic simulations, our method uses anatomical scans and B₁⁺ information to estimate localized SAR distributions. Feel free to explore the scripts and the step-by-step instructions to obtain conductivity and SAR maps. We hope this makes your SAR mapping work a little simpler and more enjoyable!

If you have any questions, please send me an email:  
[jessica.a.martinezm@gmail.com](mailto:jessica.a.martinezm@gmail.com)

---

## Citation

If you use this repository in your work, we kindly ask that you cite:

Martinez JA, Zanovello U, Arduino A, Hu HH, Moulin K, Ogier SE, Bottauscio O, Zilberti L, Keenan KE.  
**Feasibility study of subject‐specific, brain specific‐absorption‐rate maps retrieved from MRI data.**  
*Magnetic Resonance in Medicine*. 2025 May 27.  
DOI: <a href="https://doi.org/10.1002/mrm.30547" target="_blank" rel="noopener noreferrer">https://doi.org/10.1002/mrm.30547</a>


---

## Steps

### i. Making Masks
#### 1. Extract Masks using horos. 
	Keep 0 as 1 and 1 as 2 (that way they matrix size doesn’t change)
	Save the masks in a folder in the niftis file with the name of the sequence
#### 2. Run "dicomMASK2nifti.py" to save the masks in the NIFTI folder
### ii. Saving Data for EPTLib Format
#### 3. Run "ImportDicom_v0624.py" to create the H5 files for SAR EPT.
	 Here there are many median filtered phase data
### iii.  Run EPTLib
#### 4. Run "EPT_Analysis_inVivo_0626.py" for brain or  "EPT_Analysis_inVitro_0626.py" to get the conductivity values for the in vivo values, the median values are already assigned to the tissues.

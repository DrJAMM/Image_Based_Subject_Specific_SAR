# Making Masks
1. Extract Masks using horos. 
	Keep 0 as 1 and 1 as 2 (that way they matrix size doesn’t change)
	Save the masks in a folder in the niftis file with the name of the sequence
2. Run "dicomMASK2nifti.py" to save the masks in the NIFTI folder
# Saving Data for EPTLib Format
3. Run "ImportDicom_v0624.py" to create the H5 files for SAR EPT.
	 Here there are many median filtered phase data
# Run EPTLib
4. Run "EPT_Analysis_inVivo_0626.py" for brain or  "EPT_Analysis_inVitro_0626.py" to get the conductivity values 
for the in vivo values, the median values are already assigned to the tissues.

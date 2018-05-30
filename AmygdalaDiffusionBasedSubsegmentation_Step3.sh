# STEP 3 : Masks Divided and Binarized By KS Damme
# May 24 2017
# ZM Saygin 2011 Logic  - rewritten KSD
# In each subject, we calculated the connection probability (using FSL- FDT's probtrackX) from each amygdala voxel (seed) to all bilateral cortical and subcortical regions (targets), and normalized the distribution of probabilities for each seed voxel to [0,1] by dividing by the maximum probability. We then thresholded and binarized these results to exclude values below 0.1, such that every amygdaloid voxel contained a 0 or 1 for each target.
foreach subj (`cat $SUBJECTS_DIR/subjects.txt`)
foreach hemi (Left Right)
echo ${subj}
echo ${hemi}
foreach masks (seeds_to_${subj}_16_Brain-Stem.label.nii.gz seeds_to_${subj}_bilateral_Parahippocampus.nii.gz  seeds_to_${subj}_bilateral_Accumbens.nii.gz seeds_to_${subj}_bilateral_Pericalcarine.nii.gz  seeds_to_${subj}_bilateral_Caudal_Middle-Frontal.nii.gz seeds_to_${subj}_bilateral_Postcentral.nii.gz  seeds_to_${subj}_bilateral_Caudate.nii.gz seeds_to_${subj}_bilateral_Putamen.nii.gz  seeds_to_${subj}_bilateral_Cuneus.nii.gz seeds_to_${subj}_bilateral_Rostral_Anterior_Cingulate.nii.gz  seeds_to_${subj}_bilateral_Fusiform.nii.gz seeds_to_${subj}_bilateral_Striatum.nii.gz  seeds_to_${subj}_bilateral_Hippocampus.nii.gz seeds_to_${subj}_bilateral_Superior_Frontal.nii.gz  seeds_to_${subj}_bilateral_Inferior_Temporal.nii.gz seeds_to_${subj}_bilateral_﻿Superior_Parietal.nii.gz  seeds_to_${subj}_bilateral_Insula.nii.gz seeds_to_${subj}_bilateral_Superior_Temporal.nii.gz  seeds_to_${subj}_bilateral_Lateral_Occipital.nii.gz seeds_to_${subj}_bilateral_Temporal_Pole.nii.gz  seeds_to_${subj}_bilateral_Lateral_Orbitofrontal.nii.gz seeds_to_${subj}_bilateral_Thalamus_Proper.nii.gz  seeds_to_${subj}_bilateral_Lingual.nii.gz seeds_to_${subj}_bilateral_Ventral_Diencephalon.nii.gz  seeds_to_${subj}_bilateral_Medial_Orbitofrontal.nii.gz)
echo ${masks}
fslmaths ${subj}/dpath/Segmentation/${hemi}/${masks} -div `fslstats ${subj}/dpath/Segmentation/${hemi}/${masks} -R|cut -d\  -f2` -thr 0.1 -bin ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_${masks}
end
end
end

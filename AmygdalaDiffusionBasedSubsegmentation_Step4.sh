# Amygdala Diffustion Base By KS Damme
# May 24 2017
# ZM Saygin 2011 Logic  - rewritten KSD
# In each subject, we calculated the connection probability (using FSL- FDT's probtrackX) from each amygdala voxel (seed) to all bilateral cortical and subcortical regions (targets), and normalized the distribution of probabilities for each seed voxel to [0,1] by dividing by the maximum probability. We then thresholded and binarized these results to exclude values below 0.1, such that every amygdaloid voxel contained a 0 or 1 for each target.
# The Current script uses the boolean logic validated in Saygin et al., 2011 (https://www.ncbi.nlm.nih.gov/pubmed/21396459)
# Porduces a csv file with data in long form that includes subject name (fsid), hemisphere nuclei, and volume called subsegmentation_volume.csv
# 
## List of Relevant Cortical Targets FS labels ##
#---------------#
###  Parietal  ###
####  Superior Parietal
####  Postcentral
####Cuneus
###  Frontal  ###
####  Medial orbitofrontal
####  Superior frontal
####  caudal middle-frontal
####  Lateral orbitofrontal 
###  Temporal  ###
####  Superior Temporal
####  Inferior Temporal
####  Parahippocampus
####  Hippocampus
####  Temporal Pole
####  Fusiform
### Occipital ###
####  Lateral Occipital
####  Pericalcarine
####  Lingual
###  Other  ###
####  Rostral anterior cingulate
####  Insula
####  Accumbens
####  Brain Stem
####  Ventral Diencephalon
####  Thalamus Proper
####  Striatum

### Iterate through of the subjects making a bilateral mask for all of the lateral structures ALLTARGETS.txt
foreach line (`cat bilateral_labels.csv`)
set struct = `echo $line|cut -d\, -f1`
set lh = `echo $line|cut -d\, -f2`
set rh = `echo $line|cut -d\, -f3`
echo ${line}
echo ${struct}
echo ${lh}
echo ${rh}
fslmaths $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}${rh} -add $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}${lh} $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_${struct}.nii.gz
echo "$SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_${struct}.nii.gz" >>! $SUBJECTS_DIR/${subj}/ALLTARGETS.txt
end
foreach striatal_part (Caudate,11_Left-Caudate.label.nii.gz,50_Right-Caudate.label.nii.gz Putamen,12_Left-Putamen.label.nii.gz,51_Right-Putamen.label.nii.gz)
set struct = `echo $striatal_part|cut -d\, -f1`
set lh = `echo $striatal_part|cut -d\, -f2`
set rh = `echo $striatal_part|cut -d\, -f3`
fslmaths $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_${rh} -add $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_${lh} $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_${struct}.nii.gz
end
fslmaths $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_Accumbens.nii.gz -add $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_Caudate.nii.gz -add $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_Putamen.nii.gz $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_Striatum.nii.gz
echo "$SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_bilateral_Striatum.nii.gz" >>! $SUBJECTS_DIR/${subj}/ALLTARGETS.txt
echo "$SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_16_Brain-Stem.label.nii.gz" >>! $SUBJECTS_DIR/${subj}/ALLTARGETS.txt
fslmaths $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_4_Left-Lateral-Ventricle.label.nii.gz -add $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_43_Right-Lateral-Ventricle.label.nii.gz -add $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_14_3rd-Ventricle.label.nii.gz -add $SUBJECTS_DIR/${subj}/dlabel/diff/${subj}_15_4th-Ventricle.label.nii.gz $SUBJECTS_DIR/${subj}/dlabel/diff/all-ventricles.nii.gz
#### PROBABALISTIC TRACTOGRAPHY MAPS ####
##### Making a Seed (Amygdala) to Target ( 23 related structures ) divide by waytotal, threshold at .10 and apply logic (Not sure how well this works but this seems to be what she's done in her script)
mkdir ${subj}/dpath/Segmentation
mkdir ${subj}/dpath/Segmentation/Right
mkdir ${subj}/dpath/Segmentation/Left
/media/Data/fsl/bin/probtrackx2  -x /media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/${subj}_54_Right-Amygdala.label.nii.gz  -l --modeuler --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --avoid=/media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/all-ventricles.nii.gz --forcedir --opd -s /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/merged -m /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/nodif_brain_mask  --dir=/media/Data/Temple/Temple_Tracula/${subj}/dpath/Segmentation/Right --targetmasks=/media/Data/Temple/Temple_Tracula/${subj}/dpath/Segmentation/Right/targets.txt --os2t & 
/media/Data/fsl/bin/probtrackx2  -x /media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/${subj}_18_Left-Amygdala.label.nii.gz  -l --modeuler --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --avoid=/media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/all-ventricles.nii.gz --forcedir --opd -s /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/merged -m /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/nodif_brain_mask  --dir=/media/Data/Temple/Temple_Tracula/${subj}/dpath/Segmentation/Left --targetmasks=/media/Data/Temple/Temple_Tracula/${subj}/dpath/Segmentation/Left/targets.txt --os2t &

#       AMYGDALA SUBSEGMENTATION LOGIC SCRIPT          #

foreach hemi (Left Right)
foreach masks (seeds_to_${subj}_16_Brain-Stem.label.nii.gz seeds_to_${subj}_bilateral_Parahippocampus.nii.gz  seeds_to_${subj}_bilateral_Accumbens.nii.gz seeds_to_${subj}_bilateral_Pericalcarine.nii.gz  seeds_to_${subj}_bilateral_Caudal_Middle-Frontal.nii.gz seeds_to_${subj}_bilateral_Postcentral.nii.gz  seeds_to_${subj}_bilateral_Caudate.nii.gz seeds_to_${subj}_bilateral_Putamen.nii.gz  seeds_to_${subj}_bilateral_Cuneus.nii.gz seeds_to_${subj}_bilateral_Rostral_Anterior_Cingulate.nii.gz  seeds_to_${subj}_bilateral_Fusiform.nii.gz seeds_to_${subj}_bilateral_Striatum.nii.gz  seeds_to_${subj}_bilateral_Hippocampus.nii.gz seeds_to_${subj}_bilateral_Superior_Frontal.nii.gz  seeds_to_${subj}_bilateral_Inferior_Temporal.nii.gz seeds_to_${subj}_bilateral_﻿Superior_Parietal.nii.gz  seeds_to_${subj}_bilateral_Insula.nii.gz seeds_to_${subj}_bilateral_Superior_Temporal.nii.gz  seeds_to_${subj}_bilateral_Lateral_Occipital.nii.gz seeds_to_${subj}_bilateral_Temporal_Pole.nii.gz  seeds_to_${subj}_bilateral_Lateral_Orbitofrontal.nii.gz seeds_to_${subj}_bilateral_Thalamus_Proper.nii.gz  seeds_to_${subj}_bilateral_Lingual.nii.gz seeds_to_${subj}_bilateral_Ventral_Diencephalon.nii.gz  seeds_to_${subj}_bilateral_Medial_Orbitofrontal.nii.gz)
fslmaths ${subj}/dpath/Segmentation/${hemi}/${masks} -div `fslstats ${subj}/dpath/Segmentation/${hemi}/${masks} -R|cut -d\  -f2` -thr 0.1 -bin ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_${masks}
end
end

foreach hemi (Left Right)
## Central###
### BrainStem & Ventral Diecephalon & Thalamus Proper
### Overlap between BrainStem & Ventral Diecephalon & Thalamus Proper
### fslmaths BrainStem -add Ventral Diecephalon -add Thalamus Proper -thr 3 -bin Amygdala_Central_Nucleus.nii.gz

fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_16_Brain-Stem.label.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Ventral_Diencephalon.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Thalamus_Proper.nii.gz -thr 3 -bin ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}_Amygdala_Central_Nucleus.nii.gz

## Medial###
### ~(Brain Stem & Ventral Diencephalon & Thalamus Proper ) & (Ventral Diencephalon & (Striatum | Hippocampus))
### HINT NOT Central Nucleus --- Run that First
###  Ventral Diencephalon - Use the classification mask for Ventral Diencephalon output and take the overlap with the following mask
###   Striatum OR Hippocampus - Make additive mask for overlap with the ventral diencephalon and that out will be compared to central nucleus exclusion mask
###  MASK ONE (Additive Mask) -  fslmaths Striatum -add Hippocampus -bin MedialAdditive_Mask_1.nii.gz
###  MASK TWO (Overlap Mask) - flsmaths Ventral Diencephalon -add MedialAdditive_Mask_1.nii.gz -thr 2 -bin MedialOverlap_Mask_1.nii.gz
###  MASK THREE (Exclusion Applied) - flsmaths MedialOverlap_Mask_1.nii.gz -sub Amygdala_Central_Nucleus.nii.gz Amygdala_Medial_Nucleus.nii.gz


fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Striatum.nii.gz -add 31313422_10128/dpath/Segmentation/Right/DivThrBin_Right_seeds_to_31313422_10128_bilateral_Hippocampus.nii.gz -bin ${subj}/dpath/Segmentation/${hemi}/MedialAdditive_Mask_1.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Ventral_Diencephalon.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/MedialAdditive_Mask_1.nii.gz -thr 2 -bin ${subj}/dpath/Segmentation/${hemi}/MedialOverlap_Mask_1.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/MedialOverlap_Mask_1.nii.gz -sub ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}_Amygdala_Central_Nucleus.nii.gz ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}_Amygdala_Medial_Nucleus.nii.gz

## Lateral###
### ~(Superior parietal | Postcentral | Medial orbitofrontal | Lateral occipital | Pericalcarine | Cuneus) & (Temporal pole | Fusiform | Lateral orbitofrontal & (Superior temporal |Inferior Temporal )))
### NOT Superior parietal OR Postcentral OR Medial orbitofrontal OR Lateral occipital OR Pericalcarine OR Cuneus - Make an additive binary mask then use as a subrtaction from the other binarized masks
###  AND Temporal pole OR Fusiform OR Lateral orbitofrontal - Make an Additive mask that will be compared to the next additive masks and take the overlap 
###  AND Superior temporal OR Inferior Temporal - Make an Additive mask that will be compared to the next additive masks and take the overlap
### MASK ONE (Additive Exclusion Mask) : flsmaths  Superior parietal -add Postcentral -add Medial orbitofrontal -add Lateral occipital -add Pericalcarine -add Cuneus -bin EXCLUSION_MASK.nii
### MASK TWO (Additive For Overlap Mask 1) : flsmaths  Temporal Pole -add Fusiform -add Lateral orbitofrontal -bin LateralOverlap_1_MASK.nii
### MASK THREE (Additive For Overlap Mask 2) : flsmaths  Superior temporal -add Inferior Temporal -bin LateralOverlap_2_MASK.nii
### MASK FOUR (OVERLAP OF TWO and THREE FOR "AND" logic) : fslmaths LateralOverlap_1_MASK.nii -add LateralOverlap_2_MASK.nii -thr 2 -bin LateralOverlap_1_2_MASK.nii
### MASK FIVE (Apply Final Exclusion Mask) : fslmaths LateralOverlap_1_2_MASK.nii -sub EXCLUSION_MASK.nii Amygdala_Lateral_Nucleus.nii.gz

fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_﻿Superior_Parietal.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Postcentral.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Medial_Orbitofrontal.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Lateral_Occipital.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Pericalcarine.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Cuneus.nii.gz -bin ${subj}/dpath/Segmentation/${hemi}/LateralExclusion_Mask.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Temporal_Pole.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Fusiform.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Lateral_Orbitofrontal.nii.gz -bin ${subj}/dpath/Segmentation/${hemi}/LateralOverlap_1_MASK.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Superior_Temporal.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Inferior_Temporal.nii.gz -bin ${subj}/dpath/Segmentation/${hemi}/LateralOverlap_2_MASK.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/LateralOverlap_1_MASK.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/LateralOverlap_2_MASK.nii.gz -thr 2 -bin ${subj}/dpath/Segmentation/${hemi}/LateralOverlap_1_2_MASK.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/LateralOverlap_1_2_MASK.nii.gz -sub ${subj}/dpath/Segmentation/${hemi}/LateralExclusion_Mask.nii.gz ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}_Amygdala_Lateral_Nucleus.nii.gz

## Basal###
### (Parahippocampus & (Hippocampus | Rostral anterior cingulate | Lateral orbitofrontal | Medial orbitofrontal | caudal middle-frontal |Lateral occipital | Pericalcaraine | Cuneus | Lingual)) | (Insula & (Accumbens | Superior frontal))
### Parahippocampus  - Use the classification mask for Parahippocampus output and take the overlap with the following mask
###  Hippocampus OR Rostral anterior cingulate OR Lateral orbitofrontal OR Medial orbitofrontal OR caudal middle-frontal OR Lateral occipital OR Pericalcaraine OR Cuneus OR Lingual - Make an Additive mask that will be compared to the previous additive masks for overlap
### OR
### Insula  - Use the classification mask for Parahippocampus output and take the overlap with the following mask
### Accumbens OR  Superior frontal - Make an Additive mask that will be compared to the previous additive masks for overlap
###  MASK ONE (Additive Mask 1 to later overlap with parahipp): fslmaths Hippocampus -add Rostral anterior cingulate -add Lateral orbitofrontal -add Medial orbitofrontal -add caudal middle-frontal -add Lateral occipital -add Pericalcaraine -add Cuneus -add Lingual -bin Basal_Additive_1.nii.gz
###  MASK TWO (Additive Mask 2 overlap with insula): fslmaths Accumbens -add Superior frontal -bin Basal_Additive_2.nii.gz
###  MASK THREE (Additive Mask 3 overlap): fslmaths Parahippocampus -add Basal_Additive_1.nii.gz -thr 2 -bin Parahipp_Basal.nii.gz
###  MASK FOUR (Additive Mask 4 overlap): fslmaths Insula -add Basal_Additive_2.nii.gz -thr 2 -bin Insula_Basal.nii.gz
###  MASK FIVE (Additive Mask 4 overlap): fslmaths Parahipp_Basal.nii.gz -add Basal_Additive_2.nii.gz -bin Amygdala_Basal_Nucleus.nii.gz

fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Hippocampus.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Rostral_Anterior_Cingulate.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Lateral_Orbitofrontal.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Medial_Orbitofrontal.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Caudal_Middle-Frontal.nii.gz -add  ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Lateral_Occipital.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Pericalcarine.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Cuneus.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Lingual.nii.gz -bin ${subj}/dpath/Segmentation/${hemi}/Basal_Additive_1.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Accumbens.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Superior_Frontal.nii.gz -bin ${subj}/dpath/Segmentation/${hemi}/Basal_Additive_2.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Parahippocampus.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/Basal_Additive_1.nii.gz -thr 2 -bin ${subj}/dpath/Segmentation/${hemi}/Parahipp_Basal.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/DivThrBin_${hemi}_seeds_to_${subj}_bilateral_Insula.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/Basal_Additive_2.nii.gz -thr 2 -bin ${subj}/dpath/Segmentation/${hemi}/Insula_Basal.nii.gz
fslmaths ${subj}/dpath/Segmentation/${hemi}/Parahipp_Basal.nii.gz -add ${subj}/dpath/Segmentation/${hemi}/Insula_Basal.nii.gz ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}_Amygdala_Basal_Nucleus.nii.gz

mkdir ${subj}/dpath/Segmentation/${hemi}/SeedToMask
mkdir ${subj}/dpath/Segmentation/${hemi}/DivThrBin
mkdir ${subj}/dpath/Segmentation/${hemi}/FinalSubNuclei

mv ${subj}/dpath/Segmentation/${hemi}/seeds_to_*nii.gz ${subj}/dpath/Segmentation/${hemi}/SeedToMask
mv ${subj}/dpath/Segmentation/${hemi}/DivThrBin*nii.gz ${subj}/dpath/Segmentation/${hemi}/DivThrBin
mv ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}*nii.gz ${subj}/dpath/Segmentation/${hemi}/FinalSubNuclei


#         CALCULATE THE VOLUME OF THE RESULTING VOLUMES          #
echo "ID,Hemi_SubNucleus, Voxel_Count, Volume (mm^3)" >>! subsegmentation_volume.csv
foreach SubNucleus (Basal_Nucleus Lateral_Nucleus Medial_Nucleus Central_Nucleus)
echo "${subj}_${SubNucleus}"
echo "${subj},${Hemi}_${SubNucleus},`fslstats ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}_Amygdala_${SubNucleus}.nii.gz -V|cut -d\  -f1`,`fslstats ${subj}/dpath/Segmentation/${hemi}/${subj}_${hemi}_Amygdala_${SubNucleus}.nii.gz -V|cut -d\  -f2`" >> subsegmentation_volume.csv
end
end
end

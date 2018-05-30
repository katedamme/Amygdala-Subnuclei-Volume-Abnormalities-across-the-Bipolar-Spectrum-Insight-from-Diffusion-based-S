# STEP 2: Probtract Diffusion classification By Katherine SF Damme
# May 24 2017
# ZM Saygin 2011 Logic  - written KSD
# In each subject, we calculated the connection probability (using FSL- FDT's probtrackX) from each amygdala voxel (seed) to all bilateral cortical and subcortical regions (targets), and normalized the distribution of probabilities for each seed voxel to [0,1] by dividing by the maximum probability. We then thresholded and binarized these results to exclude values below 0.1, such that every amygdaloid voxel contained a 0 or 1 for each target.
foreach subj (`cat $SUBJECTS_DIR/subjects.txt`)
#### PROBABALISTIC TRACTOGRAPHY MAPS ####
##### Making a Seed (Amygdala) to Target ( 23 related structures ) divide by waytotal, threshold at .10 
mkdir ${subj}/dpath/Segmentation
mkdir ${subj}/dpath/Segmentation/Right
mkdir ${subj}/dpath/Segmentation/Left
/media/Data/fsl/bin/probtrackx2  -x /media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/${subj}_54_Right-Amygdala.label.nii.gz  -l --modeuler --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --avoid=/media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/all-ventricles.nii.gz --forcedir --opd -s /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/merged -m /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/nodif_brain_mask  --dir=/media/Data/Temple/Temple_Tracula/${subj}/dpath/Segmentation/Right --targetmasks=$SUBJECTS_DIR/${subj}/ALLTARGETS.txt --os2t  

/media/Data/fsl/bin/probtrackx2  -x /media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/${subj}_18_Left-Amygdala.label.nii.gz  -l --modeuler --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --avoid=/media/Data/Temple/Temple_Tracula/${subj}/dlabel/diff/all-ventricles.nii.gz --forcedir --opd -s /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/merged -m /media/Data/Temple/Temple_Tracula/${subj}/dmri.bedpostX/nodif_brain_mask  --dir=/media/Data/Temple/Temple_Tracula/${subj}/dpath/Segmentation/Left --targetmasks=$SUBJECTS_DIR/${subj}/ALLTARGETS.txt --os2t 

end

from algorithm import save
from nibabel import freesurfer
import numpy as np

annotpath_l = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/333330/Native/lh.aparc.DKTatlas40.annot'
label_l, ctab_l, names_l = freesurfer.read_annot(annotpath_l)


# label_oc  = np.unique(label_l)
label_oc = [11, 13, 5, 21]
occipital_label_l = np.zeros(label_l.shape)
for i in label_oc:
    occipital_label_l_index = np.where((label_l == i))
    occipital_label_l[occipital_label_l_index] = i

outpath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/333330/ROI/lh_DKTatlas_lateral_occipital_4.gii'
save.save_gifti(occipital_label_l,outpath,'CortexLeft')
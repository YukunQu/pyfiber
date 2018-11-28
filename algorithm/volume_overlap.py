import numpy as np
import copy
import nibabel as nib


def volume_overlap(img1, img2):
    """Calculate overlap between two volume files which are in same coordinate space.
    :parameter
    img1: NIfTI file pathway
    img2: NIfTI file pathway
    method : calculate method:pre;dice
             pre: precent of img2 in img1
             dice: dice index
    """
    img1_data = nib.load(img1).get_data()
    img2_data = nib.load(img2).get_data()
    overlap = img1_data*img2_data
    n = float(np.count_nonzero(overlap))
    pre = n / (np.count_nonzero(img1_data))
    return n,pre

img1 = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/l_vis_map.nii.gz"
img2 = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/l_vis_map_saclar_map_count.nii.gz"
num = volume_overlap(img1,img2)
print(num)


import os.path as op
import numpy as np
import nibabel as nib
from surfer import Brain

print(__doc__)

"""
Initialize the visualization.
"""
brain = Brain("fsaverage_sym", "lh", "inflated", background="white")

"""
Read both of the activation maps in using
surfer's io functions.
"""
sig_v1v = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/ep_v1v_101410.nii.gz").get_data().reshape(163842,1)
sig_v1d = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/ep_v1d_101410.nii.gz").get_data().reshape(163842,1)
sig_v2v = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/ep_v2v_101410.nii.gz").get_data().reshape(163842,1)
sig_v2d= nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/ep_v2d_101410.nii.gz").get_data().reshape(163842,1)
sig_v3v = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/ep_v3v_101410.nii.gz").get_data().reshape(163842,1)
sig_v3d = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/ep_v3d_101410.nii.gz").get_data().reshape(163842,1)
sig_other = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/ep_other_101410.nii.gz").get_data().reshape(163842,1)
"""
Zero out the vertices that do not meet a threshold.
"""

"""
A conjunction is defined as the minimum significance
value between the two maps at each vertex.
"""
conjunct = np.min(np.vstack((sig_v1v, sig_v1d,sig_v2v,sig_v2d,sig_v3v,sig_v3d)), axis=0)

"""
Now load the numpy array as an overlay.
Use a high saturation point so that the
blob will largely be colored with lower
values from the lookup table.
"""
brain.add_overlay(sig_v1v,min= 0, name="sig_v1v")

"""
A pointer to the overlay's color manager
gets stored in the overlays dictionary.
Change the lookup table to "Reds" and turn the
color bar itself off, as otherwise the bars
for the three maps will get confusingly stacked.
"""
brain.overlays["sig_v1v"].pos_bar.lut_mode = "Reds"
brain.overlays["sig_v1v"].pos_bar.visible = False

"""
Now load the other two maps and again change
the lookup table and turn off the color bar itself.
"""
brain.add_overlay(sig_v1d,min=2, name="sig_v1d")
brain.overlays["sig_v1d"].pos_bar.lut_mode = "Blues"
brain.overlays["sig_v1d"].pos_bar.visible = False

"""
Display the overlap as purple, which is what you
get when you mix red and blue.
"""
brain.add_overlay(conjunct,min=0, name="conjunct")
brain.overlays["conjunct"].pos_bar.lut_mode = "Purples"
brain.overlays["conjunct"].pos_bar.visible = False
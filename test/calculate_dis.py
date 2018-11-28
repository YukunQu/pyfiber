import numpy as np
import nibabel as nib
from nibabel.affines import apply_affine
from scipy.spatial.distance import cdist

#The program calculates and plots the curve of distances between fiebr endpoints and gray matter.

f_pathway = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k.tck'
r_pathway = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/ROI/ribbon_binary.nii.gz'

streamlines = nib.streamlines.load(f_pathway).streamlines
roi = nib.load(r_pathway)
affine = roi.affine
roi_data = roi.get_data()

#extract the endpoints
s = [[0,0,0],[0,0,0]]
endpoints= np.array(s)
for streamline in streamlines:
    s = [streamline[0], streamline[-1]]
    endpoints = np.concatenate((endpoints,s))
endpoints = np.delete(endpoints,[0,1],axis=0)
endpoints_1 = endpoints[::2,:]
endpoints_1.astype('float16')

#calculate the distance between fiber endpoint and gray matter
roi_coords = np.array(np.where(roi_data)).T
roi_coords_a = apply_affine(affine, roi_coords)
roi_coords_a.astype('float16')
dist = cdist(endpoints_1, roi_coords_a,'euclidean')

#plot (may be not use curve rather than histogram)
dist_r = np.rint(dist)
dist_min = np.min(dist)
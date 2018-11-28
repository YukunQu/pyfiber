import os
from scipy.spatial.distance import cdist
import nibabel as nib
from nibabel import  streamlines
from nibabel import freesurfer
import numpy as np
from dipy.tracking import streamline

# read data:fiber, geometry, index)

roi_path = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/ROI/benson_retinotopic_map/"
retinotopic_areas_list = ['lh_v1v_101410.nii.gz',  'lh_v1d_101410.nii.gz', 'lh_v2v_101410.nii.gz',
                          'lh_v2d_101410.nii.gz',  'lh_v3v_101410.nii.gz', 'lh_v3d_101410.nii.gz']

retinotopic_areas = [roi_path+area for area in retinotopic_areas_list]
fiber_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/'
out_fiber_list = [fiber_path+area.replace('.nii.gz', '.tck') for area in retinotopic_areas_list]

for i in range(len(retinotopic_areas)):
    label = nib.load(retinotopic_areas[i]).get_data()
    label_index = np.where(label)
    label_index = label_index[0]
    white_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/surf/lh.white'
    surf_coords, surf_faces = freesurfer.read_geometry(white_path)
    roi_coords = surf_coords[label_index]

    ori_fiber = nib.streamlines.TckFile.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion"
                "/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k_rh_occipital.tck")
    ori_fiber_streamlines = ori_fiber.streamlines

    #extract occipital terminus of streamlines
    s_0 = [s[0] for s in ori_fiber_streamlines]
    s_1 = [s[-1] for s in ori_fiber_streamlines]
    s_oc = []
    for s in range(len(s_0)):
        if s_0[s][1] < s_1[s][1]:
            s_oc.append(s_0[s])
        else:s_oc.append(s_1[s])
    stream_terminus_oc = np.array(s_oc)  # stream_terminus[stream_terminus[:, 0] < 0]

    # select fiber by rois
    dist = cdist(roi_coords, stream_terminus_oc)
    stream_index = np.array(len(ori_fiber_streamlines) * [False])

    for d in range(len(dist[:])):
        temp_index = np.array(dist[d] < 1.1)
        stream_index += temp_index
    rois_streamlines = ori_fiber_streamlines[stream_index]

    # save fiber
    header = ori_fiber.header
    data_per_streamline = ori_fiber.tractogram.data_per_streamline
    data_per_point = ori_fiber.tractogram.data_per_point
    affine_to_rasmm = ori_fiber.tractogram.affine_to_rasmm

    save_streamline = streamlines.array_sequence.ArraySequence(rois_streamlines)
    tractogram = streamlines.tractogram.Tractogram(streamlines=save_streamline, data_per_streamline=data_per_streamline,
                                                   data_per_point=data_per_point, affine_to_rasmm=affine_to_rasmm)
    datdat = streamlines.tck.TckFile(tractogram=tractogram, header=header)
    datdat.save(out_fiber_list[i])
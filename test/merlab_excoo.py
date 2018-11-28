from nibabel import freesurfer
from dipy.viz import window, actor
from dipy.viz.colormap import line_colors
import nibabel as nib
import numpy as np
import os
from scipy.spatial.distance import cdist
import nibabel.streamlines.array_sequence as nibas
from nibabel.spatialimages import ImageFileError
from nibabel import streamlines as ns


def _sort_streamlines(fasciculus_data):
    """Store order of streamline is from left to right."""
    fasciculus_data_sort = nibas.ArraySequence()
    for i in range(len(fasciculus_data)):
        if fasciculus_data[i][0][0] < 0:
            fasciculus_data_sort.append(fasciculus_data[i])
        else:
            fasciculus_data_sort.append(fasciculus_data[i][::-1])
    return fasciculus_data_sort

def separation_fib_to_hemi(data):
    """Separating a bundle fiber to both hemispheres"""

    streamlines = data

    streamlines = _sort_streamlines(streamlines)
    fib_lh = nibas.ArraySequence()
    fib_rh = nibas.ArraySequence()

    for i in range(len(streamlines)):
        l = streamlines[i][:, 0]
        l_ahead = list(l[:])
        a = l_ahead.pop(0)
        l_ahead.append(a)
        x_stemp = np.array([l, l_ahead])
        x_stemp_index = x_stemp.prod(axis=0)
        index0 = np.argwhere(x_stemp_index <= 0)
        if len(index0) != 0:
            index_term = np.argmin((abs(streamlines[i][index0[0][0]][0]),
                                    abs(streamlines[i][index0[0][0] + 1][0])))
            index = index0[0][0] + index_term
            fib_lh.append(streamlines[i][:index + 1])
            fib_rh.append(streamlines[i][index:])
        else:
            if streamlines[i][0][0] <= 0:
                fib_lh.append(streamlines[i])
            else:
                fib_rh.append(streamlines[i])

    return fib_lh, fib_rh


def select_by_surf_rois(streamlines_ori, surf_rois, geo_path):
    """
    Include or exclude the streamlines according to some surface ROIs
    Parameters
    ----------
    streamlines_ori: origin streamlines
    surf_rois: left and right surface rois
    geo_path: left and right surface geometry
    Return
    ------
    lh_rois_streamlines: extracted streamlines by left roi
    rh_rois_streamlines: extracted streamlines by right roi
    lrh_rois_streamlines: intersection of lh_rois_streamlines and rh_rois_streamlines
    """
    streamlines = _sort_streamlines(streamlines_ori)
    s0 = [s[0] for s in streamlines]
    s_t = [s[-1] for s in streamlines]
    # s = s0 + s_t
    # stream_terminus = np.array(s)
    stream_terminus_lh = np.array(s0)  # stream_terminus[stream_terminus[:, 0] < 0]
    stream_terminus_rh = np.array(s_t)  # stream_terminus[stream_terminus[:, 0] > 0]

    # load surf geometry
    suffix = os.path.split(geo_path[0])[1].split('.')[-1]
    if suffix in ('white', 'inflated', 'pial'):
        coords_lh, faces_lh = nib.freesurfer.read_geometry(geo_path[0])
    elif suffix == 'gii':
        gii_data = nib.load(geo_path[0]).darrays
        coords_lh, faces_lh = gii_data[0].data, gii_data[1].data
    else:
        raise ImageFileError('This file format-{} is not supported at present.'.format(suffix))

    suffix = os.path.split(geo_path[1])[1].split('.')[-1]
    if suffix in ('white', 'inflated', 'pial'):
        coords_rh, faces_rh = nib.freesurfer.read_geometry(geo_path[1])
    elif suffix == 'gii':
        gii_data = nib.load(geo_path[1]).darrays
        coords_rh, faces_rh = gii_data[0].data, gii_data[1].data
    else:
        raise ImageFileError('This file format-{} is not supported at present.'.format(suffix))

    dist_lh = cdist(coords_lh[surf_rois[0]], stream_terminus_lh)
    lh_stream_index = np.array(len(streamlines) * [False])
    for i in range(len(dist_lh[:])):
        temp_index = np.array(dist_lh[i] <= 1.1)
        lh_stream_index += temp_index

    dist_rh = cdist(coords_rh[surf_rois[1]], stream_terminus_rh)
    rh_stream_index = np.array(len(streamlines) * [False])
    for j in range(len(dist_rh[:])):
        temp_index = np.array(dist_rh[j] <= 1.1)
        rh_stream_index += temp_index

    lh_rois_streamlines = streamlines_ori[lh_stream_index]
    rh_rois_streamlines = streamlines_ori[rh_stream_index]
    lrh_rois_streamlines = streamlines_ori[np.array([lh_stream_index, rh_stream_index]).sum(axis=0) == 2]
    if len(lrh_rois_streamlines) == 0:
        print "ROI-1 to ROI-2 have no connection!"
        return lh_rois_streamlines, rh_rois_streamlines,lrh_rois_streamlines
    else:
        return lh_rois_streamlines, rh_rois_streamlines, lrh_rois_streamlines


#get label_index
annotpath_l = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Native/lh.aparc.DKTatlas40.annot'
label_l, ctab_l, names_l = freesurfer.read_annot(annotpath_l)
annotpath_r = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Native/rh.aparc.DKTatlas40.annot'
label_r, ctab_r, names_r = freesurfer.read_annot(annotpath_r)

occipital_label_l = np.where(((label_l ==11) | (label_l==13) | (label_l==5) | (label_l==21)))
occipital_label_r = np.where(((label_r ==11) | (label_r==13) | (label_r==5) | (label_r==21)))
occipital_label = list(occipital_label_l + occipital_label_r)

white_path_L = ['/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Native/101410.L.white.native.surf.gii']
white_path_R = ['/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Native/101410.R.white.native.surf.gii']
white_path = tuple(white_path_L + white_path_R)

ori_fiber = nib.streamlines.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410"
                                   "/Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k_rh.tck")

out_fiber = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410" \
            "/Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k_rh_occipital.tck"

streamlines = ori_fiber.streamlines
header = ori_fiber.header
data_per_streamline = ori_fiber.tractogram.data_per_streamline
data_per_point = ori_fiber.tractogram.data_per_point
affine_to_rasmm = ori_fiber.tractogram.affine_to_rasmm

lh_rois_streamlines, rh_rois_streamlines, lrh_rois_streamlines = select_by_surf_rois(streamlines_ori=streamlines,surf_rois=occipital_label,geo_path=white_path)
lh_rois_streamlines,lh_rois_streamlines_inrh = separation_fib_to_hemi(lh_rois_streamlines)


tractogram = ns.tractogram.Tractogram(streamlines=lh_rois_streamlines, data_per_streamline=data_per_streamline,
                                               data_per_point=data_per_point, affine_to_rasmm=affine_to_rasmm)
datdat = ns.tck.TckFile(tractogram=tractogram, header=header)
datdat.save(out_fiber)

#visualization
#
# T1w = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/T1w_acpc_dc_restore_brain1.25.nii.gz")
# img_fiber_1 = lh_rois_streamlines
# #img_fiber_2 = rh_rois_streamlines
# #img_fiber_3 = lrh_rois_streamlines
# affine = np.linalg.inv(T1w.affine)
#
# interactive = True
#
# color = line_colors(img_fiber_1)
# cc_streamlines_actor1 = actor.line(img_fiber_1, color)
# #color = line_colors(img_fiber_2)
# #cc_streamlines_actor2 = actor.line(img_fiber_2, color)
# #color = line_colors(img_fiber_3)
# #cc_streamlines_actor3 = actor.line(img_fiber_3, color)
#
# img_T1w_data = T1w.get_data()
#
# vol_actor = actor.slicer(img_T1w_data,T1w.affine)
# vol_actor.display(x=72)
# vol_actor2 = vol_actor.copy()
# vol_actor2.display(z=40)
#
# r = window.Renderer()
# r.add(vol_actor)
# #r.add(vol_actor2)
# r.add(cc_streamlines_actor1)
# # r.add(cc_streamlines_actor2)
# # r.add(cc_streamlines_actor3)
#
# if interactive:
#     window.show(r)
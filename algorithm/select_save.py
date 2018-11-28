from nibabel import streamlines
from dipy.tracking import streamline
import algorithm.volume_extract
from nibabel import freesurfer
from dipy.viz import window, actor
from dipy.viz.colormap import line_colors
import nibabel as nib
import numpy as np
import os
from scipy.spatial.distance import cdist
import nibabel.streamlines.array_sequence as nibas
from nibabel.spatialimages import ImageFileError


def save_selectfibers(Ori_fiber,Out_fiber,roi1,roi2):
    # roi1 : roi coordinate             format:nii.gz
    # roi2 : roi coordinate             format:nii.gz
    # Ori_fiber : Original fiber        format:tck
    # Out_fiber : selected fiber name   format:tck
    #load(pathname)

    # load
    roi1 = nib.load(roi1)
    roi2 = nib.load(roi2)
    Ori_fiber = streamlines.TckFile.load(Ori_fiber)

    roi_affine = roi1.affine
    roi1_data = roi1.get_data().astype(bool)
    roi2_data = roi2.get_data().astype(bool)

    stream_lines = Ori_fiber.streamlines
    header = Ori_fiber.header
    data_per_streamline = Ori_fiber.tractogram.data_per_streamline
    data_per_point = Ori_fiber.tractogram.data_per_point
    affine_to_rasmm = Ori_fiber.tractogram.affine_to_rasmm

    #select
    roi1_streamlines = streamline.select_by_rois(stream_lines,[roi1_data],[True],mode='any',affine=roi_affine)
    streamlines1 = list(roi1_streamlines)
    roi2_streamlines = streamline.select_by_rois(streamlines1,[roi2_data],[True],mode='any',affine=roi_affine)
    streamlines2 = list(roi2_streamlines)

    #save
    save_streamline = streamlines.array_sequence.ArraySequence(streamlines2)

    tractogram = streamlines.tractogram.Tractogram(streamlines=save_streamline,data_per_streamline=data_per_streamline,
                                                   data_per_point=data_per_point,affine_to_rasmm=affine_to_rasmm)
    datdat = streamlines.tck.TckFile(tractogram=tractogram,header=header)
    datdat.save(Out_fiber)


def select_by_prob(ori_fiber,prob_atlas,out_fiber,volume,th):
    # Ori_fiber : Original fiber                        format:tck
    # prob_atlas: probabilistic atlas as reference       format:tck
    # Out_fiber : selected fiber name                   format:tck
    # volume :to select special fiber in probabilistic atlas, [1,20]
    # th : to select higher probability fiber

    #load(pathway)
    prob_atlas = nib.load(prob_atlas)
    affine = prob_atlas.affine
    prob_atlas_data = prob_atlas.get_data()

    ori_fiber = streamlines.TckFile.load(ori_fiber)
    stream_lines = ori_fiber.streamlines
    header = ori_fiber.header
    data_per_streamline = ori_fiber.tractogram.data_per_streamline
    data_per_point = ori_fiber.tractogram.data_per_point
    affine_to_rasmm = ori_fiber.tractogram.affine_to_rasmm

    #extract special fiber probabilistic atlas
    special_atlas_data = algorithm.volume_extract.extract_label(prob_atlas_data, volume, th).astype(bool)

    #select fiber by probabilistic atlas
    selected_streamlines = streamline.select_by_rois(stream_lines, [special_atlas_data], [True], mode='any',
                                                     affine=affine)
    selected_fiber = list(selected_streamlines)

    #save
    save_streamline = streamlines.array_sequence.ArraySequence(selected_fiber)

    tractogram = streamlines.tractogram.Tractogram(streamlines=save_streamline, data_per_streamline=data_per_streamline,
                                                   data_per_point=data_per_point, affine_to_rasmm=affine_to_rasmm)
    datdat = streamlines.tck.TckFile(tractogram=tractogram, header=header)
    datdat.save(out_fiber)


def save_selectfiber(Ori_fiber,Out_fiber,roi):
    # roi : roi coordinate             format:nii.gz
    # Ori_fiber : Original fiber        format:tck
    # Out_fiber : selected fiber name   format:tck
    #load(pathname)

    # load
    roi = nib.load(roi)
    roi_affine = roi.affine
    roi_data = roi.get_data().astype(bool)

    Ori_fiber = streamlines.TckFile.load(Ori_fiber)
    stream_lines = Ori_fiber.streamlines
    header = Ori_fiber.header
    data_per_streamline = Ori_fiber.tractogram.data_per_streamline
    data_per_point = Ori_fiber.tractogram.data_per_point
    affine_to_rasmm = Ori_fiber.tractogram.affine_to_rasmm

    # select
    roi_streamlines = streamline.select_by_rois(stream_lines, [roi_data], [True], mode='either_end', affine=roi_affine)
    fibers = list(roi_streamlines)

    #save
    save_streamline = streamlines.array_sequence.ArraySequence(fibers)
    tractogram = streamlines.tractogram.Tractogram(streamlines=save_streamline,data_per_streamline=data_per_streamline,
                                                   data_per_point=data_per_point,affine_to_rasmm=affine_to_rasmm)
    datdat = streamlines.tck.TckFile(tractogram=tractogram,header=header)
    datdat.save(Out_fiber)
    return 'Complete Done.'

Ori_fiber = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/l_vis.tck'
Out_fiber = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/l_vis_endinribbon.tck'
roi = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/ROI/ribbon_binary.nii.gz'

save_selectfiber(Ori_fiber,Out_fiber,roi)

# prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects'
# subjects = ['597869','102008','200008']
#
# for subject in subjects:
#     Ori_fiber = os.path.join(prepath,subject,'Diffusion/tractography/Det/'
#                                              'SD_Stream_angle20_cutoff0.1_length50_250_seedFP_100k.tck')
#     Out_fiber1 = os.path.join(prepath,subject,'Diffusion/tractography/Det/'
#                                               'SD_Stream_angle20_cutoff0.1_length50_250_seedFP_roiOccipital_L_100k.tck')
#     Out_fiber2 = os.path.join(prepath,subject,'Diffusion/tractography/Det/'
#                                               'SD_Stream_angle20_cutoff0.1_length50_250_seedFP_roiOccipital_R_100k.tck')
#     roi1 = os.path.join(prepath,subject,'ROI/L_Occipital.nii.gz')
#     roi2 = os.path.join(prepath,subject,'ROI/R_Occipital.nii.gz')
#     save_selectfiber(Ori_fiber,Out_fiber1,roi1)
#     save_selectfiber(Ori_fiber,Out_fiber2,roi2)


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

def select_by_hemi_surf_rois(streamlines_ori, surf_rois, geo_path):
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

    # load surf geometry
    suffix = os.path.split(geo_path)[1].split('.')[-1]
    if suffix in ('white', 'inflated', 'pial'):
        coords, faces = nib.freesurfer.read_geometry(geo_path)
    elif suffix == 'gii':
        gii_data = nib.load(geo_path).darrays
        coords, faces = gii_data[0].data, gii_data[1].data
    else:
        raise ImageFileError('This file format-{} is not supported at present.'.format(suffix))

    rois_streamlines = streamline.select_by_rois(streamlines_ori,[coords[surf_rois]],[True],mode="either_end")

    return rois_streamlines


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
    sort_streamlines = _sort_streamlines(streamlines_ori)

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


    lh_rois_streamlines = streamline.select_by_rois(sort_streamlines.data,[coords_lh[surf_rois[0]]],[True],mode="either_end")
    rh_rois_streamlines = streamline.select_by_rois(sort_streamlines.data, [coords_rh[surf_rois[1]]], [True], mode="either_end")

    return lh_rois_streamlines, rh_rois_streamlines

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

streamlines_fiber = ori_fiber.streamlines
header = ori_fiber.header
data_per_streamline = ori_fiber.tractogram.data_per_streamline
data_per_point = ori_fiber.tractogram.data_per_point
affine_to_rasmm = ori_fiber.tractogram.affine_to_rasmm

lh_rois_streamlines, rh_rois_streamlines, lrh_rois_streamlines = select_by_surf_rois(streamlines_ori=streamlines_fiber,surf_rois=occipital_label,geo_path=white_path)
lh_rois_streamlines, lh_rois_streamlines_inrh = separation_fib_to_hemi(lh_rois_streamlines)

# #visualization
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
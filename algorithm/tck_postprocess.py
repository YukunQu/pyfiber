import os
import numpy as np
import nibabel as nib
from nibabel import streamlines
from nibabel import freesurfer
from scipy.spatial.distance import cdist
import nibabel.streamlines.array_sequence as nibas

def single_point_mid_sag(fiber):
    """
    A single time through the sagittal plane in the middle.
    The function is implemented after function endpoint_dissimilarity
    """
    fasciculus_data = fiber.streamlines
    index = len(fasciculus_data) * [False]   #dengtongyu xianwei shumu de suoyi2
    for i in range(len(fasciculus_data)):    #circulate for every streamline
        l = fasciculus_data[i][:, 0]         #l =  x coordinate of all points of the fiber
        l_ahead = list(l[:])                 #listful
        a = l_ahead.pop(0)                   # a = first x
        l_ahead.append(a)                    # add the first x to end
        x_stemp = np.array([l, l_ahead])
        x_stemp_index = x_stemp.prod(axis=0)
        if len(np.argwhere(x_stemp_index < 0)) == 2 or len(np.argwhere(x_stemp_index == 0)) == 2:
            index[i] = True
    fasciculus_data = fasciculus_data[index]
    return fasciculus_data


def _sort_streamlines(fasciculus_data):
    """Store order of streamline is from left to right."""
    fasciculus_data_sort = nibas.ArraySequence()
    for i in range(len(fasciculus_data)):
        if fasciculus_data[i][0][0] < 0:
            fasciculus_data_sort.append(fasciculus_data[i])
        else:
            fasciculus_data_sort.append(fasciculus_data[i][::-1])
    return fasciculus_data_sort


def separation_fib_to_hemi(fasciculus_data):
    """Separating a bundle fiber to both hemispheres"""

    streamlines = fasciculus_data
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


def select_by_surf_rois(streamlines_ori, surf_rois, geo_coords,threshold):
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
    s_0 = np.array([streamline[0] for streamline in streamlines_ori])
    s_1 = np.array([streamline[-1] for streamline in streamlines_ori])

    # select those fiber
    # that the distance between the either ends of the fibers and the coordinates is less than five centimeters.
    dist_0 = cdist(geo_coords[surf_rois], s_0)
    stream_index = np.array(len(streamlines_ori) * [False])
    for i in range(len(dist_0[:])):
        temp_index = np.array(dist_0[i] < threshold)
        stream_index += temp_index

    dist_1 = cdist(geo_coords[surf_rois], s_1)
    for j in range(len(dist_1[:])):
        temp_index = np.array(dist_1[j] < threshold)
        stream_index += temp_index

    rois_streamlines = streamlines_ori[stream_index]

    return rois_streamlines


def tck_postprocess(ori_fiber_path, annot_path, geo_coords_path,out_fiber_path,hemi = 'lh'):
    """Operate the fiber from the Mrtrix tractograhy to delete the non-target fibers.

     Parameters
    ----------
    ori_fiber:
    surf_roi: annot file
    geo_coords:
     """
    # load_fiber_rh,lh;white.surf.gii,surf_roi.annot
    ori_fiber = nib.streamlines.TckFile.load(ori_fiber_path)
    label, ctab, names = freesurfer.read_annot(annot_path)
    occipital_label = np.where(((label == 11) | (label == 13) | (label == 5) | (label == 21)))
    gii_data = nib.load(geo_coords_path).darrays
    coords, faces = gii_data[0].data, gii_data[1].data

    # postprocess : middle_sag, select by rois of surface, separation_fib_to_hemi)
    fiber_mid_sag = single_point_mid_sag(ori_fiber)
    rois_streamlines = select_by_surf_rois(fiber_mid_sag,occipital_label,coords)
    fiber_seg_lh, fiber_seg_rh = separation_fib_to_hemi(rois_streamlines)

    # save read from header file and save fiber_through_mid_sag_surf_roi_2hemi.
    header = ori_fiber.header
    data_per_streamline = ori_fiber.tractogram.data_per_streamline
    data_per_point = ori_fiber.tractogram.data_per_point
    affine_to_rasmm = ori_fiber.tractogram.affine_to_rasmm

    assert hemi in ['lh','rh']
    if hemi is 'lh':
        tractogram = streamlines.tractogram.Tractogram(streamlines=fiber_seg_lh, data_per_streamline=data_per_streamline,
                                                       data_per_point=data_per_point, affine_to_rasmm=affine_to_rasmm)
        datdat = streamlines.tck.TckFile(tractogram=tractogram, header=header)
        datdat.save(out_fiber_path)
    else:
        tractogram = streamlines.tractogram.Tractogram(streamlines=fiber_seg_rh, data_per_streamline=data_per_streamline,
                                                       data_per_point=data_per_point, affine_to_rasmm=affine_to_rasmm)
        datdat = streamlines.tck.TckFile(tractogram=tractogram, header=header)
        datdat.save(out_fiber_path)


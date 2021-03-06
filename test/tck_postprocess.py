import os
import numpy as np
import nibabel as nib
from nibabel import streamlines
from dipy.tracking import streamline
from nibabel import freesurfer
from scipy.spatial.distance import cdist
import nibabel.streamlines.array_sequence as nibas


def single_point_mid_sag(fiber):
    """
    A single time through the sagittal plane in the middle.
    The function is implemented after function endpoint_dissimilarity
    """
    fasciculus_data = fiber.streamlines
    index = len(fasciculus_data) * [False]  # dengtongyu xianwei shumu de suoyi2
    for i in range(len(fasciculus_data)):  # circulate for every streamline
        l = fasciculus_data[i][:, 0]  # l =  x coordinate of all points of the fiber
        l_ahead = list(l[:])  # listful
        a = l_ahead.pop(0)  # a = first x
        l_ahead.append(a)  # add the first x to end
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


def select_by_surf_rois(streamlines_ori, surf_rois, geo_coords, threshold):
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


def tck_postprocess(ori_fiber_path, annot_path, geo_coords_path, roi_path, out_fiber_path, threshold):
    """Operate the fiber from the Mrtrix tractograhy to delete the non-target fibers.

     Parameters
    ----------
    ori_fiber:
    surf_roi:  annot file
    geo_coords:
     """
    # load_fiber_rh,lh;white.surf.gii,surf_roi
    ori_fiber = nib.streamlines.TckFile.load(ori_fiber_path)
    ori_streamlines = ori_fiber.streamlines
    label, ctab, names = freesurfer.read_annot(annot_path)
    occipital_label = np.where(((label == 11) | (label == 13) | (label == 5) | (label == 21)))
    gii_data = nib.load(geo_coords_path).darrays
    coords, faces = gii_data[0].data, gii_data[1].data

    # load roi,affine
    roi_sp = nib.load(roi_path)
    roi_affine = roi_sp.affine
    roi_splenium_data= roi_sp.get_data().astype(bool)

    # postprocess :  select by rois of surface and splenium mask
    rois_oc_streamlines = select_by_surf_rois(ori_streamlines, occipital_label, coords, threshold)
    rois_sp_streamlines = streamline.select_by_rois(rois_oc_streamlines, [roi_splenium_data],[True],mode='any',affine=roi_affine)
    out_streamlines = list(rois_sp_streamlines)

    # save read from header file and save fiber
    header = ori_fiber.header
    data_per_streamline = ori_fiber.tractogram.data_per_streamline
    data_per_point = ori_fiber.tractogram.data_per_point
    affine_to_rasmm = ori_fiber.tractogram.affine_to_rasmm

    save_streamline = streamlines.array_sequence.ArraySequence(out_streamlines)
    tractogram = streamlines.tractogram.Tractogram(streamlines=save_streamline,
                                                   data_per_streamline=data_per_streamline,
                                                   data_per_point=data_per_point, affine_to_rasmm=affine_to_rasmm)
    datdat = streamlines.tck.TckFile(tractogram=tractogram, header=header)
    datdat.save(out_fiber_path)

# postprocessing batch

prepath = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects"
subjects = ['100206', '100307', '100408', '100610', '101006', '101107', '101309', '101410', '101915', '102008',
            '102109', '102311', '102513', '102614', '102715', '102816', '103010', '103111', '111211', '200008',
            '201818', '333330', '530635', '555651', '597869', '686969', '767464', '888678', '995174', '996782']

for subject in subjects:
    ori_fiber_lh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k.tck")
    ori_fiber_rh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k.tck")

    annot_lh_path = os.path.join(prepath, subject, "surf/lh.aparc.DKTatlas40.annot")
    annot_rh_path = os.path.join(prepath, subject, "surf/rh.aparc.DKTatlas40.annot")

    geo_coords_lh_path = os.path.join(prepath, subject, "surf/lh.white.surf.gii")
    geo_coords_rh_path = os.path.join(prepath, subject, "surf/rh.white.surf.gii")

    out_fiber_lh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc1.1.tck")
    out_fiber_rh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc1.1.tck")
    roi_path = os.path.join(prepath, subject, "ROI/splenium_mask.nii.gz")

    tck_postprocess(ori_fiber_lh_path, annot_lh_path, geo_coords_lh_path, roi_path, out_fiber_lh_path, 1.1)
    tck_postprocess(ori_fiber_rh_path, annot_rh_path, geo_coords_rh_path, roi_path, out_fiber_rh_path, 1.1)

prepath = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects"
subjects = ['100206', '100307', '100408', '100610', '101006', '101107', '101309', '101410', '101915', '102008',
            '102109', '102311', '102513', '102614', '102715', '102816', '103010', '103111', '111211', '200008',
            '201818', '333330', '530635', '555651', '597869', '686969', '767464', '888678', '995174', '996782']

for subject in subjects:
    ori_fiber_lh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k.tck")
    ori_fiber_rh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k.tck")

    annot_lh_path = os.path.join(prepath, subject, "surf/lh.aparc.DKTatlas40.annot")
    annot_rh_path = os.path.join(prepath, subject, "surf/rh.aparc.DKTatlas40.annot")

    geo_coords_lh_path = os.path.join(prepath, subject, "surf/lh.white.surf.gii")
    geo_coords_rh_path = os.path.join(prepath, subject, "surf/rh.white.surf.gii")

    out_fiber_lh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc3.tck")
    out_fiber_rh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc3.tck")
    roi_path = os.path.join(prepath, subject, "ROI/splenium_mask.nii.gz")

    tck_postprocess(ori_fiber_lh_path, annot_lh_path, geo_coords_lh_path, roi_path, out_fiber_lh_path, 3)
    tck_postprocess(ori_fiber_rh_path, annot_rh_path, geo_coords_rh_path, roi_path, out_fiber_rh_path, 3)

prepath = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects"
subjects = ['100206', '100307', '100408', '100610', '101006', '101107', '101309', '101410', '101915', '102008',
            '102109', '102311', '102513', '102614', '102715', '102816', '103010', '103111', '111211', '200008',
            '201818', '333330', '530635', '555651', '597869', '686969', '767464', '888678', '995174', '996782']

for subject in subjects:
    ori_fiber_lh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k.tck")
    ori_fiber_rh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k.tck")

    annot_lh_path = os.path.join(prepath, subject, "surf/lh.aparc.DKTatlas40.annot")
    annot_rh_path = os.path.join(prepath, subject, "surf/rh.aparc.DKTatlas40.annot")

    geo_coords_lh_path = os.path.join(prepath, subject, "surf/lh.white.surf.gii")
    geo_coords_rh_path = os.path.join(prepath, subject, "surf/rh.white.surf.gii")

    out_fiber_lh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc5.tck")
    out_fiber_rh_path = os.path.join(prepath, subject,
                                     "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc5.tck")
    roi_path = os.path.join(prepath, subject, "ROI/splenium_mask.nii.gz")

    tck_postprocess(ori_fiber_lh_path, annot_lh_path, geo_coords_lh_path, roi_path, out_fiber_lh_path, 5)
    tck_postprocess(ori_fiber_rh_path, annot_rh_path, geo_coords_rh_path, roi_path, out_fiber_rh_path, 5)

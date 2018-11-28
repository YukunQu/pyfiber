import os
from algorithm import save
import subprocess
import numpy as np
import nibabel as nib
from nibabel import streamlines
from dipy.tracking import streamline
from nibabel import freesurfer
from scipy.spatial.distance import cdist




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
    # that the distance between the either ends of the fibers and the coordinates is less than threshold.
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


def fiber_density_surfmap(fiber_path, white_path, output, hemi='lh', threshold=1.1):
    # aim to get a Nx3 array. N is fiber_num,3 is coordinates fiber endpoints(x,y,z)
    fiber = streamlines.TckFile.load(fiber_path)
    f_streamlines = fiber.streamlines

    s_0 = [streamline[0] for streamline in f_streamlines]
    s_1 = [streamline[-1] for streamline in f_streamlines]
    s_oc = []

    assert hemi in ['lh', 'rh']
    if hemi == 'lh':
        for n in range(len(s_0)):
            if s_0[n][0] < s_1[n][0]:
                s_oc.append(s_0[n])
            else:
                s_oc.append(s_1[n])
        s_oc = np.array(s_oc)
    else:
        for n in range(len(s_0)):
            if s_0[n][0] > s_1[n][0]:
                s_oc.append(s_0[n])
            else:
                s_oc.append(s_1[n])
        s_oc = np.array(s_oc)

    # aim to get a Mx3 array. M is num of occipital coordinates,3 is coordinates fiber endpoints(x,y,z)
    gii_data = nib.load(white_path).darrays
    white_coords = gii_data[0].data

    # aim to get a MxN matrix.M is num of occipital coordinates, N is fiber_num. scalar value is distance.
    dist = cdist(white_coords, s_oc)

    # aim to get a Mx1 matrix.M is num of occipital coordinates, N is the amount of endpoints meeting the threshold.

    dist_oc = dist < threshold
    dist_oc_num = dist_oc.sum(axis=1)  # the scalar data

    # save to a gifiti file
    if hemi == 'lh':
        save.save_gifti(data=dist_oc_num, save_name=output, hemisphere='CortexLeft')
    else:
        save.save_gifti(data=dist_oc_num, save_name=output, hemisphere='CortexRight')



def surf2surf(src_subject_name,src_val_path,target_subject_name,target_vol_path,hemi):
    #
    subprocess.call('mri_surf2surf --srcsubject {} --sval {} --trgsubject {} --trgsurfval {} --hemi {}'.format(
        src_subject_name,src_val_path,target_subject_name, target_vol_path,hemi),shell= True)



def cal_normal(input_name = None, output_name=None):
    prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/density_map/5/'
    density_map_list = os.listdir(prepath)  #acquire the density map of each subjects
    temp = []
    for density_map_name in density_map_list:
        density_map_path = os.path.join(prepath,density_map_name)
        density_map = nib.load(density_map_path).darrays
        density_map_data = density_map[0].data
        dmdmax = np.max(density_map_data)
        dmdmin = np.min(density_map_data)
        normal_density_map_data = [(element - dmdmin) / (dmdmax - dmdmin) for element in density_map_data]
        normal_density_map_data = np.array(normal_density_map_data)
        temp.append(normal_density_map_data)

    temp = np.array(temp)
    density_map_data_sum = temp.sum(axis=0)
    save_normal_dmap = os.path.join(prepath,
     'normal_iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc5_map_sym.gii')
    save.save_gifti(density_map_data_sum, save_normal_dmap, hemisphere='CortexLeft')


def cal_zscore(input_name = None, output_name=None):
    prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/density_map/5/'
    density_map_list = os.listdir(prepath)  #acquire the density map of each subjects
    temp = []
    for density_map_name in density_map_list:
        density_map_path = os.path.join(prepath, density_map_name)
        density_map = nib.load(density_map_path).darrays
        density_map_data = density_map[0].data
        dmd_std = np.average(density_map_data)
        dmd_average = np.std(density_map_data)
        normal_density_map_data = [((element - dmd_average) / dmd_std) for element in density_map_data]
        normal_density_map_data = np.array(normal_density_map_data)
        temp.append(normal_density_map_data)

    temp = np.array(temp)
    density_map_data_zscore_sum = temp.sum(axis=0)
    save_normal_dmap = os.path.join(prepath,
     'average_iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc5_map_sym.gii')
    save.save_gifti(density_map_data_zscore_sum, save_normal_dmap, hemisphere='CortexLeft')
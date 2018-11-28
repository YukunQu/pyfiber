import os
import numpy as np
import nibabel as nib
from ep_density_plot import density_distribution_plot
from algorithm import save


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
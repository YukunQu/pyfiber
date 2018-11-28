import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np


def single_np(arr, target):
    arr = np.array(arr)
    mask = (arr == target)
    arr_new = arr[mask]
    return arr_new.size


def all_np(arr):
    # calculate the nums of different elements
    arr = np.array(arr)
    key = np.unique(arr)
    result = {}
    for k in key:
        mask = (arr == k)
        arr_new = arr[mask]
        v = arr_new.size
        result[k] = v
    return result


def all_voxel_np(arr):
    key = np.unique(arr)
    result = {}
    for k in key:
        if k == 0:
            continue
        else:
            mask = (arr == k)
            arr_new = arr[mask]
            v = arr_new.size
            result[k] = v
    return result

def density_distribution_plot(data):
    density = all_voxel_np(data)
    print(density)
    plt.figure(figsize=(8, 6), dpi=80)
    plt.subplot(1, 1, 1)
    plt.bar(density.keys(), density.values(), width=0.6, label="fiber density", color='#87CEFA')
    plt.xlabel("density")
    plt.ylabel("number")
    plt.show()

#fiber  density-num Histogram
fiber_ep = nib.load('/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/lh_ep_101410.mgz')
ep_data = fiber_ep.get_data()


#
#
# f_streamlines = nib.streamlines.TckFile.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/888678/Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_surfocc5.tck").streamlines
#
# s_0 = [streamline[0] for streamline in f_streamlines]
# s_1 = [streamline[-1] for streamline in f_streamlines]
# s_oc = []
# for n in range(len(s_0)):
#     if s_0[n][1] > s_1[n][1]:
#         s_oc.append(s_0[n])
#     else:s_oc.append(s_1[n])
# s_oc = np.array(s_oc)
# s_oc_y = np.rint(s_oc[:,1])
# s_oc_y_density = all_voxel_np(s_oc_y)
# average_value = np.average(s_oc_y)
# average_density = np.average(s_oc_y_density.values())
# print(average_density)
# s_oc_y_density[0] = average_density
# plt.figure(figsize=(8, 6), dpi=80)
# plt.subplot(1, 1, 1)
# plt.bar(s_oc_y_density.keys(), s_oc_y_density.values(), width=0.6, label="fiber density", color='#87CEFA')
# plt.xlabel("coordinate")
# plt.ylabel("fiber density")
# plt.show()
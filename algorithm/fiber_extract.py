import numpy as np
import nibabel as nib


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
#
# def extract_hemisphere_fiber(streamlines,hemispere):
#     if hemispere == 'rh':
#         streamlines_rh =  streamlines[np.argwhere(streamlines[:,0] > 0,:]
#         return streamlines_rh
#     elif hemispere == 'lh':
#         streamlines_lh = streamlines[streamlines[:,0] < 0,:]
#         return streamlines_lh

fiber = nib.streamlines.load("/nfs/h1/workingshop/quyukun/fiberdata/subjects/101410/Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k.tck")
Out_fiber = '/nfs/h1/workingshop/quyukun/fiberdata/subjects/101410/Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k_rh_ex.tck'
stream_lines = fiber.streamlines
header = fiber.header
data_per_streamline = fiber.tractogram.data_per_streamline
data_per_point = fiber.tractogram.data_per_point
affine_to_rasmm = fiber.tractogram.affine_to_rasmm

fiber_rh = single_point_mid_sag(fiber)

save_streamline = nib.streamlines.array_sequence.ArraySequence(fiber_rh)
tractogram = nib.streamlines.tractogram.Tractogram(streamlines=save_streamline,data_per_streamline=data_per_streamline,
                                                   data_per_point=data_per_point,affine_to_rasmm=affine_to_rasmm)
datdat = nib.streamlines.tck.TckFile(tractogram=tractogram,header=header)
datdat.save(Out_fiber)
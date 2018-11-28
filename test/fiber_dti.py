import numpy as np
import os
import nibabel as nib
from scipy.spatial.distance import cdist
from nibabel import  streamlines
from algorithm import save


def fiber_density_surfmap(fiber_path,white_path,output,hemi='lh',threshold = 1.1):

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
            else:s_oc.append(s_1[n])
        s_oc = np.array(s_oc)
    else:
        for n in range(len(s_0)):
            if s_0[n][0] > s_1[n][0]:
                s_oc.append(s_0[n])
            else:s_oc.append(s_1[n])
        s_oc = np.array(s_oc)

    # aim to get a Mx3 array. M is num of occipital coordinates,3 is coordinates fiber endpoints(x,y,z)
    gii_data = nib.load(white_path).darrays
    white_coords = gii_data[0].data

    #aim to get a MxN matrix.M is num of occipital coordinates, N is fiber_num. scalar value is distance.
    dist = cdist(white_coords, s_oc)

    #aim to get a Mx1 matrix.M is num of occipital coordinates, N is the amount of endpoints meeting the threshold.

    dist_oc = dist < threshold
    dist_oc_num = dist_oc.sum(axis=1)  #the scalar data

    #save to a gifiti file
    if hemi == 'lh':
        save.save_gifti(data=dist_oc_num, save_name=output, hemisphere='CortexLeft')
    else:
        save.save_gifti(data=dist_oc_num, save_name=output, hemisphere='CortexRight')


prepath = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects"

subjects = ['100206', '100307', '100408', '100610', '101006', '101107', '101309', '101410', '101915', '102008',
            '102109', '102311', '102513', '102614', '102715', '102816', '103010', '103111', '111211', '200008',
            '201818', '333330', '530635', '555651', '597869', '686969', '767464', '888678', '995174', '996782']
for subject in subjects:
    fiber_lh_path = os.path.join(prepath, subject,
        "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc1.1.tck")
    fiber_rh_path = os.path.join(prepath, subject,
        "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc1.1.tck")

    white_path = os.path.join(prepath, subject, "surf/lh.white.surf.gii")
    output_lh = os.path.join(prepath, subject, "Diffusion/tractography/endpoints/" 
              "iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc1.1_map.gii")
    output_rh = os.path.join(prepath, subject, "Diffusion/tractography/endpoints/"
              "iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc1.1_map.gii")

    fiber_density_surfmap(fiber_lh_path, white_path, output_lh, 'lh', 1.1)
    fiber_density_surfmap(fiber_rh_path, white_path, output_rh, 'rh', 1.1)


    fiber_lh_path = os.path.join(prepath, subject,
                                 "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc3.tck")
    fiber_rh_path = os.path.join(prepath, subject,
                                 "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc3.tck")

    white_path = os.path.join(prepath, subject, "surf/lh.white.surf.gii")
    output_lh = os.path.join(prepath, subject, "Diffusion/tractography/endpoints/"
                                               "iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc3_map.gii")
    output_rh = os.path.join(prepath, subject, "Diffusion/tractography/endpoints/"
                                               "iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc3_map.gii")

    fiber_density_surfmap(fiber_lh_path, white_path, output_lh, 'lh', 3)
    fiber_density_surfmap(fiber_rh_path, white_path, output_rh, 'rh', 3)


    fiber_lh_path = os.path.join(prepath, subject,
                                 "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc5.tck")
    fiber_rh_path = os.path.join(prepath, subject,
                                 "Diffusion/tractography/Prob/iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc5.tck")

    white_path = os.path.join(prepath, subject, "surf/lh.white.surf.gii")
    output_lh = os.path.join(prepath, subject, "Diffusion/tractography/endpoints/"
                                               "iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc5_map.gii")
    output_rh = os.path.join(prepath, subject, "Diffusion/tractography/endpoints/"
                                               "iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc5_map.gii")

    fiber_density_surfmap(fiber_lh_path, white_path, output_lh, 'lh', 5)
    fiber_density_surfmap(fiber_rh_path, white_path, output_rh, 'rh', 5)

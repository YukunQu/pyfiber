import os
import subprocess

def antsRegistration(fixed_image,moving_image,output_name):
    #  fixed image to which we register the moving image,absolute path,type:string;
    #  moving image to be mapped to fixed space,absolute path,type:string;
    subprocess.call('antsRegistrationSyNQuick.sh -d 3 -f {} -m {} -o {}'.format(fixed_image,moving_image,output_name),
                    shell = True)
    print('Affine and Warp files acquired')


def antsApplyTransforms_e0(fixed_image,moving_image,affine,warp,outputname):
    # affine : path of GenericAffine.mat
    # warp : path of Warp.nii.gz'
    subprocess.call('antsApplyTransforms -d 3 -e 0 -i {} -r {} -t {} -t {} -o {}'.format(moving_image,
                    fixed_image,affine,warp,outputname),shell = True)


def antsApplyTransforms_e3(fixed_image,moving_image,affine,warp,outputname):
    # t1 : path of GenericAffine.mat
    # t2 : path of Warp.nii.gz'
    subprocess.call('antsApplyTransforms -d 3 -e 3 -i {} -r {} -t {} -t {} -o {}'.format(moving_image,
                    fixed_image,affine,warp,outputname),shell = True)

def surf2surf(src_subject_name,src_val_path,target_subject_name,target_vol_path,hemi):
    #
    subprocess.call('mri_surf2surf --srcsubject {} --sval {} --trgsubject {} --trgsurfval {} --hemi {}'.format(
        src_subject_name,src_val_path,target_subject_name, target_vol_path,hemi),shell= True)



prepath = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects"

subjects = ['100206', '100307', '100408', '100610', '101006', '101107', '101309', '101410', '101915', '102008',
            '102109', '102311', '102513', '102614', '102715', '102816', '103010', '103111', '111211', '200008',
            '201818', '333330', '530635', '555651', '597869', '686969', '767464', '888678', '995174', '996782']

for subject in subjects:
    src_subject_name = 'fsaverage_sym'
    target_subject_name = subject
    src_val_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/benson_retinotopic_map/areas-template-2.5.sym.mgh'

    target_vol_path = os.path.join(prepath, subject,'ROI/benson_retinotopic_map/areas-template-2.5.sym.mgh')
    surf2surf(src_subject_name, src_val_path, target_subject_name, target_vol_path, hemi='lh')


# for subject in subjects:
#     prepath= '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects'
#     src_subject_name = subject
#     target_subject_name = 'fsaverage_sym'
#
#
#     src_val_path = os.path.join(prepath,subject,'Diffusion/tractography/endpoints/'
#                         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc1.1_map.gii')
#     target_vol_path = os.path.join(prepath, 'fsaverage_sym/density_map',
#          'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc1.1_map_{}sym.gii'.format(subject))
#     surf2surf(src_subject_name, src_val_path, target_subject_name, target_vol_path, hemi='lh')
#
#
#     src_val_path = os.path.join(prepath,subject,'Diffusion/tractography/endpoints/'
#                         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc3_map.gii')
#     target_vol_path = os.path.join(prepath, 'fsaverage_sym/density_map',
#         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc3_map_{}sym.gii'.format(subject))
#     surf2surf(src_subject_name, src_val_path, target_subject_name, target_vol_path, hemi='lh')
#
#
#     src_val_path = os.path.join(prepath,subject,'Diffusion/tractography/endpoints/'
#                         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc5_map.gii')
#     target_vol_path = os.path.join(prepath, 'fsaverage_sym/density_map',
#         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_L_100k_spsurfoc5_map_{}sym.gii'.format(subject))
#     surf2surf(src_subject_name, src_val_path, target_subject_name, target_vol_path, hemi='lh')
#
#
#     src_val_path = os.path.join(prepath, subject, 'Diffusion/tractography/endpoints/'
#                         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc1.1_map.gii')
#     target_vol_path = os.path.join(prepath, 'fsaverage_sym/density_map',
#         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc1.1_map_{}sym.gii'.format(subject))
#     surf2surf(src_subject_name, src_val_path, target_subject_name, target_vol_path, hemi='rh')
#
#
#     src_val_path = os.path.join(prepath, subject, 'Diffusion/tractography/endpoints/'
#                         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc3_map.gii')
#     target_vol_path = os.path.join(prepath, 'fsaverage_sym/density_map',
#         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc3_map_{}sym.gii'.format(subject))
#     surf2surf(src_subject_name, src_val_path, target_subject_name, target_vol_path, hemi='rh')
#
#     src_val_path = os.path.join(prepath, subject, 'Diffusion/tractography/endpoints/'
#                         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc5_map.gii')
#     target_vol_path = os.path.join(prepath, 'fsaverage_sym/density_map',
#         'iFOD2_angle20_cutoff0.05_length50_250_seedFP_roiOccipital_R_100k_spsurfoc5_map_{}sym.gii'.format(subject))
#     surf2surf(src_subject_name, src_val_path, target_subject_name, target_vol_path, hemi='rh')
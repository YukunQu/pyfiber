# antsRegistration
#Apply transform for ROI from standard space  to native space.
#Apply transform for tract_prob from standard space to native space.
#extract special fiber probability atlas from tract_prob in native space
#!!attention: The structure of directroy must be coincident.
import os
import nibabel as nib
import registration
import volume_extract as ex

path_MNI = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/standard/MNI152_T1_2mm_brain.nii.gz'
path_JHU_tract_prob = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/standard/JHU-ICBM-tracts-prob-2mm.nii.gz'
path_ROI = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/tract_ROI'

prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects'

pospath_T1w = 'T1w/T1w_acpc_dc_restore_brain1.25.nii.gz'
pospath_transform = 'regis_ROI/transform'
pospath_regis_ROI = 'regis_ROI'
pospath_regis_atlas = 'prob_atlas'

subjects = [101107]

for sj_id in subjects:
    path_T1w = os.path.join(prepath, str(sj_id), pospath_T1w)
    path_output = os.path.join(prepath, str(sj_id), pospath_transform, str(sj_id))

    # antsRegistration
    registration.antsRegistration(path_T1w,path_MNI,path_output)

    # Apply transform for ROI from standard space  to native space.
    t1 = str(sj_id) + '0GenericAffine.mat'
    t2 = str(sj_id) + '1Warp.nii.gz'
    path_affine = os.path.join(prepath, str(sj_id), pospath_transform, t1)
    path_warp = os.path.join(prepath, str(sj_id), pospath_transform, t2)

    rois = ['L_Occipital.nii.gz','R_Occipital.nii.gz']
    for roi in rois:
        roi_path = os.path.join(path_ROI,roi)
        roi_output_name = os.path.join(prepath, str(sj_id), pospath_regis_ROI,roi)
        registration.antsApplyTransforms_e0(path_T1w,roi_path,path_affine,path_warp,roi_output_name)

    #Apply transform for tract_prob from standard space to native space.
    t3 = str(sj_id) + '_' + 'JHU_tract_prob.nii.gz'
    tract_outputname = os.path.join(prepath,str(sj_id),pospath_regis_atlas,t3)
    registration.antsApplyTransforms_e3(path_T1w,path_JHU_tract_prob,path_affine,path_warp,tract_outputname)

    # extract special fiber probability atlas from tract_prob in native space
    pro_atlas_path = tract_outputname
    t4 = 'FP_prob.nii.gz'
    FP_output_name = os.path.join(prepath, str(sj_id), pospath_regis_atlas, t4)

    proatlas = nib.load(pro_atlas_path)
    affine = proatlas.affine

    proatlas_data = proatlas.get_data()
    extract_atlas = ex.extract_label(proatlas_data, 8)

    # save as Nifti
    dm_img = nib.Nifti1Image(extract_atlas, affine)
    dm_img.to_filename(FP_output_name)
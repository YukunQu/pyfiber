import nibabel as nib
import os
import algorithm.volume_extract as volex


def map_f2t(pro_file,out_file):
    # 4D to 3D
    profile = nib.load(pro_file)
    affine = profile.affine
    profile_data = profile.get_data()

    extract_map = volex.extract_label(profile_data,0)
    dm_img = nib.Nifti1Image(extract_map.astype("float32"), affine)
    dm_img.to_filename(out_file)


def map_sp_val(pro_file,out_file,sp):
    #extract special map data below a certain threshold.
    profile = nib.load(pro_file)
    affine = profile.affine
    profile_data = profile.get_data()

    extract_map = volex.extract_sp_val(profile_data,sp)
    dm_img = nib.Nifti1Image(extract_map.astype("float32"), affine)
    dm_img.to_filename(out_file)


def map_endpoints(pro_file,out_file):
    #extract fiber endpoints map near occipital,reject the fiber endpoints near corpus callosum
    profile = nib.load(pro_file)
    affine = profile.affine
    file_data = profile.get_data()

    file_data[72].flat = 0

    dm_img = nib.Nifti1Image(file_data.astype("float32"), affine)
    dm_img.to_filename(out_file)

def map_gt_th(pro_file, out_file):
    # extract map that values greater threshold in a batch run.
    sessid = os.listdir('/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects')
    for subjects_id in sessid:
        proatlas = nib.load(pro_file.format(subjects_id=subjects_id))
        affine = proatlas.affine
        proatlas_data = proatlas.get_data()

        extract_map = volex.extract_threshold(proatlas_data, 0.45)
        dm_img = nib.Nifti1Image(extract_map.astype("float32"), affine)
        dm_img.to_filename(out_file.format(subjects_id=subjects_id))


pro_file = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/' \
           'iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k_rh_occipital_map.nii.gz'
out_file = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/' \
           'iFOD2_angle20_cutoff0.05_length50_250_seedFP__roiOccipital_L_100k_rh_occipital_map_rrs.nii.gz'

map_endpoints(pro_file,out_file)
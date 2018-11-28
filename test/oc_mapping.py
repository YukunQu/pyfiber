import nibabel as nib
import subprocess
import numpy as np
import os

label = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/areas-template-2.5.sym.mgh")
label_data = label.get_data().reshape(163842, 1, 1)
affine = label.affine

# classify to seven files
outputpath = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/benson_retinotopic_map"
output = { 1: "v1v_sym.nii.gz", -1: "v1d_sym.nii.gz",  2: "v2v_sym.nii.gz",
          -2: "v2d_sym.nii.gz",  3: "v3v_sym.nii.gz", -3: "v3d_sym.nii.gz"}

prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects'
subjects = os.listdir(prepath)
subjects.remove('fsaverage_sym')
subjects = ['101410']

for value,outputname in output.items():
    ep_vison_data = label_data.copy()
    ep_vison_data[ep_vison_data != value] = 0
    dm_img = nib.Nifti1Image(ep_vison_data.astype("float32"),affine)
    dm_img.to_filename(os.path.join(outputpath,outputname))

    for subject in subjects:

        sval = os.path.join(prepath,'fsaverage_sym/benson_retinotopic_map',outputname)
        trgsubject = subject
        subject_out_lh = 'lh_' + outputname.replace("sym",subject)
        subject_out_rh = 'rh_' + outputname.replace("sym", subject)
        tval_lh = os.path.join(prepath, subject, 'ROI/benson_retinotopic_map', subject_out_lh)
        tval_rh = os.path.join(prepath, subject, 'ROI/benson_retinotopic_map', subject_out_rh)

        subprocess.call('mri_surf2surf --srcsubject fsaverage_sym --sval {} --trgsubject {} --tval {} '
                        '--hemi lh'.format(sval,trgsubject,tval_lh),shell=True)
        subprocess.call('mri_surf2surf --srcsubject fsaverage_sym --sval {} --trgsubject {} --tval {} '
                        '--hemi rh'.format(sval,trgsubject,tval_rh),shell=True)

vison_areas = []
for vison_area in vison_areas:
    index = np.where(vison_area)
    selcet_tck = selcet_by_roi(ori_tck,white_geo,index)

    def select_by_sufroi(ori_tck,white_geo,index):
        selcet_tck = select_by_roi(ori.tck,white_geo[index],affine)
        return selcet_tck

    save_to_tck = (vison_area_filename)
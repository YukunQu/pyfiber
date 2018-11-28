import nibabel as nib
import numpy as np
import os

endpoints = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/lh_ep_101410.mgz")
endpoints_data = endpoints.get_data()
affine = endpoints.affine
header = endpoints.header

label = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/areas-template-2.5.sym.mgh")
label_data = label.get_data()

ep_index = np.argwhere(endpoints_data)[:, 0]

# classify to seven files
# outputpath = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym"
# output = {1:"ep_v1v_101410.nii.gz",-1:"ep_v1d_101410.nii.gz",2:"ep_v2v_101410.nii.gz",-2:"ep_v2d_101410.nii.gz",
#           3:"ep_v3v_101410.nii.gz",-3:"ep_v3d_101410.nii.gz",0:"ep_other_101410.nii.gz"}
#
# for value, outputname in output.items():
#     ep_vison_index = [ep for ep in ep_index if label_data[0][0][ep] == value]
#     ep_vison_data = np.zeros(endpoints_data.shape)
#     for i in ep_vison_index:
#         ep_vison_data[i][0][0] = endpoints_data[i][0][0]
#     dm_img = nib.Nifti1Image(ep_vison_data.astype("float32"), affine)
#     dm_img.to_filename(os.path.join(outputpath,outputname))

#classify to seven label in one file
outputname = "/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/fsaverage_sym/lh_ep_101410_classify_label.nii"
output = [1,-1,2,-2,3,-3,0]
ep_vison_data = np.zeros(endpoints_data.shape)
for value in output:
    ep_vison_index = [ep for ep in ep_index if label_data[0][0][ep] == value]
    for i in ep_vison_index:
        ep_vison_data[i][0] = value
        # ep_vison_data[i][0][0] = endpoints_data[i][0][0]

dm_img = nib.Nifti1Image(ep_vison_data.astype("float32"), affine)
dm_img.to_filename(outputname)

# endpoints = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/lh_ep_101410.mgz")
# endpoints_data = endpoints.get_data()
#
# label = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/fsaverage_sym/areas-template-2.5.sym.mgh")
# label_data = label.get_data()
#
# label_value = np.unique(label_data)
# label_index = np.array([(np.argwhere(label_data == v))[:,2] for v in label_value])
# val_index = dict(zip(label_value, label_index))
#
# for label_v, index in val_index.items():
#     mask = (label_data == label_v)
#
# yong zuihou yizu zuobiao zuo wei suoyin ,xunhuan shengcheng xin de wenjian he
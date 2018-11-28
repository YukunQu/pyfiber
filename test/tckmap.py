import subprocess

retinotopic_areas_list = ['lh_v1v_101410.nii.gz',  'lh_v1d_101410.nii.gz', 'lh_v2v_101410.nii.gz',
                          'lh_v2d_101410.nii.gz',  'lh_v3v_101410.nii.gz', 'lh_v3d_101410.nii.gz']
fiber_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/Prob/'
fiber_list = [fiber_path+area.replace('.nii.gz', '.tck') for area in retinotopic_areas_list]

T1w = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/T1w_acpc_dc_restore_brain1.25.nii.gz'
out_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Diffusion/tractography/endpoints/map_'
out_list = [out_path+area for area in retinotopic_areas_list]

for i in range(len(fiber_list)):
    subprocess.call('tckmap {} -template {} -ends_only {}'.format(fiber_list[i],T1w,out_list[i]))



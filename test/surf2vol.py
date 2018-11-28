import os
import numpy as np
import nibabel as nib
import numpy.linalg as npl
from nibabel.affines import apply_affine
from nibabel.spatialimages import ImageFileError
from dipy.viz import window, actor


def roi_surf2vol(roi_surf_path, geo_path, vol_path):
    """
    Transform rois from surface to volume.
    Parameters
    ----------
    roi_surf_path: surface roi path
    geo_path: surface geometry path
    vol_path: target volume
    Return
    ------
    volume rois
    """
    # load parcelltion data
    vertices, colortable, label = nib.freesurfer.read_annot(roi_surf_path)

    # load surf geometry
    suffix = os.path.split(geo_path)[1].split('.')[-1]
    if suffix in ('white', 'inflated', 'pial'):
        coords, faces = nib.freesurfer.read_geometry(geo_path[0])
    elif suffix == 'gii':
        gii_data = nib.load(geo_path).darrays
        coords, faces = gii_data[0].data, gii_data[1].data
    else:
        raise ImageFileError('This file format-{} is not supported at present.'.format(suffix))

    # load volume data
    img = nib.load(vol_path)
    data_temp_1 = np.zeros(img.shape)
    data_temp_2 = np.zeros(img.shape)
    data_temp_3 = np.zeros(img.shape)
    data_temp_4 = np.zeros(img.shape)

    coords_affine = apply_affine(npl.inv(img.affine), coords)

    # parcellation of surf to volume
    for c in coords_affine[vertices == 11]:
        data_temp_1[int(c[0]), int(c[1]), int(c[2])] = 11
    for c in coords_affine[vertices == 13]:
        data_temp_2[int(c[0]), int(c[1]), int(c[2])] = 12
    for c in coords_affine[vertices == 5]:
        data_temp_3[int(c[0]), int(c[1]), int(c[2])] = 5
    for c in coords_affine[vertices == 21]:
        data_temp_4[int(c[0]), int(c[1]), int(c[2])] =21

    return data_temp_1,data_temp_2,data_temp_3,data_temp_4

roi_surf_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Native/lh.aparc.DKTatlas40.annot'
geo_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/Native/101410.L.white.native.surf.gii'
vol_path = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/T1w_acpc_dc_restore_brain1.25.nii.gz'
later_occipital_vol,lingual_vol,cuneus_vol, pericalarine_vol = roi_surf2vol(roi_surf_path, geo_path, vol_path)

T1w = nib.load("/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/101410/T1w_acpc_dc_restore_brain1.25.nii.gz")
T1w_data = T1w.get_data()

interactive = True
vol_actor = actor.slicer(T1w_data)
vol_actor.display(x=72)
vol_actor2 = vol_actor.copy()
vol_actor2.display(z=40)

# Make display objects
cc_ROI_actor1 = actor.contour_from_roi(later_occipital_vol,color=(1, 0, 0),opacity=0.8)
cc_ROI_actor2 = actor.contour_from_roi(lingual_vol,color=(0, 1, 0),opacity=0.8)
cc_ROI_actor3 = actor.contour_from_roi(cuneus_vol,color=(0, 0, 1),opacity=0.8)
cc_ROI_actor4 = actor.contour_from_roi(pericalarine_vol,color=(1, 1, 1),opacity=0.8)


# Add display objects to canvas
r = window.Renderer()
r.add(vol_actor)
r.add(vol_actor2)
r.add(cc_ROI_actor1)
r.add(cc_ROI_actor2)
r.add(cc_ROI_actor3)
r.add(cc_ROI_actor4)

if interactive:
    window.show(r)

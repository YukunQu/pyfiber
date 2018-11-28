# AFQ pipeline:
# 1.ROIs to be registrate into native space
# 2.select fiber through ROIs
# 3.eliminate unnecessary fiber by probilistic atlas
# 4.delete unreal fiber by bundle coherence method or multivariate gaussian
# 5.fibers are resampled to 100 points
# 6.calculate FA value of every point

#Registration

#select
#eliminate unnecessary fiber

import nibabel as nib
from nibabel import streamlines
from dipy.tracking import streamline
from dipy.viz import window, actor
from dipy.viz.colormap import line_colors

import algorithm.volume_extract

fiber_IFO_L_100206 = streamlines.TckFile.load('/nfs/h1/workingshop/quyukun/fiberdata/test/100206_IFO_L.tck')
prob_JHU_IFO_100206 = nib.load('/nfs/h1/workingshop/quyukun/fiberdata/test/JHU100206.nii.gz')

fiber_data = fiber_IFO_L_100206.streamlines
atlas_affine = prob_JHU_IFO_100206.affine
atlas_data = prob_JHU_IFO_100206.get_data()
IFO_atlas_data = algorithm.volume_extract.extract_volume(atlas_data, 10, threshold=25).astype(bool)

selected_streamlines = streamline.select_by_rois(fiber_data,[IFO_atlas_data],[True],mode='any',affine=atlas_affine)

selected_fiber = list(selected_streamlines)


#visulaziton------------------------------------------------------------------------------------------------------
img_T1w = nib.load('/nfs/s2/userhome/quyukun/workingdir/fiberdata/100206/'
                   'Structural image/T1w_acpc_dc_restore_brain1.25.nii.gz')
img_roi1 = nib.load('/nfs/s2/userhome/quyukun/workingdir/fiberdata/100206/regis_ROI/100206_IFO_roi1_L.nii.gz')
img_roi2 = nib.load('/nfs/s2/userhome/quyukun/workingdir/fiberdata/100206/regis_ROI/100206_IFO_roi2_L.nii.gz')

img_fiber = selected_fiber
#
roi_affine = img_roi1.affine
img_roi1_data = img_roi1.get_data()
img_roi2_data = img_roi2.get_data()
img_atlas_data = IFO_atlas_data
# Enables/disables interactive visualization
interactive = True

# Make display objects
color = line_colors(img_fiber)
cc_streamlines_actor = actor.line(img_fiber, color)
cc_ROI_actor1 = actor.contour_from_roi(img_roi1_data,affine=roi_affine,color=(1.0,0.5,0.2),
                                      opacity=0.5)
cc_ROI_actor2 = actor.contour_from_roi(img_roi2_data,roi_affine,color=(1.0, 0.5, 0.2),
                                      opacity=0.5)
cc_atlas_actor3 = actor.contour_from_roi(IFO_atlas_data,roi_affine,color=(1.0,0.5,0.2),opacity=0.5)

img_T1w_data = img_T1w.get_data()
vol_actor = actor.slicer(img_T1w_data,roi_affine)

vol_actor.display(x=72)
vol_actor2 = vol_actor.copy()
vol_actor2.display(z=50)

# Add display objects to canvas
r = window.Renderer()
r.add(vol_actor)
r.add(vol_actor2)
r.add(cc_streamlines_actor)
r.add(cc_ROI_actor1)
r.add(cc_ROI_actor2)
r.add(cc_atlas_actor3)


# Save figures

if interactive:
    window.show(r)
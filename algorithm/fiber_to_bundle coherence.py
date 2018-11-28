from dipy.tracking.fbcmeasures import FBCMeasures
from dipy.denoise.enhancement_kernel import EnhancementKernel
import nibabel.streamlines.tck as nibtck
import nibabel as nib

T1w = nib.load('/nfs/s2/userhome/quyukun/workingdir/fiberdata/100206'
                   '/Structural image/T1w_acpc_dc_restore_brain1.25.nii.gz')
fiber = nibtck.TckFile.load('/nfs/s2/userhome/quyukun/workingdir/fiberdata/100206/selectedfiber/'
                                '100206_IFO_L_prob25.tck')
fiber_data = fiber.streamlines
affine = T1w.affine
img_T1w_data = T1w.get_data()

D33 = 1.0
D44 = 0.02
t = 1
k = EnhancementKernel(D33, D44, t)

# Apply FBC measures

fbc = FBCMeasures(fiber_data, k)

#Calculate LFBC for original fibers
fbc_sl_orig, clrs_orig, rfbc_orig = fbc.get_points_rfbc_thresholded(0, emphasis=0.01)

# Apply a threshold on the RFBC to remove spurious fibers
fbc_sl_thres, clrs_thres, rfbc_thres = fbc.get_points_rfbc_thresholded(0.125, emphasis=0.01)

print("The process is already running here.")

# Visualize the results
from dipy.viz import window, actor

# Enables/disables interactive visualization
interactive = True

# Create renderer
ren = window.Renderer()
lineactor = actor.line(fiber_data,linewidth=0.2)
ren.add(lineactor)

# Horizontal (axial) slice of T1 data
vol_actor1 = actor.slicer(img_T1w_data,affine=affine)
vol_actor1.display(z=50)
ren.add(vol_actor1)

# Vertical (sagittal) slice of T1 data
vol_actor2 = actor.slicer(img_T1w_data,affine=affine)
vol_actor2.display(x=70)
ren.add(vol_actor2)

# Show original fibers
if interactive:
    window.show(ren)

ren.rm(lineactor)
ren.add(actor.line(fbc_sl_thres, clrs_thres, linewidth=0.2))
if interactive:
    window.show(ren)
from surfer import Brain

"""
Bring up the visualization window.
"""
brain = Brain("fsaverage_sym", "lh", "inflated")

"""
Get a path to the overlay file.
"""
overlay_file = "lh.sig.nii.gz"

"""
Display the overlay on the surface using the defaults to control thresholding
and colorbar saturation.  These can be set through your config file.
"""
brain.add_overlay(overlay_file)

"""
You can then turn the overlay off.
"""

"""
Now add the overlay again, but this time with set threshold and showing only
the positive activations.
"""


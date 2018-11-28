import numpy as np
from nibabel import streamlines
from dipy.segment.clustering import QuickBundles
from dipy.viz import window, actor

img_fiber = streamlines.TckFile.load('/nfs/s2/userhome/quyukun/workingdir/fiberdata/100206/selectedfiber/'
                                '100206_IFO_L_prob25.tck')
streamline = img_fiber.streamlines

qb_streamline = QuickBundles(threshold=10.0)
clusters = qb_streamline.cluster(streamline)
print(clusters.centroids)


# Enables/disables interactive visualization
interactive = True

ren = window.Renderer()
colormap = actor.create_colormap(np.arange(len(clusters)))
#colormap = [0.050383,0.029803,0.527975]

#window.clear(ren)
ren.SetBackground(1, 1, 1)
#ren.add(actor.streamtube(streamline, window.colors.red, opacity=0.2))
ren.add(actor.streamtube(clusters.centroids, colormap, linewidth=0.5))
#window.record(ren, out_path='fornix_centroids.png', size=(1000, 1000))
if interactive:
    window.show(ren)
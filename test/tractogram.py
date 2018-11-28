from nibabel.streamlines import tck
class Tractogram(object):
    def __init__(self, lines=None, data=None, space=None):
        """
        Parameters
        ----------
        lines: a Lines object
        data: the scalar image data
        space: str, native, mni152
        """
        self.lines = lines
        self.data = data
        self.space = space

    def load_lines(self, filename):
        """ Load tractography from a tractography file, include tck, trk, vtk(tck file will be accepted  temporarily  )

        Parameters
        ----------
        filename: str
            Pathstr to a tractography file

        Returns
        -------
        self: a Lines object
        """
        self.lines = tck.TckFile.load(filename)
        self.data = self.lines.streamlines

filename = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects/100408/Diffusion/tractography/Det/SD_Stream_angle20_cutoff0.1_length20_80_seedAC0.5_100k.tck'
tractogram = Tractogram()
tractogram.load_lines(filename)
print(type(tractogram.data))
print(len(tractogram.data))
print(tractogram.data)



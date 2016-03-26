from astropy.io.fits import HDUList, PrimaryHDU
from .ytcube import ytCube
from distutils.version import StrictVersion

def to_yt(self, spectral_factor=1.0, nprocs=None, **kwargs):
    """
    Convert a spectral cube to a yt object that can be further analyzed in
    yt.

    Parameters
    ----------
    spectral_factor : float, optional
        Factor by which to stretch the spectral axis. If set to 1, one pixel
        in spectral coordinates is equivalent to one pixel in spatial
        coordinates.

    If using yt 3.0 or later, additional keyword arguments will be passed
    onto yt's ``FITSDataset`` constructor. See the yt documentation
    (http://yt-project.org/docs/3.0/examining/loading_data.html?#fits-data)
    for details on options for reading FITS data.
    """

    import yt

    if (('dev' in yt.__version__ or
         (StrictVersion(yt.__version__) >= StrictVersion('3.0')))):

        from yt.frontends.fits.api import FITSDataset
        from yt.units.unit_object import UnitParseError

        hdu = PrimaryHDU(self._get_filled_data(fill=0.),
                         header=self.wcs.to_header())

        units = str(self.unit.to_string())

        hdu.header["BUNIT"] = units
        hdu.header["BTYPE"] = "flux"

        ds = FITSDataset(hdu, nprocs=nprocs,
                         spectral_factor=spectral_factor, **kwargs)

        # Check to make sure the units are legit

        try:
            ds.quan(1.0,units)
        except UnitParseError:
            raise RuntimeError("The unit %s was not parsed by yt. " % units+
                               "Check to make sure it is correct.")

    else:

        from yt.mods import load_uniform_grid

        data = {'flux': self._get_filled_data(fill=0.).transpose()}

        nz, ny, nx = self.shape

        if nprocs is None:
            nprocs = 1

        bbox = np.array([[0.5,float(nx)+0.5],
                         [0.5,float(ny)+0.5],
                         [0.5,spectral_factor*float(nz)+0.5]])

        ds = load_uniform_grid(data, [nx,ny,nz], 1., bbox=bbox,
                               nprocs=nprocs, periodicity=(False, False,
                                                           False))

    return ytCube(self, ds, spectral_factor=spectral_factor)

def to_glue(self, name=None, glue_app=None, dataset=None, start_gui=True):
    """
    Send data to a new or existing Glue application

    Parameters
    ----------
    name : str or None
        The name of the dataset within Glue.  If None, defaults to
        'SpectralCube'.  If a dataset with the given name already exists,
        a new dataset with "_" appended will be added instead.
    glue_app : GlueApplication or None
        A glue application to send the data to.  If this is not specified,
        a new glue application will be started if one does not already
        exist for this cube.  Otherwise, the data will be sent to the
        existing glue application, `self._glue_app`.
    dataset : glue.core.Data or None
        An existing Data object to add the cube to.  This is a good way
        to compare cubes with the same dimensions.  Supercedes ``glue_app``
    start_gui : bool
        Start the GUI when this is run.  Set to False for testing.
    """
    if name is None:
        name = 'SpectralCube'

    from glue.qt.glue_application import GlueApplication
    from glue.core import DataCollection, Data
    from glue.core.coordinates import coordinates_from_header
    from glue.qt.widgets import ImageWidget

    if dataset is not None:
        if name in [d.label for d in dataset.components]:
            name = name+"_"
        dataset[name] = self

    else:
        result = Data(label=name)
        result.coords = coordinates_from_header(self.header)

        result.add_component(self, name)

        if glue_app is None:
            if hasattr(self,'_glue_app'):
                glue_app = self._glue_app
            else:
                # Start a new glue session.  This will quit when done.
                # I don't think the return statement is ever reached, based on
                # past attempts [@ChrisBeaumont - chime in here if you'd like]
                dc = DataCollection([result])

                #start Glue
                ga = self._glue_app = GlueApplication(dc)
                self._glue_viewer = ga.new_data_viewer(ImageWidget,
                                                       data=result)

                if start_gui:
                    self._glue_app.start()

                return self._glue_app

        glue_app.add_datasets(self._glue_app.data_collection, result)


def to_pvextractor(self):
    """
    Open the cube in a quick viewer written in matplotlib that allows you
    to create PV extractions within the GUI
    """
    from pvextractor.gui import PVSlicer

    return PVSlicer(self)

def to_ds9(self, ds9id=None, newframe=False):
    """
    Send the data to ds9 (this will create a copy in memory)

    Parameters
    ----------
    ds9id: None or string
        The DS9 session ID.  If 'None', a new one will be created.
        To find your ds9 session ID, open the ds9 menu option
        File:XPA:Information and look for the XPA_METHOD string, e.g.
        ``XPA_METHOD:  86ab2314:60063``.  You would then calll this
        function as ``cube.to_ds9('86ab2314:60063')``
    newframe: bool
        Send the cube to a new frame or to the current frame?
    """
    try:
        import ds9
    except ImportError:
        import pyds9 as ds9

    if ds9id is None:
        dd = ds9.ds9(start=True)
    else:
        dd = ds9.ds9(target=ds9id, start=False)

    if newframe:
        dd.set('frame new')

    dd.set_pyfits(HDUList(self.hdu))

    return dd

def to_ginga(self, base_url="http://localhost:9909/app"):
    """
    Open the cube in a quick viewer made from Ginga pieces
    """
    from ginga.web.pgw import Widgets
    from ginga.misc import log
    from .ginga_viewer import FITSViewer

    logger = log.get_logger("spectral_cube")

    # establish our widget application
    app = Widgets.Application(logger=logger, host='localhost', port=9909)

    #  create top level window
    window = app.make_window("Spectral Cube Ginga Webapp Viewer")

    # our own viewer object, customized with methods (see above)
    viewer = FITSViewer(logger, window)

    return viewer

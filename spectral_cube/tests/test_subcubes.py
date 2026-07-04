import pytest

from astropy import units as u
from astropy import wcs
from astropy.io import fits
import numpy as np

from . import path
from .helpers import assert_allclose, assert_array_equal
from .test_spectral_cube import cube_and_raw
from ..spectral_axis import doppler_gamma, doppler_beta, doppler_z, get_rest_value_from_wcs
from ..spectral_cube import SpectralCube
from ..utils import WCSWarning

try:
    import regions
    regionsOK = True
except ImportError:
    regionsOK = False

try:
    import scipy
    scipyOK = True
except ImportError:
    scipyOK = False


def test_subcube(data_advs, use_dask):

    cube, data = cube_and_raw(data_advs, use_dask=use_dask)

    sc1x = cube.subcube(xlo=1, xhi=3)
    sc2x = cube.subcube(xlo=24.06269*u.deg, xhi=24.06206*u.deg)
    sc2b = cube.subcube(xlo=24.06206*u.deg, xhi=24.06269*u.deg)
    # Mixed should be equivalent to above
    sc3x = cube.subcube(xlo=24.06269*u.deg, xhi=3)
    sc4x = cube.subcube(xlo=1, xhi=24.06206*u.deg)

    assert sc1x.shape == (2,3,2)
    assert sc2x.shape == (2,3,2)
    assert sc2b.shape == (2,3,2)
    assert sc3x.shape == (2,3,2)
    assert sc4x.shape == (2,3,2)
    assert sc1x.wcs.wcs.compare(sc2x.wcs.wcs)
    assert sc1x.wcs.wcs.compare(sc2b.wcs.wcs)
    assert sc1x.wcs.wcs.compare(sc3x.wcs.wcs)
    assert sc1x.wcs.wcs.compare(sc4x.wcs.wcs)

    sc1y = cube.subcube(ylo=1, yhi=3)
    sc2y = cube.subcube(ylo=29.93464 * u.deg,
                        yhi=29.93522 * u.deg)
    sc3y = cube.subcube(ylo=1, yhi=29.93522 * u.deg)
    sc4y = cube.subcube(ylo=29.93464 * u.deg, yhi=3)

    assert sc1y.shape == (2, 2, 4)
    assert sc2y.shape == (2, 2, 4)
    assert sc3y.shape == (2, 2, 4)
    assert sc4y.shape == (2, 2, 4)
    assert sc1y.wcs.wcs.compare(sc2y.wcs.wcs)
    assert sc1y.wcs.wcs.compare(sc3y.wcs.wcs)
    assert sc1y.wcs.wcs.compare(sc4y.wcs.wcs)


    # Test mixed slicing in both spatial directions
    sc1xy = cube.subcube(xlo=1, xhi=3, ylo=1, yhi=3)
    sc2xy = cube.subcube(xlo=24.06269*u.deg, xhi=3,
                         ylo=1,yhi=29.93522 * u.deg)

    sc3xy = cube.subcube(xlo=1, xhi=24.06206*u.deg,
                         ylo=29.93464 * u.deg, yhi=3)


    assert sc1xy.shape == (2, 2, 2)
    assert sc2xy.shape == (2, 2, 2)
    assert sc3xy.shape == (2, 2, 2)
    assert sc1xy.wcs.wcs.compare(sc2xy.wcs.wcs)
    assert sc1xy.wcs.wcs.compare(sc3xy.wcs.wcs)


    sc1z = cube.subcube(zlo=1, zhi=2)
    sc2z = cube.subcube(zlo=-320*u.km/u.s, zhi=-319*u.km/u.s)
    sc3z = cube.subcube(zlo=1, zhi=-319 * u.km / u.s)
    sc4z = cube.subcube(zlo=-320*u.km/u.s, zhi=2)
    assert sc1z.shape == (1, 3, 4)
    assert sc2z.shape == (1, 3, 4)
    assert sc3z.shape == (1, 3, 4)
    assert sc4z.shape == (1, 3, 4)
    assert sc1z.wcs.wcs.compare(sc2z.wcs.wcs)
    assert sc1z.wcs.wcs.compare(sc3z.wcs.wcs)
    assert sc1z.wcs.wcs.compare(sc4z.wcs.wcs)

    sc5 = cube.subcube()

    assert sc5.shape == cube.shape
    assert sc5.wcs.wcs.compare(cube.wcs.wcs)
    assert np.all(sc5._data == cube._data)


def make_allsky_car_hdu(nx=49, ny=25, nz=5):
    """
    An all-sky plate carree (CAR) cube like the HI4PI one in issue #934:
    the pixel grid extends slightly past +-90 deg latitude and +-180 deg
    longitude, so the edge pixels fall outside the valid region of the
    projection and have NaN world coordinates.
    """
    header = fits.Header()
    header['NAXIS'] = 3
    header['NAXIS1'] = nx
    header['NAXIS2'] = ny
    header['NAXIS3'] = nz
    header['CTYPE1'] = 'GLON-CAR'
    header['CRVAL1'] = 0.0
    header['CDELT1'] = -360.2 / (nx - 1)  # spans a bit more than 360 deg
    header['CRPIX1'] = (nx + 1) / 2.
    header['CUNIT1'] = 'deg'
    header['CTYPE2'] = 'GLAT-CAR'
    header['CRVAL2'] = 0.0
    header['CDELT2'] = 180.2 / (ny - 1)  # spans a bit more than 180 deg
    header['CRPIX2'] = (ny + 1) / 2.
    header['CUNIT2'] = 'deg'
    header['CTYPE3'] = 'VRAD'
    header['CRVAL3'] = 0.0
    header['CDELT3'] = 1000.0
    header['CRPIX3'] = (nz + 1) / 2.
    header['CUNIT3'] = 'm/s'
    header['BUNIT'] = 'K'

    data = np.arange(nz * ny * nx, dtype=float).reshape((nz, ny, nx))

    return fits.PrimaryHDU(data=data, header=header)


def test_subcube_allsky(use_dask):
    # regression test for issue #934: edge pixels of all-sky cubes can fall
    # outside of the valid region of the projection, so they have NaN world
    # coordinates.  These should be ignored when computing the world extrema
    # and when converting world limits to pixel limits in subcube.

    hdu = make_allsky_car_hdu()
    nz, ny, nx = hdu.data.shape

    cube = SpectralCube.read(hdu, use_dask=use_dask)

    # corner pixels are outside of the valid region of the projection
    assert all(np.isnan(coord.value) for coord in cube.world[0, 0, 0][1:])

    # ...but the world extrema should ignore the invalid edge pixels
    assert not np.any(np.isnan(cube.longitude_extrema))
    assert not np.any(np.isnan(cube.latitude_extrema))

    # the reported extrema should be those of the valid interior pixels
    dlat = 180.2 / (ny - 1)
    max_lat = np.floor(90. / dlat) * dlat
    assert_allclose(cube.latitude_extrema.value, [-max_lat, max_lat])

    # the repr exercises the world extrema, too
    assert 'nan' not in repr(cube)

    # extracting a latitude strip should keep the full longitude range
    subcube = cube.subcube(ylo=-30 * u.deg, yhi=30 * u.deg)

    lat_strip = (np.arange(ny) - (ny - 1) / 2) * dlat
    nsel = ((lat_strip >= -30) & (lat_strip <= 30)).sum()

    # end-inclusive world slicing rounds to the nearest pixel on each side
    assert subcube.shape in [(nz, nsel, nx), (nz, nsel + 2, nx)]
    assert 'nan' not in repr(subcube)

    # pixel-based slicing should be unaffected
    assert cube.subcube(xlo=1, xhi=4).shape == (nz, ny, 3)


def test_subcube_fully_invalid_wcs(use_dask):
    # companion to test_subcube_allsky: if *all* pixels have invalid world
    # coordinates, the world extrema are NaN and a warning is raised
    hdu = make_allsky_car_hdu()
    # push the entire pixel grid past the pole so no pixel is valid
    hdu.header['CRVAL2'] = 0.0
    hdu.header['CRPIX2'] = hdu.header['NAXIS2'] * 100

    cube = SpectralCube.read(hdu, use_dask=use_dask)

    with pytest.warns(WCSWarning, match='the world extrema cannot'):
        extrema = cube.world_extrema

    assert np.all(np.isnan(extrema))

    with pytest.raises(ValueError, match='cannot convert the world limit'):
        cube.subcube(ylo=-30 * u.deg, yhi=30 * u.deg)


@pytest.mark.skipif('not scipyOK', reason='Could not import scipy')
@pytest.mark.skipif('not regionsOK', reason='Could not import regions')
@pytest.mark.parametrize('regfile',
                         ('255-fk5.reg', '255-pixel.reg'),
                        )
def test_ds9region_255(regfile, data_255, use_dask):
    # specific test for correctness
    cube, data = cube_and_raw(data_255, use_dask=use_dask)

    shapelist = regions.Regions.read(path(regfile))

    subcube = cube.subcube_from_regions(shapelist)
    assert_array_equal(subcube[0, :, :].value,
                           np.array([11, 12, 16, 17]).reshape((2, 2)))


@pytest.mark.skipif('not scipyOK', reason='Could not import scipy')
@pytest.mark.skipif('not regionsOK', reason='Could not import regions')
@pytest.mark.parametrize(('regfile', 'result'),
                             (('fk5.reg', (slice(None), 1, slice(None))),
                              ('fk5_twoboxes.reg', (slice(None), 1, slice(None))),
                              ('image.reg', (slice(None), 1, slice(None))),
                              (
                              'partial_overlap_image.reg', (slice(None), 1, 1)),
                              ('no_overlap_image.reg', ValueError),
                              ('partial_overlap_fk5.reg', (slice(None), 1, 1)),
                              ('no_overlap_fk5.reg', ValueError),
                              ))
def test_ds9region_new(regfile, result, data_adv, use_dask):
    cube, data = cube_and_raw(data_adv, use_dask=use_dask)

    regionlist = regions.Regions.read(path(regfile))

    if isinstance(result, type) and issubclass(result, Exception):
        with pytest.raises(result):
            sc = cube.subcube_from_regions(regionlist)
    else:
        sc = cube.subcube_from_regions(regionlist)

        # Shapes and size should be the same.
        # squeeze on the cube is b/c is retains dimensions of size 1
        assert sc.size == data[result].size
        assert sc.filled_data[:].squeeze().shape == data[result].shape

        # If sizes are the same, values should then be the same.
        assert (sc.unitless_filled_data[:].squeeze() == data[result]).all()

        scsum = sc.sum()
        dsum = data[result].sum()
        assert_allclose(scsum, dsum)

    #region = 'fk5\ncircle(29.9346557, 24.0623827, 0.11111)'
    #subcube = cube.subcube_from_ds9region(region)
    # THIS TEST FAILS!
    # I think the coordinate transformation in ds9 is wrong;
    # it uses kapteyn?

    #region = 'circle(2,2,2)'
    #subcube = cube.subcube_from_ds9region(region)


@pytest.mark.skipif('not scipyOK', reason='Could not import scipy')
@pytest.mark.skipif('not regionsOK', reason='Could not import regions')
def test_regions_spectral(data_adv, use_dask):
    cube, data = cube_and_raw(data_adv, use_dask=use_dask)
    rf_cube = get_rest_value_from_wcs(cube.wcs).to("GHz",
                                                         equivalencies=u.spectral())

    # content of image.reg
    regpix = regions.RectanglePixelRegion(regions.PixCoord(0.5, 1), width=4, height=2)

    # Velocity range in doppler_optical same as that of the cube.
    vel_range_optical = u.Quantity([-318 * u.km/u.s, -320 * u.km/u.s])
    regpix.meta['range'] = list(vel_range_optical)
    sc1 = cube.subcube_from_regions([regpix])
    scsum1 = sc1.sum()

    freq_range = vel_range_optical.to("GHz",
                                      equivalencies=u.doppler_optical(rf_cube))
    regpix.meta['range'] = list(freq_range)
    sc2 = cube.subcube_from_regions([regpix])
    scsum2 = sc2.sum()

    regpix.meta['restfreq'] = rf_cube
    vel_range_gamma = freq_range.to("km/s", equivalencies=doppler_gamma(rf_cube))
    regpix.meta['range'] = list(vel_range_gamma)
    regpix.meta['veltype'] = 'GAMMA'
    sc3 = cube.subcube_from_regions([regpix])
    scsum3 = sc3.sum()

    vel_range_beta = freq_range.to("km/s",
                                    equivalencies=doppler_beta(rf_cube))
    regpix.meta['range'] = list(vel_range_beta)
    regpix.meta['veltype'] = 'BETA'
    sc4 = cube.subcube_from_regions([regpix])
    scsum4 = sc4.sum()

    vel_range_z = freq_range.to("km/s",
                                    equivalencies=doppler_z(rf_cube))
    regpix.meta['range'] = list(vel_range_z)
    regpix.meta['veltype'] = 'Z'
    sc5 = cube.subcube_from_regions([regpix])
    scsum5 = sc5.sum()

    dsum = data[1:-1, 1, :].sum()
    assert_allclose(scsum1, dsum)
    # Proves that the vel/freq conversion works
    assert_allclose(scsum1, scsum2)
    assert_allclose(scsum2, scsum3)
    assert_allclose(scsum3, scsum4)
    assert_allclose(scsum4, scsum5)

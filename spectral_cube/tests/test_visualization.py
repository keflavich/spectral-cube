import pytest

from .test_spectral_cube import cube_and_raw


@pytest.fixture(autouse=True)
def close_figures():
    # ensure that figures don't accumulate (and leak state) across tests
    yield
    try:
        from matplotlib import pyplot
    except ImportError:
        return
    pyplot.close('all')


def test_projvis(data_vda_jybeam_lower, use_dask):
    pytest.importorskip('matplotlib')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    mom0 = cube.moment0()
    mom0.quicklook()


def test_proj_imshow(data_vda_jybeam_lower, use_dask):
    plt = pytest.importorskip('matplotlib.pyplot')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    mom0 = cube.moment0()
    plt.imshow(mom0)


def test_projvis_saved(tmp_path, data_vda_jybeam_lower, use_dask):
    pytest.importorskip('matplotlib')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    mom0 = cube.moment0()
    mom0.quicklook(filename=str(tmp_path / 'test.png'))


def test_projvis_aplpy_deprecated(data_vda_jybeam_lower, use_dask):
    pytest.importorskip('matplotlib')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    mom0 = cube.moment0()
    with pytest.warns(DeprecationWarning, match='use_aplpy'):
        mom0.quicklook(use_aplpy=True)
    with pytest.warns(DeprecationWarning, match='use_aplpy'):
        mom0.quicklook(use_aplpy=False)


def test_mask_quicklook(data_vda_jybeam_lower, use_dask):
    pytest.importorskip('matplotlib')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    cube.mask.quicklook(view=(0, slice(None), slice(None)))


def test_mask_quicklook_wcs(tmp_path, data_vda_jybeam_lower, use_dask):
    pytest.importorskip('matplotlib')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    cube.mask.quicklook(view=(0, slice(None), slice(None)), wcs=cube.wcs,
                        filename=str(tmp_path / 'test.png'))


def test_mask_quicklook_aplpy_deprecated(data_vda_jybeam_lower, use_dask):
    pytest.importorskip('matplotlib')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    with pytest.warns(DeprecationWarning, match='use_aplpy'):
        cube.mask.quicklook(view=(0, slice(None), slice(None)),
                            use_aplpy=True)


def test_to_glue(data_vda_jybeam_lower, use_dask):
    pytest.importorskip('glue')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    cube.to_glue(start_gui=False)


def test_to_pvextractor(data_vda_jybeam_lower, use_dask):
    pytest.importorskip('pvextractor')
    cube, data = cube_and_raw(data_vda_jybeam_lower, use_dask=use_dask)
    cube.to_pvextractor()

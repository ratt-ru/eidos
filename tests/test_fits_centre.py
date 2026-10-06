"""The FITS header must put the beam centre where the beam actually is.

eidos reconstructs the beam on the grid ``np.indices((n, n)) - n / 2`` (``Zernike.unit_disk``), so the
beam centre is at 0-based pixel ``n / 2`` of the full reconstruction, i.e. 1-based ``n / 2 + 1``.  The
writers used to record ``n / 2 + 0.5`` (write_fits_cube / write_fits_single, i.e. the ``-o8`` Jones
cubes) or ``n / 2`` (write_fits), placing the beam half a pixel to a pixel away from where any
consumer of the FITS file (e.g. DDFacet / killMS with Beam.Model FITS) believes it is.

The test generates beams through the CLI entry point and checks that the sub-pixel peak of |xx|
along the first (x) axis -- along which the MeerKAT beam is symmetric -- sits at the header's
CRPIX1 (0-based CRPIX1 - 1), for even and odd pixel counts, full-size and cropped cubes.
"""
import numpy as np
import pytest
from astropy.io import fits

from eidos.create_beam import main

TOL_PX = 0.05


def peak_x(plane):
    """Sub-pixel x (last-axis) position of the maximum of |plane|, by a 3-point parabola."""
    a = np.abs(plane)
    j, i = np.unravel_index(np.argmax(a), a.shape)
    left, mid, right = a[j, i - 1], a[j, i], a[j, i + 1]
    return i + 0.5 * (left - right) / (left - 2 * mid + right)


@pytest.mark.parametrize(
    "pixels, diameter",
    [(64, 10.0), (65, 10.0), (32, 5.0), (33, 5.0)],
    ids=["even-full", "odd-full", "even-cropped", "odd-cropped"],
)
def test_eight_file_cube_centre(tmp_path, pixels, diameter):
    prefix = str(tmp_path / "beam")
    main(["-p", str(pixels), "-d", str(diameter), "-f", "1300", "-P", prefix, "-o8"])
    with fits.open(prefix + "_xx_re.fits") as h:
        hdr, data = h[0].header, h[0].data
    with fits.open(prefix + "_xx_im.fits") as h:
        data = data + 1j * h[0].data
    plane = data.reshape(-1, data.shape[-2], data.shape[-1])[0]
    assert plane.shape == (pixels, pixels)
    offset = peak_x(plane) - (hdr["CRPIX1"] - 1)
    assert abs(offset) < TOL_PX, f"beam peak is {offset:+.3f} px from CRPIX1"


@pytest.mark.parametrize("pixels", [64, 65], ids=["even", "odd"])
def test_single_file_cube_centre(tmp_path, pixels):
    """write_fits (the non -o8 output): axes are (FREQ, H, V, py, px) reversed in FITS order."""
    prefix = str(tmp_path / "beam")
    main(["-p", str(pixels), "-d", "10", "-f", "1300", "-P", prefix])
    with fits.open(prefix + "_re.fits") as h:
        hdr, data = h[0].header, h[0].data
    with fits.open(prefix + "_im.fits") as h:
        data = data + 1j * h[0].data
    plane = data.reshape(-1, data.shape[-2], data.shape[-1])[0]   # first Jones element, first channel
    offset = peak_x(plane) - (hdr["CRPIX1"] - 1)
    assert abs(offset) < TOL_PX, f"beam peak is {offset:+.3f} px from CRPIX1"

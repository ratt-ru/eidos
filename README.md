# Eidos
Primary beam modelling of radio telescope antennae.

For now, only supports **MeerKAT L-band out to 3° radius**.

## Dependencies
scipy, numpy, lmfit, astropy

## Installation
* Directly from github as `pip install git+https://github.com/kmbasad/eidos`.
* If you want to install a local version that you can play around with, first, download using `git clone https://github.com/kmbasad/eidos` and then install the downloaded package as `pip install -e /local/path/to/eidos`. This way you can change anything you want and run your own tests.
* The help file can be seen by typing `eidos -h`.

## Reconstructing beam
* To create a primary beam Jones matrix of MeerKAT for any frequency of L-band run the following command `eidos -p 257 -f 1400` where 257 is the number of pixels on each side, and 1350 MHz is the frequency. If you want to see the beam without normalisation, maybe in order to get the bandpass from the beams directly, then add the option `-N False`. In the output fits file, the first axis correspond to the real (0) and imaginary (1) parts, and the second and third axes the elements of the Jones matrix where 0 represents the Horizontal feed and 1 the Vertical feed.
* If you want Stokes I beam instead of the Jones matrix, run `eidos -p 257 -f 1350 -S I`. You can put `I,Q,U,V` and any of their combinations as inputs fot the `-S` option. The `-S IQ` will give you the leakage map from Q to I and so on for all the other combinations.
* If you want the full Mueller matrix, which will be real for a single dish, because it is the auto-correlation of a complex Jones matrix, then run `eidos -p 257 -f 1350 -S M`.
* To create beam cube for a list of frequencies use `eidos -p 257 -f 1300 1400 5` where 1300 MHz is the start frequency, 1400 MHz is the end, and 5 MHz is the frequency resolution.
* Normally the beams are given as a hypercube containing all the elemnts of the Jones matrix at all frequencies. However, you can get 8 separate files for the real and imaginary parts of the 4 Jones elements by adding the `-o8` option.

## C A U T I O N

* The diameter of the beams can only be **6 degrees**, i. e. the beam can only be created within 6 degrees for now. If you want smaller field of view, please cut the necessary portion after creating a 6 deg beam.
* The model is not available outside a circle of 6 deg diameter; so, if you need a square beam for imaging purposes, please cut a square portion from within the 6 deg circle.
* The beam model is less accurate at higher frequencies, and performs much better below 1400 MHz.
* The diagonal terms of the model beam Jones matrix are much better known than the off-diagonal terms.

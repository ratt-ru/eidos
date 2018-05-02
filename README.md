# Eidos
Primary beam modelling of radio telescope antennae.

For now, only supports **MeerKAT L-band out to 3° radius**.

## Dependencies
scipy, numpy, lmfit, astropy

## Reconstructing beam
* Download the package using: `pip install git+https://github.com/kmbasad/eidos`
* To see the help file: `eidos -h`
* To create a primary beam Jones matrix of MeerKAT for any frequency of L-band run the following command `eidos -p 257 -f 1400`.
* Or to create beam cube for a list of frequencies use `eidos -p 257 -f 1300 1400 5` where 1300 MHz is the start frequency, 1400 MHz is the end, and 5 MHz is the frequency resolution.

## C A U T I O N

* The diameter of the beams can only be **6 degrees**, i. e. the beam can only be created within 6 degrees for now. If you want smaller field of view, please cut the necessary portion after creating a 6 deg beam.
* The model is not available outside a circle of 6 deg diameter; so, if you need a square beam for imaging purposes, please cut a square portion from within the 6 deg circle.
* The beam model is less accurate at higher frequencies, and performs much better below 1400 MHz.
* The diagonal terms of the model beam Jones matrix are much better known than the off-diagonal terms.

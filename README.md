# Eidos
Primary beam modelling of radio telescope antennae.

For now, only supports **MeerKAT L-band out to 10° diameter**. These beams can be calculated from both Holographic (AH) observation and EM simulation. The beam datasets are calculated from a specific AH and EM dataset and should not be considered standard for MeerKAT.

## Dependencies
scipy, numpy, astropy

## Installation
Directly from github as `pip install git+https://github.com/ratt-ru/eidos`.

If you want to install a local version that you can play around with, first, download using 

`git clone https://github.com/ratt-ru/eidos`  

and then install the downloaded package as 

`pip install -e /local/path/to/eidos`. 

This way you can change anything you want and run your own tests.

The help file can be seen by typing `eidos -h`.

## Reconstructing beam
To create a primary beam Jones matrix of MeerKAT for any frequency of L-band run

`eidos -p 256 -d 10 -f 1070` 

where 256 is the number of pixels on each side, 10 is the diameter of the beam in degrees, and 1070 MHz is the frequency. 

You can add `-c me` option to get the EM beam models. By default, it will produce models from holography: `-c mh` 

If you want Stokes I beam instead of the Jones matrix, run 

`eidos -p 256 -d 10 -f 1070 -S I`.

You can put `I,Q,U,V` and any of their permutations as inputs fot the `-S` option. The `-S IQ` will give you the leakage map from Q to I and so on for all the other permutations.

If you want the full Mueller matrix, which will be real for a single dish, because it is the auto-correlation of a complex Jones matrix, then run 

`eidos -p 256 -d 10 -f 1070 -S M`.

To create beam cube for a list of frequencies use 

`eidos -p 256 -d 10 -f 1300 1400 5` 

where 1300 MHz is the start frequency, 1400 MHz is the end, and 5 MHz is the frequency resolution. This operation is parallelized by default.

Normally the beams are given as a hypercube containing all the elemnts of the Jones matrix at all frequencies. However, you can get 8 separate files for the real and imaginary parts of the 4 Jones elements by adding the `-o8` option.

## C A U T I O N

* The model is not available outside a circle of 10 deg diameter; so, if you need a square beam for imaging purposes, please cut a square portion from within the 10 deg circle.
* The beam model is less accurate at higher frequencies, and performs much better below 1400 MHz.
* The diagonal terms of the model beam Jones matrix are much better known than the off-diagonal terms.
* Finally, these beams are as good as the given AH and EM datasets. As we get more accurate AH models (from M. de Villiers at SARAO) and EM simulations (from Robert Lehmensiek at EMSS Antennas), this pipeline can be used to create more accurate sparse representation of primary beams using Zernike polynomials.

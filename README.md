# Eidos
Primary beam modelling of radio astronomy antennas:  
Paper I: [Modelling the Karl G. Jansky Very Large Array (VLA) L-band beam using holography](https://academic.oup.com/mnras/article/485/3/4107/5374534)  
Paper II: [Modelling the MeerKAT L-band beam](https://arxiv.org/abs/1904.07155)  

The current version can be used to create only MeerKAT L-band beams from both holographic (AH) observations and EM simulations within a maximum diameter of 10 degrees. L-band AH beam models for JVLA and the UHF-band models for MeerKAT will be added soon.  

## Dependencies
scipy, numpy, astropy

## Installation
`pip install eidos`

To create a local developer version:  
`git clone https://github.com/ratt-ru/eidos`  
`pip install -e /local/path/to/eidos`  

Help available via `eidos -h`

## Creating beam models
To create a primary beam Jones matrix of MeerKAT for any frequency of L-band run  
`eidos -p 256 -d 10 -f 1070`  
where `256` is the number of pixels on each side, `10` is the diameter of the beam in degrees, and `1070` MHz is the frequency. 

One can add `-c me` (`m` for MeerKAT, `e` for EM) option to get the EM models. The default is `-c mh` for AH.

To create Stokes I beam instead of the Jones matrix:  
`eidos -p 256 -d 10 -f 1070 -S I`

One can use `I,Q,U,V` and any of their permutations as inputs fot the `-S` option.  
`-S IQ` will produce the leakage map from Q to I and so on for all the other permutations.  

To create full Mueller matrix, which will be real for a single dish because it is the auto-correlation of a complex Jones matrix:  
`eidos -p 256 -d 10 -f 1070 -S M`

To create beam cube for a list of frequencies:  
`eidos -p 256 -d 10 -f 1300 1400 5`  
where 1300 MHz is the start frequency, 1400 MHz is the end, and 5 MHz is the frequency resolution. This operation is parallelized by default.

Normally the beams are given as a hypercube containing all the elemnts of the Jones matrix at all frequencies. However, one can get 8 separate files for the real and imaginary parts of the 4 Jones elements by adding the `-o8` option.  

## C A U T I O N

* The model is not available outside a circle of 10 deg diameter; so, if you need a square beam for imaging purposes, please cut a square portion from within the 10 deg circle.
* The beam model is less accurate at higher frequencies, and performs much better below 1400 MHz.
* The diagonal terms of the model beam Jones matrix are much better known than the off-diagonal terms.
* Finally, these beams are as good as the given AH and EM datasets. As we get more accurate AH models (from M. de Villiers at SARAO) and EM simulations (from Robert Lehmensiek at EMSS Antennas), this pipeline can be used to create more accurate sparse representation of primary beams using Zernike polynomials.

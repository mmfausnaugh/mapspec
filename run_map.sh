#!/usr/bin/bash

#An example of how to run mapspec.  It assumes the mapspec directory
#lives in your home directory.

#ref.smooth.txt is a 3 column ascii file with the reference spectrum
#that you want to calibrate to

#oiii.window is a file with wavelengths for the emission line used as
#a reference (usually oxygen 5007 angstroms for optical spectra).  It
#is 2 columns, 3 lines, like this:

#lineblue linered
#c1blue   c1red
#c2blue   c2red

#lineblue/linered and the bluest/reddest wavelengths of the emission
#line.  c1 and c2 are the continuum windows on either side of the
#line---it doesn't matter if c1 or c2 is bluer, but on a given line,
#you must have blue then red wavelengths.  See EmissionLine class in
#spectrum.py for more info.

#linear is the type of interpolation.  Other options exists, but
#linear is the only option that rigorously does the errors

#speclist_use is a 1 column ascii file, corresponding to the list of
#spectra that will be matched to the reference.

python ~/python/mapspec/do_map.py ref.smooth.txt oiii.window  linear speclist_use
#d="$(pwd)"
#echo "mapspec finsished in $d" | mail -s 'process is done' your_email@address.edu

#see do_map.py for more---as a default, it will output rescaled
#spectra, covariance matrices, MCMC chains, and a file listing the
#chi^2, acceptance fractions, and parameters for each spectrum.  Note
#that it does some crude model comparisons between smoothing with a
#Gaussian, Gauss-Hermite polynomial, and delta-function kernels.
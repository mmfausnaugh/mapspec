#!/usr/bin/bash

#python ./do_map.py ref.smooth.txt oiii.window linear speclist_use mapspec.params no_covar chains
python ../../scripts/do_map.py ref.smooth.txt oiii.window linear speclist_use mapspec.params no_covar chains
#d="$(pwd)"
#echo "mapspec finsished in $d" | mail -s 'process is done' your_email@address.edu


#An example of how to run mapspec.  It assumes you will run this from
#the mapspec/examples/mapspec_test/ directory---in general, you should
#change the relative path to suite your own needs.

#ref.smooth.txt is a 3 column ascii file with the reference spectrum
#that you want to calibrate to

#oiii.window is a file with wavelengths for the emission line used as
#a reference (usually oxygen 5007 angstroms for optical spectra).  It
#is 2 columns, 3 lines, like this:

#lineblue linered
#c1blue   c1red
#c2blue   c2red

#lineblue/linered are the bluest/reddest wavelengths of the emission
#line.  c1 and c2 are the continuum windows on either side of the
#line---it doesn't matter if c1 or c2 is bluer, but on a given window
#you must specify blue before red wavelengths.  See EmissionLine class
#in spectrum.py for more info.

#linear is the type of interpolation.  Other options exists, but
#linear is the only option that rigorously does the errors

#speclist_use is a 1 column ascii file, corresponding to the list of
#spectra that will be matched to the reference.

#mapspec.params is an output file that lists the chi^2, best fit
#parameters, and acceptence fraction for each spectrum.

#no_covar is a keyword meaning 'don't write out the covariance
#matrix.'  These are pretty big and eat up memory fast, so only save
#them if you need them

#chains is a keyword meaning "write out the MCMC chains"--these give
#you posterior distributions for the rescaling parameters.  To turn
#off, set 'no_chains'

#see do_map.py for more.  As a default, it will output rescaled
#spectra, MCMC chains, and a summary file (mapspec.params). Note that
#do_map.py does some crude model comparisons between smoothing with
#Gaussians, Gauss-Hermite polynomials, and delta-functions.

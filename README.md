# mapspec #

mapspec stands for **M**cmc **A**lgorithm for **P**arameters of **SPEC**tra.  It was originally designed to rescale astronomical spectra to a standard flux scale, following the formalism of [van Groningen & Wanders 1992](http://adsabs.harvard.edu/abs/1992PASP..104..700V).  Along the way, it has become a general purpose package for analyzing spectra, particularly those with bright emission lines (such as AGN).

* * *

# Installation and Dependencies#

mapspec is a fairly simple python package---just download it and put it in your python installation, your PYTHONPATH, or working directory.

You will also need the following python packages installed:

* numpy
* scipy
* matplotlib (for using built-in plotting functions)
* astropy (for reading and writting fits files)

* * *

# Outline #

mapsec is object-oriented and is desigend to be modular and highly extensible.  It makes heavy use of numpy and scipy libraries, and was designed to easily integrate with these packages.

There are two main parts of the program.

First, there are general-usage spectrum utilities, which are bundled in `spectrum.py`.  Second is the rescaling procedure, which is kept in `mapspec.py`.  We begin by describing (with illustrative examples) the spectrum utilities.

* * *

## Spectrum ##

The primary object is the `Spectrum` class.  `Spectrum` has 3 fundamental attributes:  wavelength (`Spectrum.wv`), flux  (`Spectrum.f`), and flux error (or measurement error, `Spectrum.ef`).  Each attribute is a numpy array, so all the functionality of numpy/scipy is preserved.  For example, to get the mean flux-level:

    from mapspec.spectrum import Spectrum
    my_spec = Spectrum()
    print my_spec.f.mean()

or to get the number of pixels:

    print my_spec.wv.size

Spectrum also has many useful methods to interpolate, rebin, extinction-correct, and smooth (or log-smooth) the data.  For example:

    import numpy as np
    dx = my_spec.wv[1] - my_spec.wv[0] #assumes even wavelength spacing throughout
    xnew = np.arange(my_spec.wv.min(), my_spec.wv.max() + 2.*dx, 2.*dx)
    my_spec.rebin(xnew)

will rebin the spectrum to twice the original wavelength spacing.

* * *

## Emission Lines ##

The next object is the `EmissionLine` class.  `EmissionLine` inherits from `Spectrum`, so it can be used in approximately the same way.  To instantiate an `EmissionLine`, you need to specify the line window and surrounding continuum:

    from mapspec.spectrum import EmissionLine
    oiii = [ [5065 5115],     #OIII 5007, at redshift 1.018
             [5059 5065],     #blue continuum
             [5115 5129]  ]   # red continuum
    my_line = EmissionLine(my_spec, oiii[0], [oiii[1],oiii[2]] )

At instantization, `EmissionLine` fits a linear model to the local continuum (`oiii[1]` and `oiii[2]`) and subtracts the model underneath the emission line.  That is, `my_line.wv` is the wavelengths of `my_spec` between `oiii[0,0]` and `oiii[0,1]` and `my_line.f` is the corresponding section of `my_spec.f` minus the continuum fit.

`EmissionLine` also has other useful methods, such as `my_line.integrate_line()`, which will use `scipy.integrate.simps` to return the total line flux, or `my_line.fwhm()`, which will return the full-width-at-half-maximum of the line.

* * *

## Line Model ##

The last class is `LineModel`.  At instantization, `LineModel` fits a function to an input `EmissionLine`:

    from mapspec.spectrum import LineModel
    my_lmodel = LineModel(my_line, func='gauss-hermite')

The _raison d'être_ of `LineModel` is to generate flux predictions at arbitrary wavelengths.  This is accomplished by overloading the class's call method, so you can think of `LineModel` as a kind of function:

    xeval = np.linspace( my_line.wv.min(), my_line.wv.max(), 500 )
    fpredict = my_lmodel(xeval)

You could use these predictions to, for example, subtract the line from the original spectrum:

    mask = ( my_spec.wv > oiii[0,0])*(my_spec.wv < oiii[0,1] )
    fsubtract = my_spec.f[mask] - my_lmodel(my_spec.wv[mask])

* * *
## Reading Data ##

Input from a file is handled with inheritance.  If you know what format your data are in, the following classes know how to read these formats and return fully functioning `Spectrum` objects:

    from mapspec.spectrum import TextSpec, TextSpec_2c, FitsSpec
    
    s1 = TextSpec('mockdata.txt')  #for 3 column format
    s2 = TestSpec_2c('mockdata_2col.txt') #for 2 columns, assumes no errors
    s3 = FitsSpec('mockdata.fits')

`FitsSpec` contains many keywords in order to deal with the vast diversity of .fits formats:

    s3 = FitsSpec('mockdata.fits', extension = 1, data_axis = 1, error_axis = 2,
                  x1key = 'CRVAL1',
                  dxkey = 'CD1_1',
                  p1key = 'CRPIX1'
                  )


* * *

## Rescaling Spectra ##

The tools for rescaling spectra are in `mapspec.py`.  Most of the work is done by the `RescaleModel` class:

    from mapspec.spectrum import *
    from mapspec.mapspec import RescaleModel
    
    sref = TextSpec('ref.txt')
    lref = EmissionLine(sref, oiii[0], [ oiii[1], oiii[2] ] )

    my_model = RescaleModel(lref, kernel = 'Gauss')

`my_model` will take an input emission line and try to align it with `lref`.  The model consists of a wavelength shift, rescaling factor, and smoothing---the form of the smoothing kernal is specified by the `kernel` keyword (in the above example, a simple Gaussian).

The parameters of the model are stored in a dictionary, so they are easy to access.  If you know what parameters you want, it is easy to transform the line:

    #shift by 0.25 x-units, multiply by 2, and smooth with kernel of width 1 x-units
    my_model.p = {'shift': 0.25, 'scale': 2.0, 'width': 1.0} 
    sout, lost_pix = my_model.output(my_spec)

`lost_pix` is a boolean array that tells you how much of the original spectrum is lost---you loose at least one pixel on either side for any non-zero shift.  If you don't care about the original wavelength grid, it is suggested to change the wavelengths in place (i.e., `my_spec.wv -= shift`).

Usually, you will want to optimize these parameters so that they match the reference.  To do this, a function `metro_hast` is included:

    from mapspec.mapspec import metro_hast
    chi2, p_best, frac_accept = metro_hast(5000, my_line, my_model)

This will run an MCMC for 5000 steps, and try to optimize `my_model.p` so that `my_line` matches `lref` (which was specified when instantiating `my_model`).  `chi2` is the best (minimum) chi^2, p_best are the corresponding parameters, and frac_accept is the acceptance fraction of the chain.

You then apply the model as above:

    my_model.p = p_best
    sout,lost_pix =  my_model.output(my_spec)

`metro_hast` has some other useful options:

    chi2, p_best, frac_accept, p_chain = metro_hast(5000, my_line, my_model, keep = True)

This will return the MCMC chain of parameters in `p_chain`.  `p_chain` is actually its own object (`mapspec.mapspec.Chain`), but it is little more than a glorified NxM array that knows how to plot itself (it also stors the log-likelihood for each set of parameters in a separate attribute).


# Details #

## Important Warning ##

Error propagation indicates that for a given operation on the flux array, the same operation should be carried out on the *square* of the error array.  Unfortunately, the squares are somewhat less stable numerically, which can result in strange behavior---for example, if you interpolate the array with b-splines, it is possible that the interpolated variances (square-errors) will go negative.

Since this depends both on the input data and the operation, mapspec does nothing more than print a warning, and then happily take the square-root of the negative variances as if nothing was wrong.  This procedure serves to flag the troublesome points with complex values in the errors.

It is up to the user to decide how to deal with this issues.  For example, one might change the operation (maybe linear interpolation instead of b-splines), or modifying the data (either changing the errors before hand, or masking/editing the flagged fluxes/uncertainties).

## Calculation of the Likelihood ##

mapspec aligns the data and the model by minimzing: (data - model)^2/error^2, where data is the input line, model is the reference value chosen at instantization, and error^2 = (reference error)^2 + (data error)^2.  In other words, mapspec takes a maximum likelihood approach, assuming normally distributed residuals and independent uncertainties on the data and reference.

In order to keep the likelihood well-defined, mapsepc raises an exception if the variance array goes negative during the fit.  This can be a problem for spectra with sharp features and models with complicated smoothing kernels, and is an exception to the above note (where the negative variances are ignored, except for being flagged by complex errors).

Because of wavelength shifts and edge-effects from convolution, there can be problems aligning the edges of the fitting region---in order to censor these points, it would be necessary to change the amount of data, and therefore the number of degrees-of-freedom, during the fit.  To prevent this while mitigating the edge-effects, mapspec ignores 10% of the data (5% on each end) when calculating the likelihood.  Although this is not an optimal solution, it does produces reasonable results and is an acceptable compromise for the sake of simplicity  and practicability.

It is therefore suggested to choose the fitting region (line wavelengths) slightly larger than would naively be expected for isolating the emission line flux.  It is also suggested to remove large shifts (greater than 1 pixel) before performing the fit.  A convenience function `mapspec.mapspec.get_cc` is provided to cross correlate two spectra, and calculate the relative shift to the nearest pixel (but reported in wavelength units).  See `do_map.py` for an example of how to apply this.
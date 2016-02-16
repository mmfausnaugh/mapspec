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

`Spectrum` also has many useful methods to interpolate, rebin, extinction-correct, and smooth (or log-smooth) the data.  For example:

    import numpy as np
    dx = my_spec.wv[1] - my_spec.wv[0] #assumes even wavelength spacing throughout
    xnew = np.arange(my_spec.wv.min(), my_spec.wv.max() + 2.*dx, 2.*dx)
    my_spec.rebin(xnew)

will rebin the spectrum to twice the original wavelength spacing.

* * *

## Emission Lines ##

The next object is the `EmissionLine` class.  `EmissionLine` inherits from `Spectrum`, so it can be used in approximately the same way.  To instantiate an `EmissionLine`, you need to specify the line window and surrounding continuum:

    from mapspec.spectrum import EmissionLine
    oiii = [ [5065, 5115],     #OIII 5007, at redshift 1.018
             [5059, 5065],     #blue continuum
             [5115, 5129]  ]   # red continuum
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
    my_spec.f[mask] -= my_lmodel(my_spec.wv[mask])

If the parameters of the fit are of interest, you can access then with `my_lmodel.p` (and the covariance matrix, for errors, is stored in `my_lmodel.covar`), but you will have to check the individual functions to see which parameter is which.

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

This will run an MCMC for 5000 steps, and try to optimize `my_model.p` so that `my_line` matches `lref` (which was specified when instantiating `my_model`).  `chi2` is the best (minimum) chi^2, `p_best` are the corresponding parameters, and `frac_accept` is the acceptance fraction of the chain.

You then apply the model as above:

    my_model.p = p_best
    sout,lost_pix =  my_model.output(my_spec)

`metro_hast` has some other useful options:

    chi2, p_best, frac_accept, p_chain = metro_hast(5000, my_line, my_model, keep = True)

This will return the MCMC chain of parameters in `p_chain`.  `p_chain` is actually its own object (`mapspec.mapspec.Chain`), but it is little more than a glorified NxM array that knows how to plot itself (it also stores the log-likelihood for each set of parameters in a separate attribute).


# Details #

## Important Warning ##

Error propagation indicates that for a given operation on the flux array, the same operation should be carried out on the **square** of the error array.  Unfortunately, the squares are numerically somewhat less stable, which can result in strange behavior---for example, if you interpolate the array with b-splines, it is possible that the interpolated variances (square-errors) will go negative.

Since this depends both on the input data and the operation, mapspec does nothing out of the normal except print a warning---it then happily take the square-root of the negative variances as if nothing were wrong.  This procedure serves to flag the troublesome points with complex values in the errors.

It is up to the user to decide how to deal with this issues.  For example, one might change the operation (maybe linear interpolation instead of b-splines), or modifying the data (either changing the errors beforehand, or masking/editing the flagged fluxes and uncertainties).

## Calculation of the Likelihood and Fitting Procedure##

mapspec aligns the data and the model by minimzing: (data - model)^2/error^2, where data is the input line, model is the reference value chosen at instantization, and error^2 = (reference error)^2 + (data error)^2.  In other words, mapspec takes a maximum likelihood approach, assuming normally distributed residuals and independent uncertainties on the data and reference.

In order to keep the likelihood well-defined, mapsepc raises an error if the variance array goes negative during the fit.  This can be a problem for spectra with sharp features and models with complicated smoothing kernels.  This treatment is an exception to the above note (where the negative variances are largely ignored).

Because of wavelength shifts and edge-effects from convolution, there can be problems aligning the edges of the fitting region---in order to censor these points, it would be necessary to change the amount of data, and therefore the number of degrees-of-freedom, during the fit.  To prevent this while mitigating the edge-effects, mapspec ignores 10% of the data (5% on each end) when calculating the likelihood.  Although this is not an optimal solution, it does produces reasonable results and appears to be a reasonable compromise between correctness, simplicity,  and practicability.

It is therefore suggested to choose the fitting region (line wavelengths) slightly larger than would naively be expected for isolating the emission line flux.  It is also suggested to remove large shifts (greater than 1 pixel) before performing the fit.  A convenience function `mapspec.mapspec.get_cc` is provided to cross correlate two spectra---see `do_map.py` for an example.

In the Bayesian framework, priors are *always* added into the likelihood calculation.  mapspec assumes uniform priors on all parameters, so the priors do not explicitly enter into the calculation of the log-likelihood.  Note that a flat prior *is informative* for the case of the scaling parameter---however, because the rescaling is a factor of a few (not many orders of magnitude), this situation makes little difference.

When using a Gauss-Hermite Kernel, h3 and h4 are constrained to be > -0.3 and < 0.3; experimentation shows that this range is sufficient to produce a wide range of emission line shapes.  The finite interval is enforced by setting the prior probability to 0 outside of this interval.  Similarly, strange behavior will occur if the kernel width drops below the spectral resolution (a kind of undersampling/aliasing effect).  A minimum kernel width is imposed at 1/2 of the wavelength spacing of the reference data.

`RescaleModel` objects have the option to add priors, either through a function provided by the user or from posterior distributions stored in `Chain` objects.  See `do_map.py` for an example.

* * *
# Methods #

## Spectrum ##

| Object | Method | args | Description |
|--------|--------|------|-------------|
|`Spectrum`| `__init__`| `style = 'sinc'`| Change the interpolation style at instantization (implicitly calls `set_interp`).|
|`Spectrum`|`restore`|    | Restore attributes (`Spectrum.wv`,`Spectrum.f`, and `Spectrum.ef`) to values when spectrum was instantiated.|
|`Spectrum`|`interp` | `xnew`| Returns interpolated flux and flux errors at location of array `xnew`.|
|`Spectrum`|`set_interp`| kwargs| Change interpolation style (see `style` below).  Keywords for custom b-splines and sinc interpolation are also provided here.|
|`Spectrum`|`rebin`| `xnew` | Rebins the spectrum to wavelengths `xnew`.  Modifies `Spectrum.wv`,`Spectrum.f`, and `Spectrum.ef`.|
|`Spectrum`|`extinction_correct`| `E_BV`, `RV`| De-extinguish the spectrum for a [CCM89](http://adsabs.harvard.edu/abs/1989ApJ...345..245C) extinction law, for a given E(B-V) (`E_BV`) and R_V (`RV`).  Modifies `Spectrum.f`, and `Spectrum.ef`.|
|`Spectrum`|`smooth`| `width`,`name`| Smooths the spectrum with a kernel of width `width` specified with `name` (calculated with `scipy.signal.get_window`).  Modifies `Spectrum.f`, and `Spectrum.ef`.  Edge-effects are treated by replacing with the original spectrum (a number of pixels equal to `width` on both edges).|
|`Spectrum`|`velocity_smooth`| `v_width`| Smooths the spectrum with a velocity dispersion of width `v_wdith` (in km/s) by rebinning evenly in log wavelength.  Modifies `Spectrum.f`, and `Spectrum.ef`.  Converts back to original wavelengths, but cannot support sinc interpolation (unevenly spaced wavelengths).|
|`Spectrum`| `wv`| Attribute| numpy array of wavelengths.|
|`Spectrum`| `f`| Attribute| numpy array of fluxes.|
|`Spectrum`| `ef`| Attribute| numpy array of flux errors.|
|`Spectrum`| `style`| Attribute|  String with instructions for `set_interp`.  Options are anything supported by `scipy.interpolate.interp1d`.  Also available is a custom sinc interpolator (`style = 'sinc'`), which uses convolution in real space and a Lanczos window to prevent ringing (very slow).  Bsplines (`style= 'bspline'`) are also available (using a custom class to interface with `splrep` and `splev` from `scipy.interpolate` ).


## EmissionLine ##

| Object | Method | args | Description |
|--------|--------|------|-------------|
|`EmissionLine`|`__init__`|`Spec`, `window`, `cwindow`| `Spec` is a `Spectrum` object from which to extract the `EmissionLine`.  `window` gives the wavelength region of the line (inclusive, blue edge first).  `cwindow` is a 2x2 list/array that specifies the continuum.  Each row is a continuum region (exclusive), which must have the blue edge first, but the blue and red regions may be in any order.|
|`EmissionLine`|`integrate_line`|        | Returns total line flux using `scipy.integrate.simps`.|
|`EmissionLine`|`equivalent_width`|      | Returns integrated flux divided by mean continuum.|
|`EmissionLine`|`wv_mean` |  | Returns first moment of the line.|
|`EmissionLine`|`wv_median`|  | Returns median wavelength (wavelength that splits the integrated flux in 2).|
|`EmissionLine`|`dispersion`| | Returns second moment of the line.|
|`EmissionLine`|`mad`| | Returns Mean Absolute Deviation.|
|`EmissionLine`|`fwhm`| `center`|  Returns a full-width-at-half-maximum, following the procedure in [Peterson et al. 2004](http://adsabs.harvard.edu/abs/2004ApJ...613..682P); `center` is the central wavelength, to check if the line is double-peaked.  Actually returns a tuple with fwhm, width on blue side, width on red-side, and a boolean for which `True` indicates double-peaked.|
|`EmissionLine`|`ip_wv`| `frac`| Returns inter-percentile region in wavelength units.  The inter-percentile region contains `frac` of the total flux, where `frac` is between 0 and 1.  Actually returns a tuple with the inter-percentile region, the corresponding blue wavelength, and the corresponding red wavelength.|
|`EmissionLine`|`ipv`| `frac`,`center`| Returns the inter-percentile region in velocity units.  `frac` is defined as above, `center` is used to convert to velocities.  Actually returns a tuple with the inter-pecentile region, blue width, and red width.  Assumes km/s.|
|`EmissionLine`| `wv`| Attribute| numpy array of wavelengths (subset of `Spec.wv` within `window`).|
|`EmissionLine`| `f`| Attribute| numpy array of fluxes (matching subset of `Spec.f` minus continuum flux).|
|`EmissionLine`| `ef`| Attribute| numpy array of flux errors (as `EmissionLine.f`, but error propagation from the continuum fit is included).|
|`EmissionLine`| `cf`| Attribute| Continuum flux underneath the line.|
|`EmissionLine`| `ecf`| Attribute| Error on continuum flux.  This is added in quadrature to the original measurement errors in the line wavelength region.|

## LineModel ##

| Object | Method | args | Description |
|--------|--------|------|-------------|
|`LineModel`|`__init__`| `EmLine`,`func='gaussian'`,`window=None`,`nline = 1`,`floating=False`,`linedata=None`|  `EmLine` is an `EmissionLine` to which the model is fit.  `func` is a string that specifies the functional form to fit to the line; can be 'gaussian', 'gauss-hermite', 'lorentzian', 'approx-voig' (sum of lorentzian and gaussian), or 'data' (empirical template).  `window` is a boolean array that masks the input `EmissionLine` (`True` indicates that the pixel is included in the fit).  `nline` indicates the number of lines (components) used in the fit.  `floating=True` adds a parameter to the model that allows a constant offset or 'pedestal' flux.  `linedata` is an `EmissionLine` that serves as a template if `func = 'data'`.  See examples in `line_fitting.py`.|
|`LineModel`|`__call__`|`x`| Returns fitted function evaluated at input `x` locations (i.e., flux-predictions, `x` is an array).|
|`LineModel`|`chi2`| | Returns the value of chi^2 from the initial fit.
|`LineModel`|`plot_components`| `x`,`ax`,`cuse`| Plots individual components (if `nline` is > 1) on axis `ax` in color `cuse` at abcissas `x`.|
|`LineModel`|`p`|Attribute| Parameters of the fitted model.|
|`LineModel`|`covar`|Attribute| Covariance matrix of fitted parameters.|
|`LineModel`| `wv`| Attribute|  Saves `EmLine.wv`.|
|`LineModel`| `lf`| Attribute|  Saves `EmLine.f`.|
|`LineModel`| `elf`| Attribute|  Saves `EmLine.ef`.|
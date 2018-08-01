# meas_extensions_convolved

This is a measurement plugin for the LSST stack that convolves images before measuring aperture fluxes.
The convolution provides a rough match to the PSF in other observations, and so aperture measurements on the convolved image are useful for:

1. Measuring the flux that would be obtained with a fiber spectrograph, by using an aperture equal to the size of the spectrograph fiber and measuring with multiple target seeings to allow interpolation of the flux at the seeing of the spectroscopic observation.
2. Measuring the color of galaxies in variable seeing, by using a common aperture and target seeing in each band.

This is similar to the usual degrade-to-a-common-seeing-then-measure (PSF-matched coadds) technique, except that the matching is done for each individual source, rather than over the entire image at once.

In the interests of simplicity and speed, we assume that all the PSFs are Gaussian and explicitly do not allow for deconvolution.
This does not provide a precise match to the desired target PSF, but should be close enough for most purposes.

The aperture measurements that are made are:

1. Multiple pre-defined circular apertures.
2. Elliptical Kron aperture, where the shape and Kron radius are provided from the catalog.


## Key configuration parameters

* `seeing`: a list of target FWHM seeings, in pixels.
* `aperture.radii`: a list of aperture radii, in pixels.
* `maxSincRadius` and `aperture.maxSincRadius`: maximum radius for which to use the 'sinc' aperture measurement algorithm ([Bickerton & Luption, 2013](2013MNRAS.431.1275B)), which provides greater precision in the fluxes for small radii apertures on properly-sampled images.
* `kronRadiusName`: column name with the Kron radius.
* `kronRadiusForFlux`: number of Kron radii to use for the Kron aperture.


## Output columns

* `<algName>_flag`: general failure flag; if `True`, the measurement algorithm failed catastrophically (e.g., an unexpected exception was thrown).
* `<algName>_seeing`: PSF determinant radius at source position, in pixels.
* `<algName>_<seeingIndex>_deconv`: deconvolution flag; if `True`, the PSF at the position is wider than the target seeing, so PSF matching would require deconvolution. No aperture measurements will have been made with this target seeing.
* `<algName>_<seeingIndex>_<aperName>_flux`: flux within circular aperture after convolution to target seeing.
* `<algName>_<seeingIndex>_<aperName>_fluxErr`: estimated error for circular aperture flux measurement after convolution to target seeing. See [this note](#errors) for a warning about their reliability.
* `<algName>_<seeingIndex>_<aperName>_flag`: aperture flux flag; if `True`, the corresponding aperture flux measurement is unreliable. This may be set, e.g., if the aperture extends off the edge of the image.
* `<algName>_<seeingIndex>_kron_flux`: flux within Kron aperture after convolution to target seeing.
* `<algName>_<seeingIndex>_kron_fluxErr`: estimated error for Kron flux measurement after convolution to target seeing. See [this note](#errors) for a warning about their reliability.
* `<algName>_<seeingIndex>_kron_flag`: Kron flux flag; if `True`, the Kron flux measurement is unreliable. This may be set, e.g., if the Kron aperture extends off the edge of the image.

In the above list of output columns:

* `<algName>` stands for the algorithm name, typically `"ext_convolved_ConvolvedFlux"`.
* `<seeingIndex>` stands for the index of the target seeing of interest (starting from `0`).
* `<aperName>` stands for the aperture name, which is the radius in pixels to one decimal precision and the decimal point replaced by an underscore (e.g., for a radius of 1.2 pixels, this would be `1_2`).

There are some useful methods on `ConvolvedFluxConfig` that help put the column names together.  Note that these require that the configuration (especially the seeing values) is set correctly.

* `getBaseNameForSeeing(seeing)`: base name for measurement, given seeing (i.e., `<algName>_<seeingIndex>`).
* `getApertureResultName(seeing, radius)`: name for aperture measurement result (i.e., `<algName>_<seeingIndex>_aperName>`).
* `getKronResultName(seeing)`: name for Kron measurement result (i.e., `<algName>_<seeingIndex>_kron`).
* `getAllApertureResultNames()`: list of all aperture result names (the results of `getApertureResultName` called for all seeings and aperture radii).
* `getAllKronResultNames()`: list of all Kron result names (the results of `getKronResultName` called for all seeings).
* `getAllResultNames()`: list of all aperture and Kron result names (the results of `getAllApertureResultNames` and `getAllKronResultNames`).


## Errors
<a name="errors"></a>

Note that no effort is made to account for covariance when estimating errors, and hence the reported errors are likely underestimated.  This is especially true for coadds and other images that already have covariance.


## Configuration for use with `SingleFrameMeasurementTask`

This is a plugin, so all you need to do is load, activate and configure the plugin, and the `SingleFrameMeasurementTask` will run it.
To do so, use something like the following in a configuration file:

    # Load plugin
    import lsst.meas.extensions.convolved
    # Activate plugin
    config.plugins.names.add("ext_convolved_ConvolvedFlux")
    # Configure plugin
    config.algorithms["ext_convolved_ConvolvedFlux"].seeing = [2.5, 4.0, 6.0]
    config.algorithms["ext_convolved_ConvolvedFlux"].aperture.radii = [3.0, 5.0, 7.0]
    
If you haven't done so elsewhere, you will also need to load and activate the Kron measurement plugin in order to get Kron radii for the convolved measurement:

    # Load and activate Kron
    import lsst.meas.extensions.photometryKron
    config.plugins.names.add("ext_photometryKron_KronFlux")
    
In the above snippets, `config` is a `SingleFrameMeasurementConfig`.  I recommend putting the above configuration in a new file (e.g., `/path/to/convolvedConfig.py`) and loading that file from a top-level task configuration file.  For example, with `processCcd.py`:

    config.charImage.measurement.load("/path/to/convolvedConfig.py")  # to calculate aperture corrections
    config.calibrate.measurement.load("/path/to/convolvedConfig.py")


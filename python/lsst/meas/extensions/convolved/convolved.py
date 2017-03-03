#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import math
import numpy as np

from lsst.pex.config import Config, Field, ListField, ConfigField, makeConfigClass
from lsst.pipe.base import Struct
from lsst.meas.extensions.photometryKron import KronAperture, KronFluxPlugin
from lsst.meas.base.wrappers import WrappedSingleFramePlugin, WrappedForcedPlugin

import lsst.meas.base
import lsst.afw.math
import lsst.afw.geom
import lsst.afw.image

__all__ = ("SingleFrameConvolvedFluxPlugin", "SingleFrameConvolvedFluxConfig",
           "ForcedConvolvedFluxPlugin", "ForcedConvolvedFluxConfig",)


SIGMA_TO_FWHM = 2.0*math.sqrt(2.0*(math.log(2.0)))  # Multiply sigma by this to get FWHM
PLUGIN_NAME = "ext_convolved_ConvolvedFlux"  # Usual name for plugin


class DeconvolutionError(RuntimeError):
    """Convolving to the target seeing would require deconvolution"""
    pass


ApertureFluxConfig = makeConfigClass(lsst.meas.base.ApertureFluxControl)


class ConvolvedFluxData(Struct):
    """A `lsst.pipe.base.Struct` for convolved fluxes

    Attributes
    ----------
    deconvKey : `lsst.afw.table.Key_Flag`
        Key to set flag indicating no measurement was made due to the need to deconvolve
    aperture : `lsst.meas.base.CircularApertureFluxAlgorithm`
        Measurement algorithm to perform aperture flux measurements
    kronKeys : `lsst.pipe.base.Struct`
        Container for Kron results or `None` if no Kron radius is available; when set,
        includes `result` (`lsst.meas.base.FluxResultKey`: keys to set results from Kron
        flux measurement) and `flag` (`lsst.afw.table.Key_Flag`: key to set failure flag
        for Kron measurement).
    """

    def __init__(self, name, schema, seeing, config, metadata):
        deconvKey = schema.addField(name + "_deconv", type="Flag",
                                    doc="deconvolution required for seeing %f; no measurement made" %
                                    (seeing,))
        aperture = lsst.meas.base.CircularApertureFluxAlgorithm(config.aperture.makeControl(), name,
                                                                schema, metadata)
        kronKeys = Struct(
            result=lsst.meas.base.FluxResultKey.addFields(schema, name + "_kron",
                                                          doc="convolved Kron flux: seeing %f" % (seeing,)),
            flag = schema.addField(name + "_kron_flag", type="Flag",
                                   doc="convolved Kron flux failed: seeing %f" % (seeing,)),
        )
        Struct.__init__(self, deconvKey=deconvKey, aperture=aperture, kronKeys=kronKeys)


class BaseConvolvedFluxConfig(Config):
    # convolution
    seeing = ListField(dtype=float, default=[3.5, 5.0, 6.5], doc="list of target seeings (FWHM, pixels)")
    kernelScale = Field(dtype=float, default=4.0, doc="scaling factor of kernel sigma for kernel size")
    # aperture flux
    aperture = ConfigField(dtype=ApertureFluxConfig, doc="Aperture photometry parameters")
    # Kron flux
    kronRadiusName = Field(dtype=str, default="ext_photometryKron_KronFlux_radius",
                           doc="name of Kron radius field in reference")
    maxSincRadius = Field(dtype=float, default=10.0,
                          doc="Largest aperture for which to use the sinc aperture code for Kron (pixels)")
    kronRadiusForFlux = Field(dtype=float, default=2.5, doc="Number of Kron radii for Kron flux")

    def setDefaults(self):
        Config.setDefaults(self)
        # Don't need the full set of apertures because the larger ones aren't affected by the convolution
        self.aperture.radii = [3.3, 4.5, 6.0]

    def getBaseNameForSeeing(self, seeing, name=PLUGIN_NAME):
        """Return base name for measurement, given seeing

        Parameters
        ----------
        seeing : `float`
            The seeing value; it is required that the `ConvolvedFluxConfig.seeing` list
            include this value.
        name : `str`, optional
            The name of the plugin.

        Returns
        -------
        baseName : `str`
            Base name for measurement with nominated seeing.
        """
        indices = [ii for ii, target in enumerate(self.seeing) if seeing == target]
        if len(indices) != 1:
            raise RuntimeError("Unable to uniquely identify index for seeing %f: %s" % (seeing, indices))
        return name + "_%d" % (indices[0],)

    def getApertureResultName(self, seeing, radius, name=PLUGIN_NAME):
        """Return name for aperture measurement result

        Parameters
        ----------
        seeing : `float`
            The seeing value; it is required that the `ConvolvedFluxConfig.seeing` list
            include this value.
        radius : `float`
            The aperture radius. If this doesn't correspond to a value in the
            `ConvolvedFluxConfig.aperture.radii` then the returned name may not be useful.
        name : `str`, optional
            The name of the plugin.

        Returns
        -------
        resultName : `str`
            Result name for aperture measurement with nominated seeing and radius.
        """
        baseName = self.getBaseNameForSeeing(seeing, name=name)
        return lsst.meas.base.CircularApertureFluxAlgorithm.makeFieldPrefix(baseName, radius)

    def getKronResultName(self, seeing, name=PLUGIN_NAME):
        """Return name for Kron measurement result

        Parameters
        ----------
        seeing : `float`
            The seeing value; it is required that the `ConvolvedFluxConfig.seeing` list
            include this value.
        name : `str`, optional
            The name of the plugin.

        Returns
        -------
        resultName : `str`
            Result name for Kron measurement with nominated seeing.
        """
        return self.getBaseNameForSeeing(seeing, name=name) + "_kron"

    def getAllApertureResultNames(self, name=PLUGIN_NAME):
        """Return all names for aperture measurements

        Parameters
        ----------
        name : `str`, optional
            The name of the plugin.

        Returns
        -------
        results : `list` of `str`
            List of names for aperture measurements (for all seeings)
        """
        return [lsst.meas.base.CircularApertureFluxAlgorithm.makeFieldPrefix(seeingName, radius) for
                seeingName in [name + "_%d" % (ii,) for ii in range(len(self.seeing))] for
                radius in self.aperture.radii]

    def getAllKronResultNames(self, name=PLUGIN_NAME):
        """Return all names for Kron measurements

        Parameters
        ----------
        name : `str`, optional
            The name of the plugin.

        Returns
        -------
        results : `list` of `str`
            List of names for Kron measurements (for all seeings)
        """
        return [name + "_%d_kron" % (ii,) for ii in range(len(self.seeing))]

    def getAllResultNames(self, name=PLUGIN_NAME):
        """Return all names for measurements

        Parameters
        ----------
        name : `str`, optional
            The name of the plugin.

        Returns
        -------
        results : `list` of `str`
            List of names for measurements (for all seeings and apertures and Kron)
        """
        return self.getAllApertureResultNames(name=name) + self.getAllKronResultNames(name=name)


class BaseConvolvedFluxPlugin(lsst.meas.base.BaseMeasurementPlugin):
    """Calculate aperture fluxes on images convolved to target seeing.

    This measurement plugin convolves the image to match a target seeing
    and measures fluxes within circular apertures and within the Kron
    aperture (defined as a multiple of the Kron radius which is already
    available in the catalog).

    Throughout, we assume a Gaussian PSF to simplify and optimise the
    convolution for speed. The results are therefore not exact, but should
    be good enough to be useful.

    The measurements produced by this plugin are useful for:
    * Fiber mags: the flux within a circular aperture in a particular seeing
      can be used to calibrate fiber-fed spectroscopic observations.
    * Galaxy photometry: the flux within an aperture in common seeing can
      be used to measure good colors for an object without assuming a model.

    The error handling is a bit different from most measurement plugins (which
    are content to fail anywhere and have the entire algorithm flagged as having
    failed), because we have multiple components (circular apertures and Kron)
    and we don't want the whole to fail because one component failed. Therefore,
    there's a few more try/except blocks than might be otherwise expected.
    """

    @classmethod
    def getExecutionOrder(cls):
        return KronFluxPlugin.getExecutionOrder() + 0.1  # Should run after Kron because we need the radius

    def __init__(self, config, name, schema, metadata):
        """Ctor

        Parameters
        ----------
        config : `ConvolvedFluxConfig`
            Configuration for plugin.
        name : `str`
            Name of plugin (used as prefix for columns in schema).
        schema : `lsst.afw.table.Schema`
            Catalog schema.
        metadata : `lsst.daf.base.PropertyList`
            Algorithm metadata to be recorded in the catalog header.
        """
        lsst.meas.base.BaseMeasurementPlugin.__init__(self, config, name)
        self.seeingKey = schema.addField(name + "_seeing", type="F",
                                         doc="original seeing (Gaussian sigma) at position",
                                         units="pixel")
        self.data = [ConvolvedFluxData(self.config.getBaseNameForSeeing(seeing), schema, seeing,
                                       self.config, metadata) for seeing in self.config.seeing]

        flagDefs = lsst.meas.base.FlagDefinitionList()
        flagDefs.addFailureFlag("error in running ConvolvedFluxPlugin")
        self.flagHandler = lsst.meas.base.FlagHandler.addFields(schema, name, flagDefs)
        # Trigger aperture corrections for all flux measurements
        for apName in self.config.getAllApertureResultNames(name):
            lsst.meas.base.addApCorrName(apName)
        for kronName in self.config.getAllKronResultNames(name):
            lsst.meas.base.addApCorrName(kronName)

        self.centroidExtractor = lsst.meas.base.SafeCentroidExtractor(schema, name)

    def measure(self, measRecord, exposure):
        """Measure source on image

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source to be measured.
        exposure : `lsst.afw.image.Exposure`
            Image to be measured.
        """
        return self.measureForced(measRecord, exposure, measRecord, None)

    def measureForced(self, measRecord, exposure, refRecord, refWcs):
        """Measure source on image in forced mode

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source to be measured.
        exposure : `lsst.afw.image.Exposure`
            Image to be measured.
        refRecord : `lsst.afw.table.SourceRecord`
            Record providing reference position and aperture.
        refWcs : `lsst.afw.image.Wcs` or `None`
            Astrometric solution for reference, or `None` for no conversion
            from reference to measurement frame.
        """
        psf = exposure.getPsf()
        if psf is None:
            raise lsst.meas.base.MeasurementError("No PSF in exposure")

        refCenter = self.centroidExtractor(refRecord, self.flagHandler)

        if refWcs is not None:
            measWcs = exposure.getWcs()
            if measWcs is None:
                raise lsst.meas.base.MeasurementError("No WCS in exposure")
            transform = lsst.afw.image.XYTransformFromWcsPair(measWcs, refWcs)
            transform = transform.linearizeForwardTransform(refCenter)
        else:
            transform = lsst.afw.geom.AffineTransform()

        kron = self.getKronAperture(refRecord, transform)

        center = refCenter if transform is None else transform(refCenter)
        seeing = psf.computeShape(center).getDeterminantRadius()
        measRecord.set(self.seeingKey, seeing)

        maxRadius = self.getMaxRadius(kron)
        for ii, target in enumerate(self.config.seeing):
            try:
                convolved = self.convolve(exposure, seeing, target/SIGMA_TO_FWHM, measRecord.getFootprint(),
                                          maxRadius)
            except (DeconvolutionError, RuntimeError):
                # Record the problem, but allow the measurement to run in case it's useful
                measRecord.set(self.data[ii].deconvKey, True)
                convolved = exposure
            self.measureAperture(measRecord, convolved, self.data[ii].aperture)
            if kron is not None:
                self.measureForcedKron(measRecord, self.data[ii].kronKeys, convolved.getMaskedImage(), kron)

    def fail(self, measRecord, error=None):
        """Record failure

        Called by the measurement framework when it catches an exception.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source on which measurement failed.
        error : `Exception`, optional
            Error that occurred, or None.
        """
        self.flagHandler.handleFailure(measRecord)

    def getKronAperture(self, refRecord, transform):
        """Determine the Kron radius

        Because we need to know the size of the area beforehand (we don't want to convolve
        the entire image just for this source), we are not measuring an independent Kron
        radius, but using the Kron radius that's already provided in the `refRecord` as
        `ConvolvedFluxConfig.kronRadiusName`.

        Parameters
        ----------
        refRecord : `lsst.afw.table.SourceRecord`
            Record for source defining Kron aperture.
        transform : `lsst.afw.geom.AffineTransform`
            Transformation to apply to reference aperture.

        Returns
        -------
        aperture : `lsst.meas.extensions.photometryKron.KronAperture`
            Kron aperture.
        """
        try:
            radius = refRecord.get(self.config.kronRadiusName)
        except:
            return None
        if not np.isfinite(radius):
            return None
        return KronAperture(refRecord, transform, radius)

    def getMaxRadius(self, kron):
        """Determine the maximum radius we care about

        Because we don't want to convolve the entire image just for this source,
        we determine the maximum radius we care about for this source and will
        convolve only that.

        Parameters
        ----------
        kron : `lsst.meas.extensions.photometryKron.KronAperture` or `None`
            Kron aperture, or `None` if unavailable.

        Returns
        -------
        maxRadius : `int`
            Maximum radius of interest.
        """
        kronRadius = kron.getAxes().getDeterminantRadius() if kron is not None else 0.0
        return int(max(max(self.config.aperture.radii), self.config.kronRadiusForFlux*kronRadius) + 0.5)

    def convolve(self, exposure, seeing, target, footprint, maxRadius):
        """Convolve image around source to target seeing

        We also record the original seeing at the source position.

        Because we don't want to convolve the entire image just for this source,
        we cut out an area corresponding to the source's footprint, grown by the
        radius provided by `maxRadius`.

        We assume a Gaussian PSF to simplify and speed the convolution.
        The `seeing` and `target` may be either Gaussian sigma or FWHM, so long
        as they are the same.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Image to convolve.
        seeing : `float`
            Current seeing, pixels.
        target : `float`
            Desired target seeing, pixels.
        footprint : `lsst.afw.detection.Footprint`
            Footprint for source.
        maxRadius : `int`
            Maximum radius required by measurement algorithms.

        Returns
        -------
        convExp : `lsst.afw.image.Exposure`
            Sub-image containing the source, convolved to the target seeing.

        Raises
        ------
        DeconvolutionError
            If the target seeing requires deconvolution.
        RuntimeError
            If the bounding box is too small after clipping.
        """

        if target < seeing:
            raise DeconvolutionError("Target seeing requires deconvolution")
        kernelSigma = math.sqrt(target*target - seeing*seeing)
        kernelRadius = int(self.config.kernelScale*kernelSigma + 0.5)
        kernelWidth = 2*kernelRadius + 1
        gauss = lsst.afw.math.GaussianFunction1D(kernelSigma)
        kernel = lsst.afw.math.SeparableKernel(kernelWidth, kernelWidth, gauss, gauss)

        bbox = footprint.getBBox()
        bbox.grow(kernelRadius + maxRadius)  # add an extra buffer?
        bbox.clip(exposure.getBBox())
        if bbox.getWidth() < kernelWidth or bbox.getHeight() < kernelWidth:
            raise RuntimeError("Bounding box is too small following clipping")

        image = exposure.getMaskedImage()
        subImage = image.Factory(image, bbox)
        convolved = image.Factory(bbox)
        lsst.afw.math.convolve(convolved, subImage, kernel, lsst.afw.math.ConvolutionControl(True, True))

        # This is ugly, but necessary; should be resolved following RFC-217, DM-5503
        convExp = lsst.afw.image.makeExposure(convolved)
        convInfo = convExp.getInfo()
        origInfo = exposure.getInfo()
        for method in dir(origInfo):
            if not method.startswith("get"):
                continue
            setter = "s" + method[1:]
            if not hasattr(convInfo, setter):
                continue
            getattr(convInfo, setter)(getattr(origInfo, method)())

        return convExp

    def measureAperture(self, measRecord, exposure, aperturePhot):
        """Perform aperture photometry

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source to be measured.
        exposure : `lsst.afw.image.Exposure`
            Image to be measured.
        aperturePhot : `lsst.meas.base.CircularApertureFluxAlgorithm`
            Measurement plugin that will do the measurement.
        """
        try:
            aperturePhot.measure(measRecord, exposure)
        except:
            aperturePhot.fail(measRecord)

    def measureForcedKron(self, measRecord, keys, image, aperture):
        """Measure forced Kron

        Because we need to know the size of the area beforehand (we don't want to convolve
        the entire image just for this source), we are doing forced measurement using the
        Kron radius previously determined.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source to be measured.
        keys : `lsst.pipe.base.Struct`
            Struct containing `result` (`lsst.meas.base.FluxResult`) and
            `flag` (`lsst.afw.table.Key_Flag`); provided by
            `ConvolvedFluxData.kronKeys`.
        image : `lsst.afw.image.MaskedImage`
            Image to be measured.
        aperture : `lsst.meas.extensions.photometryKron.KronAperture`
            Kron aperture to measure.
        """
        measRecord.set(keys.flag, True)  # failed unless we survive to switch this back
        if aperture is None:
            return  # We've already flagged it, so just bail
        try:
            flux = aperture.measureFlux(image, self.config.kronRadiusForFlux, self.config.maxSincRadius)
        except:
            return  # We've already flagged it, so just bail
        measRecord.set(keys.result.getFlux(), flux[0])
        measRecord.set(keys.result.getFluxSigma(), flux[1])
        measRecord.setFlag(keys.flag, bool(np.any(~np.isfinite(flux))))


def wrapPlugin(Base, PluginClass=BaseConvolvedFluxPlugin, ConfigClass=BaseConvolvedFluxConfig,
               name=PLUGIN_NAME, factory=BaseConvolvedFluxPlugin):
    """Wrap plugin for use

    A plugin has to inherit from a specific base class in order to be used
    in a particular context (e.g., single frame vs forced measurement).

    Parameters
    ----------
    Base : `type`
        Base class to give the plugin.
    PluginClass : `type`
        Plugin class to wrap.
    ConfigClass : `type`
        Configuration class; should subclass `lsst.pex.config.Config`.
    name : `str`
        Name of plugin.
    factory : callable
        Callable to create an instance of the `PluginClass`.

    Returns
    -------
    WrappedPlugin : `type`
        The wrapped plugin class (subclass of `Base`).
    WrappedConfig : `type`
        The wrapped plugin configuration (subclass of `Base.ConfigClass`).
    """
    WrappedConfig = type("ConvolvedFlux" + Base.ConfigClass.__name__, (Base.ConfigClass, ConfigClass), {})
    typeDict = dict(AlgClass=PluginClass, ConfigClass=WrappedConfig, factory=factory,
                    getExecutionOrder=PluginClass.getExecutionOrder)
    WrappedPlugin = type("ConvolvedFlux" + Base.__name__, (Base,), typeDict)
    Base.registry.register(name, WrappedPlugin)
    return WrappedPlugin, WrappedConfig

def wrapPluginForced(Base, PluginClass=BaseConvolvedFluxPlugin, ConfigClass=BaseConvolvedFluxConfig,
                     name=PLUGIN_NAME, factory=BaseConvolvedFluxPlugin):
    """Wrap plugin for use in forced measurement

    A version of `wrapPlugin` that generates a `factory` suitable for
    forced measurement. This is important because the required signature
    for the factory in forced measurement includes a 'schemaMapper' instead
    of a 'schema'.

    Parameters
    ----------
    Base : `type`
        Base class to give the plugin.
    PluginClass : `type`
        Plugin class to wrap.
    ConfigClass : `type`
        Configuration class; should subclass `lsst.pex.config.Config`.
    name : `str`
        Name of plugin.
    factory : callable
        Callable to create an instance of the `PluginClass`.

    Returns
    -------
    WrappedPlugin : `type`
        The wrapped plugin class (subclass of `Base`).
    WrappedConfig : `type`
        The wrapped plugin configuration (subclass of `Base.ConfigClass`).
    """

    def forcedPluginFactory(name, config, schemaMapper, metadata):
        return factory(name, config, schemaMapper.editOutputSchema(), metadata)
    return wrapPlugin(Base, PluginClass=PluginClass, ConfigClass=ConfigClass, name=name,
                      factory=staticmethod(forcedPluginFactory))

SingleFrameConvolvedFluxPlugin, SingleFrameConvolvedFluxConfig = wrapPlugin(WrappedSingleFramePlugin)
ForcedConvolvedFluxPlugin, ForcedConvolvedFluxConfig = wrapPluginForced(WrappedForcedPlugin)

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

from lsst.pex.config import Field, ListField, ConfigField, makeConfigClass
from lsst.pipe.base import Struct
from lsst.meas.extensions.photometryKron import KronAperture, KronFluxPlugin
from lsst.meas.base.flagDecorator import addFlagHandler

import lsst.meas.base
import lsst.afw.math
import lsst.afw.geom
import lsst.afw.image

__all__ = ("ConvolvedFluxConfig", "ConvolvedFluxPlugin")


SIGMA_TO_FWHM = 2.0*math.sqrt(2.0*(math.log(2.0)))  # Multiply sigma by this to get FWHM
PLUGIN_NAME = "ext_convolved_ConvolvedFlux"  # Usual name for plugin


class DeconvolutionError(RuntimeError):
    """Convolving to the target seeing would require deconvolution"""
    pass


class NoKronError(RuntimeError):
    """There's no Kron radius available"""
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

    def __init__(self, name, schema, seeing, config, metadata, doKron):
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
        ) if doKron else None
        Struct.__init__(self, deconvKey=deconvKey, aperture=aperture, kronKeys=kronKeys)


class ConvolvedFluxConfig(lsst.meas.base.SingleFramePluginConfig):
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
        lsst.meas.base.SingleFramePluginConfig.setDefaults(self)
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


@lsst.meas.base.register(PLUGIN_NAME)
@addFlagHandler(("flag", "error in running ConvolvedFluxPlugin"),)
class ConvolvedFluxPlugin(lsst.meas.base.SingleFramePlugin):
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

    ConfigClass = ConvolvedFluxConfig

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
        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.seeingKey = schema.addField(name + "_seeing", type="F",
                                         doc="original seeing (Gaussian sigma) at position",
                                         units="pixel")
        try:
            self.kronRadiusKey = schema.find(config.kronRadiusName).key
        except:
            self.haveKron = False
        else:
            self.haveKron = True

        self.data = [ConvolvedFluxData(self.config.getBaseNameForSeeing(seeing), schema, seeing,
                                       self.config, metadata, self.haveKron) for seeing in self.config.seeing]

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
        psf = exposure.getPsf()
        if psf is None:
            raise lsst.meas.base.MeasurementError("No PSF in exposure")
        center = self.centroidExtractor(measRecord, self.flagHandler)
        seeing = psf.computeShape(center).getDeterminantRadius()
        measRecord.set(self.seeingKey, seeing)

        maxRadius = self.getMaxRadius(measRecord)
        for ii, target in enumerate(self.config.seeing):
            try:
                convolved = self.convolve(exposure, seeing, target/SIGMA_TO_FWHM, measRecord.getFootprint(),
                                          maxRadius)
            except (DeconvolutionError, RuntimeError):
                measRecord.set(self.data[ii].deconvKey, True)
                continue
            self.measureAperture(measRecord, convolved, self.data[ii].aperture)
            if self.haveKron:
                self.measureForcedKron(measRecord, self.data[ii].kronKeys, convolved.getMaskedImage(), center)

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

    def getKronRadius(self, measRecord):
        """Determine the Kron radius

        Because we need to know the size of the area beforehand (we don't want to convolve
        the entire image just for this source), we are not measuring an independent Kron
        radius, but using the Kron radius that's already provided in the `measRecord` as
        `ConvolvedFluxConfig.kronRadiusName`.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source to be measured.

        Returns
        -------
        radius : `float`
            Kron radius.

        Raises
        ------
        `NoKronError`
            If the Kron radius could not be found or is non-finite.
        """
        try:
            radius = measRecord.get(self.kronRadiusKey)
        except:
            # recast the exception to something we can expect
            raise NoKronError("Unable to find Kron radius")
        if not np.isfinite(radius):
            raise NoKronError("Bad Kron radius")
        return radius

    def getMaxRadius(self, measRecord):
        """Determine the maximum radius we care about

        Because we don't want to convolve the entire image just for this source,
        we determine the maximum radius we care about for this source and will
        convolve only that.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source to measured.

        Returns
        -------
        maxRadius : `int`
            Maximum radius of interest.
        """
        try:
            kronRadius = self.getKronRadius(measRecord)
        except NoKronError:
            kronRadius = 0.0
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
        convolved : `lsst.afw.image.Exposure`
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

    def measureForcedKron(self, measRecord, keys, image, center):
        """Measure forced Kron

        Because we need to know the size of the area beforehand (we don't want to convolve
        the entire image just for this source), we are doing forced measurement using the
        Kron radius previously determined.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record for source to be measured.
        exposure : `lsst.afw.image.MaskedImage`
            Image to be measured.
        keys : `lsst.pipe.base.Struct`
            Struct containing `result` (`lsst.meas.base.FluxResult`) and
            `flag` (`lsst.afw.table.Key_Flag`); provided by
            `ConvolvedFluxData.kronKeys`.
        center : `lsst.afw.geom.Point2D`
            Center for Kron aperture.
        """
        measRecord.set(keys.flag, True)  # failed unless we survive to switch this back
        try:
            radius = self.getKronRadius(measRecord)
        except NoKronError:
            return  # We've already flagged it, so just bail
        aperture = KronAperture(measRecord, lsst.afw.geom.AffineTransform(), radius)
        try:
            flux = aperture.measureFlux(image, self.config.kronRadiusForFlux, self.config.maxSincRadius)
        except:
            return  # We've already flagged it, so just bail
        measRecord.set(keys.result.getFlux(), flux[0])
        measRecord.set(keys.result.getFluxSigma(), flux[1])
        measRecord.setFlag(keys.flag, bool(np.any(~np.isfinite(flux))))

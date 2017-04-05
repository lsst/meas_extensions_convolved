#
# LSST Data Management System
# Copyright 2017 LSST/AURA.
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
from __future__ import absolute_import, division, print_function

import math
import unittest
import lsst.utils.tests
import lsst.daf.base as dafBase
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
import lsst.meas.extensions.convolved  # Load flux.convolved algorithm

import lsst.afw.display as afwDisplay

try:
    type(display)
except NameError:
    display = False
    frame = 1

SIGMA_TO_FWHM = 2.0*math.sqrt(2.0*math.log(2.0))


def makeExposure(bbox, scale, psfFwhm, flux):
    """Make a fake exposure

    Parameters
    ----------

    bbox : `lsst.afw.geom.Box2I`
        Bounding box for image.
    scale : `lsst.afw.geom.Angle`
        Pixel scale.
    psfFwhm : `float`
        PSF FWHM (arcseconds)
    flux : `float`
        PSF flux (ADU)

    Returns
    -------
    exposure : `lsst.afw.image.ExposureF`
        Fake exposure.
    center : `lsst.afw.geom.Point2D`
        Position of fake source.
    """
    image = afwImage.ImageF(bbox)
    image.set(0)
    center = afwGeom.Box2D(bbox).getCenter()
    psfSigma = psfFwhm/SIGMA_TO_FWHM/scale.asArcseconds()
    psfWidth = 2*int(4.0*psfSigma) + 1
    psf = afwDetection.GaussianPsf(psfWidth, psfWidth, psfSigma)
    psfImage = psf.computeImage(center).convertF()
    psfFlux = psfImage.getArray().sum()
    psfImage *= flux/psfFlux

    subImage = afwImage.ImageF(image, psfImage.getBBox(afwImage.PARENT), afwImage.PARENT)
    subImage += psfImage

    exp = afwImage.makeExposure(afwImage.makeMaskedImage(image))
    exp.setPsf(psf)
    exp.getMaskedImage().getVariance().set(1.0)
    exp.getMaskedImage().getMask().set(0)

    exp.setWcs(afwImage.makeWcs(afwCoord.Coord(0.0*afwGeom.degrees, 0.0*afwGeom.degrees),
                                center, scale.asDegrees(), 0.0, 0.0, scale.asDegrees()))
    return exp, center


class ConvolvedFluxTestCase(lsst.utils.tests.TestCase):
    """A test case for measuring convolved fluxes"""

    def checkSchema(self, schema, names):
        """Check that the schema includes flux, fluxSigma and flag elements for each measurement

        Parameters
        ----------
        schema : `lsst.afw.table.Schema`
            Schema to check.
        names : `list` of `str`
            List of measurement algorithm names
        """
        for name in names:
            self.assertIn(name + "_flux", schema)
            self.assertIn(name + "_fluxSigma", schema)
            self.assertIn(name + "_flag", schema)

    def check(self, psfFwhm=0.5, flux=1000.0, forced=False):
        """Check that we can measure convolved fluxes

        We create an image with a Gaussian PSF and a single point source.
        Measurements of the point source should match expectations for a
        Gaussian of the known sigma and known aperture radius.

        Parameters
        ----------
        psfFwhm : `float`
            PSF FWHM (arcsec)
        flux : `float`
            Source flux (ADU)
        forced : `bool`
            Forced measurement?
        """
        bbox = afwGeom.Box2I(afwGeom.Point2I(12345, 6789), afwGeom.Extent2I(200, 300))

        # We'll only achieve the target accuracy if the pixel scale is rather smaller than Gaussians
        # involved. Otherwise it's important to consider the convolution with the pixel grid, and we're
        # not doing that here.
        scale = 0.1*afwGeom.arcseconds

        TaskClass = measBase.ForcedMeasurementTask if forced else measBase.SingleFrameMeasurementTask

        exposure, center = makeExposure(bbox, scale, psfFwhm, flux)
        measConfig = TaskClass.ConfigClass()
        algName = "ext_convolved_ConvolvedFlux"
        measConfig.plugins.names.add(algName)
        if not forced:
            measConfig.plugins.names.add("ext_photometryKron_KronFlux")
        else:
            measConfig.copyColumns = {"id": "objectId", "parent": "parentObjectId"}
        values = [ii/scale.asArcseconds() for ii in (0.6, 0.8, 1.0, 1.2)]
        algConfig = measConfig.plugins[algName]
        algConfig.seeing = values
        algConfig.aperture.radii = values
        algConfig.aperture.maxSincRadius = max(values) + 1  # Get as exact as we can

        if forced:
            offset = lsst.afw.geom.Extent2D(-12.3, 45.6)
            kronRadiusName = "my_Kron_Radius"
            kronRadius = 12.345
            refWcs = exposure.getWcs().clone()
            refWcs.shiftReferencePixel(offset)
            measConfig.plugins[algName].kronRadiusName = kronRadiusName
            refSchema = afwTable.SourceTable.makeMinimalSchema()
            centroidKey = afwTable.Point2DKey.addFields(refSchema, "my_centroid", doc="centroid",
                                                        unit="pixel")
            shapeKey = afwTable.QuadrupoleKey.addFields(refSchema, "my_shape", "shape")
            refSchema.getAliasMap().set("slot_Centroid", "my_centroid")
            refSchema.getAliasMap().set("slot_Shape", "my_shape")
            refSchema.addField("my_centroid_flag", type="Flag", doc="centroid flag")
            refSchema.addField("my_shape_flag", type="Flag", doc="shape flag")
            refSchema.addField(kronRadiusName, type=float, doc="my custom kron radius", units="pixel")
            refCat = afwTable.SourceCatalog(refSchema)
            refSource = refCat.addNew()
            refSource.set(centroidKey, center + offset)
            refSource.set(shapeKey, afwEll.Quadrupole(afwEll.Axes(kronRadius, kronRadius, 0)))
            refSource.set(kronRadiusName, kronRadius)
            refSource.setCoord(refWcs.pixelToSky(refSource.get(centroidKey)))
            taskInitArgs = (refSchema,)
            taskRunArgs = (refCat, refWcs)
        else:
            taskInitArgs = (afwTable.SourceTable.makeMinimalSchema(),)
            taskRunArgs = ()

        algMetadata = dafBase.PropertyList()
        task = TaskClass(*taskInitArgs, config=measConfig, algMetadata=algMetadata)

        schema = task.schema
        measCat = afwTable.SourceCatalog(schema)
        source = measCat.addNew()
        source.getTable().setMetadata(algMetadata)
        ss = afwDetection.FootprintSet(exposure.getMaskedImage(), afwDetection.Threshold(0.1))
        fp = ss.getFootprints()[0]
        source.setFootprint(fp)

        task.run(measCat, exposure, *taskRunArgs)

        disp = afwDisplay.Display(frame)
        disp.mtv(exposure)
        disp.dot("x", *center, origin=afwImage.PARENT, title="psfFwhm=%f" % (psfFwhm,))

        self.checkSchema(schema, algConfig.getAllApertureResultNames())
        self.checkSchema(schema, algConfig.getAllKronResultNames())
        self.checkSchema(schema, algConfig.getAllResultNames())

        if not forced:
            kronRadius = source.get("ext_photometryKron_KronFlux_radius")

        self.assertFalse(source.get(algName + "_flag"))  # algorithm succeeded
        originalSeeing = psfFwhm/scale.asArcseconds()
        for ii, targetSeeing in enumerate(algConfig.seeing):
            deconvolve = targetSeeing < originalSeeing
            self.assertEqual(source.get(algName + "_%d_deconv" % ii), deconvolve)
            seeing = originalSeeing if deconvolve else targetSeeing

            def expected(radius, sigma=seeing/SIGMA_TO_FWHM):
                """Return expected flux for 2D Gaussian with nominated sigma"""
                return flux*(1.0 - math.exp(-0.5*(radius/sigma)**2))

            # Kron succeeded and match expectation
            if not forced:
                kronName = algConfig.getKronResultName(targetSeeing)
                kronApRadius = algConfig.kronRadiusForFlux*kronRadius
                self.assertFloatsAlmostEqual(source.get(kronName + "_flux"),
                                             expected(kronApRadius), rtol=1.0e-3)
                self.assertGreater(source.get(kronName + "_fluxSigma"), 0)
                self.assertFalse(source.get(kronName + "_flag"))

            # Aperture measurements suceeded and match expectation
            for jj, radius in enumerate(measConfig.algorithms[algName].aperture.radii):
                name = algConfig.getApertureResultName(targetSeeing, radius)
                self.assertFloatsAlmostEqual(source.get(name + "_flux"), expected(radius), rtol=1.0e-3)
                self.assertFalse(source.get(name + "_flag"))
                self.assertGreater(source.get(name + "_fluxSigma"), 0)

    def testConvolvedFlux(self):
        for forced in (True, False):
            for psfFwhm in (0.5,  # Smaller than all target seeings
                            0.9,  # Larger than half the target seeings
                            1.3,  # Larger than all the target seeings
                            ):
                self.check(psfFwhm=psfFwhm, forced=forced)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module, backend="virtualDevice"):
    lsst.utils.tests.init()
    try:
        afwDisplay.setDefaultBackend(backend)
    except:
        print("Unable to configure display backend: %s" % backend)


if __name__ == "__main__":
    import sys

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--backend', type=str, default="virtualDevice",
                        help="The backend to use, e.g. 'ds9'. Be sure to 'setup display_<backend>'")
    args = parser.parse_args()

    setup_module(sys.modules[__name__], backend=args.backend)
    unittest.main()

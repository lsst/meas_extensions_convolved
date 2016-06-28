#!/usr/bin/env python

import math
import unittest
import lsst.utils.tests as tests
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.extensions.convolved  # Load flux.convolved algorithm

try:
    type(display)
except NameError:
    display = False
    frame = 1

import lsst.afw.display.ds9 as ds9

SIGMA_TO_FWHM = 2.0*math.sqrt(2.0*math.log(2.0))


def makeExposure(bbox, scale, psfFwhm, flux):
    """Make a fake exposure

    @param bbox: Bounding box for image (Box2I)
    @param scale: Pixel scale (Angle)
    @param psfFwhm: PSF FWHM (arcseconds)
    @param flux: PSF flux (ADU)
    @return Exposure, source center
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


class ConvolvedFluxTestCase(tests.TestCase):
    """A test case for measuring convolved fluxes"""

    def check(self, psfFwhm=0.5, flux=1000.0):
        """Check that we can measure convolved fluxes
    
        We create an image with a Gaussian PSF and a single point source.
        Measurements of the point source should match expectations for a
        Gaussian of the known sigma and known aperture radius.

        @param psfFwhm: PSF FWHM in arcsec
        @param flux: source flux in ADU
        """
        bbox = afwGeom.Box2I(afwGeom.Point2I(12345, 6789), afwGeom.Extent2I(200, 300))

        # We'll only achieve the target accuracy if the pixel scale is rather smaller than Gaussians
        # involved. Otherwise it's important to consider the convolution with the pixel grid, and we're
        # not doing that here.
        scale = 0.1*afwGeom.arcseconds

        exposure, center = makeExposure(bbox, scale, psfFwhm, flux)
        msConfig = measAlg.SourceMeasurementConfig()
        msConfig.algorithms.names.add("flux.convolved")
        values = [ii/scale.asArcseconds() for ii in (0.6, 0.8, 1.0, 1.2)]
        msConfig.algorithms["flux.convolved"].seeing = values
        msConfig.algorithms["flux.convolved"].radius = values

        schema = afwTable.SourceTable.makeMinimalSchema()
        ms = msConfig.makeMeasureSources(schema)

        table = afwTable.SourceTable.make(schema)
        msConfig.slots.setupTable(table)
        source = table.makeRecord()

        ss = afwDetection.FootprintSet(exposure.getMaskedImage(), afwDetection.Threshold(0.1))
        fp = ss.getFootprints()[0]
        source.setFootprint(fp)

        ms.apply(source, exposure, center)

        if display:
            ds9.mtv(exposure, frame=frame)
            ds9.dot("x", center.getX() - exposure.getX0(), center.getY() - exposure.getY0(), frame=frame)
            import pdb;pdb.set_trace()

        self.assertFalse(source.get("flux.convolved.flag"))  # algorithm succeeded
        for ii, seeing in enumerate(msConfig.algorithms["flux.convolved"].seeing):
            deconvolve = seeing < psfFwhm/scale.asArcseconds()
            self.assertTrue(source.get("flux.convolved.%d.deconv" % ii) == deconvolve)
            if deconvolve:
                # Not worth checking anything else
                continue
            for jj, radius in enumerate(msConfig.algorithms["flux.convolved"].radius):
                sigma = seeing/SIGMA_TO_FWHM
                expected = flux*(1.0 - math.exp(-0.5*(radius/sigma)**2))
                name = "flux.convolved.%d.%d" % (ii, jj)
                self.assertClose(expected, source.get(name), rtol=1.0e-3)
                self.assertFalse(source.get(name + ".flags"))
                self.assertGreater(source.get(name + ".err"), 0)

    def testConvolvedFlux(self):
        self.check(psfFwhm=0.5) # Smaller than all target seeings
        self.check(psfFwhm=0.9) # Larger than half the target seeings
        self.check(psfFwhm=1.3) # Larger than all the target seeings


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(ConvolvedFluxTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit = False):
    """Run the tests"""
    tests.run(suite(), exit)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--display", default=False, action="store_true", help="Activate display?")
    args = parser.parse_args()
    display = args.display
    run(True)

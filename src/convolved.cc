// -*- LSST-C++ -*-

#include <cmath>
#include <algorithm>

#include "boost/make_shared.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/table.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/math/ConvolveImage.h"
#include "lsst/meas/algorithms/Algorithm.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/extensions/photometryKron.h"

#include "lsst/meas/extensions/convolved/convolved.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace convolved {
namespace {

double const SIGMA_TO_FWHM = 2.0*std::sqrt(2.0*(std::log(2.0)));

LSST_EXCEPTION_TYPE(DeconvolutionException, pex::exceptions::RuntimeErrorException,
                    lsst::meas::extensions::convolved::DeconvolutionException);
LSST_EXCEPTION_TYPE(NoKronException, pex::exceptions::RuntimeErrorException,
                    lsst::meas::extensions::convolved::NoKronException);

class ConvolvedFlux : public algorithms::Algorithm {
public:

    ConvolvedFlux(ConvolvedFluxControl const & ctrl, afw::table::Schema & schema) :
        algorithms::Algorithm(ctrl),
        _failKey(schema.addField<afw::table::Flag>(ctrl.name + ".flag", "convolved algorithm failed")),
        _seeingKey(schema.addField<double>(ctrl.name + ".seeing", "original seeing at position"))
    {
        _deconvKeys.reserve(ctrl.seeing.size());
        _apFluxKeys.reserve(ctrl.seeing.size());
        _kronKeys.reserve(ctrl.seeing.size());
        for (std::size_t ii = 0; ii < ctrl.seeing.size(); ++ii) {
            _deconvKeys.push_back(
                schema.addField<afw::table::Flag>(
                    (boost::format("%s.%d.deconv") % ctrl.name % ii).str(),
                    (boost::format("deconvolution required for seeing %f; no measurement made") %
                        ctrl.seeing[ii]).str()
                )
            );
            _kronKeys.push_back(
                afw::table::addFluxFields(
                    schema,
                    (boost::format("%s.%d.kron") % ctrl.name % ii).str(),
                    "Forced Kron flux"
                )
            );
            std::vector<afw::table::KeyTuple<afw::table::Flux> > apFluxKeys;
            apFluxKeys.reserve(ctrl.radius.size());
            for (std::size_t jj = 0; jj < ctrl.radius.size(); ++jj) {
                apFluxKeys.push_back(
                    afw::table::addFluxFields(
                        schema,
                        (boost::format("%s.%d.%d") % ctrl.name % ii % jj).str(),
                        (boost::format("convolved flux: seeing %f radius %f") %
                         ctrl.seeing[ii] % ctrl.radius[jj]).str()
                    )
                );
            }
            _apFluxKeys.push_back(apFluxKeys);
        }
    }

private:
    std::vector<std::vector<afw::table::KeyTuple<afw::table::Flux> > > _apFluxKeys; // [seeing][radius]
    std::vector<afw::table::Key<afw::table::Flag> > _deconvKeys; // [seeing]
    afw::table::Key<afw::table::Flag> _failKey;
    std::vector<afw::table::KeyTuple<afw::table::Flux> > _kronKeys; // [seeing]
    afw::table::Key<double> _seeingKey;

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    template <typename PixelT>
    void _applyForced(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas
    ) const;

    template <typename PixelT>
    afw::image::MaskedImage<PixelT> _convolve(
        double seeingFwhm,
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center,
        double maxRadius
    ) const;

    template <typename PixelT>
    void _applyAperture(
        afw::table::SourceRecord & source,
        std::vector<afw::table::KeyTuple<afw::table::Flux> > const & keys,
        afw::image::MaskedImage<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    template <typename PixelT>
    void _applyForcedKron(
        afw::table::SourceRecord & source,
        afw::table::KeyTuple<afw::table::Flux> const & keys,
        afw::image::MaskedImage<PixelT> const & image,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas=afw::geom::AffineTransform()
    ) const;

    float _getKronRadius(
        afw::table::SourceRecord const & reference
    ) const;

    double _getMaxRadius(
        afw::table::SourceRecord const & reference // reference source, for Kron radius
    ) const;

    void _setFailFlags(
        afw::table::SourceRecord & source
    ) const;

    ConvolvedFluxControl const& getControl() const {
        return static_cast<ConvolvedFluxControl const &>(algorithms::Algorithm::getControl());
    }

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(ConvolvedFlux);
};


/************************************************************************************************************/

float ConvolvedFlux::_getKronRadius(afw::table::SourceRecord const & reference) const
{
    ConvolvedFluxControl const & ctrl = getControl();
    double radius;
    try {
        radius = reference.get(reference.getSchema().find<float>(ctrl.kronName + ".radius").key);
    } catch (...) {
        // recast the exception to something we can expect
        throw LSST_EXCEPT(NoKronException, "Unable to find Kron radius");
    }
    if (!utils::isfinite(radius)) {
        throw LSST_EXCEPT(NoKronException, "Bad Kron radius");
    }
    return radius;
}

double ConvolvedFlux::_getMaxRadius(afw::table::SourceRecord const & reference) const
{
    float kronRadius = 0.0;
    try {
        kronRadius = _getKronRadius(reference);
    } catch (NoKronException const&) {} // don't care
    ConvolvedFluxControl const & ctrl = getControl();
    double const maxAperture = *std::max_element(ctrl.radius.begin(), ctrl.radius.end());
    return std::max(ctrl.kronRadiusForFlux*kronRadius, maxAperture);
}

void ConvolvedFlux::_setFailFlags(afw::table::SourceRecord & source) const
{
    source.set(_failKey, true);
    ConvolvedFluxControl const & ctrl = getControl();
    for (std::size_t ii = 0; ii < ctrl.seeing.size(); ++ii) {
        for (std::size_t jj = 0; jj < ctrl.radius.size(); ++jj) {
            source.set(_apFluxKeys[ii][jj].flag, true);
        }
        source.set(_kronKeys[ii].flag, true);
    }
}

template <typename PixelT>
afw::image::MaskedImage<PixelT> ConvolvedFlux::_convolve(
    double seeingFwhm,
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center,
    double maxRadius
    ) const
{
    if (!exposure.getPsf()) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "No PSF in exposure");
    }
    // Check whether our target seeing requires deconvolution
    ConvolvedFluxControl const& ctrl = getControl();
    double const seeing = exposure.getPsf()->computeShape(center).getDeterminantRadius();
    source.set(_seeingKey, seeing);
    double const target = seeingFwhm/SIGMA_TO_FWHM;
    if (target < seeing) {
        throw LSST_EXCEPT(DeconvolutionException, "Target seeing requires deconvolution");
    }

    // Convolve to our target seeing
    double const kernelSigma = std::sqrt(target*target - seeing*seeing);
    std::size_t const kernelRadius = ctrl.kernelScale*kernelSigma;
    std::size_t const kernelWidth = 2.0*kernelRadius + 1;

    afw::math::GaussianFunction1<double> const gauss(kernelSigma);
    afw::math::SeparableKernel const kernel(kernelWidth, kernelWidth, gauss, gauss);

    afw::geom::Box2I bbox = source.getFootprint()->getBBox();
    bbox.grow(kernelRadius + maxRadius); // XXX extra buffer?
    bbox.clip(exposure.getBBox(afw::image::PARENT));

    afw::image::MaskedImage<PixelT> subImage(exposure.getMaskedImage(), bbox, afw::image::PARENT);
    afw::image::MaskedImage<PixelT> convolved(subImage, true);
    afw::math::convolve(convolved, subImage, kernel, afw::math::ConvolutionControl(true, true));
    return convolved;
}

template <typename PixelT>
void ConvolvedFlux::_applyAperture(
    afw::table::SourceRecord & source,
    std::vector<afw::table::KeyTuple<afw::table::Flux> > const & keys,
    afw::image::MaskedImage<PixelT> const & image,
    afw::geom::Point2D const & center
    ) const
{
    ConvolvedFluxControl const& ctrl = getControl();
    for (std::size_t ii = 0; ii < ctrl.radius.size(); ++ii) {
        double const radius = ctrl.radius[ii];
        afw::geom::ellipses::Ellipse aperture(afw::geom::ellipses::Axes(radius, radius), center);
        std::pair<double, double> flux;
        try {
            flux = algorithms::photometry::calculateSincApertureFlux(image, aperture);
        } catch (pex::exceptions::Exception const&) {
            continue;
        }
        source.set(keys[ii].meas, flux.first);
        source.set(keys[ii].err, flux.second);
        source.set(keys[ii].flag, utils::isnan(flux.first) || utils::isnan(flux.second));
    }
}

template <typename PixelT>
void ConvolvedFlux::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    _setFailFlags(source);
    double const maxRadius = _getMaxRadius(source);
    ConvolvedFluxControl const & ctrl = getControl();
    for (std::size_t ii = 0; ii < ctrl.seeing.size(); ++ii) {
        afw::image::MaskedImage<PixelT> convolved;
        try {
            convolved = _convolve(ctrl.seeing[ii], source, exposure, center, maxRadius);
        } catch (DeconvolutionException const&) {
            source.set(_deconvKeys[ii], true);
            continue;
        }
        _applyAperture(source, _apFluxKeys[ii], convolved, center);
        _applyForcedKron(source, _kronKeys[ii], convolved, center, source);
    }
    source.set(_failKey, false);
}

template <typename PixelT>
void ConvolvedFlux::_applyForced(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center,
    afw::table::SourceRecord const & reference,
    afw::geom::AffineTransform const & refToMeas
    ) const
{
    _setFailFlags(source);
    double const maxRadius = _getMaxRadius(reference);
    ConvolvedFluxControl const & ctrl = getControl();
    for (std::size_t ii = 0; ii < ctrl.seeing.size(); ++ii) {
        afw::image::MaskedImage<PixelT> convolved;
        try {
            convolved = _convolve(ctrl.seeing[ii], source, exposure, refToMeas(center), maxRadius);
        } catch (DeconvolutionException const&) {
            source.set(_deconvKeys[ii], true);
            continue;
        }
        _applyAperture(source, _apFluxKeys[ii], convolved, refToMeas(center));
        _applyForcedKron(source, _kronKeys[ii], convolved, refToMeas(center), reference);
    }
    source.set(_failKey, false);
}

template <typename PixelT>
void ConvolvedFlux::_applyForcedKron(
        afw::table::SourceRecord & source,
        afw::table::KeyTuple<afw::table::Flux> const & keys,
        afw::image::MaskedImage<PixelT> const & image,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas
    ) const
{
    source.set(keys.flag, true);  // failed unless we survive to switch this back
    ConvolvedFluxControl const& ctrl = getControl();
    float radius;
    try {
        radius = _getKronRadius(reference);
    } catch (NoKronException const&) {
        // No Kron, we've already flagged it, so bail
        return;
    }
    photometryKron::KronAperture const aperture(reference, refToMeas, radius);
    std::pair<double, double> result;
    try {
        result = aperture.measureFlux(image, ctrl.kronRadiusForFlux, ctrl.maxSincRadius);
    } catch (pex::exceptions::LengthErrorException const& e) {
        // We hit the edge of the image; there's no reasonable fallback or recovery
        source.set(keys.flag, true);
        return;
    }

    source.set(keys.meas, result.first);
    source.set(keys.err, result.second);
    source.set(keys.flag, utils::isnan(result.first) || utils::isnan(result.second));
}



LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(ConvolvedFlux);

} // anonymous namespace

PTR(algorithms::AlgorithmControl) ConvolvedFluxControl::_clone() const {
    return boost::make_shared<ConvolvedFluxControl>(*this);
}

PTR(algorithms::Algorithm) ConvolvedFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<ConvolvedFlux>(*this, boost::ref(schema));
}


}}}} // namespace lsst::meas::extensions::convolved

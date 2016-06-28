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

#include "lsst/meas/extensions/convolved/convolved.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace convolved {
namespace {

double const SIGMA_TO_FWHM = 2.0*std::sqrt(2.0*(std::log(2.0)));

class ConvolvedFlux : public algorithms::Algorithm {
public:

    ConvolvedFlux(ConvolvedFluxControl const & ctrl, afw::table::Schema & schema) :
        algorithms::Algorithm(ctrl),
        _failKey(schema.addField<afw::table::Flag>(ctrl.name + ".flag", "convolved algorithm failed"))
    {
        for (std::size_t ii = 0; ii < ctrl.seeing.size(); ++ii) {
            _deconvKeys.push_back(
                schema.addField<afw::table::Flag>(
                    (boost::format("%s.%d.deconv") % ctrl.name % ii).str(),
                    (boost::format("deconvolution required for seeing %f; no measurement made") %
                        ctrl.seeing[ii]).str()
                )
            );
            std::vector<afw::table::KeyTuple<afw::table::Flux> > apFluxKeys;
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

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

#if 0
    template <typename PixelT>
    void _applyForced(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas
    ) const;
#endif

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(ConvolvedFlux);
};


/************************************************************************************************************/

template <typename PixelT>
void ConvolvedFlux::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    // bad unless we get all the way to success at the end
    source.set(_failKey, true);
    ConvolvedFluxControl const & ctrl = static_cast<ConvolvedFluxControl const &>(this->getControl());
    for (std::size_t ii = 0; ii < ctrl.seeing.size(); ++ii) {
        for (std::size_t jj = 0; jj < ctrl.radius.size(); ++jj) {
            source.set(_apFluxKeys[ii][jj].flag, true);
        }
        source.set(_deconvKeys[ii], true);
    }

    if (!exposure.getWcs()) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "No WCS in exposure");
    }
    if (!exposure.getPsf()) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "No PSF in exposure");
    }
    double const seeing = exposure.getPsf()->computeShape(center).getDeterminantRadius();
    double const maxRadius = *std::max_element(ctrl.radius.begin(), ctrl.radius.end());

    afw::image::MaskedImage<PixelT> const mimage = exposure.getMaskedImage();

    for (std::size_t ii = 0; ii < ctrl.seeing.size(); ++ii) {
        // Check whether our target seeing requires deconvolution
        double const target = ctrl.seeing[ii]/SIGMA_TO_FWHM;
        if (target < seeing) {
            continue;
        }
        source.set(_deconvKeys[ii], false);

        // Convolve to our target seeing
        double const kernelSigma = std::sqrt(target*target - seeing*seeing);
        std::size_t const kernelRadius = ctrl.kernelScale*kernelSigma;
        std::size_t const kernelWidth = 2.0*kernelRadius + 1;

        afw::math::GaussianFunction1<double> const gauss(kernelSigma);
        afw::math::SeparableKernel const kernel(kernelWidth, kernelWidth, gauss, gauss);

        afw::geom::Box2I bbox = source.getFootprint()->getBBox();
        bbox.grow(kernelRadius + maxRadius); // XXX extra buffer?
        bbox.clip(mimage.getBBox(afw::image::PARENT));

        afw::image::MaskedImage<PixelT> subImage(mimage, bbox, afw::image::PARENT);
        afw::image::MaskedImage<PixelT> convolved(subImage, true);
        afw::math::convolve(convolved, subImage, kernel, afw::math::ConvolutionControl(true, true));

        // Perform aperture photometry
        for (std::size_t jj = 0; jj < ctrl.radius.size(); ++jj) {
            double const radius = ctrl.radius[jj];
            afw::geom::ellipses::Ellipse aperture(afw::geom::ellipses::Axes(radius, radius), center);
            std::pair<double, double> flux;
            try {
                flux = algorithms::photometry::calculateSincApertureFlux(convolved, aperture);
            } catch(pex::exceptions::Exception const&) {
                continue;
            }
            source.set(_apFluxKeys[ii][jj].meas, flux.first);
            source.set(_apFluxKeys[ii][jj].err, flux.second);
            source.set(_apFluxKeys[ii][jj].flag, false);
        }
    }

    source.set(_failKey, false);
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

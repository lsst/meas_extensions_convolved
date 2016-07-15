// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_CONVOLVED_H
#define LSST_MEAS_EXTENSIONS_CONVOLVED_H

#include "lsst/meas/algorithms/FluxControl.h"

namespace lsst { namespace meas { namespace extensions { namespace convolved {

/**
 *  @brief C++ control object for convolved fluxes.
 */
class ConvolvedFluxControl : public algorithms::FluxControl {
public:
    // convolution
    LSST_CONTROL_FIELD(seeing, std::vector<double>, "list of target seeings (FWHM, pixels)");
    LSST_CONTROL_FIELD(kernelScale, double, "scaling factor of kernel sigma for kernel size");

    // aperture flux
    LSST_CONTROL_FIELD(radius, std::vector<double>, "list of radii for aperture flux (pixels)");

    // Kron flux
    LSST_CONTROL_FIELD(kronName, std::string, "base name of Kron fields in reference");
    LSST_CONTROL_FIELD(maxSincRadius, double,
                       "Largest aperture for which to use the sinc aperture code for Kron");
    LSST_CONTROL_FIELD(kronRadiusForFlux, double, "Number of Kron radii for Kron flux");

    ConvolvedFluxControl() :
        algorithms::FluxControl("flux.convolved"),
        kernelScale(4.0),
        kronName("flux.kron"),
        maxSincRadius(10.0),
        kronRadiusForFlux(2.5)
    {
        seeing.push_back(3.5);
        seeing.push_back(5.0);
        seeing.push_back(6.5);
        radius.push_back(3.3);  // for 1.1 arcsec diameter, corresponding to PFS fiber size
        radius.push_back(4.5);  // 1.5 arcsec diameter for HSC
        radius.push_back(6.0);  // 2.0 arcsec diameter for HSC
    }

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::convolved

#endif // !LSST_MEAS_EXTENSIONS_CONVOLVED_H

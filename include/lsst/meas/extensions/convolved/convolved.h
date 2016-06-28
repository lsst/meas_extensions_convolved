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
    LSST_CONTROL_FIELD(seeing, std::vector<double>, "list of target seeings (FWHM, pixels)");
    LSST_CONTROL_FIELD(radius, std::vector<double>, "list of radii for aperture flux (pixels)");
    LSST_CONTROL_FIELD(kernelScale, double, "scaling factor of kernel sigma for kernel size");

    ConvolvedFluxControl() :
        algorithms::FluxControl("flux.convolved"),
        kernelScale(4.0)
    {
        seeing.push_back(3.5);
        seeing.push_back(5.0);
        seeing.push_back(6.5);
        radius.push_back(6.5);  // for 1.1 arcsec diameter, corresponding to PFS fiber size
    }

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::convolved

#endif // !LSST_MEAS_EXTENSIONS_CONVOLVED_H

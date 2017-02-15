#ifndef NURBS_NORM_H
#define NURBS_NORM_H

#include "MultiForest.h"
#include <vector>
#include <complex>

namespace nurbs {

    /// Compute the H_{-1/2}(div) norm
    double hdivNorm(const MultiForest& f,
                    const std::vector<std::complex<double>>& phi);
    
    /// Comptue the L2 graph norm for the Mie sphere scattering problem
    /// assuming wave is travelling in +z direction and polarised in the
    /// x-direction.
    double L2graphNormMieSphere(const MultiForest& f,
                                const double k,
                                const std::vector<std::complex<double>>& soln);
}

#endif
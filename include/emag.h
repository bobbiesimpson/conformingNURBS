#ifndef EMAG_H
#define EMAG_H

#include <complex>
#include "Point3D.h"
#include "base.h"

namespace nurbs 
{

    /// The Green's function for electromagnetics
    std::complex<double> emagKernel(const double r, const double k);
    
    /// Get the complex components of an electromagnetic plane wave
    /// given the wavevector (k), polarisation vector (pvec) and
    /// physical coordinate (x)
    ComplexDoubleVec planeWave(const Point3D& k,
                               const Point3D& polarvec,
                               const Point3D& x);
    
    /// Get plane wave component tangent to surface
    ComplexDoubleVec planeWaveTangent(const Point3D& k,
                                      const Point3D& polarvec,
                                      const Point3D& x,
                                      const Point3D& n);
}

#endif

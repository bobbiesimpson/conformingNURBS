#include "emag.h"
#include <cmath>

namespace nurbs 
{
	std::complex<double> emagKernel(const double r, const double k)
    {
        std::complex<double> power(0.0, k * r);
        return std::exp(power) / (4.0 * PI * r);
    }
    
    ComplexDoubleVec planeWave(const Point3D& k,
                               const Point3D& polarvec,
                               const Point3D& x)
    {
        ComplexDoubleVec result{ polarvec[0], polarvec[1], polarvec[2] };
        const auto wave = std::exp(-std::complex<double>(0.0, dot(k, x)));
        for(auto& r : result)
            r *= wave;
        return result;
    }
    
    ComplexDoubleVec planeWaveTangent(const Point3D& k,
                                      const Point3D& polarvec,
                                      const Point3D& x,
                                      const Point3D& n)
    {
        const auto pw = planeWave(k, polarvec, x);
        Point3D pw_r, pw_i;
        for(uint i = 0; i < 3; ++i) {
            pw_r[i] = pw[i].real(); pw_i[i] = pw[i].imag();
        }
        const auto pwt_r = cross(n, cross(pw_r, n));
        const auto pwt_i = cross(n, cross(pw_i, n));
        ComplexDoubleVec result(3);
        for(uint i = 0; i < 3; ++i)
            result[i] = std::complex<double>(pwt_r[i], pwt_i[i]);
        return result;
        
    }
}

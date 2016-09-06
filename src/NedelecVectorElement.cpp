#include "NedelecVectorElement.h"
#include "BezierVectorElement.h"
#include "Point4D.h"
#include "Geometry.h"
#include "NURBSCache.h"

#include <stdexcept>

namespace nurbs
{
    Point3D NedelecVectorElement::eval(const double u, const double v) const
    {
        
        // We get a reference to the parent element, convert to its local
        // parent coordinates and evaluate using the basis that belongs to the
        // geometry (often much coarse then the analysis basis)

            const auto pel = parent();
            const auto p_gpt = transformToParentElParentCoord(GPt2D(u,v));
            const auto& b = pel->basis(p_gpt.s, p_gpt.t);
            const auto& gvec = pel->globalBasisFuncI();
            Point4D x;
            for(uint i = 0; i < gvec.size(); ++i)
                x += geometry().controlPt(gvec[i]) * b[i];
            Point3D p = x.asCartesian();

            return p;

    }
    
    Point3D NedelecVectorElement::tangent(const double u,
                                         const double v,
                                         const ParamDir dir) const
    {
        // efficiency is crucial here for fast BE analysis

            const auto pel = parent();
            const auto ppt = transformToParentElParentCoord(GPt2D(u,v));
            const auto gvec = pel->globalBasisFuncI();
            const auto& basis = pel->basis(ppt.s, ppt.t);
            
            double w_interp = 0.0;
            for(uint i = 0; i < gvec.size(); ++i)
                w_interp += geometry().controlPt(gvec[i]).getWeight() * basis[i];
            
            // get non-rational bspline basis function parametric derivatives
            DerivType dtype = (ParamDir::S == dir) ? DS : DT;
            const auto& basis_der = pel->localBasisDers(ppt.s, ppt.t, dtype);
            
            // calculate weight function derivative
            double w_der = 0.0;
            for(uint i = 0; i < gvec.size(); ++i)
                w_der += geometry().controlPt(gvec[i]).getWeight() * basis_der[i];
            
            Point3D result;
            for(uint i = 0; i < gvec.size(); ++i) {
                const auto term = geometry().controlPt(gvec[i]).getWeight() *
                (1.0 / w_interp * basis_der[i] -
                 1.0 / (w_interp * w_interp) * w_der * basis[i]);
                result += geometry().controlPt(gvec[i]).asCartesian() * term;
            }

            return result;

    }
    
    std::vector<double> LegendreBasis1D(const double xi,
                                        const uint p)
    {
        switch(p)
        {
            case 0:
                return {1.0};
                
            case 1:
                return
            {
                0.5 * (1.0 - xi),
                0.5 * (1.0 + xi)
            };
                
            case 2:
                return
            {
                0.5 * xi * (xi - 1.0),
                0.5 * (1.0 - xi * xi),
                0.5 * xi * (xi + 1.0)
            };
                
            case 3:
                return
            {
                -9.0 / 16.0 * ( xi + 1.0 / 3.0 ) * ( xi - 1.0 / 3.0 ) * ( xi - 1.0 ),
                27.0 / 16.0 * ( xi + 1.0 ) * ( xi - 1.0 / 3.0 ) * ( xi - 1.0 ),
                -27.0 / 16.0 * ( xi + 1.0 ) * ( xi + 1.0 / 3.0 ) * ( xi - 1.0 ),
                9.0 / 16.0 * ( xi + 1.0 ) * ( xi + 1.0 / 3.0 ) * ( xi - 1.0 / 3.0 )
            };
                
            default:
                throw std::runtime_error("Legendre basis function of degree > 3 not yet implemented");
        }
    }
    
    std::vector<double> LegendreBasisDer1D(const double xi,
                                           const uint p)
    {
        switch(p)
        {
            case 0:
                return {0.0};
                
            case 1:
                return
            {
                -0.5,
                0.5
            };
                
            case 2:
                return
            {
                (xi - 0.5),
                -2.0 * xi,
                xi + 0.5
            };
                
            case 3:
                return
            {
                -9.0 / 16.0 * ( 3.0 * xi * xi - 2.0 * xi - 1.0 / 9.0 ),
                27.0 / 16.0 * ( 3.0 * xi * xi - 2.0 * xi / 3.0 - 1.0 ),
                -27.0 / 16.0 * ( 3.0 * xi * xi + 2.0 * xi / 3.0 - 1.0 ),
                9.0 / 16.0 * ( 3.0 * xi * xi + 2.0 * xi - 1.0 / 9.0 )
            };
                
            default:
                throw std::runtime_error("Legendre basis function of degree > 3 not yet implemented");
        }
    }
}
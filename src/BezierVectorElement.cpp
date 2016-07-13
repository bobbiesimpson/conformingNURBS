#include "BezierVectorElement.h"
#include "Point4D.h"
#include "Geometry.h"

namespace nurbs {
    
    Point3D BezierVectorElement::eval(const double u, const double v) const
    {
        
        // We get a reference to the parent element, convert to its local
        // parent coordinates and evaluate using the basis that belongs to the
        // geometry (often much coarse then the analysis basis)
        const auto pel = parent();
        const auto p_gpt = transformToParentElParentCoord(GPt2D(u,v));
        const auto& b = pel->basis(p_gpt.s, p_gpt.t);
        const auto gvec = pel->globalBasisFuncI();
        Point4D x;
        for(uint i = 0; i < gvec.size(); ++i)
            x += geometry().controlPt(gvec[i]) * b[i];
        return x.asCartesian();
    }
    
    Point3D BezierVectorElement::tangent(const double u,
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
    
}
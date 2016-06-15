#include "GeometryElement2D.h"
#include "Geometry2D.h"

namespace nurbs {
    
    /// Constructor
    GeometryElement2D::GeometryElement2D(const Geometry2D& g,
                                         const uint sp,
                                         const DoublePair& knots)
    :
    mGeom(g),
    mSpaceI(sp)
    {
        assert((knots.second - knots.first) > 0.0);
        knotInterval() = boost::icl::construct<boost::icl::continuous_interval<double>>(knots.first, knots.second);
    }
    
    /// Interpolate the physical coordinate given a parent coord xi \in [-1,1]
    Point3D GeometryElement2D::eval(const double xi) const
    {
        return geometry().eval(paramCoord(xi), spaceI());
    }
    
    /// Jaoobiain determinant that maps from parent to physical space
    /// Implementation requires a function composition from
    /// parent space -> parameter space -> physical space
    double GeometryElement2D::jacDet(const double xi) const
    {
        return geometry().jacDet(paramCoord(xi), spaceI()) * jacDetParam(xi);
    }
    
    /// Normal evaluation at parent coord xi \in [-1,1]
    Point3D GeometryElement2D::normal(const double xi) const
    {
        return geometry().normal(paramCoord(xi), spaceI());
    }
    
    /// Normal evaluation at parent coord xi \in [-1,1]
    Point3D GeometryElement2D::tangent(const double xi) const
    {
        return geometry().tangent(paramCoord(xi), spaceI());
    }
    
    /// Get geometry order
    uint GeometryElement2D::geometryDegree() const
    {
        return geometry().degree(spaceI());
    }
    
    /// Print function. May be overridden.
    void GeometryElement2D::print(std::ostream& ost) const
    {
        ost << "Geometry: " << &geometry() << "\n";
        ost << "Space index: " << spaceI() << "\n";
        ost << "Knot interval: [" << lowerBound() << "," << upperBound() << "]";
    }
    
}
#ifndef GEOMETRY_ELEMENT_2D_H
#define GEOMETRY_ELEMENT_2D_H

#include <vector>
#include <cassert>

#include "base.h"
#include "Point3D.h"

#include <boost/icl/continuous_interval.hpp>

namespace nurbs {
    
    /// Base class for an element which provides all the 'geometry'
    /// information. Basis function data is provided by the derived
    /// class 'AnalysisElement2D'.
    
    /// This is not an abstract class and can therefore be used directly.
    
    /// Forward declarations
    class Geometry2D;
    
    class GeometryElement2D {
        
    public:
        
        /// Constructor
        GeometryElement2D(const Geometry2D& g,
                          const uint sp,
                          const DoublePair& knots);
        
        /// Const geometry getter
        const Geometry2D& geometry() const { return mGeom; }
        
        /// Const space index getter
        uint spaceI() const { return mSpaceI; }
        
        /// Interpolate the physical coordinate given a parent coord xi \in [-1,1]
        Point3D eval(const double xi) const;
        
        /// Jaoobiain determinant that maps from parent to physical space
        /// Implementation requires a function composition from
        /// parent space -> parameter space -> physical space
        double jacDet(const double xi) const;
        
        /// Normal evaluation at parent coord xi \in [-1,1]
        Point3D normal(const double xi) const;
        
        /// Normal evaluation at parent coord xi \in [-1,1]
        Point3D tangent(const double xi) const;
        
        /// Get geometry order
        uint geometryDegree() const;
        
        /// Print function. May be overridden.
        virtual void print(std::ostream& ost) const;
        
    protected:
        
        /// Get lower bound of knot interval
        double lowerBound() const { return knotInterval().lower(); }
        
         /// Get upper bound of knot interval
        double upperBound() const { return knotInterval().upper(); }
        
        /// Jacobian determinant for mapping from parent space to parametric space
        double jacDetParam(const double xi) const { return 0.5 * (upperBound() - lowerBound()); }
        
        /// Get parametric coordinate given parent coordinate
        double paramCoord(const double xi) const { return xi * (upperBound() - lowerBound()) * 0.5 + (upperBound() + lowerBound()) * 0.5; }
        
    private:
        
        /// Non-const interval getter
        Interval& knotInterval() { return mKnotInterval; }
        
        /// Const interval getter
        const Interval& knotInterval() const { return mKnotInterval; }
        
        /// Reference to the geometry instance
        const Geometry2D& mGeom;
        
        /// The space index this element lies in.
        const uint mSpaceI;
        
        /// Knot interva that defines this geometry element
        Interval mKnotInterval;
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const GeometryElement2D& e) { e.print(ost); return ost; }
        
    };
}

#endif
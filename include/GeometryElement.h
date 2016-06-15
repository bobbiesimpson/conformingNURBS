#ifndef NURBS_BASE_ELEMENT_H
#define NURBS_BASE_ELEMENT_H

#include "base.h"
#include "Point3D.h"
#include "IElemIntegrate.h"

#include <boost/icl/continuous_interval.hpp>

#include <vector>
#include <cassert>
#include <algorithm>
#include <map>
#include <utility>
#include <mutex>

namespace nurbs
{
    /// The base element class. Can be instantianted (non-abstract class)
    /// Represents a geometry element providing the neccessary interface
    /// for derived classes.
    
    /// (u,v) \in [-1,1] x [-1,1]
    
    class Geometry;
    
    class GeometryElement {
        
    public:
        
        /// Construct with a geometry object, space index
        /// and knot coordinates which define the element bounds
        GeometryElement(const Geometry& g,
                        const uint sp,
                        const DoublePairVec& knots);
        
        /// Get the const geometry object
        const Geometry& geometry() const {return mGeom;}
        
        /// Space index accessor
        uint spaceI() const {return mSpaceI;}
        
        /// Evaluate the physical point given a parent coordinate.
        Point3D eval(const double u, const double v) const;
        
        /// Same as above but passing a 2d guass point
        Point3D eval(const GPt2D
                     & gp) const {return eval(gp.s, gp.t); }
        
        /// Get the physical coordinate given a vertex in the parent space
        Point3D evalVertex(const Vertex v) const;
        
        /// Get jacobian determinant from parent to physical space
        double jacDet(const double u, const double v) const;
        
        /// Same as above but passing a 2d guass point
        double jacDet(const GPt2D& gp) const {return jacDet(gp.s, gp.t); }
        
        /// Get jacobian
        DoubleVecVec jacob(const double u, const double v) const;
        
        /// Same as above but passing a 2d guass point
        DoubleVecVec jacob(const GPt2D& gp) const {return jacob(gp.s, gp.t); }
        
        /// Get the normal
        Point3D normal(const double u, const double v) const;
        
        /// Same as above but passing a 2d guass point
        Point3D normal(const GPt2D& gp) const {return normal(gp.s, gp.t); }
        
        /// Get the tangent
        Point3D tangent(const double u, const double v, const ParamDir d) const;
        
        /// Same as above but passing a 2d guass point
        Point3D tangent(const GPt2D& gp, const ParamDir d) const {return tangent(gp.s, gp.t, d); }
        
        /// Get the degree for all parametric directions.
        UIntVec geometryDegree() const;
        
        /// Does this element contain the given parametric point? If so, return true
        /// with local parent coordinate
        std::pair<bool, GPt2D> containsParamPt(const double s, const double t) const;
        
        /// Get approximate element size
        double size() const
        {
            return std::max(dist(eval(-1.0, -1.0), eval(1.0, 1.0)),
                            dist(eval(1.0, -1.0), eval(-1.0, 1.0)));
        }
        
        /// Print function that may be overridden as required.
        virtual void print(std::ostream& ost) const;
        
    protected:
        
        DoubleVecVec jacobParam(const double u, const double v) const
        {
            return { { (upperBound(S) - lowerBound(S)) / 2.0, 0.0},
                { 0.0, (upperBound(T) - lowerBound(T)) / 2.0}};
        }
        /// Jacobian determinant for mapping from parent to parametric
        /// space
        double jacDetParam(const double u, const double v) const
        {
            return (upperBound(S) - lowerBound(S)) / 2.0 *
                   (upperBound(T) - lowerBound(T)) / 2.0;
        }
                                            
        /// Transform parent coordinate to parametric coordinate
        ParamCoord paramCoord(const double u, const double v) const
        {
            ParamCoord p;
            p.s = u * (upperBound(S) - lowerBound(S)) / 2.0 + (upperBound(S) + lowerBound(S)) / 2.0;
            p.t = v * (upperBound(T) - lowerBound(T)) / 2.0 + (upperBound(T) + lowerBound(T)) / 2.0;
            return p;
        }
        
        /// Get the parent coordinate given a paramatric coordinate
        GPt2D parentCoord(const double s, const double t) const
        {
            return GPt2D((2.0 * s - (upperBound(S) + lowerBound(S))) / (upperBound(S) - lowerBound(S)),
                         (2.0 * t - (upperBound(T) + lowerBound(T))) / (upperBound(T) - lowerBound(T)));
        }
        
        /// Get lower bound (knot coordinate) for specified parametric direction
        double lowerBound(const ParamDir d) const
        {
            assert(d < mKnotIntervals.size());
            return mKnotIntervals[d].lower();
        }
        
        /// Get upper bound (knot coordinate) for specified parametric direction
        double upperBound(const ParamDir d) const
        {
            assert(d < mKnotIntervals.size());
            return mKnotIntervals[d].upper();
        }
        
        /// Get the parametric interval for this element for the given parametric direction
        const boost::icl::continuous_interval<double>& boostInterval(const ParamDir dir) const
        {
            assert(mKnotIntervals.size() > 1);
            return mKnotIntervals[dir];
        }
        
    private:
        
        /// Reference to the forest where this element belongs. Non-owning pointer
        const Geometry& mGeom;
        
        /// The space that this element lies in
        const uint mSpaceI;
        
        /// Mutex to protected cached data
        std::shared_ptr<std::mutex> mMutex;
        
        void insertCachedPoint(const std::pair<double, double>& p, const Point3D& pt) const
        {
            std::lock_guard<std::mutex> lock(*mMutex);
            mPointCache.insert(std::make_pair(p,pt));
        }
        
        void insertCachedJDet(const std::pair<double, double>& p, const double jdet) const
        {
            std::lock_guard<std::mutex> lock(*mMutex);
            mJDetCache.insert(std::make_pair(p,jdet));
        }
        
        void insertCachedJacob(const std::pair<double, double>& p, const DoubleVecVec& jacob) const
        {
            std::lock_guard<std::mutex> lock(*mMutex);
            mJacobCache.insert(std::make_pair(p,jacob));
        }
        
        /// The non-zero knot interval that defines this element
        std::vector<boost::icl::continuous_interval<double>> mKnotIntervals;
        
        /// Cache evaluation of physical points
        mutable std::map<std::pair<double, double>, Point3D> mPointCache;
        
        /// Cache jacobian determinant
        mutable std::map<std::pair<double, double>, double> mJDetCache;
        
        /// Cache jacobian matrix
        mutable std::map<std::pair<double, double>, DoubleVecVec> mJacobCache;
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const GeometryElement& e)
        { e.print(ost); return ost; }
        
    };
    
    /// Get distance between two elements
    double dist(const GeometryElement& e1, const GeometryElement& e);
	
}

#endif

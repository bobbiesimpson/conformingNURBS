//
//  BSplineSpace2D.h
//  nurbslib
//
//  Created by Robert Simpson on 05/12/2014.
//
//

#ifndef BSPLINESPACE2D_H
#define BSPLINESPACE2D_H

#include <iostream>
#include <vector>

#include "NURBSCommon.h"
#include "base.h"

#include <boost/icl/continuous_interval.hpp>
#include <boost/icl/closed_interval.hpp>

namespace nurbs {
    
    /// A class that represents a 2d b-spline space
    /// containing the following data structres:
    ///     - knot vector
    ///     - degree
    ///
    /// Other data structures are functions of these primary
    /// data structures, including:
    ///     - unique knot vector (defining non-zero parametric intervals)
    ///     - parametric range.
    
    class BSplineSpace2D {
        
    public:
    
        /// Default constructor
        BSplineSpace2D() { clear(); }
        
        /// Construct with knot vector, degree and identifier string
        BSplineSpace2D(const DoubleVec& knotvec,
                       const uint d,
                       const std::string& n = "no_name")
        :
        mKnotVec(knotvec),
        mDegree(d),
        mName(n)
        { init(); }
        
        /// Return the total number of basis functions
        uint basisFuncN() const { return mNumBasis; }
        
        /// Get a vector of degrees
        uint degree() const { return mDegree; }
        
        /// Get knot vector for specified parametric direction
        inline DoubleVec knotVec() const { return mKnotVec; }
        
        /// Get unique knot vectr
        inline DoubleVec uniqueKnotVec() const { return mUniqueKnotVec; }
        
        /// knot coordinate getter
        inline double knot(const uint i) const
        {
            assert(i < mKnotVec.size());
            return mKnotVec[i];
        }
        
        /// unique knot coordinate getter
        inline double uniqueKnot(const uint i) const
        {
            assert(i < mUniqueKnotVec.size());
            return mUniqueKnotVec[i];
        }
        
        /// unique knot number getter
        inline uint uniqueKnotN() const { return mUniqueKnotVec.size(); }
        
        /// Total number of non-zero knot spans (i.e. elements)
        inline uint nonzeroKnotSpanN() const { return uniqueKnotN() - 1; }
        
        /// For now, we are dealing with parametric spaces in R^1
        inline uint paramDimN() const { return 1; }
        
        /// Lower limit of parametric domain
        inline double lowerParamLimit() const { return mInterval.lower(); }
        
        /// Upper limit of parametric domain
        inline double upperParamLimit() const { return mInterval.upper(); }
        
        /// Get bspline basis given parametric coordinaten
        DoubleVec basis(const double s) const
        {return basis(s,span(s));}
        
        /// Get basis bspline basis with given knot span
        DoubleVec basis(const double s, const uint span) const
        {
            return nurbshelper::getBsplineBasis(s, span, knotVec(), degree());
        }
        
        /// Get basis function derivatives with tensor product multiplied out
        DoubleVec basisDers(const double s, const DerivOrder der = D1) const
        { return basisDers(s, span(s), der); }
        
        /// Get basis function derivatives with tensor product multiplied out
        DoubleVec basisDers(const double s, const uint span,
                            const DerivOrder der = D1) const
        {
            return nurbshelper::getBsplineBasisDers(s, span, knotVec(), degree(), der).at(der);
        }
        
        
        /// Return the non-zero local basis function indices for this parametric coordinate
        /// Returns a set of vectors corresponding to the indices in each parametric direction
        UIntVec localBasisFuncI(const double s) const
        {
            return nurbshelper::getBasisFnIndices(s, knotVec(), degree());
        }
        
        /// Get the knot indices for given parametric coordinate
        uint span(const double s) const
        {
            assert(validCoord(s));
            return nurbshelper::getKnotSpan(s, knotVec(), degree());
        }
        
        /// Const name accessor
        std::string name() const { return mName; }
        
        /// Non-const name accessor
        std::string& name() { return mName; }
        
        /// Load from an input stream
        void load(std::istream& ist);
        
    private:
        
        /// name setter
        void setName(const std::string& s) { mName = s; }
        
        /// degree setter
        void setDegree(const uint d) { mDegree = d; }
        
        /// Knot vector setter
        void setKnotVec(const DoubleVec& kv) { mKnotVec = kv; }
        
        /// Unique knot vector
        void setUniqueKnotVec(const DoubleVec& kv) { mUniqueKnotVec = kv; }
        
        /// recalculate unique knots and intervals
        void init()
        {
            if(essentiallyEqual(knotVec().front(), knotVec().back(), TOL))
                error("Cannot prescribe a b-spline space with non-zero parametric area.");
            mInterval = boost::icl::construct<Interval>
                      (knotVec().front(), knotVec().back(),
                       boost::icl::interval_bounds::closed());
            mNumBasis = knotVec().size() - degree() - 1;
            // construct unique knot vectors
            DoubleVec kv_copy(knotVec());
            auto last = std::unique(kv_copy.begin(), kv_copy.end());
            setUniqueKnotVec(DoubleVec(kv_copy.begin(), last));
        }
        
        /// Clear all data
        void clear()
        {
            mKnotVec.clear();
            mUniqueKnotVec.clear();
            mDegree = 0;
            mNumBasis = 0;
            mInterval = boost::icl::construct<Interval>
            (0, 0, boost::icl::interval_bounds::closed());
        }
        
        /// Is this a valid parametric coordinate? I.e. does it lie in the parametric space?
        bool validCoord(const double s) const
        {
            if(!boost::icl::contains(mInterval, s))
                    return false;
            return true;
        }
        
    private:
        
        /// Member variables
        
        /// Knot vector
        DoubleVec mKnotVec;
        
        /// Degree
        uint mDegree;
        
        /// Name of bspline space
        std::string mName;
        
        /// Unique knot vector
        DoubleVec mUniqueKnotVec;
        
        /// Basis func. number
        uint mNumBasis;
        
        /// Range
        Interval mInterval;
        
        /// Private member functions
        
        /// Print member data
        void printData(std::ostream& ost) const;
        
        /// overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const BSplineSpace2D& s)
        {
            s.printData(ost);
            return ost;
        }
        
        /// Overload input operator
        friend std::istream& operator>>(std::istream& ist, BSplineSpace2D& s)
        {
            s.load(ist);
            return ist;
        }
        
    };
    
}

#endif

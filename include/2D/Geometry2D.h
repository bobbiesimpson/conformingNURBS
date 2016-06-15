//
//  Geometry2D.h
//  nurbslib
//
//  Created by Robert Simpson on 05/12/2014.
//
//

#ifndef GEOMETRY_2D_H
#define GEOMETRY_2D_H

#include <stdio.h>
#include <cassert>

#include "Forest2D.h"
#include "Point4D.h"

namespace nurbs {
    
    /// A class that represent the geometry of a 2d surface using nurbs
    /// interpolation. Basically a wrapper around a forest and
    /// set of control points
    
    class Geometry2D {
        
    public:
        
        Geometry2D()
        :
        mPrimalForest(nullptr) {}
        
        /// Construct with a forest and set of control points
//        Geometry2D(const Forest2D& f,
//                   const std::vector<Point4D>& cpts)
//        :
//        mPrimalForest(f),
//        mCPts(cpts) {}
        
        Geometry2D(const Geometry2D& g)
        {
            mPrimalForest = g.mPrimalForest;
            mPrimalForest.setGeometry(this);
            mCPts = g.mCPts;
        }
        
        Geometry2D& operator=(const Geometry2D& g)
        {
            if(this == &g)
                return *this;
            mPrimalForest = g.mPrimalForest;
            mPrimalForest.setGeometry(this);
            mCPts = g.mCPts;
            return *this;
        }
        
        /// clear member data
        void clear()
        {
            primalForest().clear();
            cPtVec().empty();
        }
        
        /// primal forest getter
        inline const Forest2D& primalForest() const { return mPrimalForest; }
        
        /// const control point set getter
        const std::vector<Point4D>& cPtVec() const { return mCPts; }
        
        /// Get a control point at a specified index
        const Point4D& controlPt(const uint i) const
        {
            assert(i < cPtVec().size());
            return mCPts[i];
        }
        
        /// Number of B-spline spaces
        uint spaceN() const
        { return primalForest().spaceN(); }
        
        /// Degree of given space index
        uint degree(const uint sp) const
        {
            assert(sp < spaceN());
            return primalForest().space(sp).degree();
        }
        
        /// Get the jacobian det. at given parametric coord and space index
        double jacDet(const double s, const uint sp) const;
        
        /// Get normal at given parametric coord and space index
        Point3D normal(const double s, const uint sp) const;
        
        /// Get tangent at given parametric coord and space index
        Point3D tangent(const double s, const uint sp) const;
        
        /// Interpolate geometry at given parametric coord and space index
        Point3D eval(const double s, const uint sp) const;
        
        /// Write to vtu file.
        void writeVTPOutput(const std::string& file,
                            const uint nsample = DEFAULT_NGRID_PTS) const;
        
        /// Load from a file stream
        bool load(std::istream& ist);
        
        /// Print to a file stream
        void print(std::ostream& ost) const;
        
        
    protected:
        
        /// Get a control point given a space index and local node index.
        const Point4D& point(const uint sp, const uint i) const { return controlPt(primalForest().globalI(sp, i)); }
        
    private:
        
        /// non-const primal forest getter
        inline Forest2D& primalForest() { return mPrimalForest; }
        
        /// Primal forest setter
//        inline void setPrimalForest(const Forest2D& f) { mPrimalForest = f; }
        
        /// non-const control point set getter
        std::vector<Point4D>& cPtVec() { return mCPts; }
        
        /// control point vector setter
        inline void setCPtVec(const std::vector<Point4D>& cp_vec) { mCPts = cp_vec; }
        
        /// The primal forest (i.e. geometry basis)
        Forest2D mPrimalForest;
        
        /// Control point set
        std::vector<Point4D> mCPts;
        
        /// Input operator
        friend std::istream& operator>>(std::istream& ist, Geometry2D& g);
        
        /// Output operator
        friend std::ostream& operator<<(std::ostream& ost, const Geometry2D& g);
        
        
    };
}

#endif

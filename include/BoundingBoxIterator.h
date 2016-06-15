#ifndef NURBS_BOUNDING_BOX_H
#define NURBS_BOUNDING_BOX_H

/* Written by Robert N. Simpson
   University of Glasgow
   15th June 2016
*/

#include "base.h"

namespace nurbs {
    
    /// Forward declarations
    class Forest;
    class Point3D;

    /// A class that represents an iterator for computing the set of points and
    /// bounding boxes for the set of of basis functions defined by the given forest.
    ///
    /// The points are based on Greville abscissa and the bounding boxes are computed
    /// approximated from the span of each of the NURBS basis functions.
    ///
    /// The main application of this class is for fast kernel approximation using H-matrices
    /// with the HLibPro library (hlibpro.com)
    
    class BoundingBoxIterator {
        
    public:
        
        /// Construct with a given forest
        BoundingBoxIterator(const Forest& f);
        
        /// Increment iterator to next basis function.
        BoundingBoxIterator& operator++();
        
        /// Have we finished iterating over all global basis functions for this forest?
        bool isDone() const;
        
        /// Reset the counter
        void reset()
        {
            mCurrentIndex = 0;
        }
        
        /// Current index getter
        uint currentIndex() const
        {
            return mCurrentIndex;
        }
        
        /// Get the current point that defines the 'anchor' for the current global basis function.
        Point3D currentPt() const;
        
        /// Get the limits that define the current bounding box for the current global basis function.
        std::pair<Point3D, Point3D> currentBoundingBox() const;
        
        /// Get lower bound
        Point3D currentLowerBound() const;
        
        /// Get current upper bound
        Point3D currentUpperBound() const;
        
        /// Forest getter
        const Forest& forest() const { return mForest; }
        
    private:
        
        /// Current basis function indeix
        uint mCurrentIndex;
        
        /// Reference to (non-null) forest
        const Forest& mForest;
        
        /// Vector of greville point data
        std::vector<Point3D> mPtData;
        
        /// Vector of bounding box data
        std::vector<std::pair<Point3D, Point3D>> mBBData;
        
    };
    
}

#endif
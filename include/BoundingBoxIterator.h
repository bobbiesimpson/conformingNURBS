#ifndef NURBS_BOUNDING_BOX_H
#define NURBS_BOUNDING_BOX_H

/* Written by Robert N. Simpson
   University of Glasgow
   15th June 2016
*/

#include "Point3D.h"

namespace nurbs {
    
    /// Forward declarations
    class Forest;

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
        BoundingBoxIterator(const Forest& f)
        :
        mForest(f)
        {}
        
        /// Increment iterator to next basis function.
        BoundingBoxIterator& operator++();
        
        /// Have we finished iterating over all global basis functions for this forest?
        bool isDone() const;
        
        /// Get the current point that defines the 'anchor' for the current global basis function.
        Point3D currentPt() const;
        
        /// Get the limits that define the current bounding box for the current global basis function.
        DoublePairVec currentBoundingBox() const;
        
    private:
        
        /// Reference to (non-null) forest
        const Forest& mForest;
        
    };
    
    
}
#endif
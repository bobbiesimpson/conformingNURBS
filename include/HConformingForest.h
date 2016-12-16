#ifndef HDIV_FOREST_H
#define HDIV_FOREST_H

#include "MultiForest.h"
#include "base.h"

namespace nurbs {
    
    /// Define the HDiv MultiForest. Simply a wrapper around MultiForest.
    class HDivForest : public MultiForest {
        
    public:
        
        /// Default constructor
        HDivForest() : MultiForest() {}
        
        /// Construct with geometry. Initialise spaces with degrees:
        /// {p, p-1} x {p-1, p}
        HDivForest(const Geometry& g) : MultiForest(g)
        {
            initSpaces({0,1},{1,0});
        }
        
        /// Destructor
        virtual ~HDivForest() {}
        
        /// Override the Piola transform function
        DoubleVecVec transformBasis(const DoubleVecVec& basis,
                                    const DoubleVecVec& jacob,
                                    const double jdet) const override;
        
    private:
        
        /// Override continuity type getter
        virtual ContinuityType continuityType() const override
        { return ContinuityType::NORMAL; }
    };
    
    /// Define the HDiv MultiForest. Simply a wrapper around MultiForest.
    class HCurlForest : public MultiForest {
        
    public:
        
        /// Default constructor
        HCurlForest() : MultiForest() {}
        
        /// Construct with geometry. Initialise spaces with degrees:
        /// {p, p-1} x {p-1, p}
        HCurlForest(const Geometry& g) : MultiForest(g)
        {
            initSpaces({1,0},{0,1});
        }
        
        /// Destructor
        virtual ~HCurlForest() {}
        
        /// Override the Piola transform function
        DoubleVecVec transformBasis(const DoubleVecVec& basis,
                                    const DoubleVecVec& jacob,
                                    const double jdet) const override;
    private:

        
        /// Override continuity type getter
        virtual ContinuityType continuityType() const override
        { return ContinuityType::TANGENT; }
        
    };
}

#endif

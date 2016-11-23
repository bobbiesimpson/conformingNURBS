#ifndef NURBS_ISUBELEM_H
#define NURBS_ISUBELEM_H

#include "base.h"
#include "IElemIntegrate.h"

namespace nurbs
{
    namespace elem
    {
        
        /// A class for computing quadrature points in subelements
        class ISubElem {
            
        public:
            
            /// Construct with number of subintervals in s- and t-
            /// parametric directions
            ISubElem(const uint ns = 1, const uint nt = 1);
            
            /// Consructor with vector
            ISubElem(const UIntVec& n)
            :
            ISubElem(n[0], n[1]) {}
            
            /// Construct subelements with specified knots to divide
            /// the subelements
            ISubElem(const nurbs::DoubleVec& sknots,
                     const nurbs::DoubleVec& tknots);
            
            
            /// get current subelement range
            DoubleVec getRange() const;
            
            /// get subelement jacobian from subelement space to parent space [-1,1] x [-1,1]
            double jacob() const;
            
            /// convert gauss point to parent space (s,t)
            GPt2D get(const GPt2D& pt) const;
            
            /// Have we finished iterating?
            inline bool isDone() const { return mCurrentIndex == mNSubEls; }
            
            /// Restart the iterator
            void restart() { mCurrentIndex = 0; }
            
            /// Increment the iterator
            ISubElem& operator++()
            {
                ++mCurrentIndex;
                return *this;
            }
            
            /// Current sub cell index getter
            uint currentIndex() const { return mCurrentIndex;}
            
        private:
            
            /// Number of subelements
            uint mNSubEls;
            
            /// Current sub element we are on
            uint mCurrentIndex;
            
            /// Vectors of subelement ranges
            DoubleVecVec mRanges;
            
        };
    }
}


#endif

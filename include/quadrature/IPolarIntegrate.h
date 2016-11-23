#ifndef IPOLAR_INTEGRATE_H
#define IPOLAR_INTEGRATE_H

#include "base.h"
#include "IElemIntegrate.h"
#include "ISubElem.h"
#include "IBaseIntegrate.h"

namespace nurbs {
    
    namespace elem {
        
        /// An iterator class for calculating integrals containing a 1/R singularity defined
        /// over a quadrilateral element in R^2.
        ///
        /// The tranformation to polar coordinates automatically cancels the 1/R singularity.
        
        class IPolarIntegrate : public IBaseIntegrate  {
            
        public:
            
            /// Construct with equal number of quadrature points in each parametric
            /// direction. nrho specifies number of sub elements in radial
            /// direction.
            IPolarIntegrate(const GPt2D& spt,
                            const uint n  = DEFAULT_NGPS,
                            const uint nang = 1,
                            const uint nradial = 1)
            :
            IPolarIntegrate(spt, {n,n}, {nang, nradial}) {}
            
            /// Construct with vector of integration orders
            IPolarIntegrate(const GPt2D& spt,
                            const UIntVec& orders,
                            const uint nang = 1,
                            const uint nradial = 1)
            :
            IPolarIntegrate(spt, orders, {nang, nradial}) {}
            
            
            /// Constructor with specified orders and offset and location
            /// of singularity in the parametric interval [-1,1]^2
            IPolarIntegrate(const GPt2D& spt,
                            const UIntVec& orders,
                            const UIntVec& nsub,
                            const uint offset = 0)
            :
            mSPt(spt),
            mCurrentSubCellI(0),
            mGLIntegrator(orders, offset),
            mSubElIntegrator(nsub) { initPolarTerms(); }
            
            /// Check if iterator is finished
            bool isDone() const { return 4 == mCurrentSubCellI; }
            
            /// Get the current quadrature point
            GPt2D get() const { return mCurrentQPt; }
            
            /// Get the quadrature point component
            double get(const uint comp) const { return get().get(comp); }
            
            /// Get the current quadrature weight
            double getWeight() const { return mCurrentWeight; }
            
            /// Restart the iterator
            void restart()
            {
                mGLIntegrator.restart();
                mCurrentSubCellI = 0;
                initPolarTerms();
            }
            
            /// Prefix increment
            IPolarIntegrate& operator++()
            {
                incrementImpl();
                return *this;
            }
            
            /// Subcell index getter
            uint currentSubCellI() const { return mCurrentSubCellI; }
            
            /// Sub cell index of subcells within current triangular domain
            uint currentSubSubCellI() const { return mSubElIntegrator.currentIndex(); }
            
            /// Source point getter
            const GPt2D& sourcePt() const { return mSPt; }
            
            /// Update the source point and restart iterator
            void updateSrcAndRestart(const GPt2D s)
            {
                mSPt = s;
                restart();
            }
            
            /// Get current inner quadrature point
            const GPt2D currentInnerPt() const { return mGLIntegrator.get(); }
            
            /// Get current inner weight
            const double currentInnerWt() const { return mGLIntegrator.getWeight(); }
            
        private:
            
            void incrementImpl()
            {
                ++mSubElIntegrator;
                if(mSubElIntegrator.isDone()) { // reached end of subcells within triangle
                    mSubElIntegrator.restart();
                    ++mGLIntegrator;
                }
                if(mGLIntegrator.isDone()) // reached end of quadrature points on triangle
                    restartNextSubCell();
                
                if(!isDone()) // prevent call to non-real triangle when finished
                    initPolarTerms();
            }
            
            /// A wrapper function for calculating the current point and weight.
            /// Takes account of the source point lieing on edges.
            void initPolarTerms();
            
            /// Compute the current quadrature point and weight and store
            void computeCurrentPtWt();
            
            /// Move to next subcell
            void restartNextSubCell()
            {
                //mSubElIntegrator.restart();
                mGLIntegrator.restart();
                ++mCurrentSubCellI;
            }
            
            /// Subcell integrator
            const ISubElem& subElem() const { return mSubElIntegrator; }
            
            /// Location of singularity
            GPt2D mSPt;
            
            /// Current subcell index 0,1,2,3
            uint mCurrentSubCellI;
            
            /// Member instance of G-L integrator used in each subcell
            IElemIntegrate mGLIntegrator;
            
            /// Member instance of subelement integrator
            ISubElem mSubElIntegrator;
            
            /// Current quadrature point. Calculated in initPolarTerms()
            GPt2D mCurrentQPt;
            
            /// Current quadrature weight. Calculated in initPolarTerms()
            double mCurrentWeight;
            
        };
    }
}
#endif
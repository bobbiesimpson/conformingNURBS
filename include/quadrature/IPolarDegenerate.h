#ifndef IPOLAR_DEGENERATE_INTEGRATE_H
#define IPOLAR_DEGENERATE_INTEGRATE_H

#include "base.h"
#include "IElemIntegrate.h"
#include "ISubElem.h"
#include "IBaseIntegrate.h"

namespace nurbs {
    
    namespace elem {
        
        /// An iterator class for calculating integrals containing a 1/R singularity defined
        /// over a degenerate quadrilateral element in R^2.
        ///
        /// A tranformation to polar coordinates automatically cancels the 1/R singularity.
        /// The implementation also splits the two triangular subcells whose edges are perpendicular
        /// to the degenerate edge into two further subcells. This is done because on elements that are
        /// degenerate results in a line singularity rather than the point singularity found in
        /// non-degenerate elements.
        
        /// Furthermore, we apply a Telles transformation to eliminate the singularity on the edges of these
        /// sub (sub!) elements.
        ///
        
        class IPolarDegenerate : public IBaseIntegrate
        {
            
        public:
            
            /// Construct with equal number of quadrature points in each parametric
            /// direction. nrho specifies number of sub elements in radial
            /// direction.
            IPolarDegenerate(const GPt2D& spt,
                             const nurbs::Edge edge,
                             const uint n  = DEFAULT_NGPS)
            :
            IPolarDegenerate(spt, edge, {n,n})
            {}
            
            /// Construct with vector of integration orders
            IPolarDegenerate(const GPt2D& spt,
                             const nurbs::Edge edge,
                             const UIntVec& orders)
            :
            IPolarIntegrate(spt, edge, orders)
            {}
            
            /// Constructor with specified orders and offset and location
            /// of singularity in the parametric interval [-1,1]^2
            IPolarDegenerate(const GPt2D& spt,
                             const nurbs::Edge edge,
                             const UIntVec& orders,
                             const uint offset = 0)
            :
            mSPt(spt),
            mDegenerateEdge(edge),
            mCurrentSubCellI(0),
            mGLIntegrator(orders, offset),
            {
                initSubCellRanges();
                initPolarTerms();
            }
            
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
            
            /// Degenerate edge getter
            const nurbs::Edge degenerateEdge() const { return mDegenerateEdge; }
            
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
            
            /// Calculate appropriate subdivsion intervals given
            /// degenerate edge
            void initSubCellRanges();
            
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
            
            /// Subcell integrator for given triangle index
            const ISubElem& subElem(const uint i) { return mSubElemVec.at(i); }
            
            /// Location of singularity
            GPt2D mSPt;
            
            /// The degenerate edge enumeration
            const nurbs::Edge mDegenerateEdge;
            
            /// Current subcell index 0,1,2,3
            uint mCurrentSubCellI;
            
            /// Member instance of G-L integrator used in each subcell
            IElemIntegrate mGLIntegrator;
            
            /// Member instance of subelement integrator
            ISubElem mSubElIntegrator;
            
            /// Vector of subelements
            std::vector<ISubElem> mSubElemVec;
            
            /// Current quadrature point. Calculated in initPolarTerms()
            GPt2D mCurrentQPt;
            
            /// Current quadrature weight. Calculated in initPolarTerms()
            double mCurrentWeight;
            
        };
    }
}
#endif
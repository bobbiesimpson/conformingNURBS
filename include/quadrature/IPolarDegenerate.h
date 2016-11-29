#ifndef IPOLAR_DEGENERATE_INTEGRATE_H
#define IPOLAR_DEGENERATE_INTEGRATE_H

#include "base.h"
#include "IElemIntegrate.h"
#include "ISubElem.h"
#include "IBaseIntegrate.h"
#include "ITellesIntegrate.h"

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
        
        /// Enumeration which specifies if subcell division applies
        /// Telles transformation to the internal edge or the external
        /// edges
        enum class DivisionType
        {
            INTERNAL,
            EXTERNAL
        };
        
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
            
            /// Constructor with specified orders and offset and location
            /// of singularity in the parametric interval [-1,1]^2
            IPolarDegenerate(const GPt2D& spt,
                             const nurbs::Edge edge,
                             const UIntVec& orders)
            :
            mSPt(spt),
            mDegenerateEdge(edge),
            mOrders({2 * orders[0], orders[1]}),
            mCurrentSubCellI(0)
            {
                init();
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
                for(auto& tellesvec : mTellesIntegratorVec)
                    for(auto& i : tellesvec)
                        i.restart();
                
                mCurrentSubCellI = 0;
                initPolarTerms();
            }
            
            /// Prefix increment
            IPolarDegenerate& operator++()
            {
                incrementImpl();
                return *this;
            }
            
            /// Subcell index getter
            uint currentSubCellI() const { return mCurrentSubCellI; }
            
            /// Sub cell index of subcells within current triangular domain
            uint currentSubSubCellI() const { return subElem().currentIndex(); }
            
            /// Number of triangular subcells (always 4)
            const uint subCellN() const { return 4; }
            
            /// Number of subsubcells in the given subcell
            const uint subsubCellN(const uint isubcell) const { return subElem(isubcell).subCellN(); }
            
            /// Source point getter
            const GPt2D& sourcePt() const { return mSPt; }
            
            /// Update the source point and restart iterator
            void updateSrcAndRestart(const GPt2D s)
            {
                mSPt = s;
                restart();
            }
            
            /// Get current inner quadrature point
            const GPt2D currentInnerPt() const { return /*mGLIntegrator.get();*/  baseIntegrator().get(); }
            
            /// Get current inner weight
            const double currentInnerWt() const { return /*mGLIntegrator.getWeight(); */ baseIntegrator().getWeight(); }
            
            /// Degenerate edge getter
            const nurbs::Edge degenerateEdge() const { return mDegenerateEdge; }
            
            /// Orders getter
            const UIntVec& orders() const { return mOrders; }
            
            /// Base integrator const getter
            const ITellesIntegrate& baseIntegrator() const { return mTellesIntegratorVec.at(currentSubCellI()).at(currentSubSubCellI()); }
            
            /// Base integrator non-const getter
            ITellesIntegrate& baseIntegrator() { return mTellesIntegratorVec[currentSubCellI()][currentSubSubCellI()]; }
            
        private:
            
            void incrementImpl()
            {
                ++baseIntegrator();
                if(baseIntegrator().isDone())
                {
                    baseIntegrator().restart();
                    ++subElem();
                }
                
                if(subElem().isDone())
                {
                    restartNextSubCell();
                    initPolarTerms();
                }
                
                if(!isDone()) // prevent call to non-real triangle when finished
                    initPolarTerms();
                
                
            }
            
            /// Initiate theta values and appropriate subcell divisions
            void init();
            
            /// A wrapper function for calculating the current point and weight.
            /// Takes account of the source point lieing on edges.
            void initPolarTerms();
            
            /// Compute the current quadrature point and weight and store
            void computeCurrentPtWt();
            
            /// Move to next subcell
            void restartNextSubCell()
            {
                //baseIntegrator().restart();
                //while(subCellHasZeroArea(currentSubCellI()))
                ++mCurrentSubCellI;
            }
            
            /// Subcell integrator
            ISubElem& subElem() { return subElem(currentSubCellI()); }
            
            /// Subcell integrator
            const ISubElem& subElem() const { return subElem(currentSubCellI()); }
            
            /// Subcell integrator for given triangle index
            ISubElem& subElem(const uint i) { return mSubElemVec[i]; }
            
            /// Subcell integrator for given triangle index
            const ISubElem& subElem(const uint i) const { return mSubElemVec.at(i); }
            

            
            // Get the range [\theta_1, \theta_2] that defines the triangular subcell
            const std::pair<double, double>& thetaRange(const uint isubcell) const
            {
                return mThetaVec.at(isubcell);
            }
            
            const double thetaMap(const double xi, const uint isubcell) const
            {
                const auto theta_pair = thetaRange(isubcell);
                const auto t1 = theta_pair.first;
                const auto t2 = theta_pair.second;
                
                return (t2 - t1) * 0.5 * xi + (t2 + t1) * 0.5;
            }
            
            // Get xi given theta in a given subcell
            const double inverseThetaMap(const double theta,
                                         const uint isubcell) const
            {
                const auto t_range = thetaRange(isubcell);
                const double t1 = t_range.first;
                const double t2 = t_range.second;
                return (theta - 0.5 * (t1 + t2)) / ((t2 - t1) * 0.5);
            }
            
            const bool subCellHasZeroArea(const uint isubcell) const
            {
                switch(isubcell)
                {
                    case 0:
                        if(approximatelyEqual(sourcePt().get(0), 1.0, TOL))
                            return true;
                        break;
                    case 1:
                        if(approximatelyEqual(sourcePt().get(1), 1.0, TOL))
                            return true;
                        break;
                    case 2:
                        if(approximatelyEqual(sourcePt().get(0), -1.0, TOL))
                            return true;
                        break;
                    case 3:
                        if(approximatelyEqual(sourcePt().get(1), -1.0, TOL))
                            return true;
                        break;
                    default:
                        throw std::runtime_error("Bad subcell in function 'subcellHasZeroArea()'");
                }
                return false;
            }
            
            /// Divide a triangular subcell such that the new boundary
            /// is perpindicular to the parent element edges
            void divideSubcell(const PolarSubCell pcell,
                               const DivisionType dtype);
            
            /// Location of singularity
            GPt2D mSPt;
            
            /// The degenerate edge enumeration
            const nurbs::Edge mDegenerateEdge;
            
            /// Quadrature orders
            const UIntVec mOrders;
            
            /// Current subcell index 0,1,2,3
            uint mCurrentSubCellI;
            
            /// Vector of subelements
            std::vector<ISubElem> mSubElemVec;
            
            /// Vector of Telles integrators (one for each sub-sub cell)
            std::vector<std::vector<ITellesIntegrate>> mTellesIntegratorVec;
            
            /// Values of theta that define the four triangles centred around the source point
            std::vector<std::pair<double, double>> mThetaVec;
            
            /// Current quadrature point. Calculated in initPolarTerms()
            GPt2D mCurrentQPt;
            
            /// Current quadrature weight. Calculated in initPolarTerms()
            double mCurrentWeight;
            
        };
    }
}
#endif
#ifndef NURBS_ITELLES_INTEGRATE_H
#define NURBS_ITELLES_INTEGRATE_H

#include "base.h"
#include "IElemIntegrate.h"

namespace nurbs
{
    namespace elem
    {
        class ITellesIntegrate : public IBaseIntegrate
        {
            
        public:
            
            /// Construct with equal orders in each direction
            ITellesIntegrate(const Edge e,
                             uint ngps = DEFAULT_NGPS)
            :
            mSingularEdge(e),
            mGLIntegrator(ngps) { init(); }
            
            /// construct with vector of orders and optional offset
            ITellesIntegrate(const Edge e,
                             const UIntVec& orders,
                             const uint offset = 0)
            :
            mSingularEdge(e),
            mGLIntegrator(orders, offset) { init(); }
            
            /// Construct with no transformation applied
            ITellesIntegrate(uint ngps = DEFAULT_NGPS)
            :
            mSingularEdge(Edge::EDGE0),
            mGLIntegrator(ngps) { mEtaVals = std::make_pair(0.0, 0.0); }
            
            /// Construct with no transformation applied
            ITellesIntegrate(const UIntVec& orders,
                             const uint offset = 0)
            :
            mSingularEdge(Edge::EDGE0),
            mGLIntegrator(orders, offset) { mEtaVals = std::make_pair(0.0, 0.0); }

            /// Is integrator finished?
            bool isDone() const { return mGLIntegrator.isDone(); }
            
            /// Get the gauss point
            GPt2D get() const;
            
            /// Get the gauss point component
            double get( uint component ) const;
            
            /// Get the gauss weight
            double getWeight() const;
            
            /// Restart the iterator to point to the first gauss point
            void restart()
            { mGLIntegrator.restart(); }
            
            /// Return the total number of gauss points
            uint pointN() const
            { return mGLIntegrator.pointN(); }
            
            /// Current gauss point index
            uint currentIndex() const
            { return mGLIntegrator.currentIndex(); }
            
            /// Increment iterator
            ITellesIntegrate& operator++() {  incrementImpl(); return *this; };
            
            
        private:
            
            /// Initialise data structures
            void init();
            
            /// Singular edge getter
            const Edge singularEdge() const { return mSingularEdge; }
            
            const double etaBar1() const { return mEtaVals.first; }
            
            const double etaBar2() const { return mEtaVals.second; }
            
            /// Increment implementation
            void incrementImpl() { ++mGLIntegrator; }
            
            /// The edge containing the singularity
            const Edge mSingularEdge;
            
            /// Eta values for each parent coordinate for Telles transformation
            /// Value of zero indicates no transformation applied
            std::pair<double, double> mEtaVals;
            
            /// Gauss-Legendre integrator
            IElemIntegrate mGLIntegrator;
            
            
            
        };
        
        double tellesTransform(const double xi, const double etabar);
        
        double tellesJacob(const double xi, const double etabar);
    }
}

#endif

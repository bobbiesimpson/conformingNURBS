#ifndef NURBS_ISUB_ELEM_INTEGRATE_H
#define NURBS_ISUB_ELEM_INTEGRATE_H

#include "base.h"
#include "IElemIntegrate.h"
#include "IBaseIntegrate.h"
#include "ISubElem.h"

namespace nurbs {
    
    class ISubElemIntegrate : public IBaseIntegrate
    {
        public:
        
        /// Constructor
        ISubElemIntegrate(const uint n = DEFAULT_NGPS,
                          const uint nsubs = 1,
                          const uint nsubt = 1)
        :
        ISubElemIntegrate({n, n}, {nsubs, nsubt}) {}
        
        /// Constructor
        ISubElemIntegrate(const UIntVec& orders,
                          const uint nsubs = 1,
                          const uint nsubt = 1)
        :
        ISubElemIntegrate(orders, {nsubs, nsubt}) {}
        
        /// Constructor
        ISubElemIntegrate(const UIntVec& orders,
                          const UIntVec& nsub,
                          const uint offset = 0)
        :
        mGLIntegrator(orders, offset),
        mSubElIntegrator(nsub) {}
        
        /// Check if iterator is finished
        bool isDone() const
        {
            return mGLIntegrator.isDone();
        }
        
        /// Get the current quadrature point
        GPt2D get() const
        {
            return mSubElIntegrator.get(currentInnerPt());
        }
        
        /// Get the current quadrature weight
        double getWeight() const
        {
            return currentInnerWt() * subElem().jacob();
        }
        
        /// Restart the iterator
        void restart()
        {
            mGLIntegrator.restart();
            mSubElIntegrator.restart();
        }
        
        /// Prefix increment
        ISubElemIntegrate& operator++()
        {
            incrementImpl();
            return *this;
        }
        
        private:
        
        void incrementImpl()
        {
            ++mSubElIntegrator;
            if(mSubElIntegrator.isDone()) { // reached end of subcells within triangle
                mSubElIntegrator.restart();
                ++mGLIntegrator;
            }
        }
        
        /// Get current inner quadrature point
        const GPt2D currentInnerPt() const { return mGLIntegrator.get(); }
        
        /// Get current inner weight
        const double currentInnerWt() const { return mGLIntegrator.getWeight(); }
        
        /// Subcell integrator
        const ISubElem& subElem() const { return mSubElIntegrator; }
        
        /// Member instance of G-L integrator used in each subcell
        IElemIntegrate mGLIntegrator;
        
        /// Member instance of subelement integrator
        ISubElem mSubElIntegrator;
        
    };
    
}

#endif
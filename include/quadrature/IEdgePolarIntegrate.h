#ifndef IEDGE_POLAR_H
#define IEDGE_POLAR_H

#include "IGalerkinIntegrate.h"
#include "IElemIntegrate.h"
#include "IPolarIntegrate.h"
#include "base.h"
#include "MultiForest.h"

namespace nurbs {
    
    /// A Galerkin quadrature iterator that uses a polar transformation
    /// to remove the singularity.
    
    class IEdgePolarIntegrate : public IGalerkinIntegrate {
      
    public:
        
        /// Construct with parent coordinate of polar point in field element,
        /// vectors of quadrature orders for source and field elements,
        /// and number of subcells in angular and radial direction resp.
        IEdgePolarIntegrate(const UIntVec& srcorders,
                             const UIntVec& fieldorders,
                             const UIntVec& ncells,
                             const Edge es,
                             const Edge ef)
        :
        IGalerkinIntegrate(),
        mSIntegrate(srcorders),
        mFIntegrate(projectPt(mSIntegrate.get(), es, ef),
                    fieldorders,
                    ncells)
        {}
        
        void restart()
        {
            mSIntegrate.restart();
            mFIntegrate.restart();
            setCurrentSubCellI(0);
        }
        
    private:
        
        /// Increment implementation
        void incrementImpl();
        
        /// For the current set of quadrature indices and sub element,
        /// return the corresponding 4-d quadrature point
        GPt4D getImpl() const;
        
        /// For the current 4-d quadrature point return the quadrature weight
        /// including any necessary jacobian determinants.
        double getWeightImpl() const;
        
        /// 6 subcells for this transformation (see p.318 of Sauter and Schwab)
        uint subCellN() const { return 1; };
        
        /// Source element integrator
        IElemIntegrate mSIntegrate;
        
        /// Field element integrator
        IPolarIntegrate mFIntegrate;
    };
}

#endif
#include "IEdgePolarIntegrate.h"

namespace nurbs {
    
    void IEdgePolarIntegrate::incrementImpl()
    {
        ++mFIntegrate;
        if(mFIntegrate.isDone()) {
            ++mSIntegrate;
            mFIntegrate.updateSrcAndRestart(mSIntegrate.get());
            if(mSIntegrate.isDone())
                incrementSubCellI();
        }
    }
    
    GPt4D IEdgePolarIntegrate::getImpl() const
    {
        return GPt4D(mSIntegrate.get(), mFIntegrate.get());
    }
    
    double IEdgePolarIntegrate::getWeightImpl() const
    {
        return mSIntegrate.getWeight() * mFIntegrate.getWeight();
    }
}
#include "ITellesIntegrate.h"

namespace nurbs
{
    namespace elem
    {
        
        void ITellesIntegrate::init()
        {
            switch(singularEdge())
            {
                case Edge::EDGE0:
                    mEtaVals = std::make_pair(0.0, -1.0);
                    break;
                case Edge::EDGE1:
                    mEtaVals = std::make_pair(0.0, 1.0);
                    break;
                case Edge::EDGE2:
                    mEtaVals = std::make_pair(-1.0, 0.0);
                    break;
                case Edge::EDGE3:
                    mEtaVals = std::make_pair(0.0, 1.0);
                    break;
            }
        }
        
        
        GPt2D ITellesIntegrate::get() const
        {
            const auto base_pt = mGLIntegrator.get();
            return GPt2D(tellesTransform(base_pt.s, etaBar1()),
                         tellesTransform(base_pt.t, etaBar2()));
            
        }
        
        double ITellesIntegrate::get( uint component ) const
        {
            const GPt2D gpt = get();
            return gpt.get(component);
        }
        
        double ITellesIntegrate::getWeight() const
        {
            const double baseweight = mGLIntegrator.getWeight();
            const auto base_pt = mGLIntegrator.get();
            return baseweight * tellesJacob(base_pt.s, etaBar1()) * tellesJacob(base_pt.t, etaBar2());
        }
        
        double tellesTransform(const double xi, const double etabar)
        {
            return (1.0 - xi * xi) * etabar * 0.5 + xi;
        }
        
        double tellesJacob(const double xi, const double etabar)
        {
            return 1.0 - xi * etabar;
        }
        
    }
}
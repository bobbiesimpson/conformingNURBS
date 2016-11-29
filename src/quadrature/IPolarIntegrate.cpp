#include "IPolarIntegrate.h"
#include "base.h"

#include <stdexcept>
#include <cmath>

namespace nurbs {
    
    namespace elem {
        
        void IPolarIntegrate::initPolarTerms()
        {
            if(approximatelyEqual(sourcePt().get(0), 1.0, TOL) && 0 == currentSubCellI())
                restartNextSubCell();
            if(approximatelyEqual(sourcePt().get(1), 1.0, TOL) && 1 == currentSubCellI())
                restartNextSubCell();
            if(approximatelyEqual(sourcePt().get(0), -1.0, TOL) && 2 == currentSubCellI())
                restartNextSubCell();
            if(approximatelyEqual( sourcePt().get(1), -1.0, TOL) && 3 == currentSubCellI())
                restartNextSubCell();
            computeCurrentPtWt();
        }
        
        void IPolarIntegrate::computeCurrentPtWt()
        {
            if(isDone())
                return;
            double theta1, theta2, theta, rhohat;
            
            const GPt2D qpt_old = subElem().get(currentInnerPt());
            GPt2D qpt = qpt_old;
            const GPt2D old_inner_pt = currentInnerPt();
            
            double tweight = 1.0;
            
            /// Variable Telles transformation
//            if(3 == currentSubCellI() || currentSubCellI() == 1)
//            {
//                const double src = 0.0;
//                const double s_d = src * src - 1.0;
//                const double y_d = std::cbrt( src * s_d + std::abs( s_d ) )
//                + std::cbrt( src * s_d - std::abs( s_d ) ) + src;
//                tweight = 3.0 * (( qpt_old.s - y_d ) * ( qpt_old.s - y_d ) ) / ( 1.0 + 3.0 * y_d * y_d );
//                qpt.s = ( std::pow( ( qpt_old.s - y_d ), 3 ) + y_d * ( y_d * y_d + 3.0 ) ) /
//                ( 1.0 + 3.0 * y_d * y_d );
//                
//            }
//            
            
            // Simple Telles transformation
//            if(3 == currentSubCellI())
//            {
//                if(currentSubSubCellI() == 1)
//                {
//                    qpt = subElem().get(GPt2D((1.0 - old_inner_pt.s * old_inner_pt.s) * -0.5  + old_inner_pt.s,
//                                              old_inner_pt.t));
//                    tweight = (1.0 + old_inner_pt.s);
//                }
//                else if(currentSubSubCellI() == 0)
//                {
//                    qpt = subElem().get(GPt2D((1.0 - old_inner_pt.s * old_inner_pt.s) * 0.5  + old_inner_pt.s,
//                                old_inner_pt.t));
//                    tweight = (1.0 - old_inner_pt.s);
//                }
//            }
//            
//            else if(1 == currentSubCellI())
//            {
//                if(currentSubSubCellI() == 1)
//                {
//                    qpt = subElem().get(GPt2D((1.0 - old_inner_pt.s * old_inner_pt.s) * -0.5  + old_inner_pt.s,
//                                old_inner_pt.t));
//                    tweight = (1.0 + old_inner_pt.s);
//                }
//                
//                else if(currentSubSubCellI() == 0)
//                {
//                    qpt = subElem().get(GPt2D((1.0 - old_inner_pt.s * old_inner_pt.s) * 0.5  + old_inner_pt.s,
//                                old_inner_pt.t));
//                    tweight = (1.0 - old_inner_pt.s);
//                }
//            }
            
        
            switch(currentSubCellI())
            {
                    // first triangle
                case 0:
                    // If on edge xi=+1, omit this triangle
                    theta1 = -std::atan((1.0 + sourcePt().get(1)) / (1.0 - sourcePt().get(0)));
                    theta2 = std::atan((1.0 - sourcePt().get(1)) / (1.0 - sourcePt().get(0)));
                    theta = 0.5 * (theta2 - theta1) * qpt.get(0) + 0.5 * (theta2 + theta1);
                    rhohat = (1.0 - sourcePt().get(0)) / std::cos(theta);
                    break;
                    // second triangle
                case 1:
                    // If on edge eta = +1, omit this triangle
                    theta1 = PI / 2.0  - std::atan((1.0 - sourcePt().get(0)) / (1.0 - sourcePt().get(1)));
                    theta2 = PI / 2.0 + std::atan((sourcePt().get(0) + 1.0) / (1.0 - sourcePt().get(1)));
                    theta = 0.5 * (theta2 - theta1) * qpt.get(0) + 0.5 * (theta2 + theta1);
                    rhohat = (1.0 - sourcePt().get(1)) / std::sin(theta);
                    break;
                    // third triangle
                case 2:
                    theta1 = PI - std::atan((1.0 - sourcePt().get(1)) / (sourcePt().get(0) + 1.0));
                    theta2
                    = PI + std::atan((sourcePt().get(1) + 1.0) / (sourcePt().get(0) + 1.0));
                    theta = 0.5 * (theta2 - theta1) * qpt.get(0) + 0.5 * (theta2 + theta1);
                    rhohat = -(1.0 + sourcePt().get(0)) / std::cos(theta);
                    break;
                    // fourth triangle
                case 3:
                    // If on edge eta = -1, omit this triangle
                    theta1 = 3.0 * PI / 2.0 - std::atan((sourcePt().get(0) + 1.0) / (sourcePt().get(1) + 1.0));
                    theta2 = 3.0 * PI / 2.0 + std::atan((1.0 - sourcePt().get(0)) / (sourcePt().get(1) + 1.0));
                    theta = 0.5 * (theta2 - theta1) * qpt.get(0) + 0.5 * (theta2 + theta1);
                    rhohat = -(1.0 + sourcePt().get(1)) / std::sin(theta);
                    break;
                default:
                    throw std::runtime_error("Bad subcell");
            }
            const double rho = 0.5 * (qpt.get(1) + 1.0) * rhohat;
            const double polarjacob = 0.25 * (theta2 - theta1) * rhohat;
            // get coordinates in parent coordinate system
            mCurrentQPt.s = sourcePt().get(0) + rho * std::cos(theta);
            mCurrentQPt.t = sourcePt().get(1) + rho * std::sin(theta);
            mCurrentWeight = rho * polarjacob * /*currentInnerWt() * */subElem().jacob() * tweight;
            
        }
        
    }
}
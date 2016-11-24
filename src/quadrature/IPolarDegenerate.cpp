#include "IPolarDegenerate.h"
#include "base.h"

#include <stdexcept>
#include <cmath>

namespace nurbs
{
    namespace elem
    {
        
        void IPolarDegenerate::init()
        {
            // initate theta values
            
            // Triangle 0 \equiv 'East' triangle
            mThetaVec.push_back(std::make_pair(-std::atan((1.0 + sourcePt().get(1)) / (1.0 - sourcePt().get(0))),
                                               std::atan((1.0 - sourcePt().get(1)) / (1.0 - sourcePt().get(0)))));
            
            // Triangle 1 \equiv 'North' triangle
            mThetaVec.push_back(std::make_pair(PI / 2.0  - std::atan((1.0 - sourcePt().get(0)) / (1.0 - sourcePt().get(1))),
                                               PI / 2.0 + std::atan((sourcePt().get(0) + 1.0) / (1.0 - sourcePt().get(1)))));

            // Triangle 2 \equiv 'West' triangle
            mThetaVec.push_back(std::make_pair(PI - std::atan((1.0 - sourcePt().get(1)) / (sourcePt().get(0) + 1.0)),
                                               PI + std::atan((sourcePt().get(1) + 1.0) / (sourcePt().get(0) + 1.0))));

            // Triangle 3 \equiv 'South' triangle
            mThetaVec.push_back(std::make_pair(3.0 * PI / 2.0 - std::atan((sourcePt().get(0) + 1.0) / (sourcePt().get(1) + 1.0)),
                                               3.0 * PI / 2.0 + std::atan((1.0 - sourcePt().get(0)) / (sourcePt().get(1) + 1.0))));
            
            // Basic idea is to split triangles into subcells such that
            // the division lies along the line singularity created by the degenerate mapping
            
            // First fill up subcells with default 1 subcell per triangle
            mSubElemVec.clear();
            for(size_t i = 0; i < 4; ++i)
                mSubElemVec.push_back(ISubElem(1,1));
            
            // Compute
            if(Edge::EDGE0 == degenerateEdge() || Edge::EDGE1 == degenerateEdge())
            {
                // 'West' triangle
                const uint west = 2;
                if(!subCellHasZeroArea(west))
                {
                    const DoubleVec xivec_west{inverseThetaMap(PI, west)};
                    mSubElemVec[west] = ISubElem(xivec_west, {});
                }
                
                // 'east' triangle
                const uint east = 0;
                if(!subCellHasZeroArea(east))
                {
                    const DoubleVec xivec_east{inverseThetaMap(0.0, east)};
                    mSubElemVec[east] = ISubElem(xivec_east, {});
                }
            }
            else if(Edge::EDGE2 == degenerateEdge() || Edge::EDGE3 == degenerateEdge())
            {
                // 'South' triangle
                const uint south = 3;
                if(!subCellHasZeroArea(south))
                {
                    const DoubleVec xivec_south{inverseThetaMap(3*PI/2.0, south)};
                    mSubElemVec[south] = ISubElem(xivec_south, {});
                }
                
                // 'North' triangle
                const uint north = 1;
                if(!subCellHasZeroArea(north))
                {
                    const DoubleVec xivec_north{inverseThetaMap(PI/2.0, north)};
                    mSubElemVec[north] = ISubElem(xivec_north, {});
                }
            }
        }
        
        void IPolarDegenerate::initPolarTerms()
        {
            if(subCellHasZeroArea(currentSubCellI()))
                restartNextSubCell();
            computeCurrentPtWt();
        }
        
        void IPolarDegenerate::computeCurrentPtWt()
        {
            if(isDone())
                return;
            double rhohat = 0.0;
            
            const GPt2D qpt_old = subElem().get(currentInnerPt());
            GPt2D qpt = qpt_old;
            double tweight = 1.0;

            const double theta = thetaMap(qpt.get(0), currentSubCellI());
            const auto theta_range = thetaRange(currentSubCellI());
            const double theta1 = theta_range.first;
            const double theta2 = theta_range.second;
            
            // now compute Telles transformation
            // TODO
            
            switch(currentSubCellI())
            {
                    // first triangle
                case 0:
                    // If on edge xi=+1, omit this triangle
                    rhohat = (1.0 - sourcePt().get(0)) / std::cos(theta);
                    break;
                    // second triangle
                case 1:
                    // If on edge eta = +1, omit this triangle
                    rhohat = (1.0 - sourcePt().get(1)) / std::sin(theta);
                    break;
                    // third triangle
                case 2:
                    rhohat = -(1.0 + sourcePt().get(0)) / std::cos(theta);
                    break;
                    // fourth triangle
                case 3:
                    // If on edge eta = -1, omit this triangle
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
            mCurrentWeight = rho * polarjacob * currentInnerWt() * subElem().jacob() * tweight;
        }
        
    }
}
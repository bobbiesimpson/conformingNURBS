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
            for(size_t i = 0; i < subCellN(); ++i)
                mSubElemVec.push_back(ISubElem(1,1));
            
            // and fill up integrators with default number of one per triangle
            mTellesIntegratorVec.clear();
            for(size_t i = 0; i < subCellN(); ++i)
                mTellesIntegratorVec.push_back({ITellesIntegrate(orders())});
            
            // Compute
            if(Edge::EDGE0 == degenerateEdge() || Edge::EDGE1 == degenerateEdge())
            {
                // 'West' triangle
                const uint west = 2;
                if(!subCellHasZeroArea(west))
                {
                    const DoubleVec xivec_west{inverseThetaMap(PI, west)};
                    mSubElemVec[west] = ISubElem(xivec_west, {});
                    mTellesIntegratorVec[west].clear();
                    mTellesIntegratorVec[west].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
                    mTellesIntegratorVec[west].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
                    
                }
                
                // 'east' triangle
                const uint east = 0;
                if(!subCellHasZeroArea(east))
                {
                    const DoubleVec xivec_east{inverseThetaMap(0.0, east)};
                    mSubElemVec[east] = ISubElem(xivec_east, {});
                    mTellesIntegratorVec[east].clear();
                    mTellesIntegratorVec[east].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
                    mTellesIntegratorVec[east].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
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
                    mTellesIntegratorVec[south].clear();
                    mTellesIntegratorVec[south].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
                    mTellesIntegratorVec[south].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
                }
                
                // 'North' triangle
                const uint north = 1;
                if(!subCellHasZeroArea(north))
                {
                    const DoubleVec xivec_north{inverseThetaMap(PI/2.0, north)};
                    mSubElemVec[north] = ISubElem(xivec_north, {});
                    mTellesIntegratorVec[north].clear();
                    mTellesIntegratorVec[north].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
                    mTellesIntegratorVec[north].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
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
            
            const GPt2D qpt = subElem().get(currentInnerPt());
            
            const double theta = thetaMap(qpt.get(0), currentSubCellI());
            const auto theta_range = thetaRange(currentSubCellI());
            const double theta1 = theta_range.first;
            const double theta2 = theta_range.second;
            
            
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
            mCurrentWeight = rho * polarjacob * currentInnerWt() * subElem().jacob();
        }
        
    }
}
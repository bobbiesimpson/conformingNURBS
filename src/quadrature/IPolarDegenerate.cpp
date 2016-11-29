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
                divideSubcell(PolarSubCell::WEST, DivisionType::INTERNAL);
                divideSubcell(PolarSubCell::EAST, DivisionType::INTERNAL);
                
                if(Edge::EDGE0 == degenerateEdge())
                    divideSubcell(PolarSubCell::NORTH, DivisionType::EXTERNAL);
                else
                    divideSubcell(PolarSubCell::SOUTH, DivisionType::EXTERNAL);
                
//                // 'West' triangle
//                const uint west = 2;
//                if(!subCellHasZeroArea(west))
//                {
//                    const DoubleVec xivec_west{inverseThetaMap(PI, west)};
//                    mSubElemVec[west] = ISubElem(xivec_west, {});
//                    mTellesIntegratorVec[west].clear();
//                    mTellesIntegratorVec[west].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    mTellesIntegratorVec[west].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                }
//                
//                // 'east' triangle
//                const uint east = 0;
//                if(!subCellHasZeroArea(east))
//                {
//                    const DoubleVec xivec_east{inverseThetaMap(0.0, east)};
//                    mSubElemVec[east] = ISubElem(xivec_east, {});
//                    mTellesIntegratorVec[east].clear();
//                    mTellesIntegratorVec[east].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    mTellesIntegratorVec[east].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                }
//                
//                // For points very close to the edge of the domain it is necessary to split triangles
//                // into subcells that lie opposite the degenerate edge
//                if(Edge::EDGE0 == degenerateEdge())
//                {
//                    // 'North' triangle
//                    const uint north = 1;
//                    if(!subCellHasZeroArea(north))
//                    {
//                        const DoubleVec xivec_north{inverseThetaMap(PI/2.0, north)};
//                        mSubElemVec[north] = ISubElem(xivec_north, {});
//                        mTellesIntegratorVec[north].clear();
//                        mTellesIntegratorVec[north].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                        mTellesIntegratorVec[north].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    }
//                }
//                else
//                {
//                    // 'South' triangle
//                    const uint south = 3;
//                    if(!subCellHasZeroArea(south))
//                    {
//                        const DoubleVec xivec_south{inverseThetaMap(3*PI/2.0, south)};
//                        mSubElemVec[south] = ISubElem(xivec_south, {});
//                        mTellesIntegratorVec[south].clear();
//                        mTellesIntegratorVec[south].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                        mTellesIntegratorVec[south].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    }
//                }
            }
            else if(Edge::EDGE2 == degenerateEdge() || Edge::EDGE3 == degenerateEdge())
            {
                divideSubcell(PolarSubCell::SOUTH, DivisionType::INTERNAL);
                divideSubcell(PolarSubCell::NORTH, DivisionType::INTERNAL);
                
                if(Edge::EDGE3 == degenerateEdge())
                    divideSubcell(PolarSubCell::WEST, DivisionType::EXTERNAL);
                else
                    divideSubcell(PolarSubCell::EAST, DivisionType::EXTERNAL);
                
//                // 'South' triangle
//                const uint south = 3;
//                if(!subCellHasZeroArea(south))
//                {
//                    const DoubleVec xivec_south{inverseThetaMap(3*PI/2.0, south)};
//                    mSubElemVec[south] = ISubElem(xivec_south, {});
//                    mTellesIntegratorVec[south].clear();
//                    mTellesIntegratorVec[south].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    mTellesIntegratorVec[south].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                }
//                
//                // 'North' triangle
//                const uint north = 1;
//                if(!subCellHasZeroArea(north))
//                {
//                    const DoubleVec xivec_north{inverseThetaMap(PI/2.0, north)};
//                    mSubElemVec[north] = ISubElem(xivec_north, {});
//                    mTellesIntegratorVec[north].clear();
//                    mTellesIntegratorVec[north].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    mTellesIntegratorVec[north].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                }
//                
//                if(Edge::EDGE3 == degenerateEdge())
//                {
//                    // 'West' triangle
//                    const uint west = 2;
//                    if(!subCellHasZeroArea(west))
//                    {
//                        const DoubleVec xivec_west{inverseThetaMap(PI, west)};
//                        mSubElemVec[west] = ISubElem(xivec_west, {});
//                        mTellesIntegratorVec[west].clear();
//                        mTellesIntegratorVec[west].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                        mTellesIntegratorVec[west].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    }
//                }
//                else
//                {
//                    // 'east' triangle
//                    const uint east = 0;
//                    if(!subCellHasZeroArea(east))
//                    {
//                        const DoubleVec xivec_east{inverseThetaMap(0.0, east)};
//                        mSubElemVec[east] = ISubElem(xivec_east, {});
//                        mTellesIntegratorVec[east].clear();
//                        mTellesIntegratorVec[east].push_back(ITellesIntegrate(Edge::EDGE2, orders()));
//                        mTellesIntegratorVec[east].push_back(ITellesIntegrate(Edge::EDGE3, orders()));
//                    }
//                }
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
        
        void IPolarDegenerate::divideSubcell(const PolarSubCell pcell,
                                             const DivisionType dtype)
        {
            // Idea is to encapsulate the division of triangular subcells
            // within this function.
            
            const double tol = 1.0e-9;
            
            std::vector<Edge> edge_vec;
            if(DivisionType::INTERNAL == dtype)
                edge_vec = {Edge::EDGE3, Edge::EDGE2};
            else
                edge_vec = {Edge::EDGE2, Edge::EDGE3};
            
            // First naively create subcells
            double xi;
            uint ctype;
            switch(pcell)
            {
                case PolarSubCell::EAST:
                    ctype = 0;
                    xi = inverseThetaMap(0.0, ctype);
                    break;
                case PolarSubCell::NORTH:
                    ctype = 1;
                    xi = inverseThetaMap(PI/2.0, ctype);
                    break;
                case PolarSubCell::WEST:
                    ctype = 2;
                    xi = inverseThetaMap(PI, ctype);
                    break;
                case PolarSubCell::SOUTH:
                    ctype = 3;
                    xi = inverseThetaMap(3*PI/2.0, ctype);
                    break;
            }
            
            if(subCellHasZeroArea(ctype))
                return;
            
            mTellesIntegratorVec[ctype].clear();
            
            // We only need one sub-subcell
            if(essentiallyEqual(xi, 1.0, tol))
            {
                mSubElemVec[ctype] = ISubElem(1,1);
                mTellesIntegratorVec[ctype].push_back(ITellesIntegrate(edge_vec[0], orders()));
                
            }
            // We only need one sub-subcell
            else if(essentiallyEqual(xi, -1.0, tol))
            {
                mSubElemVec[ctype] = ISubElem(1,1);
                mTellesIntegratorVec[ctype].push_back(ITellesIntegrate(edge_vec[1], orders()));
                
            }
            // We need two sub-subcells
            else
            {
                mSubElemVec[ctype] = ISubElem(DoubleVec{xi}, {});
                for(size_t i = 0; i < edge_vec.size(); ++i)
                    mTellesIntegratorVec[ctype].push_back(ITellesIntegrate(edge_vec[i], orders()));
            }
        }
    }
}
#include <iostream>
#include <stdexcept>
#include <iomanip>

#include "Geometry.h"
#include "OutputVTK.h"
#include "IElemIntegrate.h"
#include "IWangQuadrature.h"
#include "IEqualQuadrature.h"
#include "HConformingForest.h"
#include "IRegularQuadrature.h"
#include "ISubElem.h"


using namespace nurbs;

int main(const int argc, const char* argv[])
{
    try {
        
        std::cout << "Trying to open hbs input file....\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            error("Cannot open file for reading\n");
        Geometry g;
        if(!g.loadHBSFile(ifs))
            error("Failed to load geometry from hbs data");
        
        Forest forest(g);
//        HDivForest mforest(g);
        
        double integral = 0.0;

//        for(uint iorder = 5; iorder < 31; iorder += 5)
//        {
//            const UIntVec o_order{iorder,iorder};
//            const UIntVec i_order{iorder - 1, iorder - 1};
//            
//            integral = 0.0;
//            
//            const uint ielem = 8;
//            
//            const auto p_sel = mforest.bezierElement(ielem);
//            const auto p_fel = p_sel;
//            
//            for(IElemIntegrate igpt_outer(o_order); !igpt_outer.isDone(); ++igpt_outer)
//            {
//                const auto sparent = igpt_outer.get();
//                const auto ws = igpt_outer.getWeight();
//                
//                const auto& t1_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
//                const auto& t2_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
//                const auto& x = p_sel->eval(sparent);
//                //            const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1_s, t2_s);
//                //            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
//                //            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
//                const auto& jdet_s = p_sel->jacDet(sparent, t1_s, t2_s);
//                //            const auto& jpiola_s = nurbs::cross(t1_s, t2_s).length();
//                
//                for(IElemIntegrate igpt_inner(i_order); !igpt_inner.isDone(); ++igpt_inner)
//                {
//                    const auto fparent = igpt_inner.get();
//                    const auto wf = igpt_outer.getWeight();
//                    
//                    // field point terms
//                    const auto& t1_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
//                    const auto& t2_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
//                    const auto y = p_fel->eval(fparent);
//                    //                const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1_f, t2_f);
//                    //                const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
//                    //                const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
//                    const double jdet_f = p_fel->jacDet(fparent, t1_f, t2_f);
//                    //                const double jpiola_f = nurbs::cross(t1_f, t2_f).length();
//                    
//                    const auto R = dist(x,y);
//                    const auto r = std::sqrt((sparent.s - fparent.s) * (sparent.s - fparent.s) +
//                                             (sparent.t - fparent.t) * (sparent.t - fparent.t));
//                    
//                    integral += r / R * jdet_s * jdet_f * ws * wf;
//                    
//                }
//                
//            }
//            
//            std::cout << std::setprecision(15) << integral << "\n";
//        }

        const auto el = forest.bezierElement(0);
        
        for(uint iorder = 2; iorder < 9; iorder += 2)
        {
            integral = 0.0;
            
            for(IEqualQuadrature igpt({iorder+1, iorder+1}, {iorder, iorder}); !igpt.isDone(); ++igpt)
            {
                const auto gpt = igpt.get();
                const auto s_p = gpt.srcPt();
                const auto x_p = gpt.fieldPt();
                const auto w = igpt.getWeight();
                
                const auto s = el->eval(s_p);
                const auto x = el->eval(x_p);
                const auto R = dist(s,x);
                const auto r = std::sqrt((s_p.s - x_p.s) * (s_p.s - x_p.s) +
                                         (s_p.t - x_p.t) * (s_p.t - x_p.t));
                integral += /*r/R*/ 1./R * el->jacDet(s_p) * el->jacDet(x_p) * w;
            }
            std::cout << "Suater singular integral " << integral << "\n";
        }

        

        
        // Wang quadrature
//        const uint worder = 2;
//        const DoubleVec scoords{0.5773502692, -0.5773502692};
//        const DoubleVec fcoords{0.5, -0.5};
//        const DoubleVec weight{3.3758170913,0.9577571422, 0.9577571422, 0.6550878207,
//                                0.9577571422,3.3758170913, 0.6550878207, 0.9577571422,
//                                0.9577571422, 0.6550878207, 3.3758170913,0.9577571422,
//                                0.6550878207, 0.9577571422, 0.9577571422,3.3758170913};
        
        const uint worder = 3;
        const DoubleVec scoords{0.7745966692, 0.0, -0.7745966692};
        const DoubleVec fcoords{0.7071067812, 0.0, -0.7071067812};
        const DoubleVec corner_vec{0.6969450722, 0.1633313861, 0.0871652173,
                                        0.1633313861, 0.1420548753, 0.0809251162,
                                        0.0871652173, 0.0809251162, 0.0649938701};
        const DoubleVec edge_vec{0.3528197212, 1.1320652889, 0.3528197212,
                                    0.2193171508, 0.2855264249, 0.2193171508,
                                    0.1275419569, 0.1434352976, 0.1275419569};
        const DoubleVec centre_vec{0.3592529166, 0.5977318272, 0.3592529166,
                                        0.5977318272, 1.8488501076, 0.5977318272,
                                        0.3592529166, 0.5977318272, 0.3592529166};
        
        const UIntVecVec corner_indices{ {0,1,2,3,4,5,6,7,8},
            {2,1,0,5,4,3,8,7,6},
            {6,7,8,3,4,5,0,1,2},
            {8,7,6,5,4,3,2,1,0}};
        
        const UIntVecVec edge_indices{ {0,1,2,3,4,5,6,7,8},
            {2,5,8,1,4,7,0,3,6},
            {6,3,0,7,4,1,8,5,2},
            {6,7,8, 3,4,5,0,1,2}};
        
        DoubleVec weight;
        
        for(const auto& v : corner_indices[0])
            weight.push_back(corner_vec[v]);
        
        for(const auto& v : edge_indices[0])
            weight.push_back(edge_vec[v]);
        
        for(const auto& v : corner_indices[1])
            weight.push_back(corner_vec[v]);
        
        for(const auto& v : edge_indices[1])
            weight.push_back(edge_vec[v]);
        
        for(const auto& v : centre_vec)
            weight.push_back(v);
        
        for(const auto& v : edge_indices[2])
            weight.push_back(edge_vec[v]);
        
        for(const auto& v : corner_indices[2])
            weight.push_back(corner_vec[v]);
        
        for(const auto& v : edge_indices[3])
            weight.push_back(edge_vec[v]);
        
        for(const auto& v : corner_indices[3])
            weight.push_back(corner_vec[v]);

//        const auto el = forest.bezierElement(20);
        const uint maxcells = 10;
        
        for(uint nsubcell = 1; nsubcell < maxcells; nsubcell+=2)
        {
            integral = 0.0;
            
            for(uint ixi = 0 ; ixi < worder; ++ixi)
                for(uint ieta = 0; ieta < worder; ++ieta)
                {
                    for(ISubElem isubelsource(nsubcell, nsubcell); !isubelsource.isDone(); ++isubelsource)
                    {
                        const GPt2D s_raw(scoords[ixi], scoords[ieta]);
                        const GPt2D s = isubelsource.get(s_raw);
                        
                        const auto s_phy = el->eval(s);
                        
                        for(uint iu = 0; iu < worder; ++iu)
                            for(uint iv = 0; iv < worder; ++iv)
                            {
                                const uint wcount = ixi * worder * worder * worder + ieta * worder * worder + iu * worder + iv;
                                for(ISubElem isubelfield(nsubcell, nsubcell); !isubelfield.isDone(); ++isubelfield)
                                {
                                    
                                    const GPt2D x_raw(fcoords[iu], fcoords[iv]);
                                    const auto x = isubelfield.get(x_raw);
                                    
                                    const auto x_phy = el->eval(x);
                                    
                                    const auto R = dist(s_phy, x_phy);
                                    
                                    const double r = sqrt((s.get(0)- x.get(0)) * (s.get(0)- x.get(0))
                                                          + (s.get(1)- x.get(1)) * (s.get(1)- x.get(1)));
                                    
                                    
                                    const double w = weight[wcount] * isubelsource.jacob() * isubelfield.jacob();
                                    
                                    if(essentiallyEqual(r,0.0, 1.0e-10))
                                        continue;
                                        
                                    integral += 1. * r / R * w * el->jacDet(s) * el->jacDet(x);
                                }
                            }
                    }
                }
            
            std::cout << "integral = " << std::setprecision(15) << integral << "\n";
        }
//        OutputVTK output("wang-test");
//        output.outputGeometry(forest);
        
    }
    catch(const std::exception& e) {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
    catch(...) {
        std::cerr << "Unknown exception occured\n";
        return EXIT_FAILURE;
    }
}

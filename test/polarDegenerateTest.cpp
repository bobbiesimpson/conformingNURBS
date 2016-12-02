#include <iostream>
#include <fstream>
#include <iomanip>

#include "Forest.h"
#include "Geometry.h"
#include "BezierNodalElement.h"
#include "IElemIntegrate.h"
#include "OutputVTK.h"
#include "HConformingForest.h"
#include "IPolarIntegrate.h"
#include "IPolarDegenerate.h"
#include "ITellesIntegrate.h"
#include "IEqualQuadrature.h"
#include "IEdgeQuadrature.h"
#include "BezierVectorElement.h"
#include "IVertexQuadrature.h"

using namespace nurbs;

int main(int argc, char* argv[])
{
    // Some constants
    const uint max_order = 12;                      // max order of quadrature that we loop unilt
    const uint max_forder = 8;
    const double k = 1.0;                         // wavenumber for emag kernel
    const std::complex<double> iconst(0.0, 1.0);    // imaginary number
    
    try
    {
        std::cout << "Running degenerate polar quadrature test.....\n";
        
        std::cout << "Trying to open hbs input file....\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            error("Cannot open file for reading\n");
        Geometry g;
        if(!g.loadHBSFile(ifs))
            error("Failed to load geometry from hbs data");
        
        Forest forest(g);
        HDivForest divforest(g);
        
        nurbs::OutputVTK output("nasa_almond");
        output.outputGeometry(forest);
        
        std::cout << "contructed hdiv forest with " << divforest.elemN() << " elements\n";
        
        
        // Element #s 8, 17 degenerate.
        // Element #7 next to #8 (non degenerate)

        const auto p_fel = divforest.bezierElement(8);
        const auto p_sel = divforest.bezierElement(26);
        
//        nurbs::VAnalysisElement* p_sel;
//        Edge e1, e2;
//        for(uint i = 0; i < divforest.elemN(); ++i)
//        {
//            p_sel = divforest.bezierElement(i);
//            const auto p = p_fel->degenerateEdge();
//            if(p.first)
//                continue;
//            if(nurbs::edgeConnected(*p_sel, *p_fel, e1, e2))
//                break;
//        }

        
        Edge e1, e2;
        if(nurbs::edgeConnected(*p_sel, *p_fel, e1, e2))
            std::cout << "Elements are edge connected\n";
        
        Vertex v1,v2;
        if(nurbs::vertexConnected(*p_sel, *p_fel, v1, v2))
            std::cout << "Elements are vertex connected\n";
        v1 = Vertex::VERTEX1;
        v2 = Vertex::VERTEX1;
            
        
        const auto degenerate_field_pair = p_fel->degenerateEdge();
        if(degenerate_field_pair.first)
            std::cout << "Field element is degenerate\n";
        
        const auto degenerate_source_pair = p_sel->degenerateEdge();
        if(degenerate_source_pair.first)
            std::cout << "Source element is degenerate\n";
        
        
        // element connectivity
        const auto fconn = p_fel->signedGlobalBasisFuncI();
        const auto sconn = p_sel->signedGlobalBasisFuncI();
        
        for(uint iforder = 2; iforder < max_forder + 1; ++iforder)
        {
            const auto forder = UIntVec{iforder,iforder};//p_fel->equalIntegrationOrder();
            //                    forder[0] *= 2;
            
            std::cout << "Field element quadrature order: " << forder << "\n";
            
            std::ofstream ofs("quadrature_test_" + std::to_string(iforder) + ".dat");
            if(!ofs)
                throw std::runtime_error("Cannot open file for writing.");
            
            for(uint order = 2; order < max_order; ++order)
            {
                const UIntVec sorder{order, order};
                
                std::cout << "Outer integration order: " << sorder << "....\n";
                
                uint ngpts = 0;                                 // number of quadrature points
                
                std::vector<std::vector<std::complex<double>>>  matrix(sconn.size());
                for(size_t i = 0; i < matrix.size(); ++i)
                    matrix[i].resize(fconn.size());
                
//                for(IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
//                {
//                    const auto sparent = igpt_s.get();
//                    const auto sw = igpt_s.getWeight();
//                    
//                    // cached tangent entries
//                    const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
//                    const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
//                    const double jdet_s = p_sel->jacDet(sparent, t1, t2);
//                    const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1, t2);
//                    const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
//                    const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
//                    const double jpiola_s = nurbs::cross(t1, t2).length();
//                    const auto x = p_sel->eval(sparent);
//                    
//                    //for(IPolarDegenerate igpt(GPt2D(1.0, 0.0)/*nurbs::projectPt(sparent, e1, e2)*/, degenerate_field_pair.second, forder); !igpt.isDone(); ++igpt)
//                    //for(IPolarIntegrate igpt(GPt2D(1.0, 0.0)/*nurbs::projectPt(sparent, e1, e2)*/, forder); !igpt.isDone(); ++igpt)
//                    for(IElemIntegrate igpt(forder); !igpt.isDone(); ++igpt)
//                    //for(IPolarIntegrate igpt(sparent, forder); !igpt.isDone(); ++igpt)
//                    //for(IPolarDegenerate igpt(sparent, degenerate_pair.second, forder); !igpt.isDone(); ++igpt)
//                    {
//                        
//                        const auto fparent = igpt.get();
//                        const Point3D y = p_fel->eval(fparent);
//                        const auto fw = igpt.getWeight();
//                        
//                        // cached tangent entries
//                        const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
//                        const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
//                        const double jdet_f = p_fel->jacDet(fparent, t1, t2);
//                        const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
//                        const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
//                        const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
//                        const double jpiola_f = nurbs::cross(t1, t2).length();
//                        
//                        // kernel entries
//                        const double r = dist(x,y);
//                        auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
//                        
//                        //if(igpt.currentSubCellI() == 0 ||igpt.currentSubCellI() == 1)
//                        //                        ekernel = std::complex<double>(1.0, 1.0);
//                        
//                        // now loop over test and trial functions
//                        for(size_t itest = 0; itest < fconn.size(); ++itest)
//                        {
//                            const auto igbasis_s = fconn[itest];
//                            if(-1 == igbasis_s)
//                                continue;
//                            
//                            const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
//                            
//                            for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
//                            {
//                                const auto igbasis_f = fconn[itrial];
//                                if(-1 == igbasis_f)
//                                    continue;
//                                
//                                const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
//                                //                                for(unsigned i = 0; i < 3; ++i)
//                                //                                    matrix[itest][itrial] += ( ekernel * basis_f[itrial][i]) * jdet_f * fw;
//                                //
//                                //                                // divergence (field)
//                                //
//                                //                                matrix[itest][itrial] -= 1.0 / (k * k) * (div_f * ekernel) * fw * jdet_f;
//                                
//                                for(unsigned i = 0; i < 3; ++i){
//                                    matrix[itest][itrial] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
//                                    
//                                }
//                                matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw * jdet_s * jdet_f;
//                                
//                            }
//                        }
//                        // increment # qpts counter
//                        ngpts += 1;
//                    }
//                }
                
                for(nurbs::IVertexQuadrature igpt(sorder, forder, v1, v2); !igpt.isDone(); ++igpt)
//                for(nurbs::IEdgeQuadrature igpt(sorder, forder, e1, e2); !igpt.isDone(); ++igpt)
                {
                    const auto gpt4d = igpt.get();
                    const auto sparent = gpt4d.srcPt();
                    const auto fparent = gpt4d.fieldPt();
                    
                    const auto w = igpt.getWeight();
                    
                    // source element parameters
                    const auto& t1_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
                    const auto& t2_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
                    const auto& x = p_sel->eval(sparent);
                    const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1_s, t2_s);
                    const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
                    const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
                    const auto& jdet_s = p_sel->jacDet(sparent, t1_s, t2_s);
                    const auto& jpiola_s = nurbs::cross(t1_s, t2_s).length();
                    
                    // field point terms
                    const auto& t1_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                    const auto& t2_f = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                    const auto y = p_fel->eval(fparent);
                    const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1_f, t2_f);
                    const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                    const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                    const double jdet_f = p_fel->jacDet(fparent, t1_f, t2_f);
                    const double jpiola_f = nurbs::cross(t1_f, t2_f).length();
                    
                    // kernel
                    const double r = dist(x, y);
                    const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
                    
                    // now loop over test and trial functions
                    for(size_t itest = 0; itest < sconn.size(); ++itest)
                    {
                        const auto igbasis_s = sconn[itest];
                        if(-1 == igbasis_s)
                            continue;
                        
                        
                        // divergence (source)
                        const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                        
                        for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                        {
                            const auto igbasis_f = fconn[itrial];
                            if(-1 == igbasis_f)
                                continue;
                            
                            // divergence (field)
                            const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                            
                            for(unsigned i = 0; i < 3; ++i)
                                matrix[itest][itrial] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * w;
                            matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * w * jdet_s * jdet_f;
                            
                        }
                    }
                    ngpts += 1;
                }
                
                std::cout << std::setprecision(15) << "integral = " << matrix[0][0].real() << "\n";
                ofs << ngpts << "\t" << std::setprecision(15) << matrix[0][0].real() << "\n";
            }
            
            // And now output the integrand using the highest order quadrature
            //        const GPt2D sourcept(0.96, -0.96);
            //        const auto x = p_fel->eval(sourcept);
            //        const UIntVec orders{max_order, max_order};
            
            //        IPolarDegenerate igpt(sourcept, degenerate_pair.second, {max_order, max_order});
            //        IPolarIntegrate igpt(sourcept, orders);
            //        std::map<std::pair<uint, uint>, std::vector<std::pair<double, GPt2D>>> data;
            //
            //        for(; !igpt.isDone(); ++igpt)
            //        {
            //            const auto ipair = std::make_pair(igpt.currentSubCellI(), igpt.currentSubSubCellI());
            //
            //            const GPt2D pt = igpt.get();
            //            const double w = igpt.getWeight();
            //            const double jdet_f = p_fel->jacDet(pt);
            //            const auto y = p_fel->eval(pt);
            //            const double r = dist(x,y);
            //            const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
            //            const double integrand = ekernel.real() * jdet_f * w / igpt.currentInnerWt();
            //
            //            const auto sample_pt = igpt.get(); //igpt.baseIntegrator().getBasePt();
            //            data[ipair].push_back(std::make_pair(integrand, sample_pt));
            //        }
            //
            //        for(const auto& p : data)
            //        {
            //            const auto& isubcell = p.first.first;
            //            const auto& isubsubcell = p.first.second;
            //            const auto& data_vec = p.second;
            //            std::vector<double> temp_dvec;
            //            std::vector<GPt2D> temp_pvec;
            //            for(const auto& p : data_vec)
            //            {
            //                temp_dvec.push_back(p.first);
            //                temp_pvec.push_back(p.second);
            //            }
            //            nurbs::OutputVTK output("subcell_integrand"  + std::to_string(isubcell) + "_" + std::to_string(isubsubcell));
            //            output.ouputQuadratureData("helmholtz integrand",
            //                                       temp_dvec,
            //                                       temp_pvec,
            //                                       orders[0],
            //                                       orders[1]);
            //        }
        }
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
        
    }
    catch(...)
    {
        std::cerr << "Unknown exception occured\n";
        return EXIT_FAILURE;
        
    }
}

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

#include "Forest.h"
#include "Geometry.h"
#include "BSplineSpace.h"
#include "BezierNodalElement.h"
#include "NodalElement.h"
#include "IElemIntegrate.h"
#include "OutputVTK.h"
#include "HConformingForest.h"
#include "ISubElem.h"
#include "IPolarIntegrate.h"

using namespace nurbs;

int main(int argc, char* argv[]) {
    
    try {
        
        std::cout << "Running degenerate geometry test.....\n";
        
        std::cout << "Trying to open hbs input file....\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            error("Cannot open file for reading\n");
        Geometry g;
        if(!g.loadHBSFile(ifs))
            error("Failed to load geometry from hbs data");
        
        Forest forest(g);
        for(uint ispace = 0; ispace < forest.spaceN(); ++ispace)
        {
            for(uint iedge = 0; iedge < NEDGES; ++iedge)
            {
                if(forest.degenerateEdge(ispace, iedge))
                    std::cout << "(Space,edge) = (" << ispace << "," << iedge << ") is degenerate\n";
            }
        }
        
        HDivForest divforest(g);
        

        
//        const uint hrefine = 0;
//        divforest.hrefine(hrefine);
        
        for(uint ielem = 0; ielem < divforest.elemN(); ++ielem)
            if(divforest.bezierElement(ielem)->degenerate())
                std::cout << "Element: " << ielem << " is degenerate\n";
        
        // edge singularity
        const auto p_sel = divforest.bezierElement(8);
        const auto p_fel = divforest.bezierElement(8);
        const double k = 100.0;
        const std::complex<double> iconst(0.0, 1.0);

        const auto& s_conn = p_sel->signedGlobalBasisFuncI();
        const auto& f_conn = p_fel->signedGlobalBasisFuncI();
        
        for(uint offset = 0; offset < 50; offset += 1)
        {

            
            std::vector<std::vector<std::complex<double>>>  matrix(s_conn.size());
            for(size_t i = 0; i < matrix.size(); ++i)
                matrix[i].resize(f_conn.size());
            
            const auto sorder = p_sel->equalIntegrationOrder(offset);
            const auto forder = p_fel->equalIntegrationOrder(offset + 1);
            
            nurbs::OutputVTK output("integrationOutput_" + std::to_string(offset));
            std::vector<double> data;
            
            const nurbs::GPt2D sparent(-0.9, -0.9);
            const auto x = p_fel->eval(sparent);
            const auto& t1_s = p_fel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
            const auto& t2_s = p_fel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
            const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
            const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
            const auto jparam_s = p_fel->jacobParam(sparent.s, sparent.t);
            const auto basis_source = p_fel->basis(sparent.s, sparent.t);
            const auto jpiola_s = nurbs::cross(t1_s, t2_s).length();
            
            
//            for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
//            {
//                const auto spt = igpt_s.get();
//                const auto swbase = igpt_s.getWeight();
//                
//                for(ISubElem isubel(1,1); !isubel.isDone(); ++isubel)
//                {
//                    const auto sparent = isubel.get(spt);
//                    
//                    const auto sw = swbase * isubel.jacob();
//                    
//                    // cached tangent entries
//                    const auto& t1_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
//                    const auto& t2_s = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
//                    
//                    // source element parameters
//                    const auto x = p_sel->eval(sparent);
//                    const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1_s, t2_s);
//                    const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
//                    const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
//                    const double jdet_s = p_sel->jacDet(sparent, t1_s, t2_s);
//                    const double jpiola_s = nurbs::cross(t1_s, t2_s).length();
//                    const auto jparam_s = p_sel->jacobParam(sparent.s, sparent.t);
            
                    // First integrate non singular term
                    for(nurbs::IElemIntegrate igpt_f(forder); !igpt_f.isDone(); ++igpt_f)
                    {
                        const auto ftemp = igpt_f.get();
                        const auto fparent = nurbs::GPt2D(0.5 * ( 1 + ftemp.s) * ( 1.0 - ftemp.t) - 1.0,
                                                          ftemp.t);
                        const double tri_jacob = (1 - ftemp.t) * 0.5;
                        
                        const auto fw = igpt_f.getWeight();
                        
                        // cached tangent entries
                        const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                        const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                        
                        // field element parameters
                        const auto y = p_fel->eval(fparent);
                        const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
                        const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                        const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                        const double jdet_f = p_fel->jacDet(fparent, t1, t2);
                        const double jpiola_f = nurbs::cross(t1, t2).length();
                        
                        // kernel
                        const double r = dist(x, y);
//                        if(r < 1.0e-9)
//                        {
////                            std::cout << "r = 0.0, continuing....\n";
//                            continue;
//                        }
                        const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
                        
                        // compute distance along tangent plane
                        const auto& rs_vec = t1_s * jparam_s[0][0] * (fparent.s - sparent.s) + t2_s * jparam_s[1][1] * (fparent.t - sparent.t);
                        const double rs = rs_vec.length();
                        
                        data.push_back(ekernel.real() * jdet_f);
                        
                        // now loop over test and trial functions
                        for(size_t itest = 0; itest < s_conn.size(); ++itest)
                        {
                            const auto igbasis_s = s_conn[itest];
                            if(-1 == igbasis_s)
                                continue;
                            
                            
                            // divergence (source)
                            const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                            
                            for(size_t itrial = 0; itrial < f_conn.size(); ++itrial)
                            {
                                const auto igbasis_f = f_conn[itrial];
                                if(-1 == igbasis_f)
                                    continue;
                                
                                // divergence (field)
                                const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                                
                                for(unsigned i = 0; i < 3; ++i)
                                    matrix[itest][itrial] += ( ekernel * basis_f[itrial][i] /*- basis_source[itrial][i] / (4.0 * nurbs::PI * rs)*/) *  jdet_f * tri_jacob * fw;
                                
//                                    matrix[itest][itrial] += (ekernel * basis_s[itest][i] * basis_f[itrial][i]
//                                                              * jdet_s * jdet_f * sw * fw;
                                
//                                matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw * jdet_s * jdet_f;
//                                matrix[itest][itrial] -= 1.0 / (k * k) * (div_f * ekernel /*- div_s / (4.0 * nurbs::PI * rs) */) * fw * jdet_f;
                            }
                        }
                    }
            
            output.ouputQuadratureData("helmholtzkernel", data, forder[0], forder[1]);
                    
                    
                    // integrate over field elements
//                    for(nurbs::IPolarIntegrate igpt_f(sparent, forder); !igpt_f.isDone(); ++igpt_f)
//                    {
//                        const auto fparent = igpt_f.get();
//                        const auto fw = igpt_f.getWeight();
//                        
//                        // cached tangent entries
//                        const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
//                        const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
//                        
//                        // field element parameters
//                        const auto y = p_fel->eval(fparent);
//                        const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
//                        const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
//                        const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
//                        const double jdet_f = p_fel->jacDet(fparent, t1, t2);
//                        const double jpiola_f = nurbs::cross(t1, t2).length();
//                        
//                        // kernel
//                        const double r = dist(x, y);
//                        const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
//                        
//                        // now loop over test and trial functions
//                        for(size_t itest = 0; itest < s_conn.size(); ++itest)
//                        {
//                            const auto igbasis_s = s_conn[itest];
//                            if(-1 == igbasis_s)
//                                continue;
//                            
//                            
//                            // divergence (source)
//                            const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
//                            
//                            for(size_t itrial = 0; itrial < f_conn.size(); ++itrial)
//                            {
//                                const auto igbasis_f = f_conn[itrial];
//                                if(-1 == igbasis_f)
//                                    continue;
//                                
//                                // divergence (field)
//                                const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
//                                
//                                for(unsigned i = 0; i < 3; ++i)
//                                    matrix[itest][itrial] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
//                                matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw * jdet_s * jdet_f;
//                                
//                            }
//                        }
//                    }
//                }
//            }
            std::cout << std::setprecision(12) << matrix[0][0].real() << "\n";
        }

        return EXIT_SUCCESS;
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

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
#include "IPolarIntegrate.h"
#include "IEqualQuadrature.h"
#include "IEdgeQuadrature.h"
#include "IEdgePolarIntegrate.h"
#include "IVertexQuadrature.h"
#include "ISubElem.h"

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
        
        const uint hrefine = 0;
        divforest.hrefine(hrefine);
        
        for(uint ielem = 0; ielem < divforest.elemN(); ++ielem)
            if(divforest.bezierElement(ielem)->degenerate())
                std::cout << "Element: " << ielem << " is degenerate\n";
        
        
        //        for(uint ielem = 0; ielem < divforest.elemN(); ++ielem)
        //        {
        //            std::cout << "Element: " << ielem << "\n";
        //            const auto el = divforest.bezierElement(ielem);
        //            if(el->degenerate())
        //            {
        //                std::cout << "Degenerate element\n";
        //
        //                const double k = 68.0;
        //                const std::complex<double> iconst(0.0, 1.0);
        //
        //                const auto& conn = el->signedGlobalBasisFuncI();
        //
        //                uint offset = 0;
        //                //        if(p_el->degenerate())
        //                //            offset += 4;
        //                const auto& forder = el->equalIntegrationOrder(offset);
        //                const auto& sorder = el->equalIntegrationOrder(offset);
        //                const nurbs::UIntVec nsubcells{5,5};
        //
        //                std::vector<std::vector<std::complex<double>>>  matrix(conn.size());
        //                for(size_t i = 0; i < matrix.size(); ++i)
        //                    matrix[i].resize(conn.size());
        //
        //                // and finally loop over all regular integrals
        ////                for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
        ////                {
        ////                    const auto sparent = igpt_s.get();
        ////                    const auto sw = igpt_s.getWeight();
        ////
        ////                    // cached tangent entries
        ////                    const auto& t1 = el->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
        ////                    const auto& t2 = el->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
        ////
        ////                    // source element parameters
        ////                    const auto x = el->eval(sparent);
        ////                    const auto& basis_s = el->basis(sparent.s, sparent.t, t1, t2);
        ////                    const auto& ds_s = el->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
        ////                    const auto& dt_s = el->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
        ////                    const double jdet_s = el->jacDet(sparent, t1, t2);
        ////                    const double jpiola_s = nurbs::cross(t1, t2).length();
        ////
        ////                    // integrate over field elements
        ////                    for(nurbs::IPolarIntegrate igpt_f(sparent, forder, nsubcells); !igpt_f.isDone(); ++igpt_f)
        ////                    {
        ////                        const auto fparent = igpt_f.get();
        ////                        const auto fw = igpt_f.getWeight();
        ////
        ////                        // cached tangent entries
        ////                        const auto& t1 = el->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
        ////                        const auto& t2 = el->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
        ////
        ////                        // field element parameters
        ////                        const auto y = el->eval(fparent);
        ////                        const auto& basis_f = el->basis(fparent.s, fparent.t, t1, t2);
        ////                        const auto& ds_f = el->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
        ////                        const auto& dt_f = el->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
        ////                        const double jdet_f = el->jacDet(fparent, t1, t2);
        ////                        const double jpiola_f = nurbs::cross(t1, t2).length();
        ////
        ////                        // kernel
        ////                        const double r = dist(x, y);
        ////                        const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
        ////
        ////                        // now loop over test and trial functions
        ////                        for(size_t itest = 0; itest < conn.size(); ++itest)
        ////                        {
        ////                            const auto igbasis_s = conn[itest];
        ////                            if(-1 == igbasis_s)
        ////                                continue;
        ////
        ////
        ////                            // divergence (source)
        ////                            const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
        ////
        ////                            for(size_t itrial = 0; itrial < conn.size(); ++itrial)
        ////                            {
        ////                                const auto igbasis_f = conn[itrial];
        ////                                if(-1 == igbasis_f)
        ////                                    continue;
        ////
        ////                                // divergence (field)
        ////                                const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
        ////
        ////                                for(unsigned i = 0; i < 3; ++i)
        ////                                    matrix[itest][itrial] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
        ////                                matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw * jdet_s * jdet_f;
        ////
        ////                            }
        ////                        }
        ////                    }
        ////                }
        ////
        ////                // output
        ////                for(auto& row : matrix)
        ////                {
        ////                    for(auto& val : row)
        ////                    {
        ////                        std::cout << val << "\t";
        ////                        val = 0.0;
        ////                    }
        ////                    std::cout << "\n";
        ////                }
        //
        //
        //                for(nurbs::IEqualQuadrature igpt(sorder, forder); !igpt.isDone(); ++igpt)
        //                {
        //                    const auto gpt4d = igpt.get();
        //                    const auto sparent = gpt4d.srcPt();
        //                    const auto fparent = gpt4d.fieldPt();
        //
        //                    const auto w = igpt.getWeight();
        //
        //                    // source element parameters
        //                    const auto& t1_s = el->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
        //                    const auto& t2_s = el->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
        //                    const auto& x = el->eval(sparent);
        //                    const auto& basis_s = el->basis(sparent.s, sparent.t, t1_s, t2_s);
        //                    const auto& ds_s = el->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
        //                    const auto& dt_s = el->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
        //                    const auto& jdet_s = el->jacDet(sparent, t1_s, t2_s);
        //                    const auto& jpiola_s = nurbs::cross(t1_s, t2_s).length();
        //
        //                    // field point terms
        //                    const auto& t1_f = el->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
        //                    const auto& t2_f = el->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
        //                    const auto y = el->eval(fparent);
        //                    const auto& basis_f = el->basis(fparent.s, fparent.t, t1_f, t2_f);
        //                    const auto& ds_f = el->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
        //                    const auto& dt_f = el->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
        //                    const double jdet_f = el->jacDet(fparent, t1_f, t2_f);
        //                    const double jpiola_f = nurbs::cross(t1_f, t2_f).length();
        //
        //                    // kernel
        //                    const double r = dist(x, y);
        //                    const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
        //
        //                    // now loop over test and trial functions
        //                    for(size_t itest = 0; itest < conn.size(); ++itest)
        //                    {
        //                        const auto igbasis_s = conn[itest];
        //                        if(-1 == igbasis_s)
        //                            continue;
        //
        //
        //                        // divergence (source)
        //                        const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
        //
        //                        for(size_t itrial = 0; itrial < conn.size(); ++itrial)
        //                        {
        //                            const auto igbasis_f = conn[itrial];
        //                            if(-1 == igbasis_f)
        //                                continue;
        //
        //                            // divergence (field)
        //                            const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
        //
        //                            for(unsigned i = 0; i < 3; ++i)
        //                                matrix[itest][itrial] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * w;
        //                            matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * w * jdet_s * jdet_f;
        //
        //                        }
        //                    }
        //                }
        //
        //                // output
        ////                for(auto& row : matrix)
        ////                {
        ////                    for(auto& val : row)
        ////                    {
        ////                        std::cout << val << "\t";
        ////                        val = 0.0;
        ////                    }
        ////                    std::cout << "\n";
        ////                }
        //
        //
        //                break;
        //            }
        //        }
        //
        
        // edge singularity
        const auto p_sel = divforest.bezierElement(8);
        const auto p_fel = divforest.bezierElement(8);
        
        const double k = 100.0;
        const std::complex<double> iconst(0.0, 1.0);
        
        const auto& s_conn = p_sel->signedGlobalBasisFuncI();
        const auto& f_conn = p_fel->signedGlobalBasisFuncI();
        
//        std::cout << p_sel->signedGlobalBasisFuncI() << "\n" << p_fel->signedGlobalBasisFuncI() << "\n";
        
        std::cout << "Default quadrature order: " << p_sel->equalIntegrationOrder() << "\n";
        
        for(uint offset = 0; offset < 15; offset += 1)
        {
//            const auto& forder = p_fel->equalIntegrationOrder(offset);
//            const auto& sorder = p_sel->equalIntegrationOrder(offset);
            
            const nurbs::UIntVec forder{2 + offset, 2 + offset};
            const nurbs::UIntVec sorder{2 + offset, 2 + offset};

            std::vector<std::vector<std::complex<double>>>  matrix(s_conn.size());
            for(size_t i = 0; i < matrix.size(); ++i)
                matrix[i].resize(f_conn.size());
            
            uint neval = 0;
            
            for(nurbs::IEqualQuadrature igpt(sorder, forder); !igpt.isDone(); ++igpt)
            {
                const auto gpt4d = igpt.get();
                const auto spt = gpt4d.srcPt();
                const auto fpt = gpt4d.fieldPt();
                
                const auto wbase = igpt.getWeight();
                
                for(ISubElem isubelsrc(1,1); !isubelsrc.isDone(); ++isubelsrc)
                {
                    const auto sparent = isubelsrc.get(spt);
                    
                    for(ISubElem isubelfield(1,1); !isubelfield.isDone(); ++isubelfield)
                    {
                        ++neval;
                        const auto fparent = isubelfield.get(fpt);
                        const auto w = wbase * isubelsrc.jacob() * isubelfield.jacob();
                        
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
                                    matrix[itest][itrial] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * w;
                                matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * w * jdet_s * jdet_f;
                                
                            }
                        }
                    }
                }
            }
            
//
//            for(nurbs::IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
//            {
//                
//                const auto spt = igpt_s.get();
//                const auto swbase = igpt_s.getWeight();
//                
//                
//                for(ISubElem isubel(1,1); !isubel.isDone(); ++isubel)
//                {
//                    const auto sparent = isubel.get(spt);
//                    
//                    const auto sw = swbase * isubel.jacob();
//                    
//                    // cached tangent entries
//                    const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
//                    const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
//                    
//                    // source element parameters
//                    const auto x = p_sel->eval(sparent);
//                    const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1, t2);
//                    const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
//                    const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
//                    const double jdet_s = p_sel->jacDet(sparent, t1, t2);
//                    const double jpiola_s = nurbs::cross(t1, t2).length();
//                    
//                    // integrate over field elements
//                    for(nurbs::IPolarIntegrate igpt_f(sparent, forder); !igpt_f.isDone(); ++igpt_f)
//                    {
//                        ++neval;
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
//            
//            
            std::cout << neval << "\t" << std::setprecision(12) << matrix[0][0].real() << "\n";
            
        }
        
        OutputVTK output("degenerate_geometry");
        //        output.outputGeometry(forest);
        
        std::vector<std::complex<double>> soln(divforest.globalDofN());
        output.outputComplexVectorField(divforest, "vectorfield", soln);
        
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

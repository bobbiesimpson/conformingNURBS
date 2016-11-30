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

using namespace nurbs;

int main(int argc, char* argv[])
{
    // Some constants
    const uint max_order = 12;                      // max order of quadrature that we loop unilt
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
        
        std::ofstream ofs("degenerate_quadrature_results.dat");
        if(!ofs)
            throw std::runtime_error("Cannot open file for writing.");
        
        Forest forest(g);
        HDivForest divforest(g);
        
        const auto p_fel = divforest.bezierElement(100);
        const auto p_sel = divforest.bezierElement(100);
        
        const auto degenerate_pair_source = p_sel->degenerateEdge();
        
        if(degenerate_pair_source.first)
            std::cout << "Source element is degenerate.\n\n";
        else
            std::cout << "non degenerate source element\n\n";
        
        const auto degenerate_pair = p_fel->degenerateEdge();
        
        if(degenerate_pair.first)
            std::cout << "Field element is degenerate.\n\n";
        else
            std::cout << "non degenerate field element\n\n";
        
        // element connectivity
        const auto fconn = p_fel->signedGlobalBasisFuncI();
        const auto sconn = p_sel->signedGlobalBasisFuncI();
        
        const auto forder = UIntVec{10,5};//p_fel->equalIntegrationOrder(2);
        
        for(uint order = 2; order < max_order; ++order)
        {
            const UIntVec sorder{order, order};
            
            std::cout << "Applying special polar quadrature with outer integration order: " << sorder << "....\n";
            
            uint ngpts = 0;                                 // number of quadrature points
            
            std::vector<std::vector<std::complex<double>>>  matrix(sconn.size());
            for(size_t i = 0; i < matrix.size(); ++i)
                matrix[i].resize(fconn.size());
            
            for(IElemIntegrate igpt_s(sorder); !igpt_s.isDone(); ++igpt_s)
            {
                const auto sparent = igpt_s.get();
                const auto sw = igpt_s.getWeight();
                
                // cached tangent entries
                const auto& t1 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::S);
                const auto& t2 = p_sel->tangent(sparent.s, sparent.t, nurbs::ParamDir::T);
                const double jdet_s = p_sel->jacDet(sparent, t1, t2);
                const auto& basis_s = p_sel->basis(sparent.s, sparent.t, t1, t2);
                const auto& ds_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DS);
                const auto& dt_s = p_sel->localBasisDers(sparent.s, sparent.t, nurbs::DerivType::DT);
                const double jpiola_s = nurbs::cross(t1, t2).length();
                const auto x = p_sel->eval(sparent);
                
                //                for(IElemIntegrate igpt(forder); !igpt.isDone(); ++igpt)
                for(IPolarIntegrate igpt(sparent, forder); !igpt.isDone(); ++igpt)
//                for(IPolarDegenerate igpt(sparent, degenerate_pair.second, forder); !igpt.isDone(); ++igpt)
                {
                    
                    const auto fparent = igpt.get();
                    const Point3D y = p_fel->eval(fparent);
                    const auto fw = igpt.getWeight();
                    
                    // cached tangent entries
                    const auto& t1 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::S);
                    const auto& t2 = p_fel->tangent(fparent.s, fparent.t, nurbs::ParamDir::T);
                    const double jdet_f = p_fel->jacDet(fparent, t1, t2);
                    const auto& basis_f = p_fel->basis(fparent.s, fparent.t, t1, t2);
                    const auto& ds_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DS);
                    const auto& dt_f = p_fel->localBasisDers(fparent.s, fparent.t, nurbs::DerivType::DT);
                    const double jpiola_f = nurbs::cross(t1, t2).length();
                    
                    // kernel entries
                    const double r = dist(x,y);
                    auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
                    
                    //if(igpt.currentSubCellI() == 0 ||igpt.currentSubCellI() == 1)
                    //                        ekernel = std::complex<double>(1.0, 1.0);
                    
                    // now loop over test and trial functions
                    for(size_t itest = 0; itest < fconn.size(); ++itest)
                    {
                        const auto igbasis_s = fconn[itest];
                        if(-1 == igbasis_s)
                            continue;
                        
                        const double div_s = 1./jpiola_s * (ds_s[itest][0] + dt_s[itest][1]);
                        
                        for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                        {
                            const auto igbasis_f = fconn[itrial];
                            if(-1 == igbasis_f)
                                continue;
                            
                            const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                            //                                for(unsigned i = 0; i < 3; ++i)
                            //                                    matrix[itest][itrial] += ( ekernel * basis_f[itrial][i]) * jdet_f * fw;
                            //
                            //                                // divergence (field)
                            //
                            //                                matrix[itest][itrial] -= 1.0 / (k * k) * (div_f * ekernel) * fw * jdet_f;
                            
                            for(unsigned i = 0; i < 3; ++i){
                                matrix[itest][itrial] += ekernel * basis_s[itest][i] * basis_f[itrial][i] * jdet_s * jdet_f * sw * fw;
                                
                            }
                            matrix[itest][itrial] -= 1.0 / (k * k) * div_s * div_f * ekernel * sw * fw * jdet_s * jdet_f;
                            
                        }
                    }
                    // increment # qpts counter
                    ngpts += 1;
                }
            }
            //                std::cout << std::setprecision(15) << "area = " << a << "\n";
            std::cout << std::setprecision(15) << "integral = " << matrix[0][0].real() << "\n";
            ofs << ngpts << "\t" << matrix[0][0].real() << "\n";
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

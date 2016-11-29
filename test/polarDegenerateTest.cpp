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
    const uint max_order = 30;                      // max order of quadrature that we loop unilt
    const double k = 100.0;                         // wavenumber for emag kernel
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
        
        const auto p_fel = divforest.bezierElement(8);
        const auto degenerate_pair = p_fel->degenerateEdge();
        
        const GPt2D sourcept(-1.0, -1.0);     // hardcoded source point in parent domain
        const Point3D x = p_fel->eval(sourcept);
        
        // element connectivity
        const auto fconn = p_fel->signedGlobalBasisFuncI();
        
        for(uint order = 1; order < max_order; ++order)
        {
            if(degenerate_pair.first)
            {
                // initiaite output variables
                nurbs::OutputVTK output("degenerate_quadrature_" + std::to_string(order));
                std::vector<double> data;
                std::vector<nurbs::GPt2D> pts;
                
                const UIntVec forder{order, order};
                std::cout << "Degenerate element. Applying special polar quadrature with order: " << forder << "....\n";
                
                double a = 0.0;                                 // element surface area (for checking)
                uint ngpts = 0;                                 // number of quadrature points

                std::vector<std::vector<std::complex<double>>>  matrix(fconn.size());
                for(size_t i = 0; i < matrix.size(); ++i)
                    matrix[i].resize(fconn.size());
                
//                for(IElemIntegrate igpt(forder); !igpt.isDone(); ++igpt)
                //for(IPolarIntegrate igpt(sourcept, forder); !igpt.isDone(); ++igpt)
                for(IPolarDegenerate igpt(sourcept, degenerate_pair.second, forder); !igpt.isDone(); ++igpt)
                {
                    const auto fparent = igpt.get();
                    const Point3D y = p_fel->eval(fparent);
                    const auto fw = igpt.getWeight();
                    
//                    std::cout << fparent << "\t" << fw << "\n";
                    
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
                    const auto ekernel = std::exp(-iconst * k * r) / (4.0 * nurbs::PI * r);
                    
                    // now loop over test and trial functions
                    for(size_t itest = 0; itest < fconn.size(); ++itest)
                    {
                        const auto igbasis_s = fconn[itest];
                        if(-1 == igbasis_s)
                            continue;
                        
                        for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
                        {
                            const auto igbasis_f = fconn[itrial];
                            if(-1 == igbasis_f)
                                continue;
                            
                            for(unsigned i = 0; i < 3; ++i)
                                matrix[itest][itrial] += ( ekernel * basis_f[itrial][i]) * jdet_f * fw;
                            
                            // divergence (field)
                            const double div_f = 1./jpiola_f * (ds_f[itrial][0] + dt_f[itrial][1]);
                            matrix[itest][itrial] -= 1.0 / (k * k) * (div_f * ekernel) * fw * jdet_f;
                        }
                    }
                    
                    if(igpt.currentSubCellI() == 3 && igpt.currentSubSubCellI() == 0)
                    {
                        data.push_back(ekernel.real() * jdet_f * fw);
//                        pts.push_back(igpt.baseIntegrator().getBasePt());
                        pts.push_back(fparent);
                    }
                    
                    // remember area calculation
                    a += jdet_f * fw;
                    
                    // increment # qpts counter
                    ngpts += 1;
                }
                std::cout << std::setprecision(15) << "area = " << a << "\n";
                std::cout << std::setprecision(15) << "integral = " << matrix[0][0].real() << "\n";
                ofs << ngpts << "\t" << matrix[0][0].real() << "\n";
                
                // If we're on the highest order create an output of the integrand
                if(max_order -1 == order)
                    output.ouputQuadratureData("polardegenerate_integrand", data, pts, forder[0], forder[1]);

            }
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

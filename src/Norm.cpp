#include "Norm.h"
#include "BezierVectorElement.h"
#include "base.h"

namespace nurbs {
    
    double hdivNorm(const MultiForest& f,
                    const std::vector<std::complex<double>>& phi)
    {
        std::cout << "Computing H_{-1/2}(div) norm....\n";
        double hdivnorm = 0.0;
            
        // loop over source elements
        for(size_t isrcel = 0; isrcel < f.elemN(); ++isrcel)
        {
            std::cout << "Source element " << isrcel + 1 << "/" << f.elemN() << "\n";
            const auto s_el = f.bezierElement(isrcel);
            const auto& sorder = s_el->equalIntegrationOrder();
            const auto& sconn = s_el->globalBasisFuncI();
            
            for(nurbs::IElemIntegrate is_gpt(sorder); !is_gpt.isDone(); ++is_gpt)
            {
                const auto& s_gpt = is_gpt.get();
                const auto& s_w = is_gpt.getWeight();
                
                // compute source element basis function terms
                const auto& t1 = s_el->tangent(s_gpt.s, s_gpt.t, nurbs::ParamDir::S);
                const auto& t2 = s_el->tangent(s_gpt.s, s_gpt.t, nurbs::ParamDir::T);
                const auto& s_basis = s_el->basis(s_gpt.s, s_gpt.t, t1, t2);
                const auto& x = s_el->eval(s_gpt);
                const auto& ds_s = s_el->localBasisDers(s_gpt.s, s_gpt.t, nurbs::DerivType::DS);
                const auto& dt_s = s_el->localBasisDers(s_gpt.s, s_gpt.t, nurbs::DerivType::DT);
                const double jdet_s = s_el->jacDet(s_gpt, t1, t2);
                
                // interpolate phi using given solution
                std::vector<std::complex<double>> s_phi;
                double s_div = 0.0;
                
                for(size_t is_basis = 0; is_basis < sconn.size(); ++is_basis)
                {
                    const auto phi_h = phi[sconn[is_basis]];

//                    for(unsigned i = 0; i < 3; ++i)
//                    {
//                        s_phi[i] +=  * s_basis[is_basis][i];
//                    }
                }
                
                // loop over field elements
                for(size_t ifieldel = 0; ifieldel < f.elemN(); ++ifieldel)
                {
                    // if singular, we'll compute this later
                    if(isrcel == ifieldel)
                        continue;
                    
                    const auto f_el = f.bezierElement(ifieldel);
                    const auto& forder = f_el->equalIntegrationOrder();
                    const auto& fconn = f_el->globalBasisFuncI();
                    
                    for(nurbs::IElemIntegrate if_gpt(forder); !if_gpt.isDone(); ++if_gpt)
                    {
                        const auto& f_gpt = if_gpt.get();
                        const auto& w = if_gpt.getWeight();
                        
                        const auto& t1 = f_el->tangent(f_gpt.s, f_gpt.t, nurbs::ParamDir::S);
                        const auto& t2 = f_el->tangent(f_gpt.s, f_gpt.t, nurbs::ParamDir::T);
                        const auto& f_basis = f_el->basis(f_gpt.s, f_gpt.t, t1, t2);
                        const auto& y = f_el->eval(f_gpt);
                        const auto& ds_f = f_el->localBasisDers(f_gpt.s, f_gpt.t, nurbs::DerivType::DS);
                        const auto& dt_f = f_el->localBasisDers(f_gpt.s, f_gpt.t, nurbs::DerivType::DT);
                        const double jdet_f = f_el->jacDet(f_gpt, t1, t2);
                        
                        const double r = dist(x,y);
                        const auto kernel = 1.0 / (4.0 * nurbs::PI * r);
                        
                        for(size_t itest = 0; itest < sconn.size(); ++itest)
                        {
//                            const double s_div = (ds_s[itest][0] + dt_s[itest][1]);
//                            const auto s_phi =
//                            for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
//                            {
//                                const double f_div = (ds_f[itrial][0] + dt_f[itrial][1]);
//                                
//                            }
                        }
                    }
                }
            }
        }
        
        // now compute all singular integrals
//        for(size_t iel = 0; iel < f.elemN(); ++iel)
//        {
//            
//        }
        
        return 0.0;
    }
}
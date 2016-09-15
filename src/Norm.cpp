#include "Norm.h"
#include "BezierVectorElement.h"
#include "IElemIntegrate.h"
#include "base.h"

namespace nurbs {
    
    double hdivNorm(const MultiForest& f,
                    const std::vector<std::complex<double>>& phi)
    {
//        std::cout << "Computing H_{-1/2}(div) norm....\n";
//        double hdivnorm = 0.0;
//            
//        // loop over source elements
//        for(size_t isrcel = 0; isrcel < f.elemN(); ++isrcel)
//        {
//            std::cout << "Source element " << isrcel + 1 << "/" << f.elemN() << "\n";
//            const auto s_el = f.bezierElement(isrcel);
//            const auto& sorder = s_el->equalIntegrationOrder();
//            const auto& sconn = s_el->globalBasisFuncI();
//            
//            for(nurbs::IElemIntegrate is_gpt(sorder); !is_gpt.isDone(); ++is_gpt)
//            {
//                const auto& s_gpt = is_gpt.get();
//                const auto& s_w = is_gpt.getWeight();
//                
//                // compute source element basis function terms
//                const auto& t1 = s_el->tangent(s_gpt.s, s_gpt.t, nurbs::ParamDir::S);
//                const auto& t2 = s_el->tangent(s_gpt.s, s_gpt.t, nurbs::ParamDir::T);
//                const auto& s_basis = s_el->basis(s_gpt.s, s_gpt.t, t1, t2);
//                const auto& x = s_el->eval(s_gpt);
//                const auto& ds_s = s_el->localBasisDers(s_gpt.s, s_gpt.t, nurbs::DerivType::DS);
//                const auto& dt_s = s_el->localBasisDers(s_gpt.s, s_gpt.t, nurbs::DerivType::DT);
//                const double jdet_s = s_el->jacDet(s_gpt, t1, t2);
//                
//                // interpolate phi using given solution
//                std::vector<std::complex<double>> s_phi;
//                double s_div = 0.0;
//                
//                for(size_t is_basis = 0; is_basis < sconn.size(); ++is_basis)
//                {
//                    const auto phi_h = phi[sconn[is_basis]];
//
////                    for(unsigned i = 0; i < 3; ++i)
////                    {
////                        s_phi[i] +=  * s_basis[is_basis][i];
////                    }
//                }
//                
//                // loop over field elements
//                for(size_t ifieldel = 0; ifieldel < f.elemN(); ++ifieldel)
//                {
//                    // if singular, we'll compute this later
//                    if(isrcel == ifieldel)
//                        continue;
//                    
//                    const auto f_el = f.bezierElement(ifieldel);
//                    const auto& forder = f_el->equalIntegrationOrder();
//                    const auto& fconn = f_el->globalBasisFuncI();
//                    
//                    for(nurbs::IElemIntegrate if_gpt(forder); !if_gpt.isDone(); ++if_gpt)
//                    {
//                        const auto& f_gpt = if_gpt.get();
//                        const auto& w = if_gpt.getWeight();
//                        
//                        const auto& t1 = f_el->tangent(f_gpt.s, f_gpt.t, nurbs::ParamDir::S);
//                        const auto& t2 = f_el->tangent(f_gpt.s, f_gpt.t, nurbs::ParamDir::T);
//                        const auto& f_basis = f_el->basis(f_gpt.s, f_gpt.t, t1, t2);
//                        const auto& y = f_el->eval(f_gpt);
//                        const auto& ds_f = f_el->localBasisDers(f_gpt.s, f_gpt.t, nurbs::DerivType::DS);
//                        const auto& dt_f = f_el->localBasisDers(f_gpt.s, f_gpt.t, nurbs::DerivType::DT);
//                        const double jdet_f = f_el->jacDet(f_gpt, t1, t2);
//                        
//                        const double r = dist(x,y);
//                        const auto kernel = 1.0 / (4.0 * nurbs::PI * r);
//                        
//                        for(size_t itest = 0; itest < sconn.size(); ++itest)
//                        {
////                            const double s_div = (ds_s[itest][0] + dt_s[itest][1]);
////                            const auto s_phi =
////                            for(size_t itrial = 0; itrial < fconn.size(); ++itrial)
////                            {
////                                const double f_div = (ds_f[itrial][0] + dt_f[itrial][1]);
////                                
////                            }
//                        }
//                    }
//                }
//            }
//        }
        
        // now compute all singular integrals
//        for(size_t iel = 0; iel < f.elemN(); ++iel)
//        {
//            
//        }
        
        return 0.0;
    }
    
    double L2graphNorm(const MultiForest& f,
                       const std::vector<std::complex<double>>& soln)
    {
        // compute the L2 grapm norm = L2 norm + L2 norm of surface divergence
        
//        assert(soln.size() == f.globalDofN());
        assert(soln.size() == f.globalNedelecDofN());
        
        // hardcoded analytical solution
        struct Function
        {
            /// override function operator
            std::vector<std::complex<double>> operator()(const nurbs::Point3D& x) const
            {
//                const double x = p[0];
//                const double y = p[1];
//                const double z = p[2];
                
                const Point3D p(0.0, 1.0, 0.0);
                const Point3D k(5.0, 0.0, 0.0);
                
                std::vector<std::complex<double>> result{p[0], p[1], p[2]};
                const auto wave = std::exp(-std::complex<double>(0.0, dot(k, x)));
                for(auto& r : result)
                    r *= wave;
                return result;
//                return
//                {
//                    x*z / std::sqrt(x*x + y*y),
//                    y*z / std::sqrt(x*x + y*y),
//                    -std::sqrt(x * x + y * y)
//                };
            }
            
//            std::vector<double> curl(const nurbs::Point3D& p) const
//            {
//                const double x = p[0];
//                const double y = p[1];
//                const double z = p[2];
//                const double const1 = x * x + y * y;
//                const double dfzdy = -y /std::sqrt(const1);
//                const double dfydz = -y /std::sqrt(const1);
//                const double dfxdz = -x /std::sqrt(const1);
//                const double dfzdx = -x /std::sqrt(const1);
//                const double dfydx = -x * y * z / (std::pow(const1, 3.0/2.0));
//                const double dfxdy = -x * y * z / (std::pow(const1, 3.0/2.0));
//                return {dfzdy - dfydz, dfxdz - dfzdx, dfydx - dfxdy};
//            }
//            
//            std::vector<double> div(const nurbs::Point3D& p) const
//            {
//                const double x = p[0];
//                const double y = p[1];
//                const double z = p[2];
//                const double const1 = x * x + y * y;
//                
//                return
//                {
//                    y * y * z / (std::pow(const1, 3.0/2.0)),
//                    x * x * z / (std::pow(const1, 3.0/2.0)),
//                    0.0
//                };
//            }
            
        } function;
        
        double norm = 0.0;
        double exact_norm = 0.0;
        
        for(uint ielem = 0; ielem < f.elemN(); ++ielem)
        {
//            const auto el = f.bezierElement(ielem);
            const auto el = f.nedelecElement(ielem);
            
            for(IElemIntegrate igpt(el->integrationOrder()); !igpt.isDone(); ++igpt)
            {
                const auto gpt = igpt.get();
                const auto w = igpt.getWeight();
                
                const auto basis = el->basis(gpt.s, gpt.t);
                const auto localbasis = el->localBasis(gpt.s, gpt.t);
                const auto jdet = el->jacDet(gpt);
                const auto x = el->eval(gpt);
                
                // Genereate 'exact' function in parametric space
                auto exact_val = function(x);
                
                const auto& t1 = el->tangent(gpt.s, gpt.t, nurbs::ParamDir::S);
                const auto& t2 = el->tangent(gpt.s, gpt.t,nurbs::ParamDir::T);
                const double jpiola = nurbs::cross(t1, t2).length();
                
                DoubleVecVec j;
                j.push_back(t1.asVec());
                j.push_back(t2.asVec());
                j.push_back(
                                {
                                    1.0/jpiola * (j[0][1] * j[1][2] - j[0][2] * j[1][1]),
                                    1.0/jpiola * (j[0][2] * j[1][0] - j[0][0] * j[1][2]),
                                    1.0/jpiola * (j[0][0] * j[1][1] - j[0][1] * j[1][0])
                                });
                auto jinv = inv3x3Mat(j);
                for(auto& row : jinv)
                    row.erase(row.begin() + 2);
                
                
                std::vector<std::complex<double>> fexact(2);
                std::vector<std::complex<double>> f_h(2);
                
                for(size_t i = 0; i < 2; ++i)
                    for(size_t j = 0; j < 3; ++j)
                        fexact[i] += jpiola * jinv[j][i] * exact_val[j];
                
//
//                const auto& ds = el->localBasisDers(gpt.s, gpt.t, nurbs::DerivType::DS);
//                const auto& dt = el->localBasisDers(gpt.s, gpt.t, nurbs::DerivType::DT);
//                const auto jacob_param = el->jacobParam(gpt.s, gpt.t);
                const auto econn = el->globalBasisFuncI();
                
                // interpolate numerical value
                std::vector<std::complex<double>> val(3,0.0); // interpolated surface current
                
                for(uint ibasis = 0; ibasis < econn.size(); ++ibasis)
                {
                    if(econn[ibasis] == -1)
                        continue;
                    
                    for(uint i = 0; i < 3; ++i)
                        val[i] += basis[ibasis][i] * soln[econn[ibasis]];
                    
                    for(size_t i = 0; i < 2; ++i)
                        f_h[i] += localbasis[ibasis][i] * soln[econn[ibasis]];
                    
//                    surface_div += 1./jpiola * (ds[ibasis][0] + dt[ibasis][1]) * soln[econn[ibasis]];
                }
                
//                std::cout << f_h << "\n";
//                std::cout << fexact << "\n";
                
                // generate 'exact' solution
                Point3D temp_pt(exact_val[0].real(),
                                exact_val[1].real(),
                                exact_val[2].real());

                // now restrict to the surface
                const auto n = el->normal(gpt);
                const auto exact_real = cross(n, cross(temp_pt, n)).asVec();
                
                temp_pt.setCoord(0, exact_val[0].imag());
                temp_pt.setCoord(1, exact_val[1].imag());
                temp_pt.setCoord(2, exact_val[2].imag());
                const auto exact_imag = cross(n, cross(temp_pt, n)).asVec();
                
                for(size_t i = 0; i < 3; ++i)
                {
                    const std::complex<double> phi(exact_real[i], exact_imag[i]);
                    const auto& phi_h = val[i];
                    norm += std::norm(phi_h - phi) * jdet * w;
                    exact_norm += std::norm(phi) * jdet * w;
                }
                
//                for(size_t i = 0; i < 2; ++i)
//                {
//                    
//                    norm += std::norm(f_h[i] - fexact[i]) * jdet * w;
//                    exact_norm += std::norm(fexact[i]) * jdet * w;
//                }

                
                // now compute analytical surface divergence
//                Point3D div_pt(function.div(x));
//                const double exact_div_val = dot(div_pt, t1) + dot(div_pt, t2);
//                
//                
//                div_pt -= dot(n, div_pt) * n;
//                const double div = div_pt[0] + div_pt[1] + div_pt[2];
//                std::cout << "Divergence (numerical/analytical): " << surface_div << "\t" << div << "\n";
//                std::cout << div/ surface_div << "\n";
            }
        }
//        std::cout << "squared error norm = " << norm << "\n";
        return std::sqrt(norm) / std::sqrt(exact_norm);
        
    }
}
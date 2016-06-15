//
//  Header.h
//  nurbslib
//
//  Created by Robert Simpson on 14/09/2015.
//
//

#ifndef nurbslib_Header_h
#define nurbslib_Header_h

#include <set>
#include <vector>
#include <thread>
#include <tuple>

#include "base.h"
#include "Forest.h"

#include "Adaptive2DPolarClenshaw.h"
#include "Adaptive2DClenshaw.h"

#include "Epetra_CrsMatrix.h"

using namespace nurbs::elem;

namespace nurbs {
    
    /// A 'work' function that is fed to a single process thread to allow parallel BEM
    /// computations on multiple cores.
    template<typename F>
    void workfunction(const Forest& f,
                      const std::vector<double>& bcdata,
                      const std::vector<std::pair<uint, uint>>& collocvec,
                      const double tol,
                      const uint lower,
                      const uint upper,
                      std::vector<double>& result)
    {
        for(uint icpt = lower; icpt < upper; ++icpt) {
            
            // Get the relevant collocation indicies
            const uint icspace = collocvec[icpt].first; // space index
            const uint icolloc = collocvec[icpt].second; // local cpt index
            const uint gindex = f.globalCollocI(icspace, icolloc);
            
            // Now construct the surface functor
            F functor(f, icspace, icolloc, bcdata);
            
            //std::cout << "Collocatiuon index : " << icpt << " (" << lower << "/" << upper << ")\n";
            // Loop over all 'field' spaces
            for(uint ifspace = 0; ifspace < f.spaceN(); ++ifspace) {
                functor.setFieldSpaceI(ifspace);
                const BSplineSpace& fspace = f.space(ifspace);
                const auto s_interval = fspace.range(S);
                const auto t_interval = fspace.range(T);
                const auto connected = f.connectedCollocPtI(icspace, icolloc, ifspace);
                
                // Singular quadrature
                if(connected.first) {
                    //std::cout << "integrating...\n";
                    const auto cpt_param = fspace.grevilleAbscissaPt(connected.second);
                    bemquad::Adaptive2DPolarClenshawIntegrator<double, F> integrator(s_interval.first,
                                                                                     s_interval.second,
                                                                                     t_interval.first,
                                                                                     t_interval.second,
                                                                                     cpt_param.s,
                                                                                     cpt_param.t,
                                                                                     functor, tol);
                    //std::cout << "inserting...\n";
                    result[gindex] += -1.0 * integrator.eval();
                }
                else {
                    //std::cout << "integrating...\n";
                    // apply standard non-singular quadrature
                    bemquad::AdaptiveClenshaw2DIntegrator<double, F> integrator(s_interval.first,
                                                                                s_interval.second,
                                                                                t_interval.first,
                                                                                t_interval.second,
                                                                                functor,
                                                                                tol);
                    //std::cout << "inserting...\n";
                    result[gindex] += -1.0 * integrator.eval();
                    
                }
            }
            std::cout << "Work item: " << icpt << " in range: [" << lower << "," << upper - 1 << "] complete\n";
        }
    }
    
    
    template<typename F>
    std::vector<double> multithreadBEMOperator(const Forest& f,
                                               const std::vector<double>& bcdata,
                                               const double tol)
    {
        assert(bcdata.size() == f.globalDofN());
        
        // This is where the magic happens.
        std::vector<std::pair<uint, uint>> interior_cpts, boundary_cpts;
        std::set<uint> boundarycpt_cache; // set of previously inserted global colloc. pt. indices
        for(uint icspace = 0; icspace < f.spaceN(); ++icspace) {
            const BSplineSpace& cspace = f.space(icspace);
            for(uint icolloc = 0; icolloc < cspace.grevilleAbscissaPtN(); ++icolloc) {
                const uint gindex = f.globalCollocI(icspace, icolloc);
                if(interiorLocalI(cspace, icolloc))
                    interior_cpts.push_back(std::make_pair(icspace, icolloc));
                else {
                    const auto find = boundarycpt_cache.find(gindex);
                    if(find != boundarycpt_cache.end())
                        continue;
                    boundary_cpts.push_back(std::make_pair(icspace, icolloc));
                    boundarycpt_cache.insert(gindex);
                }
            }
        }
        
        std::vector<double> result(f.globalDofN(), 0.0); // we store the operator [A]{x} here
        
        // Determine available threads on current hardware
        uint nthreads = std::thread::hardware_concurrency();
        std::cout << "Computing with " << nthreads << " threads\n";
        std::vector<std::thread> threads;
        
        // Split work according to number of threads
        std::vector<long int> limits = bounds(nthreads, interior_cpts.size());
        std::cout << "Processing cpts interior to parametric spaces....\n";
        for(uint i = 0; i < nthreads; ++i)
            threads.push_back(std::thread(workfunction<F>,
                                          std::ref(f),
                                          std::ref(bcdata),
                                          std::ref(interior_cpts),
                                          tol,
                                          limits[i],
                                          limits[i+1],
                                          std::ref(result)));
        // And launch the threads
        for(auto& t: threads)
            t.join();
        
        std::cout << "Now processing cpts on space edges....\n";
        threads.clear();
        limits = bounds(nthreads, boundary_cpts.size());
        for(uint i = 0; i < nthreads; ++i)
            threads.push_back(std::thread(workfunction<F>,
                                          std::ref(f),
                                          std::ref(bcdata),
                                          std::ref(boundary_cpts),
                                          tol,
                                          limits[i],
                                          limits[i+1],
                                          std::ref(result)));
        // And launch the threads
        for(auto& t: threads)
            t.join();
        
        return result;
    }
    
    
    /// This function is operated on my a single thread. It contributes to the calculation
    /// of the vector 'result'.
    template<typename F>
    void preconditionerWorkFunction(const Forest& f,
                                    const std::vector<std::tuple<uint, GPt2D, uint>>& cptdata,
                                    const uint lower,
                                    const uint upper,
                                    const double tol,
                                    Epetra_CrsMatrix& pmatrix)
    {
        for(uint icpt = lower; icpt < upper; ++icpt) { // loop over global collocation points
            const uint cindex = std::get<0>(cptdata[icpt]);
            const GPt2D cpt_parent = std::get<1>(cptdata[icpt]);
            const uint ielem = std::get<2>(cptdata[icpt]);
//            std::cout << "Generating singular terms for global cindex: " << cindex << " global element index: " << ielem << "\n";
            const auto el = f.element(ielem);
            std::vector<int> gindices;
            for(const auto& ibasis : el->globalBasisFuncI())
                gindices.push_back(static_cast<int>(ibasis));
            std::vector<double> vals(gindices.size(), 0.0);
            F functor(el, el->eval(cpt_parent));
            for(uint ibasis = 0; ibasis < el->basisFuncN(); ++ibasis) {
                functor.setLocalBasisI(ibasis);
                bemquad::Adaptive2DPolarClenshawIntegrator<double, F> integrator(-1.0,1.0,-1.0,1.0,
                                                                                 cpt_parent.s,
                                                                                 cpt_parent.t,
                                                                                 functor, tol);
                // HARDCODED -1.0 TO SHIFT SL TERMS TO LHS.
                vals[ibasis] += -1.0 * integrator.eval();
            }
            
            pmatrix.InsertGlobalValues(cindex, gindices.size(), vals.data(), gindices.data());
        }
    
    }
    
    /// A multithread implementation which generates a preconditioning matrix
    /// for a GMRES solver. This is approximated by computing the terms related to
    /// the global basis functions whose span intersects the relevant collocation point
    template<typename F>
    void multithreadGeneratePreconditioner(const Forest& f,
                                           const double tol,
                                           Epetra_CrsMatrix& pmatrix)
    {
        // First generate the set of collocation point/element pairings that we wish to evaluate
        std::vector<std::tuple<uint, GPt2D, uint>> cptdata;
        
        std::set<uint> cached_cpts; // map to determine whether cpt has been previous calculated.
        
        // Now fill the vector with the relevant data
        for(uint icspace = 0; icspace < f.spaceN(); ++icspace) {
            const BSplineSpace& cspace = f.space(icspace);
            for(uint c = 0; c < cspace.grevilleAbscissaPtN(); ++c) {
                const uint glb_cindex = f.globalCollocI(icspace, c);
                const auto search = cached_cpts.find(glb_cindex);
                if(search != cached_cpts.end()) // if point has been cached, continue
                    continue;
                for(uint ifspace = 0; ifspace < f.spaceN(); ++ifspace) { // loop over field spaces
                    const BSplineSpace& fspace = f.space(ifspace);
                    const auto connected = f.connectedCollocPtI(icspace, c, ifspace);
                    if(connected.first) { //
                        const GPt2D cpt_param = fspace.grevilleAbscissaPt(connected.second);
                        for(uint ielem = 0; ielem < f.elemN(ifspace); ++ielem) {
                            const auto el = f.element(ifspace, ielem);
                            const auto pair = el->containsParamPt(cpt_param.s, cpt_param.t);
                            if(!pair.first)
                                continue; // non-singular
                            cptdata.push_back(std::make_tuple(glb_cindex, pair.second, f.globalElI(ifspace, ielem)));
                        }
                    }
                }
                cached_cpts.insert(glb_cindex);
            }
        }

        // Determine available threads on current hardware
        const uint nthreads = std::thread::hardware_concurrency();
        std::cout << "Computing preconditioning matrix with " << nthreads << " threads\n";
        std::vector<std::thread> threads;
        
        // Split work according to number of threads
        std::vector<long int> limits = bounds(nthreads, cptdata.size());
        for(uint i = 0; i < nthreads; ++i)
            threads.push_back(std::thread(preconditionerWorkFunction<F>,
                                          std::ref(f),
                                          std::ref(cptdata),
                                          limits[i],
                                          limits[i+1],
                                          tol,
                                          std::ref(pmatrix)));
        for(auto& t : threads)
            t.join();
        
        pmatrix.FillComplete();
        
    }
    
}


#endif

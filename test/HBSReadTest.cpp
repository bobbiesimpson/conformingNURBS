#include <iostream>
#include <chrono>
#include <iomanip>
#include "Geometry.h"
#include "Forest.h"
#include "HConformingForest.h"
#include "IElemIntegrate.h"

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
        
        // Construct the necessary forests
        Forest forest(g);
        HDivForest multiforest(g);
        uint refine = 0;
        if(argc > 2)
            refine = std::atoi(argv[2]);
        forest.hrefine(refine);
        multiforest.hrefine(refine);
        
        std::cout << "Constructed multi forest with " << multiforest.elemN() << " elements\n";
        
        auto start = std::chrono::high_resolution_clock::now();
        double area = 0.0;
        for(uint ielem = 0; ielem < multiforest.elemN(); ++ielem)
        {
            const auto p_el = multiforest.element(ielem);
            const auto p_bel =  multiforest.bezierElement(ielem);
            for(IElemIntegrate igpt(p_el->equalIntegrationOrder()); !igpt.isDone(); ++igpt)
            {
                const auto& gp = igpt.get();
//                const auto& xf = p_bel->eval(gp);
                const auto& t1 = p_bel->tangent(gp, ParamDir::S);
                const auto& t2 = p_bel->tangent(gp, ParamDir::T);

//                const auto& normal = cross(t1, t2).asNormal();
                
                const auto& basis = p_bel->basis(gp.s, gp.t);
                for(const auto& b : basis)
                    std::cout << b << "\n";
                
                std::cout << "correct basis\n\n";
                
                const auto& correctbasis = p_el->basis(gp.s, gp.t);
                for(const auto& b : correctbasis)
                    std::cout << b << "\n";
                
                const auto& basisder = p_bel->localBasisDers(gp.s, gp.t, DerivType::DS);
                for(const auto& b : basisder)
                    std::cout << b << "\n";
                
                std::cout << "correct basis fn derivatives\n\n";
                
                const auto& correctbasisder = p_el->localBasisDers(gp.s, gp.t, DerivType::DS);
                for(const auto& b : correctbasisder)
                    std::cout << b << "\n";
                
                area += p_bel->jacDet(gp.s, gp.t, t1, t2) * igpt.getWeight();
            }
            
        }
        
        std::cout << "Area = " << area << "\n";
        
        auto time =  std::chrono::high_resolution_clock::now() - start;
        std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(time).count() << "(s)" << "\n";
        
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

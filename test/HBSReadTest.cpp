#include <iostream>
#include <chrono>
#include "Geometry.h"
#include "Forest.h"
#include "IPolarIntegrate.h"

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
        uint refine = 0;
        if(argc > 2)
            refine = std::atoi(argv[2]);
        forest.hrefine(refine);
        
        std::cout << "Constructed forest with " << forest.elemN() << " elements\n";
        
        std::cout << "Read successful\n";
        
        auto start = std::chrono::high_resolution_clock::now();
        for(uint ielem = 0; ielem < forest.elemN(); ++ielem) {
            const auto el = forest.element(ielem);
            for(nurbs::IPolarIntegrate igpt(nurbs::GPt2D(0.0, 0.0), el->integrationOrder()); !igpt.isDone(); ++igpt) {
                const auto gp = igpt.get();
                const auto p = el->eval(gp.s, gp.t);
                auto n = el->normal(gp.s, gp.t);
//                auto basis = el->basis(gp.s, gp.t);
                //auto jacob = el->jacDet(gp.s, gp.t);
            }
        }
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

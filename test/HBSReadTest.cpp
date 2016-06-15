#include <iostream>

#include "Geometry.h"
#include "Forest.h"

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
        std::cout << "Performing h-refinement test...\n";
        const uint nrefine = 10;
        forest.hrefine(nrefine);
        
        std::cout << "Constructed forest with " << forest.elemN() << " elements\n";
        
        std::cout << "Read successful\n";
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

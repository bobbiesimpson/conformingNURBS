#include <iostream>

#include "Geometry.h"
#include "Forest.h"
#include "BoundingBoxIterator.h"

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
        forest.hrefine(0);
        std::cout << "Constructed forest with " << forest.elemN() << " elements\n";
        
        std::cout << "Creating bounding box data for H-matrix construction\n";
        
        double max = std::numeric_limits<double>::lowest();
        double min = std::numeric_limits<double>::max();
        for(BoundingBoxIterator it(forest); !it.isDone(); ++it) {
            std::cout << "Basis function index : " << it.currentIndex() << "\n";
            const double length = dist(it.currentLowerBound(), it.currentUpperBound());
            if(length > max) max = length;
            if(length < min) min = length;
            std::cout << length << "\n";
        }
        
        std::cout << "min = " << min << " max = " << max << "\n";
        std::cout << "ratio = " << max / min << "\n";
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

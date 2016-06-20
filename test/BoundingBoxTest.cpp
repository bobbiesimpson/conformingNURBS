#include <iostream>

#include "Geometry.h"
#include "Forest.h"
#include "BoundingBoxIterator.h"
#include "OutputVTK.h"

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
        forest.hrefine(1);
        
        for(uint ielem = 0; ielem < forest.elemN(); ++ielem) {
            const auto el = forest.element(ielem);
            for(uint icpt = 0; icpt < el->collocPtN(); ++icpt){
                std::cout << el->collocParentCoord(icpt) << "\n";
                
            }
        }
        std::cout << "Constructed forest with " << forest.elemN() << " elements\n";
        
        std::cout << "Creating bounding box data for H-matrix construction\n";
        
        std::vector<std::pair<Point3D, Point3D>> bbdata;
        
        for(BoundingBoxIterator it(forest); !it.isDone(); ++it) {
            bbdata.push_back(it.currentBoundingBox());
        }
        
        // Output the bounding boxes to a vtu file
        OutputVTK output("sphere_test");
        output.outputGeometry(forest);
        output.outputBoundingBoxSet(bbdata);
        
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

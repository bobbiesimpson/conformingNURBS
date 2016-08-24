#include <iostream>
#include <fstream>
#include <chrono>
#include "Forest.h"
#include "Geometry.h"
#include "BSplineSpace.h"
#include "BezierNodalElement.h"
#include "NodalElement.h"
#include "IElemIntegrate.h"
#include "OutputVTK.h"
#include "HConformingForest.h"

using namespace nurbs;

int main(int argc, char* argv[]) {
    
    try {
        
        std::cout << "Running degenerate geometry test.....\n";

        std::cout << "Trying to open hbs input file....\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            error("Cannot open file for reading\n");
        Geometry g;
        if(!g.loadHBSFile(ifs))
            error("Failed to load geometry from hbs data");
        
        Forest forest(g);
        for(uint ispace = 0; ispace < forest.spaceN(); ++ispace)
        {
            for(uint iedge = 0; iedge < NEDGES; ++iedge)
            {
                if(forest.degenerateEdge(ispace, iedge))
                    std::cout << "(Space,edge) = (" << ispace << "," << iedge << ") is degenerate\n";
            }
        }
        
        HDivForest divforest(g);
        
        for(uint ielem = 0; ielem < divforest.elemN(); ++ielem)
        {
            const auto el = divforest.bezierElement(ielem);
            std::cout << el->globalBasisFuncI() << "\n\n";
            std::cout << el->signedGlobalBasisFuncI() << "\n\n";
        }
        
        OutputVTK output("degenerate_geometry");
        output.outputGeometry(forest);
        
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

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Geometry.h"
#include "Forest.h"
#include "HConformingForest.h"
#include "OutputVTK.h"
#include "IElemIntegrate.h"
#include "base.h"
#include "IPolarIntegrate.h"

using namespace nurbs;

int main(const int argc, const char* argv[])
{
    try
    {
        std::cout << "Running parametric quadrature test.....\n";
        
        std::cout << "Trying to open hbs input file....\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            nurbs::error("Cannot open file for reading\n");
        
        // First create a geometry object
        nurbs::Geometry g;
        if(!g.loadHBSFile(ifs))
            nurbs::error("Failed to load geometry from hbs data");
        
        // And now create a conforming forest object
        nurbs::Forest forest(g);
        nurbs::HDivForest divforest(g);
        std::cout << "contructed hdiv forest with " << divforest.elemN() << " elements\n";
        
        const auto pel = divforest.bezierElement(0);
        
        // integrate over parametric domain
        uint max_order = 50;
        
        for(uint iorder = 1; iorder < max_order; ++iorder)
        {
            double val = 0.0;
            for(nurbs::IPolarIntegrate igpt(nurbs::GPt2D(0.0), iorder); !igpt.isDone(); ++igpt)
            {
                const double jdet = pel->jacDet(igpt.get());
                val +=  jdet * igpt.getWeight();
            }
            std::cout << iorder << ":\t" << std::setprecision(15) <<  val << "\n";
        }
        

        return EXIT_SUCCESS;
        
    } catch (const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
}
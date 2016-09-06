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
#include "ISubElemIntegrate.h"

using namespace nurbs;

int main(int argc, char* argv[]) {
    
    try {
        
        std::cout << "Running Bezier extraction test.....\n";
        
        //    BSplineSpace space({{0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0},
        //                        {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0} },
        //                        {3,3});
        //
        //    BSplineSpace space2({{0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0},
        //                        {0.0, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0} },
        //                       {3,3});
        //
        
        std::cout << "Trying to open hbs input file....\n";
        std::ifstream ifs(argv[1]);
        if(!ifs)
            error("Cannot open file for reading\n");
        Geometry g;
        if(!g.loadHBSFile(ifs))
            error("Failed to load geometry from hbs data");
        
        Forest forest(g);
//        OutputVTK output("almond");
//        output.outputGeometry(forest);
        
        
//        forest.hrefine(8);
        
        std::cout   << "Running Bezier test on forest with " << forest.elemN()
                    << " elements and " << forest.globalDofN() << " dof\n";
        
        auto start = std::chrono::high_resolution_clock::now();
        double a = 0.0;
        for(uint ielem = 0; ielem < forest.elemN(); ++ielem) {
            //std::cout << "Element: " << ielem << "\n";
            //const auto el = forest.element(ielem);
            const auto b_el = forest.bezierElement(ielem);
            for(ISubElemIntegrate igpt(b_el->integrationOrder(), {10,10}); !igpt.isDone(); ++igpt) {
                const auto gp = igpt.get();
                std::cout << gp << "\n";
                const Point3D t1 = b_el->tangent(gp, S);
                const Point3D t2 = b_el->tangent(gp, T);
                const Point3D m = cross(t1,t2);
                const Point3D n = m.asNormal();
                const auto jdet = m.length() * b_el->jacDetParam(gp.s, gp.t);
                const auto x = b_el->eval(gp);

                a += jdet * igpt.getWeight();
            }
        }
        auto time =  std::chrono::high_resolution_clock::now() - start;
        std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(time).count() << "(s)" << "\n";
        std::cout << "area = " << a << "\n";
        
        
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

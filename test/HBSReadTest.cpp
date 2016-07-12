#include <iostream>
#include <chrono>
#include "Geometry.h"
#include "NodalElement.h"
#include "Forest.h"
#include "HConformingForest.h"
#include "IRegularQuadrature.h"
#include "IEdgeQuadrature.h"
#include "IVertexQuadrature.h"
#include "IEqualQuadrature.h"
#include "base.h"
#include <iomanip>

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
        HDivForest multiforest(g);
        uint refine = 0;
        if(argc > 2)
            refine = std::atoi(argv[2]);
        forest.hrefine(refine);
        multiforest.hrefine(refine);
        
        std::cout << "Constructed forest with " << forest.elemN() << " elements\n";
        
        std::cout << "Read successful\n";
        
        auto start = std::chrono::high_resolution_clock::now();

        const uint n = 5;
        const UIntVec sorder{n,n};
        const UIntVec forder{n,n};
        
        std::vector<double> data;
        
        for(uint ielem_s = 0; ielem_s < multiforest.elemN(); ++ielem_s)
        {
            const auto p_srcel = multiforest.element(ielem_s);
            
            for(uint ielem_f = 0; ielem_f < multiforest.elemN(); ++ielem_f)
            {
                double result = 0.0;
                const auto p_fieldel = multiforest.element(ielem_f);
                Edge esrc, efield;
                if(edgeConnected(*p_srcel, *p_fieldel, esrc, efield))
                {
                    for(IEdgeQuadrature iedge(sorder, forder, esrc, efield); !iedge.isDone(); ++iedge)
                    {
                        const auto& gpt = iedge.get();
                        const auto w = iedge.getWeight();
                        const auto xs = p_srcel->eval(gpt.srcPt());
                        const auto xf = p_fieldel->eval(gpt.fieldPt());
                        
                        const auto r = dist(xs, xf);
                        result += 1.0 / (4.0 * PI * r) * w * p_srcel->jacDet(gpt.srcPt()) * p_fieldel->jacDet(gpt.fieldPt());
                        
                        
                    }
                    if(esrc == Edge::EDGE2 && efield == Edge::EDGE0)
                        std::cout << result << "\n";
                    data.push_back(result);
                }
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

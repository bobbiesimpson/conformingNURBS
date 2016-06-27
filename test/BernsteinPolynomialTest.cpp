#include <iostream>
#include <fstream>

#include "NURBSCommon.h"

using namespace nurbs::nurbshelper;

int main(int argc, char* argv[]) {
    
    const size_t maxpts = 1e7;
    const size_t default_ptn = 1e3;
    const uint default_p = 4;
    
    size_t n = default_ptn;
    
    if(argc > 1) {
        const size_t input = std::atoi(argv[1]);
        if(input > 1 & input < maxpts)
            n = input;
        else
            std::cout << "Setting default to " << default_ptn << "\n";
    }
    
    size_t p = default_p;
    if(argc > 2) {
        const size_t input = std::atoi(argv[2]);
        if(input > 0)
            p = input;
        else
            std::cout << "Setting default degree to " << default_p << "\n";
    }
    
    std::ofstream ofs("bernstein_basis.dat");
    if(!ofs) {
        std::cerr << "Bad file opening\n";
        return EXIT_FAILURE;
    }
    
    std::cout << "Computing degree " << p << " Bernstein basis functions with " << n << " sample points....\n\n";
    long double xi = -1.0;
    const long double inc = 2.0 / (n - 1);
    while(xi < 1.0) {
        const auto bernsteinbasis = bernsteinPolynomial(xi, p);
        ofs << xi << "\t";
        for(const auto& b : bernsteinbasis)
            ofs << b << "\t";
        ofs << "\n";
        xi += inc;
    }
    ofs.close();
    
    std::ofstream ofs2("bspline_basis.dat");
    if(!ofs2) {
        std::cerr << "Bad file opening\n";
        return EXIT_FAILURE;
    }
    
    std::cout << "Computing degree " << p << " Bspline basis functions with " << n << " sample points....\n\n";
    long double s = 0.0;
    const long double inc2 = 1.0 / (n - 1);
    const std::vector<double> kvec{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    while(s < 1.0) {
        
        const auto basis = getBsplineBasis(s, 4, kvec, 4);
        ofs2 << s << "\t";
        for(const auto& b : basis)
            ofs2 << b << "\t";
        ofs2 << "\n";
        s += inc2;
    }
    
    return EXIT_SUCCESS;
    
}
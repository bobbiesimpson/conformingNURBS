#include <iostream>
#include "base.h"
#include "NURBSCommon.h"

#include <map>

using namespace nurbs;

int main(int argc, char* argv[])
{
    
    
    // The knot vector we are going to decompose into Bezier elements
    const DoubleVec kvec{0,0,0,1,1, 2,3,4,4, 5, 5, 5};//{0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    const uint p = 2;
    
    // create unique knot vectors
    auto uniqueknotvec = kvec;
    auto last = std::unique(uniqueknotvec.begin(), uniqueknotvec.end());
    uniqueknotvec.erase(last, uniqueknotvec.end());
    
    const size_t uniqueknotn = uniqueknotvec.size();
    const uint m = kvec.size();
    
    // now create a map of all the elements which are already in bezier
    // form and the index to resume the extraction algorithm.
    
    std::map<unsigned, unsigned> beziermap;
    
    for(unsigned i = 0; i < uniqueknotvec.size() - 1; ++i)
    {
        const auto lowerkval = uniqueknotvec[i];
        const auto upperkval = uniqueknotvec[i+1];
        const auto lowercount = std::count_if(kvec.begin(), kvec.end(), [&](double k) {return essentiallyEqual(k, lowerkval, 1.0e-5);});
        const auto upperount = std::count_if(kvec.begin(), kvec.end(), [&](double k) {return essentiallyEqual(k, upperkval, 1.0e-5);});
        
        if(lowercount >=p && upperount >= p)
        {
            if(i == uniqueknotvec.size() - 2)
                beziermap[i] = nurbshelper::getKnotSpan(uniqueknotvec[i+1], kvec, p) + p + 1 - upperount + 1;
            else
                beziermap[i] = nurbshelper::getKnotSpan(uniqueknotvec[i+1], kvec, p) - upperount + 1;
        }
    }
    
    
    uint a = p;
    uint b = a + 1;
    uint nb = 0;                // # bezier elements in this direction
    
    /// Identity matrix
    std::vector<std::vector<double>> I(p + 1, std::vector<double>(p+1, 0.0));
    for(uint i = 0; i < p + 1; ++i)
        I[i][i] = 1.0;
    
    // if linear, the extraction operator is simply the identity matrix
    // since linear Bsplines and linear Bernstein polynomials are equivalent
    if(1 == p) {
        for(uint iel = 0; iel < uniqueknotn-1; ++iel)
            std::cout << "element: " << iel << "\t" <<  I << "\n";
        return EXIT_SUCCESS;

    }
    
    // and do a check that the knot vector isn't already in bezier form
//    unsigned minrepeats = std::numeric_limits<unsigned>::max();
//    unsigned currentcount = 1;
//    for(size_t i = 0; i < kvec.size() - 1; ++i)
//    {
//        if(essentiallyEqual(kvec[i], kvec[i+1], 1.e-4))
//            ++currentcount;
//        else
//        {
//            if(currentcount < minrepeats)
//                minrepeats = currentcount;
    
//            currentcount = 1;
//        }
//    }
//    // if already in bezier form, set the identity matrix for each extraction operator.
//    if(minrepeats >= p)
//    {
//        for(uint iel = 0; iel < uniqueknotn-1; ++iel)
//            std::cout << "element: " << iel << "\t" <<  I << "\n";
//        return EXIT_SUCCESS;
//    }
    
    // code for degree > 1
    
    /// Initialise extraction matrices
    auto Ccurrent = I;
    auto Cnext = I;
    
    while(b < m) {
        
        // We're done.
        if(nb == uniqueknotn - 2)
        {
            std::cout << "element: " << nb << "\t" <<  Ccurrent << "\n";
            break;
        }
        
        // check if element is already in bezier form
//        auto search = beziermap.find(nb);
//        if(search != beziermap.end())
//        {
//            std::cout << "element: " << nb << "\t" <<  I << "\n";
//            Ccurrent = I;
//            Cnext = I;
//            b = search->second;
//            nb += 1;
//            continue;
//        }
        
        uint i = b;
        uint mult = 0;
        
        // count multiplitiy of knot at location b
        while(b < m && essentiallyEqual(kvec[b+1], kvec[b], 1.e-4))
            b += 1;
        
        mult = b - i + 1;
        
        // initialise alpha
        std::vector<double> alphas(p + 1, 0.0);
        
        if(mult < p)
        {
            const double numer = kvec[b] - kvec[a];
            for(uint j = p; j > mult; j--)
                alphas[j-mult-1] = numer / (kvec[a+j] - kvec[a]);
            
            const uint r = p - mult;
            
            // update matrix coefficients
            for(uint j = 1; j <= r; ++j) {
                const uint save = r - j;
                const uint s = mult + j;
                for(uint k = p; k >= s; k--) {
                    const double alpha = alphas[k-s];
                    for(uint irow = 0; irow < Ccurrent.size(); ++irow) {
                        Ccurrent[irow][k] = alpha * Ccurrent[irow][k]
                        + (1.0 - alpha) * Ccurrent[irow][k-1];
                    }
                }
                if(b < m )
                {
                    for(uint i = 0; i <= j; ++i)
                    {
                        const double val = Ccurrent[p-j+i][p];
                        //if(!essentiallyEqual(val, 0.0, 1.0-7))
                            Cnext[save+i][save] = val;
                    }
                }
            }
            
            // store current operator
            std::cout << "element: " << nb << "\t" <<  Ccurrent << "\n";
            
            // Set next operator equal to current and increment element counter
            Ccurrent = Cnext;
            Cnext = I;
            nb += 1;
            
            // update indices for next operator
            if(b < m) {
                a = b;
                b += 1;
            }
        }
    }

    return EXIT_SUCCESS;
}


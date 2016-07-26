#include <iostream>
#include "base.h"
#include "NURBSCommon.h"
#include <set>

#include <map>

using namespace nurbs;

int main(int argc, char* argv[])
{
    
    
    // The knot vector we are going to decompose into Bezier elements
    const DoubleVec kvec{0,0,0,0,1, 2, 3, 4,5, 5, 5, 5};//{0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    const uint p = 3;
    
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
    
    
    uint a = p + 1;
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
        
        Cnext = I;
        
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
        while(b < m && essentiallyEqual(kvec[b], kvec[b-1], 1.e-4))
            b += 1;
        
        mult = b - i + 1;
        
        // initialise alpha
        std::vector<double> alphas(p + 1, 0.0);
        
        if(mult < p)
        {
            const double numer = kvec[b-1] - kvec[a-1];
            for(uint j = p; j > mult; j--)
                alphas[j-mult-1] = numer / (kvec[a+j-1] - kvec[a-1]);
            
            const uint r = p - mult;
            
            // update matrix coefficients
            for(uint j = 1; j <= r; ++j) {
                const uint save = r - j;
                const uint s = mult + j;
                for(uint k = p+1; k > s; k--) {
                    const double alpha = alphas[k-s-1];
                    for(uint irow = 0; irow < Ccurrent.size(); ++irow)
                    {
                        Ccurrent[irow][k-1] = alpha * Ccurrent[irow][k-1]
                        + (1.0 - alpha) * Ccurrent[irow][k-2];
                    }
                }
                if(b < m)
                {
                    for(uint i = 0; i < j+1; ++i)
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
            nb += 1;
            
            // update indices for next operator
            if(b < m) {
                a = b;
                b += 1;

            }
        }
    }

    
    // let's compute the global matrix explicitly.
    
    // first figure out what knots need to be inserted to obtain the requried
    // bezier form
    
    typedef std::vector<std::vector<double>> Matrix;

    size_t currentindex = 0;
    size_t insertedknot_n = 0;
    
    Matrix global_cmatrix;
    auto kvec_copy = kvec;
    
    const auto n = kvec_copy.size() - p - 1;
    
    // keep going until we have C^0 continuity at all knots
    while(currentindex < uniqueknotvec.size())
    {
        const auto kval = uniqueknotvec[currentindex];
        const auto count = std::count_if(kvec_copy.begin(), kvec_copy.end(), [&](double k) {return essentiallyEqual(k, kval, 1.0e-5);});
        auto requiredknots = p - count;
        
        // perform required knot insertion
        while(requiredknots > 0)
        {
            const unsigned k = nurbshelper::getKnotSpan(kval, kvec_copy, p);
            std::cout << requiredknots << " knots required at " << kval << " at knot span " <<  k << "\n";
            
            // initialise Cmatrix
            Matrix Cmatrix;
            for(size_t i = 0; i < n + insertedknot_n; ++i)
                Cmatrix.push_back(std::vector<double>(n + insertedknot_n + 1, 0.0));
            
            for(unsigned a = 0; a < n + insertedknot_n + 1; ++a)
            {
                // First compute alpha
                double alpha;
                if(a <= k - p - 1)
                    alpha = 1.0;
                else if(a >= k - p && a <= k - 1)
                {
                    double numer = kval - kvec_copy[a];
                    alpha = numer / (kvec_copy[a + p] - kvec_copy[a]);
                }
                else
                    alpha = 0.0;
                
                // now populate the matrix for the current knot insertion
                if(a == 0)
                    Cmatrix[a][a] = alpha;
                else if(a == n + insertedknot_n)
                    Cmatrix[a-1][a] = 1 - alpha;
                else
                {
                    Cmatrix[a][a] = alpha;
                    Cmatrix[a-1][a] = 1- alpha;
                }
            }
            
            if(insertedknot_n == 0)
                global_cmatrix = Cmatrix;
            else
            {
                // now update the global matrix
                auto global_cmatrix_copy = global_cmatrix;
                
                global_cmatrix.clear();
                // resize global matrix
                for(size_t i = 0; i < n; ++i)
                    global_cmatrix.push_back(std::vector<double>(n + insertedknot_n + 1, 0.0));
                
                for(size_t i = 0; i < n; ++i)
                    for(size_t j = 0; j < n + insertedknot_n + 1; ++j)
                        for(size_t k = 0; k < n + insertedknot_n; ++k)
                            global_cmatrix[i][j] += global_cmatrix_copy[i][k] * Cmatrix[k][j];
            }
            // insert the new knot into the knot vector
            kvec_copy.insert(kvec_copy.begin() + k , kval);
            requiredknots--;
            ++insertedknot_n;
        }
        ++currentindex;
    }
    
    // print out cmatrix
    for(size_t i = 0; i < global_cmatrix.size(); ++i)
    {
        for(const auto& v : global_cmatrix[i])
            std::cout << v << "\t";
        std::cout << "\n";
    }
    
    return EXIT_SUCCESS;
}


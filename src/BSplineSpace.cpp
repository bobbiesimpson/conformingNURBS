#include "BSplineSpace.h"
#include "NURBSCommon.h"
#include "base.h"

#include <algorithm>
#include <stdexcept>

namespace nurbs 
{

	DoubleVecVec BSplineSpace::tensorBasis(const double s, const double t, const UIntVec& span) const
	{
		return {nurbshelper::getBsplineBasis(s, span[S], knotVec(S), degree(S)),
				nurbshelper::getBsplineBasis(t, span[T], knotVec(T), degree(T))};
	}
	
	DoubleVecVec BSplineSpace::tensorBasisDers(const double s, const double t,
											   const UIntVec& span, const DerivOrder der) const
	{
		return {nurbshelper::getBsplineBasisDers(s, span[S], knotVec(S), degree(S), der).at(der),
				nurbshelper::getBsplineBasisDers(t, span[T], knotVec(T), degree(T), der).at(der)};
	}
	
	UIntVecVec BSplineSpace::localBasisFuncI(const double s, const double t) const
	{
		return {nurbshelper::getBasisFnIndices(s, knotVec(S), degree(S)),
				nurbshelper::getBasisFnIndices(t, knotVec(T), degree(T))};
	}
	
	UIntVec BSplineSpace::span(const double s, const double t) const
	{
		return {nurbshelper::getKnotSpan(s, knotVec(S), degree(S)),
				nurbshelper::getKnotSpan(t, knotVec(T), degree(T)) };
	}
    
    void BSplineSpace::degreeReduce(const ParamDir dir)
    {
        if(0 == mDegrees[dir])
            throw std::runtime_error("Cannot degree reduce degree 0 spline\n");
        auto& knotvec = mKnotVecs[dir];
        std::vector<double> kv_reduced;
        std::copy(knotvec.begin() + 1, knotvec.end() - 1, std::back_inserter(kv_reduced));
        knotvec = kv_reduced;
        
//        auto it = knotvec.begin();
//        for(uint i = 0; i < uniqueKnotN(dir); ++i) {
//            const double knot = uniqueKnot(i, dir);
//            it = std::find(it, knotvec.end(), knot);
//            knotvec.erase(it); // only remove first and last knots
//        }
        mDegrees[dir] -= 1; // reduce degree
        init();
    }
    
    void BSplineSpace::degreeElevate(const ParamDir dir)
    {
        auto& knotvec = mKnotVecs[dir];
        for(uint i = 0; i < uniqueKnotN(dir); ++i) {
            const double knot = uniqueKnot(i, dir);
            knotvec.push_back(knot);
        }
        mDegrees[dir] += 1; // elevate degree
        std::sort(knotvec.begin(), knotvec.end());
        init();
    }
    
    void BSplineSpace::hrefine(const uint n)
    {
        for(auto& kv : mKnotVecs)
            kv = nurbshelper::uniformKnotInsertion(kv,n);
        init();
    }

	void BSplineSpace::load(std::istream& ist)
	{
		clear();
		char ch;
		if(ist >> ch && ch != '{') {
			ist.unget();
			ist.clear(std::ios_base::failbit);
			return;
		}
		std::string s;
		ist >> s;
		if(!ist )
			error("Cannot read b-spline space");
		mName = s;
		InputVec<uint> d_vec;
		if(!(ist >> d_vec))
			error("Bad degree vector");
		mDegrees = d_vec.data;
		if(ist >> ch && ch != '{') 
			error("Bad knot vector reading");
		while(true) {
			InputVec<double> k_vec;
			if(!(ist >> k_vec))
				break;
			mKnotVecs.push_back(k_vec.data);
		}
		endOfLoop(ist, '}', "Bad knot vector.");
		init(); // recalculate unique knots + intervals
	}
	
    void BSplineSpace::init()
    {
        mUniqueKnotVecs.clear(); // clear data
        mIntervalVec.clear();
        mNumBasisVec.clear();
        mSExtractionOperators.clear();
        mTExtractionOperators.clear();
        
        for(uint i = 0; i < paramDimN(); ++i) {
            mIntervalVec.push_back(boost::icl::construct<Interval>
                                   (mKnotVecs[i].front(), mKnotVecs[i].back(),
                                    boost::icl::interval_bounds::closed()));
            Interval& last = mIntervalVec.back();
            if(essentiallyEqual(last.upper(), last.lower(), TOL))
                error("Cannot prescribe a b-spline space with non-zero parametric area.");
                mNumBasisVec.push_back(mKnotVecs[i].size() - mDegrees[i] - 1);
                }
        
        // construct unique knot vectors
        for(auto kv : mKnotVecs) {
            auto last = std::unique(kv.begin(), kv.end());
            mUniqueKnotVecs.emplace_back( DoubleVec( kv.begin(), last ) );
        }
        
        //
        // construct element extraction operators from knot vectors
        //
        for(uint iparam = 0; iparam < paramDimN(); ++iparam) {
            
            const ParamDir dir = ParamDirType(iparam);
            const uint p = degree(dir);
            auto kvec = knotVec(dir);
            
            const uint m = knotVec(dir).size();
            
            uint a = p;
            uint b = a + 1;
            uint nb = 0;                // # bezier elements in this direction
            
            /// Identity matrix
            std::vector<std::vector<double>> I(p + 1, std::vector<double>(p+1, 0.0));
            for(uint i = 0; i < p + 1; ++i)
                I[i][i] = 1.0;
            
            // if linear, the extraction operator is simply the identity matrix
            // since linear Bsplines and linear Bernstein polynomials are equivalent
            if(1 == degree(dir)) {
                for(uint iel = 0; iel < uniqueKnotN(dir)-1; ++iel)
                    setExtractionOperator(iel, dir, I);
                continue;
            }
            
            // code for degree > 1
            
            /// Initialise extraction matrices
            auto Ccurrent = I;
            auto Cnext = I;
            
            while(b < m) {

                // We're done.
                if(nb == uniqueKnotN(dir) - 2) {
                    setExtractionOperator(nb, dir, Ccurrent);
//                    std::cout << Ccurrent << "\n";
                    break;
                }
                
                uint i = b;
                uint mult = 0;
                
                // count multiplitiy of knot at location b
                while(b < m && essentiallyEqual(kvec[b+1], kvec[b], 1.e-7))
                    b += 1;
                
                mult = b - i + 1;
                
                // initialise alpha
                std::vector<double> alphas(p + 1, 0.0);
                
                if(mult < p) {
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
                        if(b < m) {
                            for(uint i = 0; i <= j; ++i)
                                Cnext[save+i][save] = Ccurrent[p-j+i][p];
                        }
                    }
                    
                    // store current operator
                    setExtractionOperator(nb, dir, Ccurrent);
//                    std::cout << Ccurrent << "\n";
                    
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
        }
    }
    
	void BSplineSpace::printData(std::ostream& ost) const
	{
		ost << "--------- NURBS SPACE DATA -------------\n\n";
		for(uint i = 0; i < paramDimN(); ++i) {
			ost << "Parameteric direction: " << i << "\n";
			ost << "Degree: " << mDegrees[i] << "\n";
			ost << "Knot vector: ";
			printVector(mKnotVecs[i], ost);
			ost << "\n";
			ost << "Unique knot vector: ";
			printVector(mUniqueKnotVecs[i], ost);
			ost << "\n---------\n";
			
		}
		ost << "\n--------------------------------------\n";				
	}

}

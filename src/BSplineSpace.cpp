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

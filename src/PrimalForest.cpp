#include <ios>
#include <cassert>

#include "NURBSCommon.h"
#include "PrimalForest.h"
#include "Point4D.h"
#include "Point3D.h"

namespace nurbs 
{

	double PrimalForest::jacDet(const double s, const double t,
								const uint sp) const
	{
		return cross(tangent(s, t, sp, S), tangent(s, t, sp, T)).length();
	}
	

	DoubleVecVec PrimalForest::jacob(const double s, const double t,
								  const uint sp) const
	{
		return {tangent(s, t, sp, S).asVec(), tangent(s, t, sp, T).asVec()};
	}

	Point3D PrimalForest::normal(const double s, const double t,
								 const uint sp) const
	{
        Point3D n = cross(tangent(s, t, sp, S), tangent(s, t, sp, T));
		return n / n.length();
	}
	
	Point3D PrimalForest::tangent(const double s, const double t,
								  const uint sp, const ParamDir dir) const
	{
        const ParamDir dir2 = (dir == S) ? T : S;
		Point4D temp;
        const BSplineSpace& b_space = space(sp);
		UIntVecVec indices = b_space.localBasisFuncI(s, t);
		UIntVec span = b_space.span(s, t);
		DoubleVecVec basis = b_space.tensorBasis(s,t,span);
		DoubleVecVec ders = b_space.tensorBasisDers(s,t,span);
		Point3D a, ader;
		double w = 0.0, wder = 0.0;
		for(uint i = 0; i < b_space.degree(S) + 1; ++i ) {
			for(uint j = 0; j < b_space.degree(T) + 1; ++j ) {
				const std::size_t der_i = (dir == S) ? i : j;
				const std::size_t nder_i = (dir == S) ? j : i;
				const Point4D& p = point(sp, indices[S][i], indices[T][j]);
				a += p.asUnweighted() * basis[S][i] * basis[T][j];
				ader += p.asUnweighted() * ders[dir][der_i] * basis[dir2][nder_i];
				w += p.getWeight() * basis[S][i] * basis[T][j];
				wder += p.getWeight() * ders[dir][der_i] * basis[dir2][nder_i];
			}
		}
		return nurbshelper::getNonRationalDeriv( { a, ader }, { w, wder } );
	}

	Point3D PrimalForest::eval(const double s, const double t,
							   const uint sp) const
	{
		const BSplineSpace& b_space = space(sp);
		UIntVec indices = b_space.globalBasisFuncI(s,t);
		DoubleVec basis = b_space.basis(s,t);
		assert(basis.size() == indices.size());
		Point4D p;
		for(uint i = 0; i < basis.size(); ++i)
			p += point(sp,indices[i]) * basis[i];
		return p.asCartesian();
	}
	
	void PrimalForest::writeVTKOutput(const std::string& s, const uint nsample) const
	{
		// loop over all trunks and elements within these and sample
		// geometry
	}

	void PrimalForest::loadImpl(std::istream& is)
	{
		char ch;
		if(is >> ch && ch != '{') {
			is.unget();
			is.clear(std::ios_base::failbit);
			return;
		}
		mGeomVec.clear();
		Point4D p;
		while(is >> p) 
 			mGeomVec.push_back(p);
		endOfLoop(is, '}', "Bad control point set");
	}
	
	void PrimalForest::printImpl(std::ostream& ost) const
	{
		std::cout << mGeomVec.size() << " CPts\n";
		for(const auto& p : mGeomVec)
			ost << p << "\n";
	}
    
//    std::unique_ptr<drc::Element> PrimalForest::createElement(const uint sp, const uint iel) const
//    {
//        
//    }
	
	
}

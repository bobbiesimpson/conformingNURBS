#include <memory>
#include <iostream>

#include "NURBSSurface.h"
#include "NURBSCommon.h"
#include "IElem.h"
#include "Point.h"
#include "IElemIntegrate.h"



using namespace nurbs::nurbshelper;
using namespace nurbs::elem;

namespace nurbs
{
	
	Point3D NURBSSurface::getCoord( const double s, const double t ) const
	{
		Point4D temp;
		std::vector<UIntVec> indices = getIndices(s,t);
		UIntVec span = getSpan(s,t);
		DoubleVecVec basis = BsplineBasisImpl(s,t,span);
		for(std::size_t i = 0; i < getOrder(S) + 1; ++i)
			for(std::size_t j = 0; j < getOrder(T) + 1; ++j)
				temp += mCPts[ indices[T][j] ][ indices[S][i] ]
					* basis[S][i] * basis[T][j];
		return temp.asCartesian();
	}

	double NURBSSurface::jacDet(const double s, const double t) const
	{
		return cross(tangent(s,t,S), tangent(s,t,T)).length();
	}

	Point3D NURBSSurface::normal(const double s, const double t) const
	{
		Point3D n = cross(tangent(s,t,S), tangent(s,t,T));
		return n / n.length();
	}
	
	void NURBSSurface::printData( std::ostream& ost ) const
	{
		ost << "--------- NURBS SURFACE DATA -------------\n\n";
		ost << "order ( s direction ) = " << getOrder( S ) << "\n";
		ost << "order ( t direction ) = " << getOrder( T ) << "\n";	
		ost << "knot vector ( s direction ) = {";
		std::copy( mKnotVec[ S ].begin(), mKnotVec[ S ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mKnotVec[ S ].back() << "}\n";
		ost << "knot vector ( t direction ) = {";
		std::copy( mKnotVec[ T ].begin(), mKnotVec[ T ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mKnotVec[ T ].back() << "}\n";
		ost << "unique knot vector ( s direction ) = {";
		std::copy( mUniqueKnotVec[ S ].begin(), mUniqueKnotVec[ S ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mUniqueKnotVec[ S ].back() << "}\n";
		ost << "unique knot vector ( t direction ) = {";
		std::copy( mUniqueKnotVec[ T ].begin(), mUniqueKnotVec[ T ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mUniqueKnotVec[ T ].back() << "}\n";
		ost << "Control points:\n";
		size_t count = 1;
		for( auto row : mCPts )
		{
			for( auto p : row )
				ost << "(" << count++ << ") " << p << "\n";
		}
		ost << "\n--------------------------------------\n";		
	}

	

	
	Point3D NURBSSurface::tangent(const double s, const double t, const ParamDir dir) const
	{
		const ParamDir dir2 = (dir == S) ? T : S;
		Point4D temp;
		std::vector<UIntVec> indices = getIndices(s,t);
		UIntVec span = getSpan(s,t);
		DoubleVecVec basis = BsplineBasisImpl(s,t,span);
		DoubleVecVec ders = BsplineBasisDersImpl(s,t,span);
		Point3D a, ader;
		double w = 0.0, wder = 0.0;
		for(std::size_t i = 0; i < getOrder(S) + 1; ++i ) {
			for(std::size_t j = 0; j < getOrder(T) + 1; ++j ) {
				const std::size_t der_i = (dir == S) ? i : j;
				const std::size_t nder_i = (dir == S) ? j : i;
				const Point4D& p = mCPts[ indices[T][j] ][ indices[S][i] ];
				a += p.asUnweighted() * basis[S][i] * basis[T][j];
				ader += p.asUnweighted() * ders[dir][der_i] * basis[dir2][nder_i];
				w += p.getWeight() * basis[S][i] * basis[T][j];
				wder += p.getWeight() * ders[dir][der_i] * basis[dir2][nder_i];
			}
		}
		return getNonRationalDeriv( { a, ader }, { w, wder } );
	}

	DoubleVecVec NURBSSurface::BsplineBasisImpl(const double s, const double t,
												const UIntVec& span) const
	{
		DoubleVec s_basis = getBsplineBasis( s, span[S], mKnotVec[ 0 ], getOrder( S ) );
		DoubleVec t_basis = getBsplineBasis( t, span[T], mKnotVec[ 1 ], getOrder( T ) );
		return {s_basis, t_basis};
	}

	DoubleVecVec NURBSSurface::BsplineBasisDersImpl(const double s,
													const double t,
													const UIntVec& span,
													const DerivOrder order) const
	{
		DoubleVecVec s_basis = getBsplineBasisDers(s, span[S], mKnotVec[ 0 ], getOrder(S), order);
		DoubleVecVec t_basis = getBsplineBasisDers(t, span[T], mKnotVec[ 1 ], getOrder(T), order);
		return {s_basis.at(order), t_basis.at(order)};
	}	
	
	UIntVec NURBSSurface::getSpan(const double s, const double t) const
	{
		return {getKnotSpan(s, mKnotVec[0], getOrder(S)),
				getKnotSpan(t, mKnotVec[1], getOrder(T)) };
	}

	std::vector<UIntVec> NURBSSurface::getIndices(const double s, const double t) const
	{
		return {getBasisFnIndices(s, mKnotVec[0], getOrder(S)),
				getBasisFnIndices(t, mKnotVec[1], getOrder(T)) };
	}
	
	std::ostream& operator<<( std::ostream& ost, const NURBSSurface& s )
	{
		s.printData( ost );
		return ost;
	}
	
	
}

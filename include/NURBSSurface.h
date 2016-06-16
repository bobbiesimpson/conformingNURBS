#ifndef NURBS_SURFACE_H
#define NURBS_SURFACE_H

#include <cassert>
#include <algorithm>

#include "NURBSCommon.h"
#include "base.h"
#include "Point3D.h"
#include "Point4D.h"

namespace nurbs
{
	class NURBSSurface
	{
		public:

		/// Constructor
		NURBSSurface( const std::vector< std::vector< Point4D > >& cpts,
					  const std::vector< std::vector< double > >& knotvec,
					  const uint p,
					  const uint q )
			: mCPts( cpts ),
			  mKnotVec( knotvec ),
			  mOrderS( p ),
			  mOrderT( q )
		{
			if( knotvec.size() != 2 )
				error( "Incorrect number of knot vectors specified for NURBS surface" );
			// check number control points in each row is equal
			uint prev = cpts[ 0 ].size();
			uint npts = prev;
			for( std::size_t i = 1; i < cpts.size(); ++i )
			{
				if( cpts[ i ].size() != prev )
					error( "Incorrect control point array specified." );
				prev = cpts[ i ].size();
				npts += prev;
			}
			if( npts != ( knotvec[ 0 ].size() - p - 1 ) * ( knotvec[ 1 ].size() - q - 1 ) )
				error( "Incorrect data specified for NURBS surface" );
			mNumBasisS = knotvec[ 0 ].size() - p - 1;
			mNumBasisT = knotvec[ 1 ].size() - q - 1;
			// construct unique knot vectors
			for( auto kv : mKnotVec )
			{
				auto last = std::unique( kv.begin(), kv.end() );
				mUniqueKnotVec.push_back( DoubleVec( kv.begin(), last ) );
			}
		}

		/// get number of basis functions in each parameteric direction
		uint getBasisFuncN( const ParamDir dir ) const
		{
			return ( dir == S ) ? mNumBasisS : mNumBasisT;
		}

		/// get order of basis for each parametric direction
		uint getOrder( const ParamDir dir ) const
		{
			return ( dir == S ) ? mOrderS : mOrderT;
		}

		/// get knot coordinate at specified index and direction
		double getKnotCoord( const uint i, ParamDir dir ) const
		{
			assert( i < mKnotVec[ dir ].size() );
			return mKnotVec[ dir ][ i ];
		}

		/// get unique knot coordinate
		double getUniqueKnotCoord( const uint i, ParamDir dir ) const
		{
			assert( i < mUniqueKnotVec[ dir ].size() );
			return mUniqueKnotVec[ dir ][ i ];
		}

		/// get homogenenous form of control point
		inline Point4D getHomogeneousCPt( const uint i, const uint j ) const
		{
			assert( i < getBasisFuncN( S ) && j < getBasisFuncN( T ) );
			return mCPts[ i ][ j ];
		}

		/// get cartesian form of control point
		inline Point3D getCartesianCPt( const uint i, const uint j ) const
		{
			return getHomogeneousCPt( i, j ).asCartesian();
		}

		/// get number of knot spans ( elements ) in each direction
		inline uint getUniqueKnotN( const ParamDir dir ) const { return mUniqueKnotVec[ dir ].size(); }

		/// get total number of elements
		inline uint getElN() const { return ( ( getUniqueKnotN( S ) - 1 ) * ( getUniqueKnotN( T ) - 1 ) ); }

		/// interpolate the surface and return the cartesian coordinate
		Point3D getCoord( const double s, const double t ) const;

		/// Get the jacobian determinant
		double jacDet(const double s, const double t) const;

		/// The normal to the surface
		Point3D normal(const double s, const double t) const;

		/// Get tangent at parametric coordinate
		Point3D tangent(const double s, const double t, const ParamDir dir) const;
		

		
		protected:

		private:



        /// Get Bspline basis functions
		DoubleVecVec BsplineBasisImpl(const double s, const double t,
									  const UIntVec& span) const;

		/// Get bspline basis derivatives
		DoubleVecVec BsplineBasisDersImpl(const double s,
										  const double t,
										  const UIntVec& span,
										  const DerivOrder order = D1) const;

		/// Get vector of spans for this parmaetric coordinate
		UIntVec getSpan(const double s, const double t) const;

		/// Get vector of basis function indices for this parametric coordinate
		std::vector<UIntVec> getIndices(const double s, const double t) const;
	
		typedef std::vector< std::vector< Point4D > > Point4DVecVec;

		/// Matrix of control points size ( n x m )
		Point4DVecVec mCPts;

		/// Knot vectors
		DoubleVecVec mKnotVec;

		/// Unique knot vectors ( elements )
		DoubleVecVec mUniqueKnotVec;

		/// Basis function order in S direction
		uint mOrderS;

		/// Basis funcion order in T direction
		uint mOrderT;

		/// Number of basis functions in S directions
		uint mNumBasisS;

		/// Number of basis functions in T direction
		uint mNumBasisT;

		/// function to output current state of NURBSSurface
		void printData( std::ostream& ost ) const;
		
		/// overload output operator ( no need to be a friend of NURBSSurface )
		friend std::ostream& operator<<( std::ostream& ost, const NURBSSurface& s );
	};


}

#endif

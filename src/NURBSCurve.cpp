#include "NURBSCurve.h"
#include <algorithm>
#include "base.h"
#include "NURBSCommon.h"
#include <boost/math/special_functions/binomial.hpp>

namespace nurbs
{
	
	NURBSCurve::NURBSCurve( const std::vector< Point4D >& cpts,
							const DoubleVec& knotvec,
							const uint p,
							const std::string& n )
		: mCPts( cpts ),
		  mKnotVec( knotvec ),
		  mOrder( p ),
		  mName( n )
	{
		// check for valid input
		if( !isValid() )
			error( "Invalid data entered for NURBS curve" );
		
	}

	NURBSCurve::NURBSCurve( const std::vector< Point3D >& cpts,
							const DoubleVec& knotvec,
							const uint p,
							const std::string& n )
		: mKnotVec( knotvec ),
		  mOrder( p ),
		  mName( n )
	{
		for( auto &p : cpts )
			mCPts.push_back( Point4D( p.getCoord( 0 ), p.getCoord( 1 ), p.getCoord( 2 ), 1.0 ) );
		// check for valid input
		if( !isValid() )
			error( "Invalid data entered for NURBS curve" );
	}

	Point3D NURBSCurve::getCoord( const double s ) const
	{
		UIntVec indices = getGlobalFuncI( s );
		const uint span = getKnotSpan( s );
		return getCoordImpl( s, indices, span );
	}

	Point3D NURBSCurve::getCoord( const double s,
								  const UIntVec& indices,
								  const uint span ) const
	{
		return getCoordImpl( s, indices, span );
	}

	DoubleVec NURBSCurve::getBasisFuncs( const double s,
										 const uint span ) const
	{
		return getBasisFuncsImpl( s, span );
	}

	DoubleVec NURBSCurve::getBasisFuncs( const double s ) const
	{
		return getBasisFuncsImpl( s, getKnotSpan( s ) );
	}

	DoubleVec NURBSCurve::getBasisFuncDers( const double s,
											const uint span,
											const DerivOrder deriv ) const
	{
		return getBasisFuncDersImpl( s, span, deriv );
	}
	

	DoubleVec NURBSCurve::getBasisFuncDers( const double s,
											const DerivOrder deriv ) const
	{
		return getBasisFuncDersImpl( s, getKnotSpan( s ), deriv );
	}
	
	double NURBSCurve::getJacob( const double s ) const
	{
		return getJacobImpl( s, getGlobalFuncI( s ), getKnotSpan( s ) );
	}

	double NURBSCurve::getJacob( const double s,
									const UIntVec& indices,
									const uint span ) const
	{
		return getJacobImpl( s, indices, span );
	}
	
	Point3D NURBSCurve::getNormal( const double s ) const
	{
		return getNormalImpl( s, getGlobalFuncI( s ), getKnotSpan( s ) );
	}

	Point3D NURBSCurve::getNormal( const double s,
								   const UIntVec& indices,
								   const uint span ) const
	{
		return getNormalImpl( s, indices, span );
	}

	Point3D NURBSCurve::getTangent( const double s ) const
	{
		return getTangentImpl( s, getGlobalFuncI( s ), getKnotSpan( s ) ).asNormal();
	}

	Point3D NURBSCurve::getTangent( const double s,
								   const UIntVec& indices,
								   const uint span ) const
	{
		return getTangentImpl( s, indices, span ).asNormal();
	}

	uint NURBSCurve::getKnotSpan( const double s ) const
	{
		return nurbshelper::getKnotSpan( s, mKnotVec, getOrder() );
	}
	
	UIntVec NURBSCurve::getGlobalFuncI( const double s ) const
	{
		return nurbshelper::getBasisFnIndices( s, mKnotVec, getOrder() );
	}
	
	void NURBSCurve::outputInterpolation( std::ostream& ost, const double inc ) const
	{
		double s =  mKnotVec.front();
		while( s < mKnotVec.back() )
		{
			Point3D p = this->getCoord( s );
			ost << p.getCoord( 0 )<< "\t" << p.getCoord( 1 ) << "\t" << p.getCoord( 2 ) << "\n";
			s += inc;
		}
	}

	void NURBSCurve::hrefine( const int level )
	{
		if( level == 0 )
			return;
		mUniqueKnotVec.clear(); // clear the 'element' knot vector. Getter will recalculate.
		const int n = getBasisFuncN() - 1;
		const int p = getOrder();
		const int m = getKnotVecSize() - 1;
		DoubleVec& u = mKnotVec; // alias for current knot vector
		std::vector< double > unique; // determine unique knot coords
		std::unique_copy( u.begin(), u.end(), std::back_inserter( unique ) );
		int knot_n = 0;
		for( int i = 0; i < level; ++i )
			knot_n += std::pow( 2.0, i );  // determine number of knots to insert per interval
		std::vector< double > x; // vector of new knot coords
		for( std::size_t i = 0; i < unique.size() - 1; ++i )
			for( int j = 0; j < knot_n; ++j )
				x.push_back( ( unique[ i + 1 ] - unique[ i ] ) / ( knot_n + 1 ) * ( j + 1 )
							 + unique[ i ] ); // insert new knots
		const int r = x.size() - 1;

		// Begin Piegl and Tiller algorithm p.164
		std::vector< double > ubar( u.size() + x.size() );
		std::vector< Point4D >& pvec = mCPts;  // handy alias
		std::vector< Point4D > qvec( n + r + 2 );
		int a = nurbshelper::getKnotSpan( x[ 0 ], u, p );
		int b = nurbshelper::getKnotSpan( x[ r ], u, p );
		++b;
		for( int j = 0; j <= a - p; j++ ) qvec[ j ] = pvec[ j ];
		for( int j = b - 1; j <= n; j++ ) qvec[ j + r + 1 ] = pvec[ j ];
		for( int j = 0; j <= a; j++ ) ubar[ j ] = u[ j ];
		for( int j = b + p; j <= m; j++ ) ubar[ j + r + 1 ] = u[ j ];
		int i = b + p - 1;
		int k = b + p + r;
		for( int j = r; j >= 0; j-- )
		{
			while( x[ j ] <= u[ i ] && i > a )
			{
				qvec[ k - p - 1 ] = pvec[ i - p - 1 ];
				ubar[ k ] = u[ i ];
				--k;
				--i;
			}
			qvec[ k - p - 1 ] = qvec[ k - p ];
			for( int l = 1; l <= p; l++ )
			{
				int ind = k - p + l;
				double alfa = ubar[ k + l ] - x[ j ];
				if( essentiallyEqual( std::fabs( alfa ), 0.0, TOL ) )
					qvec[ ind - 1 ] = qvec[ ind ];
				else
				{
					alfa /= ubar[ k + l ] - u[ i - p + l ];
					Point4D& cp1 = qvec[ ind - 1 ];
					const Point4D& cp2 = qvec[ ind ];
					for( auto i = 0; i < 3; ++i )
						cp1.setCoord( i, cp1.getCoord( i ) * alfa
									  + cp2.getCoord( i ) * ( 1.0 - alfa ) );
					cp1.setWeight( cp1.getWeight() * alfa + cp2.getWeight() * ( 1.0 - alfa ) );
					
				}
			}
			ubar[ k ] = x[ j ];
			--k;
		}
		u = ubar; // set new knot vector
		pvec = qvec; // set new control points
	}

	void NURBSCurve::prefine( const uint level )
	{
		if( level == 0 )
			return;
	}

	void NURBSCurve::normalise( const double range, const double origin )
	{
		const double r = mKnotVec.back() - mKnotVec.front();
		const double shift = origin - mKnotVec.front();
		for( auto &k : mKnotVec )
		{
			k *= range / r;
			k += shift;
		}
	}

	double NURBSCurve::getMinCartesianCoord( const CartesianComponent dir ) const
	{
		double min = mCPts.front().getCartesianCoord( dir );
		for( auto it = mCPts.begin() + 1; it != mCPts.end(); ++it )
		{
			const double t = it->getCartesianCoord( dir );
			min = ( t < min ) ? t : min;
		}
		return min;
	}

	double NURBSCurve::getMaxCartesianCoord( const CartesianComponent dir ) const
	{
		double max = mCPts.front().getCartesianCoord( dir );
		for( auto it = mCPts.begin() + 1; it != mCPts.end(); ++it )
		{
			const double t = it->getCartesianCoord( dir );
			max = ( t > max ) ? t : max;
		}
		return max;
	}
	
	Point3D NURBSCurve::getCoordImpl( const double s,
									  const UIntVec& indices,
									  const uint span ) const
	{
		DoubleVec localbasis = nurbshelper::getBsplineBasis( s, span, mKnotVec, getOrder() );
		Point4D p;
		for( uint i = 0; i < getOrder() + 1; ++i )
			p += mCPts[ indices[ i ] ] * localbasis[ i ];
		return p.asCartesian();
	}

	DoubleVec NURBSCurve::getBasisFuncsImpl( const double s,
											 const uint span ) const
	{
		return nurbshelper::getBsplineBasis( s, span, mKnotVec, getOrder() );
	}

	DoubleVec NURBSCurve::getBasisFuncDersImpl( const double s,
												const uint span,
												const DerivOrder deriv ) const
	{
		return nurbshelper::getBsplineBasisDers( s, span, getKnotVec(), getOrder(), deriv ).at( deriv );
	}

	double NURBSCurve::getJacobImpl( const double s,
									 const UIntVec& indices,
									 const uint span ) const
	{
		return getTangentImpl( s, indices, span ).length();
	}

	Point3D NURBSCurve::getNormalImpl( const double s,
									   const UIntVec& indices,
									   const uint span ) const
	{
		Point3D t = getTangentImpl( s, indices, span );
		return Point3D{ t[ 1 ], -t[ 0 ], t[ 2 ] }.asNormal();
	}	

	Point3D NURBSCurve::getTangentImpl( const double s,
										const UIntVec& indices,
										const uint span ) const
	{
		DoubleVec basis = getBasisFuncs( s, span );
		DoubleVec der = getBasisFuncDers( s, span );
		Point3D a, ader;
		double w = 0.0, wder = 0.0;
		for( uint i = 0; i < getOrder() + 1; ++i )
		{
			const Point4D& p = mCPts[ indices[ i ] ];
			a += p.asUnweighted() * basis[ i ];
			ader += p.asUnweighted() * der[ i ];
			w += p.getWeight() * basis[ i ];
			wder += p.getWeight() * der[ i ];
		}
		return nurbshelper::getNonRationalDeriv( { a, ader }, { w, wder } );
	}
	
	DoubleVec NURBSCurve::getUniqueKnotVec() const
	{
		if( mUniqueKnotVec.empty() ) // if empty, we need work out unique coords and set
			std::unique_copy( mKnotVec.begin(),
							  mKnotVec.end(),
							  std::back_inserter( mUniqueKnotVec ) );
		return mUniqueKnotVec;
	}
	
	void NURBSCurve::printData( std::ostream& ost ) const
	{
		ost << "--------- NURBS CURVE DATA -------------\n\n";
		ost << "order = " << getOrder() << "\n";
		ost << "knot vector = {";
		std::copy( mKnotVec.begin(), mKnotVec.end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mKnotVec.back() << "}\n";
		ost << "Control points:\n";
		size_t count = 1;
		for( auto& p : mCPts )
			ost << "(" << count++ << ") " << p << "\n";
		ost << "\n--------------------------------------\n";
 	}

	std::ostream& operator<<( std::ostream& ost, const NURBSCurve& c )
	{
		c.printData( ost );
		return ost;
	}

	std::istream& operator>>( std::istream& ist, NURBSCurve& c )
	{
		char ch;
		if( ist >> ch && ch != '{' ) // check that we start with '{'
		{
			ist.unget();
			ist.clear( std::ios_base::failbit );
			return ist;
		}
		std::string n;
		ist >> n;
		if( !ist )
			error( "Cannot read nurbs curve name" );
		NURBSCurve::CPSet cpset;
		ist >> cpset;
		if( !ist )
			error( "Bad control point set specified" );
		char marker;
		uint p;
		ist >> marker >> p;
		if( !ist || marker != 'p' )
			error( "Bad NURBS curve order" );
		std::vector< double > kvec;
		if( ist >> marker && marker != '(' )
			error( "Bad knot vector reading" );
		while( true )
		{
			double r;
			if( !( ist >> r >> ch ) )
				break;
			if( ch == ',' || ch == ')' )
			{
				kvec.push_back( r );
				if( ch == ')' )
				{
					ist.unget();
					ist.clear( std::ios::failbit );
					break;
				}
			}
			else
				error( "Bad knot vector reading" );
		}
		endOfLoop( ist, ')', "Bad end of knot vector" );
		if( kvec.size() != cpset.pts.size() + p + 1 )
			error( "Invalid NURBS curve specified - knot vector size != n + p + 1" );
		if( ist >> ch && ch != '}' )
			error( "Bad end of NURBS curve" );
		c = NURBSCurve( cpset.pts, kvec, p, n );
		return ist;
	}

	std::istream& operator>>( std::istream& ist, NURBSCurve::CPSet& set )
	{
		char ch;
		if( ist >> ch && ch != '(' )
		{
			std::cout << "BAd char = " << ch << "\n";
			ist.unget();
			ist.clear( std::ios_base::failbit );
			return ist;
		}
		Point4D pt;
		std::vector< Point4D > pt_vec;
		while( ist >> pt )
			pt_vec.push_back( pt );
		endOfLoop( ist, ')', "Bad end of control point set" );
		set.pts = pt_vec;
		return ist;
	}

	std::ostream& operator<<(std::ostream& ost, const Element& e)
	{
		ost << "------- Element " << e.elementI() << " --------\n"
			<< "Knot span index: " << e.knotSpanI() << "\n"
			<< "Physical coordinate, midpoint: " << e.eval(0.0) << "\n"
			<< "Jacobian, midpoint: " << e.jacob(0.0);
		return ost;
	}

}

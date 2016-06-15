#include "Point4D.h"

namespace nurbs
{
	Point3D Point4D::asCartesian() const
	{
		return Point3D( getCartesianCoord( 0 ),
						getCartesianCoord( 1 ),
						getCartesianCoord( 2 ) );
	}

	Point3D Point4D::asUnweighted() const
	{
		return Point3D( getCoord( 0 ),
						getCoord( 1 ),
						getCoord( 2 ) );
	}

	void Point4D::printImpl( std::ostream& ost ) const
	{
		ost << "4D homogeneous coord: (";
	}

	std::istream& operator>>( std::istream& ist, Point4D& p )
	{
		char ch;
		if( ist >> ch && ch != '(' ) // this is not a point, return
		{
			ist.unget();
			ist.clear( std::ios_base::failbit );
			return ist;
		}
		double x, y, z, w;
		char sep;
		ist >> x >> sep >> y >> sep >> z >> sep >> w >> ch;
		if( !ist || ch != ')' )
			error( "Bad control point reading" );
		p.setCoord( 0, x * w );
		p.setCoord( 1, y * w );
		p.setCoord( 2, z * w );
		p.setCoord( 3, w );
		return ist;
	}

	Point4D operator+( const Point4D& p, const Point4D& p2 )
	{
		return Point4D( p.getCoord( 0 ) + p2.getCoord( 0 ),
						p.getCoord( 1 ) + p2.getCoord( 1 ),
						p.getCoord( 2 ) + p2.getCoord( 2 ),
						p.getCoord( 3 ) + p2.getCoord( 3 ) );
	}

	Point4D asPoint4D( const Point3D& p )
	{
		return Point4D( p.getCoord( 0 ), p.getCoord( 1 ), p.getCoord( 2 ) );
	}
	
}

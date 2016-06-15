#include "Point3D.h"

namespace nurbs
{
	void Point3D::printImpl( std::ostream& ost ) const
	{
		ost << "3D Point (";
	}
	
	double dist( const Point3D& p1, const Point3D& p2 )
	{
		double d = 0.0;
		for( uint i = 0; i < p1.getSize(); ++i )
			d += ( p1.getCoord( i ) - p2.getCoord( i ) ) *
				( p1.getCoord( i ) - p2.getCoord( i ) );
		return sqrt( d );
	}

	/// addition helper function
	Point3D operator+( const Point3D& p, const Point3D& p2 )
	{
		return Point3D( p.getCoord( 0 ) + p2.getCoord( 0 ),
						p.getCoord( 1 ) + p2.getCoord( 1 ),
						p.getCoord( 2 ) + p2.getCoord( 2 ) );
	}

	/// overload substraction operator
	Point3D operator-( const Point3D& p, const Point3D& p2 )
	{
		return Point3D( p.getCoord( 0 ) - p2.getCoord( 0 ),
						p.getCoord( 1 ) - p2.getCoord( 1 ),
						p.getCoord( 2 ) - p2.getCoord( 2 ) );
	}
	
	Point3D operator*( const double v, const Point3D& p )
	{
		Point3D temp( p );
		return temp *= v;
	}

	Point3D cross(const Point3D& p1, const Point3D& p2)
	{
		return Point3D(p1[1] * p2[2] - p1[2] * p2[1],
					   p1[2] * p2[0] - p1[0] * p2[2],
					   p1[0] * p2[1] - p1[1] * p2[0]);
	}
    
    bool operator<(const Point3D& p1, const Point3D& p2)
    {
        if(p1[0] != p2[0])
            return p1[0] < p2[0];
        if(p1[1] != p2[1])
            return p1[1] < p2[1];
        return p1[2] < p2[2];
    }
	
}

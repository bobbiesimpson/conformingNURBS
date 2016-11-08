#include "Point.h"
#include <cmath>

namespace nurbs
{
	std::ostream& operator<<( std::ostream& ost, const Point& p )
	{
		p.print( ost );
		return ost;
	}
	
    bool approximatelyEqualPoints(const Point& p1, const Point& p2, const double tol)
    {
        assert(p1.getSize() == p2.getSize());
        
        for(std::size_t i = 0; i < p1.getSize(); ++i)
            if(!essentiallyEqual(p1.getCoord(i), p2.getCoord(i), tol))
                return false;
        return true;
    }
    
	double dot( const Point& p1, const Point& p2 )
	{
		assert( p1.getSize() == p2.getSize() );
		double v = 0.0;
		for( uint i = 0; i < p1.getSize(); ++i )
			v += p1[ i ] * p2[ i ];
		return v;
	}
    
    
	
}

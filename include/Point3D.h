#ifndef POINT_3D_H
#define POINT_3D_H

#include "Point.h"
#include <cmath>

namespace nurbs
{
	/// The class which represent points in 3D cartesian space
	class Point3D : public Point
	{
		public:

		/// The one and only constructor
		explicit Point3D( double x = 0.0,
						  double y = 0.0,
						  double z = 0.0 )
			: Point( { x, y, z } )
		{}

		/// Constructor from vector
		template< typename T >
		explicit Point3D( T&& v ) 
			: Point( std::forward< T >( v ) )
		{}
		
		/// number of dimensions
		inline uint getSize() const { return 3; }

		/// Overload *= operator
		Point3D& operator*=( const double v )
		{
			Point::operator*=( v );
			return *this;
		}

		/// Overload /= operator
		Point3D& operator/=( const double v )
		{
			Point::operator/=( v );
			return *this;
		}

		/// Overload compound += operator
		Point3D& operator+=( const Point3D& p )
		{
			Point::operator+=( p );
			return *this;
		}

		/// Overload -= operator
		Point3D& operator-=( const Point3D& p )
		{
			Point::operator-=( p );
			return *this;
		}
		
		/// Overload multiplication operator
		Point3D operator*( const double v ) const
		{
			Point3D temp = *this;
			return temp *= v;
			/* for( std::size_t i = 0; i < getSize(); ++i ) */
			/* 	temp.setCoord( i, this->getCoord( i ) * v ); */
			/* return temp; */
		}

		/// overload division operator
		Point3D operator/( const double v ) const
		{
			Point3D temp( *this );
			return temp /= v;
		}
		
		/// get length
		double length() const
		{
			double result = 0.0;
			for( uint i = 0; i < getSize(); ++i )
				result += getCoord( i ) * getCoord( i );
			return sqrt( result );
		}

		/// Return unit length point
		Point3D asNormal() const
		{
			Point3D p( *this );
			return p.normalise();
		}
		
		/// normalise the point
		Point3D& normalise()
		{
			( *this ) /= length();
			return *this;
		}

		/// cast to complex double
		operator ComplexDouble() const
		{
			return std::complex< double >( getCoord( X ), getCoord( Y ) );
		}
        
        /// Cast to vector
        std::vector<double> asVec() const { return { getCoord(X), getCoord(Y), getCoord(Z)}; }
		
		private:

		/// define print implementation functions
		virtual void printImpl( std::ostream& ost ) const;

	};

	/// addition helper function
	Point3D operator+( const Point3D& p, const Point3D& p2 );

	/// subtraction helper function
	Point3D operator-( const Point3D& p, const Point3D& p2 );

	/// helper function to allow multiplication like ' 2.0 * p'
	Point3D operator*( const double v, const Point3D& p );
	
	/// non member non friend function for distance
	double dist( const Point3D& p1, const Point3D& p2 );

	/// Cross product
	Point3D cross(const Point3D& p1, const Point3D& p2);
    
    /// Less-than operator
    bool operator<(const Point3D& p1, const Point3D& p2);

}

#endif

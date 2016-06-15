#ifndef POINT4D_H
#define POINT4D_H

#include "base.h"
#include "Point.h"
#include "Point3D.h"
#include <iostream>
#include <istream>

namespace nurbs
{
	/// Point4D represent a 4-dimensional coordinate (x*w,y*w,z*w,w)
	/// where w represents the appropriate weight
	/// This allows B-spline algorithms to be used to generate NURBS curves.
	/// We simply apply the appropriate
	/// algorithms to the 4-D 'homogeneous' coordinate system ( x*w,y*w,z*w,w )
	/// and divide through by w to retrieve the cartesian coordinate (x,y,z)
	
	class Point4D : public Point
	{
		public:

		/// typedef of matrix of Point4D
		typedef std::vector< std::vector< Point4D > > Point4DVecVec;
		
		/// Constructor
		explicit Point4D( const double x = 0.0, const double y = 0.0, const double z = 0.0,
						  const double w = 0.0 )
			: Point{ x * w, y * w, z * w, w }
		{}
		
		/// Constructor with vector
		explicit Point4D( const std::vector< double >& v ) 
			: Point( v )
		{}
		
		/// getter for homogenous coordinate ( for const return type)
		inline double getCartesianCoord( const uint i ) const
		{
			return getCoord( i ) / getWeight();
		}

		/// return the weight
		inline double getWeight() const { return getCoord( 3 ); }

		/// set the weight
		inline void setWeight( const double w ) { setCoord( 3, w ); }

		/// returns size
		inline uint getSize() const { return 4;}
			
		/// getter for homogeneous coordinate
		Point3D asCartesian() const;

		/// get un-weighted coordinates
		Point3D asUnweighted() const;

		Point4D operator*( const double v ) const
		{
			Point4D temp = *this;
			for( std::size_t i = 0; i < getSize(); ++i )
				temp.setCoord( i, this->getCoord( i ) * v );
			return temp;
		}
		
		private:

		/// define print implementation functions
		virtual void printImpl( std::ostream& ost ) const;
			
		/// overload input operator
		friend std::istream& operator>>( std::istream& ist, Point4D& p );
		
	};

	/// addition helper function
	Point4D operator+( const Point4D& p, const Point4D& p2 );

	/// cast Point3D to Point4D
	Point4D asPoint4D( const Point3D& p );
	
}

#endif

#ifndef POINT_H
#define POINT_H

#include "base.h"
#include <iostream>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace nurbs
{
	class Point
	{
		public:

		/// subscript operator for non-const objects
		inline double& operator[]( const uint i ) { return mCoords[ i ]; }

		/// subscript operator for const objects
		inline const double& operator[]( const uint i ) const { return mCoords[ i ]; }

		/// getter
		inline double getCoord( uint i ) const
		{
			return mCoords[ i ];
		}

		/// get raw data
		inline double* data() { return mCoords.data(); }

		/// get const raw data
		inline const double* data() const { return mCoords.data(); }
        
        /// Convert to vector
        inline DoubleVec asVec() const {return DoubleVec(data(), data() + getSize());}
		
		/// setter
		inline void setCoord( const uint i, const double c ) { mCoords[ i ] = c; }

		/// get vector size
		virtual uint getSize() const = 0;
		
		/// overload compound *= operator
		Point& operator*=( const double v )
		{
			for( auto &c : mCoords )
				c *= v;
			return *this;
		}

		/// overload /= compound operator
		Point& operator/=( const double v )
		{
			for( auto& c : mCoords )
				c /= v;
			return *this;
		}

		/// overload compound += operator
		Point& operator+=( const Point& p )
		{
			assert( p.mCoords.size() == mCoords.size() );
			for( std::size_t i = 0; i < mCoords.size(); ++i )
				mCoords[ i ] += p.mCoords[ i ];
			return *this;
		}

		/// Overload -= operator
		Point& operator-=( const Point& p )
		{
			for( std::size_t i = 0; i < mCoords.size(); ++i )
				mCoords[ i ] -= p.mCoords[ i ];
			return *this;
		}
        
        /// Overload equality operator. Calls essentiallyEqual function
        /// defined in base.h. Users should take in its use.
        bool operator==(const Point& p) const
        {
            for(std::size_t i = 0; i < mCoords.size(); ++i)
                if(!essentiallyEqual(mCoords[i], p.mCoords[i], TOL))
                    return false;
            return true;
        }

		/// print the data for this point
		void print( std::ostream& ost ) const
		{
			printImpl( ost );
			std::copy( mCoords.begin(), mCoords.end() - 1, std::ostream_iterator< double >( ost, "," ) );
			ost << mCoords.back();
			ost << ")";
		}
		
		protected:

		/// constructor with an initializer list ( protected, since abstract )
		Point( std::initializer_list< double > l ) 
			: mCoords( l )
		{}

        Point(const std::vector<double>& v)
        : mCoords(v) {}
        
//		/// Constructor with vector
//		template< typename T >
//		Point( T&& v ) 
//			: mCoords( std::forward< T >( v ) )
//		{}
			
		private:

		/// virtual implementation function for printing data
		virtual void printImpl( std::ostream& ost ) const = 0;
		
		/// The coordinate data for the n-dimensional coordinate
		DoubleVec mCoords;

		/// make the output operator overload for Point a friend (access private members)
		friend std::ostream& operator<<( std::ostream& ost, const Point& p );
		
	};
    
    /// Helper function for determining equality of two points with specified tolerance
    bool approximatelyEqualPoints(const Point& p1, const Point& p2, const double tol);

	/// Dot product
	double dot( const Point& p1, const Point& p2 );
}

#endif

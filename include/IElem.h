#ifndef ELEMENT_ITERATOR_H
#define ELEMENT_ITERATOR_H

#include <iostream>

#include "base.h"
#include "NURBSSurface.h"

namespace nurbs
{
	namespace elem
	{
		/// simple struct to hold information about intervals for an element
		struct Element
		{
			Element( double ls, double us,
					 double lt, double ut ) 
			: lower_s( ls ),
				upper_s( us ),
				lower_t ( lt ),
				upper_t( ut ) 
				{}
		
			double lower_s;
			double upper_s;
			double lower_t;
			double upper_t;
		};

		/// An iterator class that iterates over all elements (knot spans)
		/// of a nurbs surface
		class IElem
		{
			public:

			/// constructor
			explicit IElem(const NURBSSurface& ns)
				: mSurface( ns ),
				mCurrentIndex( 0 )
				{}

			/// are we finished looping over all elements?
			bool isDone() const
			{
				return mCurrentIndex == ( mSurface.getUniqueKnotN( S ) - 1 ) *
					( mSurface.getUniqueKnotN( T ) - 1 );
			
			}

			/// return the limits of the current element to the client
			Element getCurrentEl() const
			{
				uint i = mCurrentIndex / ( mSurface.getUniqueKnotN( S ) - 1 );
				uint j = mCurrentIndex % ( mSurface.getUniqueKnotN( S ) - 1 );
				return Element( mSurface.getUniqueKnotCoord( j, S ),
								mSurface.getUniqueKnotCoord( j + 1, S ),
								mSurface.getUniqueKnotCoord( i, T ),
								mSurface.getUniqueKnotCoord( i + 1, T ) );
			}

			/// postfix increment operator
			IElem& operator++()
			{
				++mCurrentIndex;
				return *this;
			}
			
			protected:

			private:
		
			/// Reference to the nurbs surface instance		
			const NURBSSurface& mSurface;

			/// the global index counter
			uint mCurrentIndex;
		};

		std::ostream& operator<<( std::ostream& ost, const Element& e )
		{
			ost << "Element defined by intervals s[" << e.lower_s << "," << e.upper_s
				<< "]  t[" << e.lower_t << "," << e.upper_t << "]";
			return ost;
		}

		/// Simple struct to hold a parametric coordinate
		struct ParamPt
		{
			/// The constructor
			explicit ParamPt( const double sin , const double tin ) 
				: s( sin ), t( tin )
				{}

			/// S coordinate
			double s;

			/// T coordinate
			double t;
		
		};

		/// overload output operator for ParamPt
		std::ostream& operator<<( std::ostream& ost, const ParamPt& p )
		{
			ost << p.s << "\t" << p.t;
			return ost;
		}
		
		/// An iterator class that iterates over a set of evenly spaced points over
		/// an element (knot span)
		class ISamplePt
		{
			public:

			explicit ISamplePt( const Element& e, const uint npts ) 
				: mCurrentEl( e ),
				mKnotLinePts( npts ),
				mCurrentIndex( 0 )
				{
					mSamplePtN = npts * npts;
					mSIncrement = ( mCurrentEl.upper_s - mCurrentEl.lower_s ) / static_cast< double >( npts - 1 );
					mTIncrement = ( mCurrentEl.upper_t - mCurrentEl.lower_t ) / static_cast< double >( npts - 1 );
				}

            /// Default constuctors
            ISamplePt(const uint npts,
                      const double ls = -1.0,
                      const double us = 1.0,
                      const double lt = -1.0,
                      const double ut = 1.0)
            :
            ISamplePt(Element(ls,us,lt,ut), npts) {}
            
			bool isDone() const
			{
				return mCurrentIndex == mSamplePtN;
			}

			/// prefix increment operator
			ISamplePt& operator++()
			{
				++mCurrentIndex;
				return *this;
			}

			/// postfix increment operator
			ISamplePt operator++( int i )
			{
				ISamplePt temp( *this );
				++(*this);
				return temp;
			}

			ParamPt getCurrentPt() const
			{
				const uint i = mCurrentIndex / mKnotLinePts;
				const uint j = mCurrentIndex % mKnotLinePts;
				return ParamPt( mCurrentEl.lower_s + mSIncrement * j,
								mCurrentEl.lower_t + mTIncrement * i );
			}
		
			protected:

			private:

			/// The current element we are iterating over
			const Element mCurrentEl;

			/// The number of points along each knot line
			const uint mKnotLinePts;
		
			/// the total number of sample points for the element
			uint mSamplePtN;

			/// the current point index
			uint mCurrentIndex;

			/// The increment between points in the s direction
			double mSIncrement;

			/// The increment between points in the t direction
			double mTIncrement;
		
		
		};
	}
	
	
	
}

#endif

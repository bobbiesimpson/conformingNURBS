#ifndef IGaussQuadrature1D_H
#define IGaussQuadrature1D_H

#include "base.h"

namespace nurbs
{
	class IGaussQuadrature1D
	{
		public:
	
		/// Construct with integer for num gpts in both directions
		explicit IGaussQuadrature1D( uint ngps = DEFAULT_NGPS )
			: mPointN( ngps ),
			  mCurrentGP( 0 ) 
		{
			initGaussPts( ngps );
		}

		/// Check if we've reached the end
		inline bool isDone() const { return isDoneImpl(); }

		/// Get the gauss point component
		inline double get() const { return getImpl(); }
		
		/// Get the gauss weight
		double getWeight() const { return getWeightImpl(); }

		/// Prefix increment
		IGaussQuadrature1D& operator++() { incrementImpl(); return *this; };

		/// Postfix increment
		IGaussQuadrature1D operator++( int )
		{
			IGaussQuadrature1D temp( *this );
			++( *this );
			return temp;
		}

		protected:

		/// restart the iterator
		inline void reset() { mCurrentGP = 0; }
		
		/// get current gauss point implementation (can be used by derived classes)
		virtual double getImpl() const;

		/// get current weight implementation (can be used by derived classes)
		virtual double getWeightImpl() const;
		
		/// implementation of termination function. Can be overridden by derived classses
		virtual bool isDoneImpl() const
		{
			return mCurrentGP == mPointN;
		}

		/// increment implementation function
		virtual void incrementImpl()
		{
			++mCurrentGP;
		}

		private:
		
		/// 1d gauss points
		void initGaussPts( uint rule );

		/// Total number of guass points
		uint mPointN;

		/// Current gauss points
		uint mCurrentGP;

		/// gauss pts 
		DoubleVec mGaussPts;
		
		/// gauss weights
		DoubleVec mGaussWts;
		
	};
	
}

#endif

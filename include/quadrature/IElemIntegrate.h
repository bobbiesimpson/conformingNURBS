#ifndef NURBS_IElemIntegrate_H
#define NURBS_IElemIntegrate_H

#include "base.h"
#include "IBaseIntegrate.h"

namespace nurbs
{
	namespace elem 
	{
        class IElemIntegrate : public IBaseIntegrate
		{
			public:

			/// Construct with integer for num gpts in both directions
			IElemIntegrate( uint ngps = DEFAULT_NGPS ): mCurrentGP( 0 ) 
			{
				mOrders.push_back( ngps );
				mOrders.push_back( ngps );
				mPointN = ngps * ngps;
				init(mOrders);
			}

			/// Construct with a specified number of points in each direction
			IElemIntegrate(const UIntVec& orders, const uint offset = 0)
            : mCurrentGP(0)
			{
                uint n = 1;
                for(const auto& v : orders) {
                    mOrders.push_back(v + offset);
                    n *= mOrders.back();
                }
				mPointN = n;
				init(mOrders);
			}
            
		
			/// Default copy constructor will work

			/// Default assignment operator will work

			/// Check if we've reached the end
			bool isDone() const { return mCurrentGP == mPointN; }

			/// Get the gauss point
			GPt2D get() const;

			/// Get the gauss point component
			double get( uint component ) const;
		
			/// Get the gauss weight
			double getWeight() const;
            
            /// Restart the iterator to point to the first gauss point
            void restart()
            { mCurrentGP = 0; }
            
            /// Return the total number of gauss points
            uint pointN() const
            { return mPointN; }

            /// Current gauss point index
            uint currentIndex() const
            { return mCurrentGP; }
            
			/// Increment iterator
            IElemIntegrate& operator++() {  incrementImpl(); return *this; };

			/// Print status of iterator (gauss pts, weights, current index)
			void printStatus() const;
	
			private:

			/// initate gauss point vectors
			void init( const UIntVec& orders );

            /// Increment implementation
            void incrementImpl() { ++mCurrentGP; }
		
			/// Number of gauss points in each parametric direction
			UIntVec mOrders;

			/// Total number of guass points
			uint mPointN;

			/// Current gauss points
			uint mCurrentGP;

			/// gauss pts s direction
			DoubleVec mGaussPtsS;

			/// gauss pts t direction
			DoubleVec mGaussPtsT;
		
			/// gauss pts s direction
			DoubleVec mGaussWtsS;

			/// gauss pts t direction
			DoubleVec mGaussWtsT;
		
		};
        
        /// Non-member, non-friend functions
        void get1dGaussRule( uint rule, DoubleVec& pts, DoubleVec& wts );
        
        /// Wrapper class for 1D integration
        class IElemIntegrate1D {
            
        public:
            
            /// The only constructor
            IElemIntegrate1D(const uint ngps = DEFAULT_NGPS)
            :
            mPointN(ngps),
            mCurrentI(0)
            { get1dGaussRule(ngps, mGaussPts, mGaussWts); }
            
            /// Check if we've reached the end
            bool isDone() const { return currentI() == pointN(); }
            
            /// Get the gauss point component
            double get() const { return gaussPt(currentI()); }
            
            /// Get the gauss weight
            double getWeight() const { return gaussWt(currentI()); }
            
            /// Restart the iterator to point to the first gauss point
            void restart() { currentI() = 0; }
            
            /// quadrature point no. getter
            uint pointN() const { return mPointN; }
            
            /// Current index getter
            uint currentI() const { return mCurrentI; }
            
            /// Increment iterator
            IElemIntegrate1D& operator++() { ++currentI(); return *this; };
            
            /// Print status of iterator (gauss pts, weights, current index)
            void printStatus() const
            {
                std::cout << "index: " << currentI() << " / " << pointN() << "\n";
            }
            
        private:
            
            /// Gauss point getter
            double gaussPt(const uint i) const
            {
                assert(i < pointN());
                return mGaussPts[i];
            }
            
            /// Gauss weight getter
            double gaussWt(const uint i) const
            {
                assert(i < pointN());
                return mGaussWts[i];
            }
            
            /// non-const current index getter
            uint& currentI() { return mCurrentI; }

            /// Quadrature points
            DoubleVec mGaussPts;
            
            /// Quadrature weights
            DoubleVec mGaussWts;
            
            /// Number of quadrature points
            uint mPointN;
            
            /// Current quadrature index
            uint mCurrentI;
        };
	}
	
}

#endif

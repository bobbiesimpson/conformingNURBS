#include <iostream>

#include "base.h"
#include "IElemIntegrate.h"

namespace nurbs
{
	namespace elem 
	{	
		GPt2D IElemIntegrate::get() const
		{
			// find t index
			uint t_index = mCurrentGP / mGaussPtsS.size();
			// find s index (remainder)
			uint s_index = mCurrentGP % mGaussPtsS.size();
			return GPt2D(mGaussPtsS[s_index], mGaussPtsT[t_index]);
		}

		double IElemIntegrate::get( uint component ) const
		{
			if( component > 2 )
				error( "Out of bounds access of integrator gauss point: aborting" );
			if( component == 0 )
				return get().get( 0 );
			else
				return get().get( 1 );
		}

		double IElemIntegrate::getWeight() const
		{
			// find t index
			uint t_index = mCurrentGP / mGaussPtsS.size();
			// find s index (remainder)
			uint s_index = mCurrentGP % mGaussPtsS.size();
			return mGaussWtsS[ s_index ] * mGaussWtsT[ t_index ];
		}
	
		void IElemIntegrate::init( const UIntVec& orders )
		{
			// initialise gauss pts and weights in S and T directions
			get1dGaussRule( orders[ 0 ], mGaussPtsS, mGaussWtsS );
			get1dGaussRule( orders[ 1 ], mGaussPtsT, mGaussWtsT );
		}

		void IElemIntegrate::printStatus() const
		{
			std::cout << "Current gauss point index = " << mCurrentGP << " / " << mPointN << "\n\n";
		
			// print gauss points in S direction
			std::cout << "Gauss points (S direction)" << std::endl;
			for( DoubleVec::const_iterator it = mGaussPtsS.begin(); it != mGaussPtsS.end(); ++it )
				std::cout << *it << std::endl;
			std::cout << "\n";
			// print gauss points T direction
			std::cout << "Gauss points (T direction)" << std::endl;
			for( DoubleVec::const_iterator it = mGaussPtsT.begin(); it != mGaussPtsT.end(); ++it )
				std::cout << *it << std::endl;
			std::cout << "\n";
		
			// print gauss weights in S direction
			std::cout << "Gauss weights (S direction)" << std::endl;
			for( DoubleVec::const_iterator it = mGaussWtsS.begin(); it != mGaussWtsS.end(); ++it )
				std::cout << *it << std::endl;
			std::cout << "\n";
		
			// print gauss weights T direction
			std::cout << "Gauss weights (T direction)" << std::endl;
			for( DoubleVec::const_iterator it = mGaussWtsT.begin(); it != mGaussWtsT.end(); ++it )
				std::cout << *it << std::endl;		
		}
        
        void get1dGaussRule( uint rule, DoubleVec& pts, DoubleVec& wts )
        {
            if( rule == 0 ) {}
            else if( rule == 1 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 1, 0.0 );
                wts.resize( 1, 0.0 );
                
                pts[ 0 ] = 0;
                
                wts[ 0 ] = 2;
            }
            else if( rule == 2 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 2, 0.0 );
                wts.resize( 2, 0.0 );
                
                pts[ 0 ] = -1 / sqrt( 3 );
                pts[ 1 ] = 1 / sqrt( 3 );
                
                wts[ 0 ] = wts[ 1 ] = 1;
            }
            else if( rule == 3 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 3, 0.0 );
                wts.resize( 3, 0.0 );
                
                wts[0] = wts[2] = 5.0 / 9.0;
                wts[1] = 8.0 / 9.0;
                
                pts[0] = -sqrt( 3.0 / 5.0 );
                pts[1] = 0;
                pts[2] = sqrt( 3.0 / 5.0 );
            }
            else if( rule == 4 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 4, 0.0 );
                wts.resize( 4, 0.0 );
                
                wts[0] = wts[3] = ( 18.0 - sqrt( 30.0 ) ) / 36.0;
                wts[1] = wts[2] = ( 18.0 + sqrt( 30.0 ) ) / 36.0;
                
                double inner = sqrt( ( 3.0 - 2.0 * sqrt( 6.0 / 5.0 ) ) / 7.0 );
                double outer = sqrt( ( 3.0 + 2.0 * sqrt( 6.0 / 5.0 ) ) / 7.0 );
                
                pts[0] = -outer;
                pts[1] = -inner;
                pts[2] = inner;
                pts[3] = outer;
            }
            else if( rule == 5 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 5, 0.0 );
                wts.resize( 5, 0.0 );
                
                wts[0] = wts[4] = (322.0 - 13.0*sqrt(70.0))/900.0;
                wts[1] = wts[3] = (322.0 + 13.0*sqrt(70.0))/900.0;
                wts[2] = 128.0/225.0;
                
                double inner = 1.0/3.0 * sqrt( 5.0 - 2.0*sqrt(10.0/7.0));
                double outer = 1.0/3.0 * sqrt( 5.0 + 2.0*sqrt(10.0/7.0));
                
                pts[0] = -outer;
                pts[1] = -inner;
                pts[2] = 0.0;
                pts[3] = inner;
                pts[4] = outer;
            }
            else if( rule == 6 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 6, 0.0 );
                wts.resize( 6, 0.0 );
                
                wts[0] = wts[5] = 0.171324492379170;
                wts[1] = wts[4] = 0.360761573048139;
                wts[2] = wts[3] = 0.467913934572691;
                
                double inner = 0.238619186083197;
                double mid = 0.661209386466265;
                double outer = 0.932469514203152;
                
                pts[0] = -outer;
                pts[1] = -mid;
                pts[2] = -inner;
                pts[3] = inner;
                pts[4] = mid;
                pts[5] = outer;
            }
            else if( rule == 7 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 7, 0.0 );
                wts.resize( 7, 0.0 );
                
                wts[0] = wts[6] = 0.129484966168870;
                wts[1] = wts[5] = 0.279705391489277;
                wts[2] = wts[4] = 0.381830050505119;
                wts[3] = 0.417959183673469;
                
                double inner = 0.405845151377397;
                double mid = 0.741531185599394;
                double outer = 0.949107912342759;
                
                pts[0] = -outer;
                pts[1] = -mid;
                pts[2] = -inner;
                pts[3] = 0.0;
                pts[4] = inner;
                pts[5] = mid;
                pts[6] = outer;
            }
            else if( rule == 8 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 8, 0.0 );
                wts.resize( 8, 0.0 );
                
                wts[7] = wts[0] = 0.101228536290376;
                wts[6] = wts[1] = 0.222381034453374;
                wts[5] = wts[2] = 0.313706645877887;
                wts[4] = wts[3] = 0.362683783378362;
                
                pts[7] = 0.960289856497536; pts[0] = -pts[7];
                pts[6] = 0.796666477413627; pts[1] = -pts[6];
                pts[5] = 0.525532409916329; pts[2] = -pts[5];
                pts[4] = 0.183434642495650; pts[3] = -pts[4];
            }
            else if( rule == 9 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 9, 0.0 );
                wts.resize( 9, 0.0 );
                
                wts[8] = wts[0] = 0.081274388361574;
                wts[7] = wts[1] = 0.180648160694857;
                wts[6] = wts[2] = 0.260610696402935;
                wts[5] = wts[3] = 0.312347077040003;
                wts[4] = 0.330239355001260;
                
                pts[8] = 0.968160239507626; pts[0] = -pts[8];
                pts[7] = 0.836031107326636; pts[1] = -pts[7];
                pts[6] = 0.613371432700590; pts[2] = -pts[6];
                pts[5] = 0.324253423403809; pts[3] = -pts[5];
                pts[4] = 0.0;
            }
            else if( rule == 10 )
            {
                pts.clear();
                wts.clear();
                pts.resize( 10, 0.0 );
                wts.resize( 10, 0.0 );
                
                wts[9] = wts[0] = 0.066671344308688;
                wts[8] = wts[1] = 0.149451349150581;
                wts[7] = wts[2] = 0.219086362515982;
                wts[6] = wts[3] = 0.269266719309996;
                wts[5] = wts[4] = 0.295524224714753;
                
                pts[9] = 0.973906528517172; pts[0] = -pts[9];
                pts[8] = 0.865063366688985; pts[1] = -pts[8];
                pts[7] = 0.679409568299024; pts[2] = -pts[7];
                pts[6] = 0.433395394129247; pts[3] = -pts[6];
                pts[5] = 0.148874338981631; pts[4] = -pts[5];
            }
            else
            {
                // numerically determine gauss points and weights
                // slow but nice to have for higher order gauss rules.
                pts.clear();
                wts.clear();
                pts.resize( rule, 0.0 );
                wts.resize( rule, 0.0 );
                
                const double x1 = -1.0;
                const double x2 = 1.0;
                const double EPS = 1.0e-14;
                int m = ( rule + 1 ) / 2;
                double xm = 0.5 * ( x2 + x1 );
                double xl = 0.5 * ( x2 - x1 );
                for( int i = 0; i < m; ++i )
                {
                    double z = cos( PI * ( i + 0.75 ) / ( rule + 0.5 ) );
                    double z1, pp, p3;
                    do
                    {
                        double p1( 1.0 );
                        double p2( 0.0 );
                        for( uint j = 0; j < rule; ++j )
                        {
                            p3 = p2;
                            p2 = p1;
                            p1 = ( ( 2.0 * j + 1.0 ) * z * p2 - j * p3 ) / ( j + 1 );
                        }
                        pp = rule * ( z * p1 - p2 ) / ( z * z - 1.0 );
                        z1 = z;
                        z = z1 - p1 / pp;
                    } while ( fabs( z - z1 ) > EPS );
                    pts[ i ] = xm - xl * z;
                    pts[ rule - 1 - i ] = xm + xl * z;
                    wts[ i ] = 2.0 * xl / ( ( 1.0 - z * z ) * pp * pp );
                    wts[ rule - 1 - i ] = wts[ i ];
                }
            }
        }
	}
}

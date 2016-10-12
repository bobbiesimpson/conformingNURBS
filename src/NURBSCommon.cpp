#include "NURBSCommon.h"
#include "base.h"
#include "Point.h"
#include "Point4D.h"
#include "NURBSCache.h"

#include <vector>
#include <algorithm>

#include <boost/math/special_functions/binomial.hpp>

// #include "Epetra_SerialComm.h"
// #include "Epetra_Map.h"
// #include "Epetra_Vector.h"
// #include "Epetra_CrsMatrix.h"
// #include "Epetra_LinearProblem.h"
// #include "AztecOO.h"

// #include "Amesos.h"
// #include "Amesos_BaseSolver.h"

#include <cmath>

namespace nurbs
{
	namespace nurbshelper
	{
		std::vector< uint > getBasisFnIndices( const double s,
											   const DoubleVec& knotvec,
											   const uint p )
		{
			std::vector< uint > v;
			const uint i = getKnotSpan( s, knotvec, p );
			uint j = i - p;
			while( j <= i )
				v.push_back( j++ );
			return v;
		}

		uint getKnotSpan( const double s,
						  const DoubleVec& knotvec,
						  const uint p )
		{
            // First search cache
            NURBSCache& cache = NURBSCache::Instance();
            auto cp = cache.span(s, knotvec, p);
            if(cp.first)
                return cp.second;
            
			// taken from Piegl and Tiller (binary search)
			const uint n = knotvec.size() - p - 1;
			if( essentiallyEqual( s, knotvec[ n ], TOL ) )
				return n - 1;
			uint low = p;
			uint high = n + 1;
			uint mid = ( low + high ) / 2;
			while( s < knotvec[ mid ] || s >= knotvec[ mid + 1 ] )
			{
				if( s < knotvec[ mid ] )
					high = mid;
				else low = mid;
				mid = ( low + high ) / 2;
			}
            cache.cacheSpan(s, knotvec, p, mid);
			return mid;
		}

		DoubleVec getBsplineBasis( const double s,
								   const uint span,
								   const DoubleVec& knotvec,
								   const uint p )
		{
            
            // First search cache
            NURBSCache& cache = NURBSCache::Instance();
            auto cp = cache.basis(s, span, knotvec, p);
            if(cp.first)
                return cp.second;
            
			DoubleVec basis( p + 1, 0.0 );
			DoubleVec left( p + 1, 0.0 );
			DoubleVec right( p + 1, 0.0 );
			basis[ 0 ] = 1.0;
			for( uint j = 1; j <= p; ++j )
			{
				left[ j ] = s - knotvec[ span + 1 - j ];
				right[ j ] = knotvec[ span + j ] - s;
				double saved = 0.0;
				for( uint r = 0; r < j; ++r )
				{
					double temp =  basis[ r ] / ( right[ r + 1 ] + left[ j - r ] );
					basis[ r ] = saved + right[ r + 1 ] * temp;
					saved = left[ j - r ] * temp;
				}
				basis[ j ] = saved;
			}
            cache.cacheBasis(s, span, knotvec, p, basis);
			return basis;
		}


		double getBsplineBasisWithIndex( const double s,
										 const uint gindex,
										 const DoubleVec& knotvec,
										 const uint p )
		{
			const uint m = knotvec.size() - 1;
			DoubleVec N( p + 1 );
			if( ( gindex == 0 && essentiallyEqual( s, knotvec[0], TOL ) ) ||
				( gindex == ( m - p - 1 ) && essentiallyEqual( s, knotvec[ m ], TOL ) ) )
				return 1.0;
			if( s < knotvec[ gindex ] || s >= knotvec[ gindex + p + 1 ] )
				return 0.0;
			
			for( uint j = 0; j <= p; j++)
			{
				if( s >= knotvec[ gindex + j ] && s < knotvec[ gindex + j + 1 ] ) N[ j ] = 1.0;
				else N[ j ] = 0.0;
			}

			double saved, Uleft, Uright, temp;
			for( uint k = 1; k <= p; k++ )
			{
				if( essentiallyEqual( N[ 0 ], 0.0, TOL ) ) saved = 0.0;
				else saved = ( ( s - knotvec[ gindex ] ) * N[ 0 ] )
						 /  ( knotvec[ gindex + k ] - knotvec[ gindex ] );
				for( uint j = 0; j < ( p - k + 1 ); j++)
				{
					Uleft = knotvec[ gindex + j + 1 ];
					Uright = knotvec[ gindex + j + k + 1 ];
					if( essentiallyEqual( N[ j + 1 ], 0.0, TOL ) )
					{
						N[ j ] = saved; saved = 0.0;
					}
					else 
					{
						temp = N[ j + 1 ] / ( Uright - Uleft );
						N[ j ] = saved + ( Uright - s ) * temp;
						saved = ( s - Uleft ) * temp;
					}
				}
			}
			double Nip = N[ 0 ];
			return Nip;
		}		

        DoubleVecVec getBsplineBasisDers(const double s,
                                         const uint span,
                                         const DoubleVec& knotvec,
                                         const int p,
                                         const DerivOrder order )
		{
			if( order > 2 )
				error( "B-spline derivatives of order > 2 not implemented yet." );
            
            // First search cache
            NURBSCache& cache = NURBSCache::Instance();
            auto cp = cache.basisDer(s, span, knotvec, p, order);
            if(cp.first)
                return cp.second;
            
			double temp;
			DoubleVec left( p + 1 );
			DoubleVec right( p + 1 );
			DoubleVecVec ndu(  p + 1, DoubleVec( p + 1 ) );
			DoubleVecVec a( p + 1, DoubleVec( p + 1 ) );
			DoubleVecVec ders( p + 1, DoubleVec( p + 1 ) );

			ndu[ 0 ][ 0 ] = 1.0;
			for( int j = 1; j <= p; j++ )
			{
				left[ j ] = s - knotvec[ span + 1 - j ];
				right[ j ] = knotvec[ span + j ] - s;
				double saved = 0.0;
				for( int r = 0; r < j; r++ )
				{
					ndu[ j ][ r ] = right[ r + 1 ] + left[ j - r ];
					temp = ndu[ r ][ j - 1 ] / ndu[ j ][ r ];
      
					ndu[ r ][ j ] = saved + right[ r +1 ] * temp;
					saved = left[ j - r ] * temp;
				}
				ndu[ j ][ j ] = saved;
			}
			for( int j = 0; j <= p; j++ )
				ders[ 0 ][ j ] = ndu[ j ][ p ];
			if( order == 0 )
				return ders;
			for( int r = 0; r <= p; r++ )
			{
				int s1 = 0, s2 = 1;   
				a[ 0 ][ 0 ] = 1.0;

				for( int k = 1; k <= order; k++ )
				{
					double d = 0.0;
					int rk = r - k, pk= p - k;
					if( r >= k )
					{
						a[ s2 ][ 0 ] = a[ s1 ][ 0 ] / ndu[ pk + 1 ][ rk ];
						d = a[ s2 ][ 0 ] * ndu[ rk ][ pk ];
					}
					int j1 = rk >= -1 ? 1 : -rk;
					int j2 = ( r -1 <= pk ) ? k - 1 : p - r;
					for( int j = j1; j <= j2; j++ )
					{
						a[ s2 ][ j ]=( a[ s1 ][ j ] - a[ s1 ][ j - 1 ] )
							/ ndu[ pk + 1 ][ rk + j ];
						d += a[ s2 ][ j ] * ndu[ rk + j ][ pk ];
					}
					if( r <= pk )
					{
						a[ s2 ][ k ] = -a[ s1 ][ k - 1 ] / ndu[ pk + 1 ][ r ];
						d += a[ s2 ][ k ] * ndu[ r ][ pk ];
					}
					ders[ k ][ r ] = d;
					int j = s1; s1 = s2; s2 = j;  
				}
			}
			int r = p;
			for( int k = 1; k <= order; k++ )
			{
				for( int j = 0; j <= p; j++ ) 
					ders[ k ][ j ] *= r;
				r *= ( p - k );
			}
            cache.cacheBasisDer(s, span, knotvec, p, order, ders);
			return ders;
		}

		Point3D getNonRationalDeriv( const std::vector< Point3D >& aders,
									 const std::vector< double >& wders,
									 const DerivOrder order )
		{
			std::vector< Point3D > cvec( order + 1 ); // calculated non-homogeneous coords
			for( int k = 0; k <= order; ++k )
			{
				Point3D v = aders[ k ];
				for( int i = 1; i <= k; ++i )
					v -= boost::math::binomial_coefficient< double >( k, i )
						* wders[ i ] * cvec[ k - i ];
				cvec[ k ] = v / wders[ 0 ];
			}
			return cvec[ order ];
		}
		
		
		double divide( const double n, const double d )
		{
			if( fabs( d ) < TOL )  // if dividing by zero, return 0
				return 0.0;
			else
				return n / d;
		}
        
        DoubleVec uniformKnotInsertion(const DoubleVec& knotvec,
                                       const uint refine)
        {
            // assume knot vec is sorted in ascending order
            DoubleVec unique = knotvec; // first make a copy of knot vector
            auto l = std::unique(unique.begin(), unique.end());
            unique.erase(l, unique.end());
            int knot_n = 0;
            for(int i = 0; i < refine; ++i)
                knot_n += std::pow( 2.0, i);  // number of knots to insert per interval
            std::vector<double> x; // vector of new knot coords
            for(std::size_t i = 0; i < unique.size() - 1; ++i)
                for( int j = 0; j < knot_n; ++j )
                    x.push_back( ( unique[ i + 1 ] - unique[ i ] ) / ( knot_n + 1 ) * ( j + 1 )
                                + unique[ i ] ); // insert new knots
            DoubleVec refined_kv;
            auto i_new = x.begin(); // iterator to beginning of new knots
            for(const auto& k : knotvec) {
                //std::cout << *i_new << "\t" << k << "\n";
                while(k > *i_new && i_new != x.end()) {
                    refined_kv.push_back(*i_new);
                    ++i_new;
                }
                refined_kv.push_back(k);
            }
            return refined_kv;
        }
        
        DoubleVec gradedKnotInsertion(const DoubleVec& knotvec,const uint refine, const double coeff)
        {
//            DoubleVec knotvec = {0,1};
//            DoubleVec unique = {0,1};
//            int knot_n = nelem - 1;
//            double base = 0.0;
//            for (int ielem = 0; ielem < nelem; ++ielem) {
//                base += pow(coeff, ielem);
//            }
//            double knot_v = 0.0;
//            double first_element = 1.0 / base;
            DoubleVec unique = knotvec; // first make a copy of knot vector
            auto l = std::unique(unique.begin(), unique.end());
            unique.erase(l, unique.end());
            int knot_n = 0;
            for(int i = 0; i < refine; ++i)
                knot_n += std::pow( 2.0, i);  // number of knots to insert per interval

            std::vector<double> insertknots;
//            for (int iknot = 0; iknot < knot_n; ++iknot) {
//                knot_v += pow(coeff,iknot) * first_element;
//                insertknots[iknot] = knot_v;
//            }
            for(std::size_t i = 0; i < unique.size() - 1; ++i){
                if (i == 0) {
                    // graded_refine
                    int nelem = knot_n + 1;
                    double base = 0.0;
                    for (int ielem = 0; ielem < nelem; ++ielem) {
                        base += pow(coeff, ielem);
                    }
                    double initial_element = (unique[i+1]-unique[i])/base;
                    double knot_v =0.0;
                    for (int j = 0; j < knot_n; ++j) {
                        knot_v += pow(coeff,j) * initial_element;
                        insertknots.push_back(unique[i+1]-knot_v);
                    
                    }
                }
                else if (i == unique.size()-2){
                    // graded_refine
                    int nelem = knot_n + 1;
                    double base = 0.0;
                    for (int ielem = 0; ielem < nelem; ++ielem) {
                        base += pow(coeff, ielem);
                    }
                    double initial_element = (unique[i+1]-unique[i])/base;
                    double knot_v =0.0;
                    for (int j = 0; j < knot_n; ++j) {
                        knot_v += pow(coeff,j) * initial_element;
                        insertknots.push_back( unique[i] + knot_v);
                    }
                }
                else{
                    // uniform_refine
                    for( int j = 0; j < knot_n; ++j ){
                        insertknots.push_back( ( unique[ i + 1 ] - unique[ i ] ) / ( knot_n + 1 ) * ( j + 1 )
                                              + unique[ i ] ); // insert new knots
                    }
                }
            }
            DoubleVec kv;
            auto i_new = insertknots.begin();
            std::sort(insertknots.begin(),insertknots.end());
            for (const auto& k : knotvec) {
                while(k > *i_new && i_new != insertknots.end()) {
                    kv.push_back(*i_new);
                    ++i_new;
                }
                kv.push_back(k);
            }
            return kv;
        }

        
        DoubleVec bernsteinPolynomial(const double xi, const uint p)
        {
            // First search cache
            NURBSCache& cache = NURBSCache::Instance();
            auto cp = cache.bernsteinBasis(xi, p);
            if(cp.first)
                return cp.second;
            
            const double x = 0.5 * (xi + 1.0);
            
            DoubleVec rvec;
            switch(p) {
                case 0:
                    rvec = {1};
                    break;
                case 1:
                    rvec = {1.0 - x, x};
                    break;
                case 2:
                    rvec =
                    {
                        std::pow(1.0 - x, 2.0),
                        2.0 * x * (1.0 - x),
                        x * x
                    };
                    break;
                case 3:
                    rvec =
                    {
                        std::pow(1.0 - x, 3.0),
                        3.0 * x * std::pow((1.0 - x), 2),
                        3 * x * x * (1.0 - x),
                        std::pow(x, 3)
                    };
                    break;
                case 4:
                    rvec =
                    {
                        std::pow(1.0 - x, 4.0),
                        4.0 * x * std::pow((1.0 - x), 3),
                        6.0 * x * x * std::pow((1.0 - x), 2),
                        4.0 * std::pow(x,3) * (1.0 - x),
                        std::pow(x, 4)
                    };
                    break;
                default:
                    for(uint i = 0; i < p + 1; ++i)
                        rvec.push_back(boost::math::binomial_coefficient<double>(p,i) * std::pow(x, i) * std::pow(1.0 - x, p - i));
                    break;
            }
            cache.cacheBernsteinBasis(xi, p, rvec);
            return rvec;
        }
        
        DoubleVec bernsteinPolynomialDeriv(const double xi, const uint p)
        {
            NURBSCache& cache = NURBSCache::Instance();
            auto cp = cache.bernsteinBasisDeriv(xi, p);
            if(cp.first)
                return cp.second;
            
            auto basis = bernsteinPolynomial(xi, p-1);
            DoubleVec deriv;
            const double dx_dxi = 0.5; // jacobian determinant from [0,1] to [-1,1]
            for(uint k = 0; k < p+1; ++k) {
                double t1 = 0.0;
                if(k > 0)
                    t1 += p * basis[k-1];
                double t2 = 0.0;
                if(k < p)
                    t2 += p * basis[k];
                deriv.push_back(dx_dxi * (t1 - t2));
            }
            cache.cacheBernsteinBasisDeriv(xi, p, deriv);
            return deriv;
        }

	// 	/// generate a NURBS curve by interpolating a set of points
	// 	NURBSCurve interpolatePts( const std::vector< Point3D >& pts,
	// 							   const uint p,
	// 							   InterpType interptype )
	// 	{
	// 		// check order
	// 		if( pts.size() < p + 1 )
	// 			error( "Cannot interpolate curve: insufficient number of points specified." );
	// 		const uint n = pts.size();
	// 		DoubleVec ubar( n, 0.0 );
	// 		double d = 0.0;
	// 		switch( interptype )
	// 		{
	// 			case equal:
	// 				for( uint k = 0; k < n; ++k )
	// 					ubar[ k ] = static_cast< double >( k ) / ( n - 1 );
	// 				break;
	// 			case chord:
	// 				for( std::size_t i = 1; i < pts.size(); ++i )
	// 					d += dist( pts[ i ], pts[ i - 1 ] );
	// 				for( uint k = 1; k < n; ++k )
	// 					ubar[ k ] = ubar[ k - 1 ] + dist( pts[ k ], pts[ k - 1 ] ) / d;
	// 				break;
	// 			case centripetal:
	// 				for( std::size_t i = 1; i < pts.size(); ++i )
	// 					d += sqrt( dist( pts[ i ], pts[ i - 1 ] ) );
	// 				for( uint k = 1; k < n; ++k )
	// 					ubar[ k ] = ubar[ k - 1 ] + sqrt( dist( pts[ k ], pts[ k - 1 ] ) )/ d;
	// 				break;
	// 			default:
	// 				error( "Could not determine interpolation type" );
	// 		}
	// 		DoubleVec knotvec( n + p + 1, 0.0 );
	// 		for( uint j = p + 1; j < knotvec.size(); ++j )
	// 		{
	// 			if( j >= ( knotvec.size() - p - 1 ) )
	// 				knotvec[ j ] = 1.0;
	// 			else
	// 			{
	// 				for( uint i = j - p; i < j; ++i )
	// 					knotvec[ j ] += ubar[ i ];
	// 				knotvec[ j ] /= p;
	// 			}
	// 		}
	// 		//printVector( knotvec, std::cout );
	// 		const uint dim = 3;
	// 		Epetra_SerialComm comm;
	// 		Epetra_Map eq_map( -1, n * dim, 0, comm );
	// 		Epetra_Vector f( eq_map );
	// 		Epetra_Vector a( eq_map );
	// 		Epetra_CrsMatrix N( Copy, eq_map, p + 1 );

	// 		for( uint a = 0; a < ubar.size(); ++a )
	// 		{
	// 			const double u = ubar[ a ];
	// 			std::cout << "u = " << u << "\n";
	// 			for( uint i = 0; i < dim; ++i )
	// 			{
	// 				const int row = a * dim + i;
	// 				UIntVec glbindices = getBasisFnIndices( u, knotvec, p );
	// 				const uint span = getKnotSpan( u, knotvec, p );
	// 				DoubleVec basis = getBsplineBasis( u, span, knotvec, p );
	// 				printVector( glbindices, std::cout ); std::cout << "\n";
	// 				printVector( basis, std::cout );
	// 				std::cout << "\n";
	// 				for( uint b = 0; b < glbindices.size(); ++b )
	// 				{
	// 					const int glbcolI = glbindices[ b ];
	// 					for( uint j = 0; j < dim; ++j )
	// 					{
	// 						const int colindex = glbcolI * dim + j;
	// 						double val = ( i == j ) ? basis[ b ] : 0.0;
	// 						N.InsertGlobalValues( row, 1, &val, &colindex );
	// 					}
	// 				}
	// 				f.SumIntoGlobalValues( 1, &pts[ a ][ i ], &row );
	// 			}
	// 		}
	// 		N.FillComplete();
	// 		const int max_iter = 1000;
	// 		const double tol = 1.0e-15;
	// 		AztecOO solver( &N, &a, &f );
	// 		solver.SetAztecOption( AZ_precond, AZ_Jacobi );
	// 		if( 0 != solver.Iterate( max_iter, tol ) )
	// 			error( "Cannot solve system through GMRES iteration: aborting" );
	// 		solver.PrintLinearSystem( "linear_system.dat" );
	// 		std::vector< Point3D > nurbs_pts( ubar.size() );
	// 		for( std::size_t i = 0; i < ubar.size(); ++i )
	// 			nurbs_pts[ i ] = Point3D( a[ i * dim ],
	// 									  a[ i * dim + 1 ],
	// 									  a[ i * dim + 2 ] );
	// 		return NURBSCurve( nurbs_pts, knotvec, p );
	// 	}	
	}
}

#ifndef NURBS_COMMON_H
#define NURBS_COMMON_H

#include "base.h"
#include "Point3D.h"

namespace nurbs
{
	/// the namespace which holds the nurbs helper functions
	namespace nurbshelper
	{
		/// The method of knot spacing for interpolation (see p. 364-365 of Piegl and Tiller
		enum InterpType
		{
			equal = 0,
			chord,
			centripetal
		};
		
		/// return a vector of the non-zero basis function indices at parameter s
		std::vector< uint > getBasisFnIndices( const double s,
											   const DoubleVec& knotvec,
											   const uint p );

		/// determine the knot span ( range ) in which a given parameter coordinate lies
		uint getKnotSpan( const double s,
						  const DoubleVec& knotvec,
						  const uint p );

		/// get B-spline basis function
		DoubleVec getBsplineBasis( const double s,
								   const uint span,
								   const DoubleVec& knotvec,
								   const uint p );

		/// Get individual basis function 'i' at parametric coordinate 's'.
		double getBsplineBasisWithIndex( const double s,
										 const uint gindex,
										 const DoubleVec& knotvec,
										 const uint p );

		/// get the set of non-zero B-spline derivatives
		DoubleVecVec getBsplineBasisDers( const double s,
										  const uint span,
										  const DoubleVec& knotvec,
										  const int p,
										  const DerivOrder order = D1 );

		/// Calculate non-rational form of derivative
		/// hderiv: homogeneous coords of point and derivatives up to 'order'
		/// returns required non-homogeneous derivative
		Point3D getNonRationalDeriv( const std::vector< Point3D >& aders,
									 const std::vector< double >& wders,
									 const DerivOrder order = D1 );
		

		/// To ensure we don't get an error when dividing by zero
		double divide( const double n, const double d );
			
        /// Create a new knot vector by applying knot insertion
        DoubleVec uniformKnotInsertion(const DoubleVec& knotvec,
                                       const uint refine = 1);
        
        /// Return the set of Bernstein polynomails evaluate at xi \in [-1,1]
        /// for degree p
        DoubleVec bernsteinPolynomial(const double xi, const uint p);
        
		/* /// Interpolate a set of points and return the associate nurbs curve */
		/* NURBSCurve interpolatePts( const std::vector< Point3D >& pts, */
		/* 						   const uint p, */
		/* 						   InterpType interptype = chord ); */
		
	}
}

#endif

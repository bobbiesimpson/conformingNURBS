#ifndef NURBS_CURVE_H
#define NURBS_CURVE_H

#include "base.h"
#include "Point3D.h"
#include "Point4D.h"
#include <cassert>
#include <istream>
#include <ostream>

#include <boost/icl/continuous_interval.hpp>


namespace nurbs
{
	/// A NURBS curve generated through a set of control points, basis order and knot vector
	class NURBSCurve
	{
		public:
		
		/// Default constructor
		explicit NURBSCurve() = default;
		
		/// Constructor with 4-d coordinates
		explicit NURBSCurve( const std::vector< Point4D >& cpts,
							 const DoubleVec& knotvec,
							 const uint p,
							 const std::string& n = "" );

		/// Consructor with 3-d coordinates
		explicit NURBSCurve( const std::vector< Point3D >& cpts,
							 const DoubleVec& knotvec,
							 const uint p,
							 const std::string& n = "" );

		/// get number of basis functions / no. of control points
		inline uint getBasisFuncN() const { return mCPts.size(); }
		
		/// getter for basis order
		inline uint getOrder() const { return mOrder; }

		/// getter for knot vector coordinate
		inline double getKnotCoord( uint i ) const { return mKnotVec[ i ]; }
		
		/// getter for knot vector size
		inline uint getKnotVecSize() const { return mKnotVec.size(); }

		/// getter for knot vector
		inline const DoubleVec& getKnotVec() const { return mKnotVec; }

		/// getter for control points coord (in 4-d space)
		inline Point4D getHomogeneousCPt( uint i ) const
		{
			assert( i < mCPts.size() - 1 );
			return mCPts[ i ];
		}

		/// same as previous, except in cartesian coordinates
		inline Point3D getCartesianCPt( uint i ) const
		{
			assert( i < mCPts.size() );
			return mCPts[ i ].asCartesian();
		}

		/// interpolate the NURBS curve for parametric coordinate s
		Point3D getCoord( const double s ) const;

		/// Overload for more efficient call from Element class (stores indices and span)
		Point3D getCoord( const double s,
						  const UIntVec& indices,
						  const uint span ) const;

		/// Get the basis function given a parametric coord, local index and span
		DoubleVec getBasisFuncs( const double s,
								 const uint span ) const;

		/// Get basis functions given parametric coordinate
		DoubleVec getBasisFuncs( const double s ) const;
		
		/// Get the basis function derivatives
		DoubleVec getBasisFuncDers( const double s,
									const uint span,
									const DerivOrder deriv = D1 ) const;

		/// Get basis function derivatives
		DoubleVec getBasisFuncDers( const double s,
									const DerivOrder deriv = D1 ) const;
		
		
		/// Get number of unique knot coordinates
		uint getUniqueKnotN() const { return getUniqueKnotVec().size(); }

		/// Get the unique knot value at the specified index
		double getUniqueKnot( const uint i ) const { return getUniqueKnotVec().at( i ); }

		/// Get jacobian
		double getJacob( const double s ) const;

		/// Get jacobian with knot basis func indices and span
		double getJacob( const double s,
						 const UIntVec& indices,
						 const uint span ) const;

		/// Get the outward pointing normal
		Point3D getNormal( const double s ) const;

		/// Get the normal with known basis func indices and span
		Point3D getNormal( const double s,
						   const UIntVec& indices,
						   const uint span ) const;

		/// Get the tangent ( this is not normalised )
		Point3D getTangent( const double s ) const;

		/// Get the (non-normalised) tangent with known basis func indices and span
		Point3D getTangent( const double s,
						   const UIntVec& indices,
						   const uint span ) const;
		
		/// get current name
		inline std::string getName() const { return mName; }

		/// set current name
		inline void setName( const std::string& n ) { mName = n; }

		/// get the knot span for the specified parametric coordinate
		uint getKnotSpan( const double s ) const;

		/// the global basis function indices that are non-zero at this param. coord.
		UIntVec getGlobalFuncI( const double s ) const;
		
		/// interpolate the curve at increments specified by inc and output
		void outputInterpolation( std::ostream& ost, const double inc = INC ) const;

		/// apply uniform h-refinement (knot insertion) to the curve
		void hrefine( const int level = 1 );

		/// apply p-refinement (order elevation) to the curve
		void prefine( const uint level = 1 );

		/// normalise knot vector to specified range and (optionally) specify origin
		void normalise( const double range = 1.0, const double origin = 0.0 );

		/// Determin minimum cartesian coordinate in specified direction
		double getMinCartesianCoord( const CartesianComponent dir ) const;

		/// Determin maximum cartesian coordinate in specified direction
		double getMaxCartesianCoord( const CartesianComponent dir ) const;

		/// get upper bound of parametric range
		double getUB() const
		{
			return mKnotVec.back();
		}

		/// get lower bound of parametric range
		double getLB() const
		{
			return mKnotVec.front();
		}

		private:
		
		///  A struct used to help input of a NURBS curve
		struct CPSet
		{
			/// the vector of control points
			std::vector< Point4D > pts;
		};

		/// Implementation function for generating physical coordinate
		Point3D getCoordImpl( const double s,
							  const UIntVec& indices,
							  const uint span ) const;

		/// Implementation function for generating basis functions
		DoubleVec getBasisFuncsImpl( const double s,
									 const uint span ) const;

		/// Implementation function for basis function derivatives
		DoubleVec getBasisFuncDersImpl( const double s,
										const uint span,
										const DerivOrder deriv ) const;

		/// jacobian implementation
		double getJacobImpl( const double s,
							 const UIntVec& indices,
							 const uint span ) const;

		/// normal evaluation implementation
		Point3D getNormalImpl( const double s,
							   const UIntVec& indices,
							   const uint span ) const;
		
		/// Get the curve tangent
		Point3D getTangentImpl( const double s,
								const UIntVec& indices,
								const uint span ) const;
		
		/// is the current state of the object a valid NURBS curve?
		inline bool isValid() const
		{
			if( mKnotVec.size() ==  mCPts.size() + mOrder + 1 )
				return true;
			return false;
		}

		/// Get the unique knot vector. Called by public member functions
		DoubleVec getUniqueKnotVec() const;

		/// print all data for this curve
		void printData( std::ostream& ost ) const;
		
		/// Vector of control points
		std::vector< Point4D > mCPts;

		/// Knot Vector
		DoubleVec mKnotVec;

		/// Unique knot coordinates
		mutable DoubleVec mUniqueKnotVec;

		/// Order of basis 
		uint mOrder;

		/// name of curve ( assigned in input )
		std::string mName;
		
		/// overload output operator
		friend std::ostream& operator<<( std::ostream& ost, const NURBSCurve& c );

		/// overload input operator
		friend std::istream& operator>>( std::istream& ist, NURBSCurve& c );
		
		/// overload the input operator	
		friend std::istream& operator>>( std::istream& ist, CPSet& set );

		/// Make iterator class a friend
		friend class INURBSCurve;
		
	};

	/// Iterate over NURBS curve to sweep out paramterisation in physical coordinates
	class INURBSCurve
	{
		public:

		/// Default construct
		explicit INURBSCurve( const NURBSCurve& c, const uint n )
			: mCurve( c ),
			  mSampleN( n ),	
			  mSampleInc( ( mCurve.mKnotVec.back() - mCurve.mKnotVec.front() )
						/ static_cast< double >( n - 1 ) ),
			  mCurrentParamC( mCurve.mKnotVec.front() ),
			  mCurrentSampleI( 0 )
		{}

		/// get the current curve
		inline const NURBSCurve& getCurve() const { return mCurve; }

		/// get the current coordinate in physical space
		Point3D getCurrentCoord() const { return getCurve().getCoord( mCurrentParamC ); }
		
		/// get the current parametric coordinate
		inline double getCurrentParam() const { return mCurrentParamC; }

		/// get number of sample points
		inline uint getSamplePtN() const { return mSampleN; }

		/// get current sample index
		inline uint getCurrentI() const { return mCurrentSampleI; }

		/// prefix increment
		INURBSCurve& operator++()
		{
			++mCurrentSampleI;
			mCurrentParamC += mSampleInc;
			return *this;
		}

		/// post fix increment
		INURBSCurve operator++( int )
		{
			INURBSCurve copy( *this );
			++( *this );
			return copy;
		}

		/// are we done iterating?
		bool isDone() const
		{
			return mCurrentSampleI == mSampleN;
		}
		
		private:
		
		/// The currrent curve
		const NURBSCurve& mCurve;

		/// Number of sample points
		const uint mSampleN;
		
		/// number of sample points
		const double mSampleInc;
		
		/// The current parametric coordinate
		double mCurrentParamC;

		/// Current sample index
		uint mCurrentSampleI;
		
	};

	/// Object which provides all the necessary functionality to build
	/// element stiffness/mass/.. matrices for FE analysis
	class Element {

		public:

		/// Constructor
		Element(const NURBSCurve& c, const uint iel)
			: mCurve(c),
			mElementI(iel),
        mElInterval(boost::icl::construct<boost::icl::continuous_interval<double> >(c.getUniqueKnot(iel),
                                                                                    c.getUniqueKnot(iel + 1),
                                                                                    boost::icl::interval_bounds::closed())),
			mKnotSpanI(c.getKnotSpan(mElInterval.lower()))
		{}
			
		/// element index getter
		uint elementI() const {return mElementI;}

		/// knot span index
		uint knotSpanI() const { return mKnotSpanI; }

		/// curve getter
		const NURBSCurve& curve() const{ return mCurve; }

		/// element interval lower bound
		double lowerbound() const { return mElInterval.lower();}

		/// element upper bound
		double upperbound() const { return mElInterval.upper(); }
		
		/// Indices of global non-zero basis functions
		std::vector<uint> globalBasisI() const { return curve().getGlobalFuncI(lowerbound()); }

		/// interpolate the physical coordinate given a local coordinate
		Point3D eval(const double xi) const { return curve().getCoord(paramCoord(xi)); }

		/// evaluate basis function given local basis index and local coordinate
		double evalBasis(const double a, const double xi) const 
		{
			assert(validLocalIndex(a));
			const DoubleVec basis = curve().getBasisFuncs(paramCoord(xi), knotSpanI()); 
			return basis.at(a);
		}
		
		/// evaluate the derivative on a particular basis function
		/// w.r.t. the local coordinate
		double evalBasisDeriv(const uint a, const double xi) const
		{
			assert(validLocalIndex(a));
			const DoubleVec basis = curve().getBasisFuncDers(paramCoord(xi), knotSpanI()); 
			return basis.at(a) * jacobLocalTransform();
		}

		/// Evaluate basis function derivatives w.r.t. physical coordinates
		double evalGlobalBasisDeriv(const uint a, const double xi) const
		{
			assert(validLocalIndex(a));
			const double s = paramCoord(xi);
			const DoubleVec basis = curve().getBasisFuncDers(s, knotSpanI());
			return basis.at(a) * 1.0 /  jacobParamTransform(s);
		}	
			
		/// Get the jacobian from parent space to physical space
		double jacob(const double xi) const
		{
			return jacobParamTransform(paramCoord(xi)) * jacobLocalTransform();
		}
		
		private:
		
		/// Is this local index valid. That is, is i <= p + 1 ?	
		bool validLocalIndex(const unsigned int i) const { return i < curve().getOrder() + 1;}

		/// Transform parent coordinate to parametric coordinate
		double paramCoord (const double xi) const
		{
			return (upperbound() - lowerbound()) * 0.5 * (xi + 1.0) + lowerbound();
		}

		/// Jacobian from parameter space to physical space
		double jacobParamTransform(const double s) const
		{
			return curve().getJacob(s);
		}

		/// Jacobian from local parent space to parameter space
		double jacobLocalTransform() const
		{
			return 0.5 * (upperbound() - lowerbound());
		}
		
		/// Reference to nurbs curve
		const NURBSCurve& mCurve;

		/// Element index. i.e. unique knot span index
		const uint mElementI;
		
		/// Boost element interval used for range checking
		const boost::icl::continuous_interval<double> mElInterval;

        /// Knot span index
		const uint mKnotSpanI;

		/// overload output operator
		friend std::ostream& operator<<(std::ostream& ost, const Element& e);	
	};


	/// Iterator that loops over all non-zero knot spans (i.e. elements) of a NURBS curve
	class IAllEls {

		public:

		/// Construtor
		explicit IAllEls(const NURBSCurve& c)
			: mCurve(c), mCurrentElI(0) {}

		/// curve getter
		const NURBSCurve& curve() const{ return mCurve; }

		/// element index getter
		uint index() const { return mCurrentElI; }
		
		/// return element object
		Element currentEl() const
		{
			return Element(curve(), index());
		}
		
		/// prefix increment
		IAllEls& operator++() { ++mCurrentElI; return *this; }

		/// postfix increment
		IAllEls operator++(int)
		{
			IAllEls copy(*this);
			++(*this);
			return copy;
		}

		/// have we finished iterating over all non-zero knot spans?
		bool isDone() const { return index() >= curve().getUniqueKnotN() - 1; }
		
		private:

		/// Reference to nurbs curve we are operating over
		const NURBSCurve& mCurve;

		/// Current knot span index
		uint mCurrentElI;
	};
	
}

#endif

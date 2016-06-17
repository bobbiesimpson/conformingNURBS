#ifndef NURBS_NURBS_SPACE_H
#define NURBS_NURBS_SPACE_H

#include <vector>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <boost/icl/continuous_interval.hpp>
#include <boost/icl/closed_interval.hpp>

#include "base.h"
#include "InputDataStructures.h"

namespace nurbs
{
	/// A representation of the set of basis functions that
	/// constitute a b-spline space. We deliberately separate this
	/// from the geometry information since the spaces which
	/// make up the De Rham sequence require no knowledge of the
	/// geometry.
	class BSplineSpace {

		public:

		/// Default constructor
		BSplineSpace() = default;
		
		/// Constructor
		BSplineSpace(const DoubleVecVec& knotvec,
					 const UIntVec& d,
					 std::string n = "no_name")
			: mKnotVecs(knotvec),
			  mDegrees(d),
			  mName(n)
		{
			if( knotvec.size() != 2 )
				error( "Currently only R^2 nurbs spaces are implemented." );
			init();
		}

		/// Alternative constructor. Simply delegates work to former constructor
		BSplineSpace(const DoubleVecVec& knotvec,
					 const uint p,
					 const uint q,
					 std::string n = "no_name")
			: BSplineSpace(knotvec,{p,q}, n) {}

		/// Return the total number of basis functions
		uint basisFuncN() const
		{
			uint n = 1;
			for(const auto& v : mNumBasisVec)
				n *= v;
			return n;
		}
		
		/// basis function number getter
		uint basisFuncN(const ParamDir dir) const
		{
			assert(dir < mNumBasisVec.size());
			return mNumBasisVec[dir];
		}
        
        /// Get a vector of degrees
        UIntVec degree() const { return mDegrees; }
		
		/// degree getter
		inline uint degree(const ParamDir dir) const 
		{
			assert(dir < mDegrees.size());
			return mDegrees[dir];
		}
        
        /// Const name accessor
        std::string name() const { return mName; }
        
        /// Non-const name accessor
        std::string& name() { return mName; }

		/// Get knot vector for specified parametric direction
		inline DoubleVec knotVec(const ParamDir dir) const
		{
			assert(dir < paramDimN());
			return mKnotVecs[dir];
		}
		
		/// knot coordinate getter
		inline double knot(const uint i, const ParamDir dir) const
		{
			assert(dir < mKnotVecs.size());
			assert(i < mKnotVecs[dir].size());
			return mKnotVecs[dir][i];
		}
        
        /// Unique knot vector getter
        std::vector<double> uniqueKnotVec(const ParamDir dir) const
        {
            return mUniqueKnotVecs[dir];
        }

		/// unique knot coordinate getter
		inline double uniqueKnot(const uint i, const ParamDir dir) const
		{	
			assert(dir < mUniqueKnotVecs.size());
			assert(i < mUniqueKnotVecs[dir].size());
			return mUniqueKnotVecs[dir][i];
		}

		/// unique knot number getter
		inline uint uniqueKnotN(const ParamDir dir) const
		{
			assert(dir < mUniqueKnotVecs.size());
			return mUniqueKnotVecs[dir].size();
		}
        
        /// return number of times knot value is repeated
        inline uint knotRepeats(const uint i, const ParamDir dir) const
        {
            assert(i < uniqueKnotN(dir));
            const auto kv = knotVec(dir);
            return std::count(kv.begin(), kv.end(), uniqueKnot(i, dir));
        }
        
        /// Get the range for the given parametric direction
        std::pair<double, double> range(const ParamDir dir) const
        {
            return std::make_pair(mUniqueKnotVecs[dir].front(), mUniqueKnotVecs[dir].back());
        }
        
        /// Get a 'global' greville abscissa index given parametric indices
        uint globalGrevillePtI(const uint i,
                               const uint j) const
        {
            assert(i < grevilleAbscissaPtN(ParamDir::S));
            assert(j < grevilleAbscissaPtN(ParamDir::T));
            return j * grevilleAbscissaPtN(ParamDir::S) + i;
            
        }
        /// Number of greville abscissa points in given direction
        inline uint grevilleAbscissaPtN(const ParamDir dir) const { return basisFuncN(dir); }
        
        /// Number of greville abscissa points on this space
        inline uint grevilleAbscissaPtN() const { return basisFuncN(); }
        
        /// Greville abscissa point getter
        double grevilleAbscissaPt(const uint i, const ParamDir dir) const
        {
            assert(i < grevilleAbscissaPtN(dir));
            double xi = 0.0;
            for(uint k = i + 1; k < i + degree(dir) + 1; ++k)
                xi += knot(k, dir) / degree(dir);
            return xi;
        }
        
        /// Return the greville abscissa vector for a given parametric direction
        std::vector<double> grevilleAbscissa(const ParamDir dir) const
        {
            std::vector<double> v;
            for(uint i = 0; i < grevilleAbscissaPtN(dir); ++i)
                v.push_back(grevilleAbscissaPt(i, dir));
            return v;
        }
        
        /// Get a two dimensionl Greville point
        GPt2D grevilleAbscissaPt(const uint i) const
        {
            assert(i < grevilleAbscissaPtN());
            const uint sindex = i % grevilleAbscissaPtN(S);
            const uint tindex = i / grevilleAbscissaPtN(S);
            return GPt2D(grevilleAbscissaPt(sindex, S),
                         grevilleAbscissaPt(tindex, T));
        }

        /// REturn the Greville abscissa  points along the given edge
        std::vector<GPt2D> grevilleAbscissaPts(const Edge e) const
        {
            std::vector<GPt2D> v;
            switch (e) {
                case Edge::EDGE0:
                    for(uint i = 0; i < grevilleAbscissaPtN(S); ++i)
                        v.push_back(GPt2D(grevilleAbscissaPt(i, S), range(T).first));
                    break;
                case Edge::EDGE1:
                    for(uint i = 0; i < grevilleAbscissaPtN(S); ++i)
                        v.push_back(GPt2D(grevilleAbscissaPt(i, S), range(T).second));
                    break;
                case Edge::EDGE2:
                    for(uint j = 0; j < grevilleAbscissaPtN(T); ++j)
                        v.push_back(GPt2D(range(S).first, grevilleAbscissaPt(j, T)));
                    break;
                case Edge::EDGE3:
                    for(uint j = 0; j < grevilleAbscissaPtN(T); ++j)
                        v.push_back(GPt2D(range(S).second, grevilleAbscissaPt(j, T)));
                    break;
            }
            return v;
        }
        
        // Number of non-zero knot spans in a given direction
        inline uint nonzeroKnotSpanN(const ParamDir dir) const
        {
            return mUniqueKnotVecs[dir].size() - 1;
        }
        
		/// Total number of non-zero knot spans (i.e. elements)
		inline uint nonzeroKnotSpanN() const
		{
			uint c = 1;
			for(const auto& v : mUniqueKnotVecs)
				c *= v.size() - 1;
			return c;
		}

		/// For now, we are dealing with spaces in R^2
		inline uint paramDimN() const { return 2; }
		
		/// Get b-spline basis function at given parametric coordinate.
		/// Prefer to use other member function however for efficiency.
		inline DoubleVecVec tensorBasis(const double s, const double t) const
		{return tensorBasis(s,t,span(s,t));}
		
		/// Get b-spline basis with precalculated knot span index
		DoubleVecVec tensorBasis(const double s, const double t, const UIntVec& span) const;

		DoubleVec basis(const double s, const double t) const
		{return basis(s,t,span(s,t));}
		
		/// Get basis functions (not in tensor product format)
		DoubleVec basis(const double s, const double t, const UIntVec& span) const
		{
			DoubleVecVec bvec = tensorBasis(s,t,span);
			DoubleVec v;
			for(const auto& b_j : bvec[T])
				for(const auto& b_i : bvec[S])
					v.emplace_back(b_i * b_j);
			return v;
		}
		
		/// Get b-spline basis function derivatives at given parametric coordinate.
		/// Prefer to use other member function however for efficiency.
		DoubleVecVec tensorBasisDers(const double s, const double t, const DerivOrder der = D1) const
		{return tensorBasisDers(s,t,span(s,t),der);}

		/// Get b-spline basis with precalculate knot span index
		DoubleVecVec tensorBasisDers(const double s, const double t,
									 const UIntVec& span, const DerivOrder der = D1) const;

		/// Get basis function derivatives with tensor product multiplied out
		DoubleVec basisDers(const double s, const double t, const DerivOrder der = D1) const
		{return basisDers(s,t,span(s,t),der);}

		/// Get basis function derivatives with tensor product multiplied out
		DoubleVec basisDers(const double s, const double t, const UIntVec& span,
							const DerivOrder der = D1) const
		{
			DoubleVecVec bvec = tensorBasisDers(s,t,span,der);
			DoubleVec v;
			for(const auto& b_j : bvec[T])
				for(const auto& b_i : bvec[S])
					v.emplace_back(b_i * b_j);
			return v;			
		}
        
        /// Get basis function derivatives with specified derivative
        DoubleVec basisDers(const double s,
                            const double t,
                            const DerivType dtype) const
        { return basisDers(s,t,span(s,t), dtype); }
        
        /// Get basis function derivatives
        DoubleVec basisDers(const double s,
                            const double t,
                            const UIntVec& span,
                            const DerivType dtype) const
        {
            const DoubleVecVec tbasis = tensorBasis(s, t, span);
            const DoubleVecVec tbasis_d = tensorBasisDers(s, t, span);
            DoubleVec dvec;
            if(DS == dtype) {
                for(const auto& b_j : tbasis[T])
                    for(const auto& bd_i : tbasis_d[S])
                        dvec.push_back(bd_i * b_j);
            }
            else if(DT == dtype) {
                for(const auto& b_j : tbasis_d[T])
                    for(const auto& bd_i : tbasis[S])
                        dvec.push_back(bd_i * b_j);
            }
            else
                throw std::runtime_error("Bad derivative type specified.");
            return dvec;
        }
		
		/// Return the non-zero local basis function indices for this parametric coordinate
		/// Returns a set of vectors corresponding to the indices in each parametric direction
		UIntVecVec localBasisFuncI(const double s, const double t) const;

		/// Get the non-zero basis function indices using a row major numbering system.
		/// These are global in the sense that we return a single vector of indices rather
		/// than a tensor product format.
		UIntVec globalBasisFuncI(const double s, const double t) const
		{
			UIntVec gb_vec;
			const uint nb_s = basisFuncN(S);
			auto v = localBasisFuncI(s,t);
			for(const auto& t_index : v[T])
				for(const auto& s_index : v[S])
					gb_vec.push_back(t_index * nb_s + s_index);
			return gb_vec;
		}

		/// Get the knot indices for given parametric coordinate
		UIntVec span(const double s, const double t) const;

        /// Degree reduce in given parametric direction
        void degreeReduce(const ParamDir dir);
        
        /// Degree reduce space in both parametric directions
        void degreeReduce()
        {
            degreeReduce(S);
            degreeReduce(T);
        }
        
        /// Degree elevate in given parametric direction
        void degreeElevate(const ParamDir dir);
        
        /// Apply h-refinement (knot insertion) n times
        void hrefine(const uint n = 1);
		
		/// Load from an input stream
		void load(std::istream& ist);
		
		private:

		/// recalculate unique knots and intervals
		void init()
		{
            mUniqueKnotVecs.clear(); // clear data
            mIntervalVec.clear();
            mNumBasisVec.clear();
            
			for(uint i = 0; i < paramDimN(); ++i) {
				mIntervalVec.push_back(boost::icl::construct<Interval>
									   (mKnotVecs[i].front(), mKnotVecs[i].back(),
										boost::icl::interval_bounds::closed()));
				Interval& last = mIntervalVec.back();
				if(essentiallyEqual(last.upper(), last.lower(), TOL))
					error("Cannot prescribe a b-spline space with non-zero parametric area.");
				mNumBasisVec.push_back(mKnotVecs[i].size() - mDegrees[i] - 1);
			}
			// construct unique knot vectors
			for(auto kv : mKnotVecs) {
				auto last = std::unique(kv.begin(), kv.end());
				mUniqueKnotVecs.emplace_back( DoubleVec( kv.begin(), last ) );
			}			
		}
		
		/// Clear all data
		void clear()
		{
			mKnotVecs.clear();
			mUniqueKnotVecs.clear();
			mDegrees.clear();
			mNumBasisVec.clear();
			mIntervalVec.clear();
		}
		
		/// Is this a valid parametric coordinate? I.e. does it lie in the parametric space?
		bool validCoord(const double s, const double t) const
		{
			DoubleVec p{s,t};
			for(uint i = 0; i < paramDimN(); ++i) 
				if(!boost::icl::contains(mIntervalVec[i], p[i]))
					return false;
			return true;
		}
		
		/// Global knot vectors
		DoubleVecVec mKnotVecs;

		/// Unique global knot vectors (define elements with non-zero
		/// parametric area
		DoubleVecVec mUniqueKnotVecs;

		/// Vector of degrees for each direction
		UIntVec mDegrees;

		/// Vector containing number of basis function in each direction
		UIntVec mNumBasisVec;

		/// Vector of interval that specify bounds of this space
		std::vector<Interval> mIntervalVec;

		/// The name of this space
		std::string mName;
		
		/// Print function
		void printData(std::ostream& ost) const;

		/// overload output operator
		friend std::ostream& operator<<(std::ostream& ost, const BSplineSpace& s)
		{
			s.printData(ost);
			return ost;
		}

		/// Overload input operator
		friend std::istream& operator>>(std::istream& ist, BSplineSpace& s)
		{
			s.load(ist);
			return ist;
		}
		
	};

}

#endif

//
//  IGalerkinIntegrate.h
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#ifndef NURBS_IGALERKIN_INTEGRATE
#define NURBS_IGALERKIN_INTEGRATE

#include <iostream>
#include <vector>
#include <cassert>
#include <utility>
#include <algorithm>

#include "Point3D.h"
#include "IElemIntegrate.h"

namespace nurbs {
    
    using namespace elem;
    
    /// A representation of a 4D gauss point used for Galerkin BEM
    
    struct GPt4D {
        
        /// Simple constructor
        GPt4D(const double c1 = 0.0,
              const double c2 = 0.0,
              const double c3 = 0.0,
              const double c4 = 0.0) : data{c1, c2, c3, c4}
        {}
        
        /// Construct with two 2D gauss points
        GPt4D(const GPt2D& p1,
              const GPt2D& p2)
        : data{p1.s, p1.t, p2.s, p2.t}
        {}
        
        /// Construct with vector
        GPt4D(const std::vector<double>& v) : data(v)
        { assert(4 == data.size()); }
        
        /// Construct with two Point3d instances
        GPt4D(const Point3D& p1, const Point3D& p2)
        : data{ p1.getCoord(0), p1.getCoord(1), p2.getCoord(0), p2.getCoord(1) }
        {}
        
        /// Overload *= operator
        template<typename T>
        GPt4D& operator*=(const T& s)
        {
            std::for_each(data.begin(), data.end(), [=](double& c) { c *= static_cast<T>(s); });
            return *this;
        }
        
        /// overload const subscript operator
        const double& operator[](const int i) const
        {
            assert(i <= 3 && i >= 0);
            return data[i];
        }
        
        /// overload non-const subscript operator
        double& operator[](const int i)
        {
            assert(i <= 3 && i >= 0);
            return data[i];
        }
        
        /// Get source point as 2D gauss point
        GPt2D srcPt() const
        { return GPt2D(data[0], data[1]); }
        
        /// Get field point as 2D gauss point
        GPt2D fieldPt() const
        { return GPt2D(data[2], data[3]); }
        
        /// Rotate gauss point with rotation matrix
        GPt4D rotate(const DoubleVecVec& smat,
                     const DoubleVecVec& fmat) const;
        
        /// Vector while holds 4-d coordinates
        std::vector<double> data;
        
        /// overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const GPt4D& p);
    };
    
    template<typename T>
    GPt4D operator*(GPt4D p, const T& v)
    {
        return p *= v;
    }
    
    /// An iterator class that produces integration points and
    /// weights for a Galerkin 3D BEM integral.
    
    /// Points are returned as a 4D gauss pt (see GPt above) in
    /// [-1,1]^4
    
	class IGalerkinIntegrate {
        
    public:
        
		/// get the current quadrature point (in 4D space)
		GPt4D get() const { return getImpl(); }
        
		/// get the current quadrature 'weight'
		/// This silently also incorporates a jacobian
		/// for the transformation that removes the singularity
		double getWeight() const { return getWeightImpl(); }
        
		/// Have we iterated through all quadrature points?
		bool isDone() const { return subCellI() >= subCellN(); }
        
		/// Restart the quadrature iteration
		virtual void restart()
		{
			std::for_each(mBaseQuadratureVec.begin(), mBaseQuadratureVec.end(),
						  [](IElemIntegrate& i){ i.restart(); });
			mSubCellI = 0;
            mDidUpdateSrcI = true;
		}
		
		/// Prefix increment
        IGalerkinIntegrate& operator++() { incrementImpl(); return *this; }
        
        /// Did we increment the source point index? Used for efficiency
        /// gains in Galerkin BEM.
        bool didIncrementSrcI() const { return mDidUpdateSrcI; }
        
		/* /// Postfix increment */
		/* IGalerkinIntegrate operator++(int)  */
		/* { */
		/* 	IGalerkinIntegrate t(*this); */
		/* 	++(*this); */
		/* 	return t; */
		/* } */
        
		/// Return the total number of quadrature points
		/// Assumes idential quadrature routine is used in each sub cell
		uint pointN() const
		{
			uint s = 0;
			std::for_each(mBaseQuadratureVec.begin(), mBaseQuadratureVec.end(),
						  [&s](const elem::IElemIntegrate& i){ s += i.pointN(); });
			return s * subCellN();
		}
        
        /// Output all quadrature points to dat file.
        void outputDatFile(const std::string& sfile,
                           const std::string& ffile);
		
    protected:
        
		/// Enumeration for point type
		enum PointLocation { SOURCE = 0, FIELD };
		
        /// Construct empty iterator
        IGalerkinIntegrate()
        :
        mSubCellI(0),
        mDidUpdateSrcI(true) {}

        /// Constructor
		IGalerkinIntegrate(const UIntVec& sourceorder,
						   const UIntVec& fieldorder)
        : mBaseQuadratureVec{IElemIntegrate(sourceorder), IElemIntegrate(fieldorder)},
          mSubCellI(0),
          mDidUpdateSrcI(true)
		{}
        
		/// Get current base quadrature point for specified location
		GPt2D getBasePt(PointLocation p) const
		{
			return mBaseQuadratureVec.at(p).get();
		}
        
		/// Get current base quadrature weight for specified location
		double getBaseWeight(PointLocation p) const
		{
			return mBaseQuadratureVec.at(p).getWeight();
		}
        
		/// Get the current subcell index
		uint subCellI() const { return mSubCellI; }
        
		/// base class print implementation
		void print(std::ostream& ost) const;
        
		/// Get the current integration point in the space [0,1]^4
		GPt4D getCurrentUnitIntervalPt() const;
        
        /// Get the current integration point in the space [-1,1]^4
        GPt4D getCurrentBiUnitIntervalPt() const;
        
        /// Increment the sub cell index
        void incrementSubCellI() { ++mSubCellI; }
        
        /// Subcell index setter
        void setCurrentSubCellI(const uint i) { mSubCellI = i; }
        
    private:
        
        /// increment implementation. Can be overridden if necessary.
        virtual void incrementImpl();
        
		/// Non-const getter for base quadrature
		elem::IElemIntegrate& quadrature(PointLocation p)
		{
			return mBaseQuadratureVec.at(p);
		}
        
		/// Const-getter for base quadrature
		const elem::IElemIntegrate& quadrature(PointLocation p) const
		{
			return mBaseQuadratureVec.at(p);
		}
		
		/// Implementation function of get.
		virtual GPt4D getImpl() const = 0;
        
		/// Implementation function of getWeight.
		virtual double getWeightImpl() const = 0;
        
		/// Number of subcells. Must be implemented
		virtual uint subCellN() const = 0;
        
		/// Print implementation for derived classes. Optional implementation.
		virtual void printImpl(std::ostream& ost) const
		{
			ost << "Derived class print function not implemented.\n";
		}
		
		/// Vector holding quadrature for source and field
		/// quadrature points
		std::vector<elem::IElemIntegrate> mBaseQuadratureVec;
        
		/// Index of current subcell
		uint mSubCellI;
        
        /// Did we update the source index?
        bool mDidUpdateSrcI;
        
		/// Overload outpout operator
		friend std::ostream& operator<<(std::ostream& ost, const IGalerkinIntegrate& g);
	};
    
	/// non-member, non-friend helper functions
    
	/// Typedef for interval
	typedef std::pair<double, double> Range;
    
	/// Interval types
	enum IntervalType { PARENT_SPACE_INTERVAL, UNIT_PARENT_SPACE_INTERVAL };
    
	/// Get an appropriate interval given its type
	Range make_interval(IntervalType type);
	
	/// Convert a 4-d point in the interval defined by i1 to a 4d
	/// point in the interval defined by i2.
	GPt4D convert_interval(const GPt4D& p,
                           const Range& from_range,
                           const Range& to_range);
    
}

#endif

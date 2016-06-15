#ifndef NURBS_DOF_MGR
#define NURBS_DOF_MGR

#include <map>
#include <ostream>
#include <complex>
#include <vector>

#include "base.h"

namespace nurbs 
{
	/// This class stores coefficients that are interpolated using a
	/// B-spline space.

	/// The class is parameterised by the following parameters:
	/// D: dimension of each dof.  E.g. a weighted B-spline 3D geometry, D = 4.
	/// T: the datatype of each coefficient. In most cases this will be either
	/// double or complex<double>.
	
	template<uint D, typename T>
	class DOFManager {
		
		public:
		
		/// Typedef will allows the datatype to be determined
		typedef T datatype;

		/// Typedef for the vector type
		typedef std::vector<T> vectype;
		
		DOFManager(const uint dofn, T v)
			: mCoeffVec(dofn * D, v) {}
		
		/// The dimension of each dof
		inline uint componentN() const {return D;}

		/// Total number of degrees of freedom
		inline uint dofN() const {return mCoeffVec.size() / componentN();}

		inline T get(const uint i, const uint c) const
		{
			assert(i < dofN() && c < componentN());
			return mCoeffVec[map(i,c)];
		}

		inline uint getIndex(const uint i, const uint c) const
		{
			if(!(i < dofN() && c < componentN()))
				error("Bad index requested in DOFManager");
			return map(i,c);
		}
		
		/// The total number of coefficients stored.
		inline uint size() const {return dofN() * componentN();}

		/// Const access operator
		inline T operator()(const int i, const int c) const
		{ return mCoeffVec[getIndex(i,c)]; }

		/// Non-const accessor
		inline T& operator()(const int i, const int c)
		{ return mCoeffVec[getIndex(i,c)]; }

		/// Return a vector of all the components for the specified dof.
		std::vector<T> getVec(const uint i) const
		{
			assert(i < dofN());
			std::vector<T> v(componentN());
			for(uint c = 0; c < componentN(); ++c)
				v[c] = get(i,c);
			return v;
		}
		
		private:

		/// A mapping that returns a global dof given a coefficient index
		/// and local dof index.
		inline uint map(const uint i, const uint c) const {return componentN() * i + c;}

		/// Print output for debugging
		void print(std::ostream& ost) const
		{
			ost << "-\nDOFManager\n-\n";
			for(uint i = 0; i < dofN(); ++i) {
				ost << "dof = " << i << "\t{";
				for(uint c = 0; c < componentN(); ++c) {
					ost << get(i,c);
					if(c == componentN() - 1)
						ost << "}";
					else ost <<  ",";
				}
				ost << "\n";
			}
		}
		
		/// Map which holds all coefficients
		std::vector<T> mCoeffVec;

		friend std::ostream& operator<<(std::ostream& ost, const DOFManager& m)
		{ m.print(ost); return ost; }
	};

	/// Some handy typedefs
	typedef DOFManager<4,double> DOFManager4D;
	typedef DOFManager<3,double> DOFManager3D;
	typedef DOFManager<1,double> DOFManager1D;
	typedef DOFManager<1, std::complex<double>> DOFManager1C;
	typedef DOFManager<3, std::complex<double>> DOFManager3C;
}



#endif

#ifndef NURBS_BASE_H
#define NURBS_BASE_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <limits>
#include <iterator>
#include <memory>
#include <complex>

#include <boost/icl/continuous_interval.hpp>

namespace nurbs
{
    
	/// error function (wraps a runtime exception)
	void error( const std::string& msg );
	
	/// Useful typedefs
	typedef unsigned int uint;
	typedef std::vector< double > DoubleVec;
	typedef std::vector< std::vector< double > > DoubleVecVec;
    typedef std::vector<int> IntVec;
	typedef std::vector< uint > UIntVec;
	typedef std::vector<UIntVec> UIntVecVec;
	typedef std::complex< double > ComplexDouble;
    typedef std::vector<std::complex<double>> ComplexDoubleVec;
    typedef std::pair<double, double> DoublePair;
    typedef std::vector<std::pair<double, double>> DoublePairVec;

	typedef boost::icl::continuous_interval<double> Interval;

	/// Tolerance that we use throughout the code
	const double TOL = std::numeric_limits< double >::epsilon();

	/// Default increment for NURBS curve interpolation
	const double INC = 0.1;

	/// Default number of grid points for vtk output
	const uint DEFAULT_NGRID_PTS = 12;

	/// Default number of guass points
	const uint DEFAULT_NGPS = 4;

	/// Default number of elements per quadtree cell
	const uint DEFAULT_QT_ELN = 1;

	/// Default quadtree tolerance
	const double DEFAULT_QT_TOL = 1.0e-3;

	/// Double error
	const double DOUBLE_ERROR = -9999999999.99;

	/// Uint error
	const uint INVALID_UINT = std::numeric_limits<uint>::max();

	/// pi = 3.141...
	const double PI = atan( 1.0 ) * 4.0;
    
    /// A simple struct for holding a gauss point
    struct GPt2D
    {
        GPt2D(const double s_in = 0.0,
              const double t_in = 0.0)
        :
        s(s_in),
        t(t_in) {}
        
        /// Access component values
        double get(const uint i) const
        {
            if(i == 0)
                return s;
            else if(i == 1)
                return t;
            else {
                error("Bad index specified for GaussPt struct");
                return DOUBLE_ERROR;
            }
        }
        
        /// Components
        double s;
        double t;
        
    };
    
    /// Overload comparison operator
    bool operator<(const GPt2D& g1, const GPt2D& g2);
    
    /// Overload output operator
    std::ostream& operator<<(std::ostream& ost, const GPt2D& gpt);
    
	/// Specifies the order of a B-spline derivative
	enum DerivOrder
	{
		D1 = 1,
		D2
	};
    
    /// Enumeration of derivative type
    enum DerivType {
        DS = 1,
        DT = 2
    };

	/// Specifies a parametric direction
	enum ParamDir
	{
		S = 0,
		T
	};
    
    /// Cast uint to ParamDir type
    ParamDir ParamDirType(const uint d);

	/// Representing each of the carteisan components
	enum CartesianComponent
	{
		X = 0,
		Y = 1,
		Z = 2
	};

    /// Boundary condition type enumeration
	enum BCType
	{
		DIRICHLET = 0,
		NEUMANN
	};
    
    /// A simple struct for representing a parametric coordinate
    typedef struct {
        double s;
        double t;
    } ParamCoord;

    /// Vertex enumeration
    /// 2 ----- 3
    /// |       |
    /// |       |
    /// |       |
    /// 0 ------1
    enum class Vertex {
        VERTEX0 = 0,
        VERTEX1 = 1,
        VERTEX2 = 2,
        VERTEX3 = 3
    };
    
    /// Overload vertex output operator
    std::ostream& operator<<(std::ostream& ost, Vertex v);
    
    /// Get vertex type given an unsigned integer
    Vertex vertexType(const uint v);
    
    /// Edge enumeration
    ///  ---1---
    /// |       |
    /// 2       3
    /// |       |
    ///  --0----
    enum class Edge {
        EDGE0 = 0,
        EDGE1 = 1,
        EDGE2 = 2,
        EDGE3 = 3
    };
    
    /// Overload output operator
    std::ostream& operator<<(std::ostream& ost, Edge e);
    
    /// Get edge type given an unsigned integer
    Edge edgeType(const uint e);
    
    /// Sign enumeration. Used for edges and vector basis 
    enum class Sign {
        POSITIVE,
        NEGATIVE
    };
    
    /// Cast Sign to double
    double asDouble(Sign s);
    
    /// Overload output operator for Sign
    std::ostream& operator<<(std::ostream& ost, Sign s);
    
    /// Vector basis direction
    enum class ContinuityType {
        NORMAL,
        TANGENT
    };
    
    /// Enumeration of singularity types.
    enum class AdjacencyType {
        NONE, /* i.e. no singularity */
        EDGE0,
        EDGE1,
        EDGE2,
        EDGE3,
        VERTEX0,
        VERTEX1,
        VERTEX2,
        VERTEX3,
        EQUAL /* i.e. the knot spans are identical */
    };
    
    /// Total number of edges for tree assuming tensor-product structure
    const uint NEDGES = 4;
    
    /// Total number of vertices for tree assuming tensor-product structure
    const uint NVERTICES = 4;

	/// overload output operator for BCType
	std::ostream& operator<<( std::ostream& ost, const BCType& bc );
	
	/// helper functions
	template< class T >
	bool approximatelyEqual( T a, T b, T epsilon )
	{
		return fabs( a - b ) <= ( ( fabs( a ) < fabs( b ) ? fabs( b ) : fabs( a ) ) * epsilon );
	}

	template< class T >
	bool essentiallyEqual( T a, T b, T epsilon )
	{
		return fabs( a - b ) <= ( ( fabs( a ) > fabs( b ) ? fabs( b ) : fabs( a ) ) * epsilon );
	}

	template< class T >
		bool lessThanOrEqual( T a, T b, T epsilon )
	{
		return ( a < b ) || essentiallyEqual( a, b, epsilon );
	}

	template< class T >
		bool greaterThanOrEqual( T a, T b, T epsilon )
	{
		return ( a > b ) || essentiallyEqual( a, b, epsilon );
	}

	/// print a vector
	template< typename T >
	void printVector( const std::vector< T >& vec, std::ostream& ost )
	{
		ost << "{";
        if(vec.size() != 0) {
            std::copy( vec.begin(), vec.end() - 1, std::ostream_iterator< T >( ost, "," ) );
            ost << vec.back();
        }
		ost <<  "}";
	}
    
    /// overload output operator for a generic vector
    template<typename T>
    std::ostream& operator<<(std::ostream& ost, const std::vector<T>& v)
    {
        ost << "{";
        for(typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
            ost << *it;
            if(it != v.end() - 1)
                ost << ", ";
        }
        ost << "}";
        return ost;
    }

	template< class T >
		void printVal( const std::string& n, T val, std::ostream& ost )
	{
		ost << n << " = " << val << "\n";
	}
	
	/// the factory function for creating unique_ptr instances (due to be implemented in c++14)
	/// Use as e.g. std::unique_ptr< Widget > pw = make_unique< Widget >();
	template< typename T, typename ...Args >
	std::unique_ptr< T > make_unique( Args&& ...args )
	{
		return std::unique_ptr< T >( new T( std::forward< Args >(args)... ) );
	}

	/// check the input stream if we have reached the specified end char. Set the stream
	/// to good if equal, otherwise throw an exception
	void endOfLoop( std::istream& ist, const char delim, const std::string& msg );
    
    /// Does the given parent coordinate lie in the standard [-1,1] x [-1,1]
    /// interval?
    bool validParentCoord(const double u, const double v);
    
    /// N Evenly spaced values defined in the interval [a,b]
    DoubleVec range(const double a, const double b, const uint N);
    
    /// Return a vector that splits 'mem' into equal parts.
    /// If a remainder exists, the last term is modified accordingly.
    std::vector<long int> bounds(long int parts, long int mem);
    
    
}

#endif

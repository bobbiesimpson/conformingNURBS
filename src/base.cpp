#include "base.h"
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <gsl/gsl_sf_legendre.h>
#include <boost/math/special_functions.hpp>

namespace nurbs
{
	void error( const std::string& msg )
	{
		throw std::runtime_error( "Error: " + msg );
	}
    
    bool operator<(const GPt2D& g1, const GPt2D& g2)
    {
        if(g1.s != g2.s)
            return g1.s < g2.s;
        return g1.t < g2.t;
    }
    
    std::ostream& operator<<(std::ostream& ost, const GPt2D& gpt)
    {
        ost << gpt.s << "\t" << gpt.t;
        return ost;
    }
    
    ParamDir ParamDirType(const uint d)
    {
        assert(d < 2);
        if(d == 0)
            return S;
        else
            return T;
    }

    Vertex vertexType(const uint v)
    {
        switch(v)
        {
            case 0: return Vertex::VERTEX0;
            case 1: return Vertex::VERTEX1;
            case 2: return Vertex::VERTEX2;
            case 3: return Vertex::VERTEX3;
            default: throw std::runtime_error("Bad vertex specified");
        }
    }

    Edge edgeType(const uint e)
    {
        switch (e) {
            case 0: return Edge::EDGE0;
            case 1: return Edge::EDGE1;
            case 2: return Edge::EDGE2;
            case 3: return Edge::EDGE3;
            default: throw std::runtime_error("Bad edge specified");
        }
    }
    
    double asDouble(Sign s)
    {
        switch(s)
        {
            case Sign::POSITIVE: return 1.0;
            case Sign::NEGATIVE: return -1.0;
			default: throw std::runtime_error("Bad sign specified");
        }
    }
    
    std::ostream& operator<<(std::ostream& ost, Vertex v)
    {
        switch (v) {
            case Vertex::VERTEX0:
                ost << "VERTEX0";
                break;
            case Vertex::VERTEX1:
                ost << "VERTEX1";
                break;
            case Vertex::VERTEX2:
                ost << "VERTEX2";
                break;
            case Vertex::VERTEX3:
                ost << "VERTEX3";
                break;
        }
        return ost;
    }
    
    std::ostream& operator<<(std::ostream& ost, Edge e)
    {
        switch (e) {
            case Edge::EDGE0:
                ost << "EDGE0";
                break;
            case Edge::EDGE1:
                ost << "EDGE1";
                break;
            case Edge::EDGE2:
                ost << "EDGE2";
                break;
            case Edge::EDGE3:
                ost << "EDGE3";
                break;
        }
        return ost;
    }
    
    std::ostream& operator<<(std::ostream& ost, Sign s)
    {
        if(Sign::POSITIVE == s) ost << "+";
        else ost << "-";
        return ost;
    }
    
    
	std::ostream& operator<<( std::ostream& ost, const BCType& bc )
	{
		if( bc == DIRICHLET ) ost << "Dirichlet";
		else ost << "Neumann";
		return ost;
	}
	
	void endOfLoop( std::istream& ist, const char delim, const std::string& msg )
	{
		if( ist.fail() )
		{
			ist.clear();
			char ch;
			if( ist >> ch && ch == delim )
				return;
			error( msg );
		}
	}
    
    bool validParentCoord(const double u, const double v)
    {
        Interval i = boost::icl::construct<boost::icl::continuous_interval<double>>(-1.0,1.0);
        if(boost::icl::contains(i,u) && boost::icl::contains(i, v))
            return true;
        return false;
    }
    
    DoubleVec range(const double a, const double b, const uint N)
    {
        assert(N > 0);
        if (1 == N) {
            return {(a + b) / 2.0};
        }
        else {
            DoubleVec rvec;
            const double h = (b - a) / (N - 1);
            for(uint i = 0; i < N; ++i)
                rvec.push_back(i * h + a);
            return rvec;
        }
    }
    
    double LegendrePDeriv(const uint n, const double x)
    {
        if(0 == n)
            return 0.0;
        else if(1 == n)
            return x;
        else
            //return (n * legendre_p(n,x) + LegendrePDeriv(n - 1, x) ) / x;
            return (2.0 * (n - 1.0) + 1.0) * boost::math::legendre_p(n-1,x) +  LegendrePDeriv(n - 2, x);
    }
    
    ComplexDoubleVec mieSurfaceCurrent(const double a,
                                       const double k,
                                       const double theta,
                                       const double phi)
    {
        double error = 100.0;
        const double tol = 1.0e-10;
        std::complex<double> jtheta, jphi;
        double curr_mag = 0.0;
        double prev_mag;
        uint n = 1;
        const double ct = std::cos(theta);
        const double st = std::sin(theta);
        
        while(error > tol) {
            const auto i = std::complex<double>(0.0, 1.0);
            const auto an = 1.0 / k * std::pow(i, -n) * (2.0 * n + 1.0) / (n * (n + 1.0));
            const auto h = boost::math::sph_hankel_2(n, k);
            const auto Hbar = k * h;
            const auto h_d = n / k * boost::math::sph_hankel_2(n,k) - boost::math::sph_hankel_2(n + 1, k);
            const auto Hbar_d = k * h_d + h;
            
            
            if(std::abs(theta) < 1.0e-7) {
                const double tri_n = 0.5 * (n * n + n);
                jtheta += std::cos(phi) * an * ( tri_n / Hbar_d
                                                + (i * -tri_n) / Hbar );
                jphi += std::sin(phi) * an * (-tri_n /  Hbar_d
                                              - tri_n / (i * Hbar) );
            }
            else if(std::abs(theta - PI) < 1.0e-7) {
                const double tri_n = 0.5 * (n * n + n) * std::pow(-1.0, n);
                jtheta += std::cos(phi) * an * ( tri_n / Hbar_d
                                                + (i * tri_n) / Hbar );
                jphi += std::sin(phi) * an * (tri_n /  Hbar_d
                                              - tri_n / (i * Hbar) );
            }
            else {
                double l_d[n+1];
                double l[n+1];
                gsl_sf_legendre_Plm_deriv_array(n, 1, ct, l, l_d);
                const auto P = l[n-1];
                const auto Pd = l_d[n-1];

                
                jtheta += std::cos(phi) * an * ( (st * Pd) / Hbar_d
                                                + (i * P) / (st * Hbar) );
                jphi += std::sin(phi) * an * (P / (st * Hbar_d)
                                              - (st * Pd) / (i * Hbar) );
            }
            prev_mag = curr_mag;
            curr_mag = std::sqrt( std::abs(jtheta) * std::abs(jtheta)
                                  + std::abs(jphi) * std::abs(jphi) );
            error = std::abs(prev_mag - curr_mag);
            ++n;
        }
        const double cp = std::cos(phi);
        const double sp = std::sin(phi);
        return {ct * cp * jtheta - sp * jphi,
                ct * sp * jtheta + cp * jphi,
                -st * jtheta};
    }
    
    std::vector<long int> bounds(long int parts, long int mem)
    {
        std::vector<long int>bnd;
        long int delta = mem / parts;
        long int reminder = mem % parts;
        long int N1 = 0, N2 = 0;
        bnd.push_back(N1);
        for (long int i = 0; i < parts; ++i) {
            N2 = N1 + delta;
            if (i == parts - 1)
                N2 += reminder;
            bnd.push_back(N2);
            N1 = N2;
        }
        return bnd;
    }
}

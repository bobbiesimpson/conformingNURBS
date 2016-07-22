//
//  Geometry.cpp
//  nurbslib
//
//  Created by Robert Simpson on 14/07/2014.
//
//

#include <ios>
#include <memory>
#include <string>
#include <numeric>

#include "Geometry.h"
#include "BSplineSpace.h"
#include "Forest.h"
#include "NURBSCommon.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"

namespace nurbs
{
    UIntVec Geometry::degree(const uint sp) const
    {
        return primalForest().space(sp).degree();
    }
    
    double Geometry::jacDet(const double s, const double t, const uint sp) const
    {
        return cross(tangent(s,t,sp,S), tangent(s,t,sp,T)).length();
    }
    
    /// Jacobian
    DoubleVecVec Geometry::jacob(const double s, const double t, const uint sp) const
    {
        return {tangent(s,t,sp,S).asVec(), tangent(s,t,sp,T).asVec()};
    }
    
    /// Normal
    Point3D Geometry::normal(const double s, const double t, const uint sp) const
    {
        Point3D n = cross(tangent(s,t,sp,S), tangent(s,t,sp,T));
        if(normalsFlipped())
            n *= -1;
        return n / n.length();
    }
    
    /// Tangent
    Point3D Geometry::tangent(const double s,
                              const double t,
                              const uint sp,
                              const ParamDir dir) const
    {
        const ParamDir dir2 = (dir == S) ? T : S;
		Point4D temp;
        const BSplineSpace& b_space = primalForest().space(sp);
		UIntVecVec indices = b_space.localBasisFuncI(s, t);
		UIntVec span = b_space.span(s, t);
		DoubleVecVec basis = b_space.tensorBasis(s,t,span);
		DoubleVecVec ders = b_space.tensorBasisDers(s,t,span);

        std::vector<double> result;
        for(uint j = 0; j < basis[T].size(); ++j)
            for(uint i = 0; i < ders[S].size(); ++i)
                result.push_back(basis[T][j] * ders[S][i]);
//        std::cout << result << "\n";
        
		Point3D a, ader;
		double w = 0.0, wder = 0.0;
		for(uint i = 0; i < b_space.degree(S) + 1; ++i ) {
			for(uint j = 0; j < b_space.degree(T) + 1; ++j ) {
				const std::size_t der_i = (dir == S) ? i : j;
				const std::size_t nder_i = (dir == S) ? j : i;
				const Point4D& p = point(sp, indices[S][i], indices[T][j]);
				a += p.asUnweighted() * basis[S][i] * basis[T][j];
				ader += p.asUnweighted() * ders[dir][der_i] * basis[dir2][nder_i];
				w += p.getWeight() * basis[S][i] * basis[T][j];
				wder += p.getWeight() * ders[dir][der_i] * basis[dir2][nder_i];
			}
		}
//        std::cout << "w = " << w << "\n\n";
//        std::cout << "wder = " << wder << "\n\n";
        
		return nurbshelper::getNonRationalDeriv( { a, ader }, { w, wder } );
    }
    
    /// Interpolate the surface
    Point3D Geometry::eval(const double s, const double t, const uint sp) const
    {
        const BSplineSpace& b_space = primalForest().space(sp);
		UIntVec indices = b_space.globalBasisFuncI(s,t);
		DoubleVec basis = b_space.basis(s,t);
		assert(basis.size() == indices.size());
		Point4D p;
		for(uint i = 0; i < basis.size(); ++i)
			p += point(sp,indices[i]) * basis[i];
		return p.asCartesian();
    }
    
    /// Write to vtu file.
    void Geometry::writeVTKOutput(const std::string& file,
                                  const uint nsample) const
    {
        /// TODO
    }
    
    bool Geometry::load(std::istream& ist)
    {
        char ch;
		if(ist >> ch && ch != '{') {
			ist.unget();
			ist.clear(std::ios_base::failbit);
			return false;
		}
        clear();
        Forest f_load;
        if(!(ist >> f_load))
            error("Cannot load forest into geometry data structure");
        mPrimalForest = f_load;
        
		if(ist >> ch && ch != '{') { // We might be reading something else
			ist.unget();
			ist.clear(std::ios_base::failbit);
			return false;
		}
		mCPts.clear();
		Point4D p;
		while(ist >> p)
 			mCPts.push_back(p);
		endOfLoop(ist, '}', "Bad control point set");
        return true;
    }
    
    bool Geometry::loadHBSFile(std::istream& ist)
    {
        if(!ist)
            error("Cannot load HBS file: bad istream");
        clear();
        
        uint d; // get the dimension
        std::string s;
        std::getline(ist, s);
        std::istringstream ifs(s);
        ifs >> d;
        if(d != 2)
            error("Currently only parametric dimension = 2 is implemented ");
        
        uint npts;
        std::getline(ist, s);
        ifs.clear(); ifs.str(s);
        ifs >> npts;
        
        uint ntrunk;
        std::getline(ist, s);
        ifs.clear(); ifs.str(s);
        ifs >> ntrunk;
        
        std::vector<Point4D> cp_vec;
        for(uint i = 0; i < npts; ++i) {
            std::getline(ist, s);
            ifs.clear(); ifs.str(s);
            double x,y,z,w;
            ifs >> x >> y >> z >> w;
            cp_vec.emplace_back(x,y,z,w);
        }
        
        std::vector<BSplineSpace> s_vec;
        std::vector<UIntVec> conn_vec;
        for(uint i = 0; i < ntrunk; ++i) {
            UIntVec degree;
            std::getline(ist, s);
            ifs.clear(); ifs.str(s);
            for(uint c = 0; c < d; ++c) {
                uint v;
                ifs >> v;
                degree.push_back(v);
            }
            
            UIntVec point_vec;
            std::getline(ist, s);
            ifs.clear(); ifs.str(s);
            for(uint c = 0; c < d; ++c) {
                uint v;
                ifs >> v;
                point_vec.push_back(v);
            }
            uint n = 1;
            std::for_each(point_vec.begin(), point_vec.end(), [&n](uint& val){ n *= val; });
            
            std::vector<DoubleVec> knot_vec(degree.size());
            for(uint c = 0; c < point_vec.size(); ++c) {
                std::getline(ist, s);
                ifs.clear(); ifs.str(s);
                double kval;
                while(ifs >> kval)
                    knot_vec[c].push_back(kval);
            }
            s_vec.emplace_back(knot_vec, degree);
            
            UIntVec conn;
            for(uint c = 0; c < n; ++c) {
                std::getline(ist, s);
                ifs.clear(); ifs.str(s);
                uint val;
                ifs >> val;
                conn.push_back(val);
            }
            conn_vec.push_back(conn);
        }
        mPrimalForest = Forest(s_vec, conn_vec);
        mPrimalForest.setGeometry(this);
        mCPts = cp_vec;
        return true;
    }
    
    
    void Geometry::print(std::ostream& ost) const
    {
        ost << "- Geometry object -\n\n";
        ost << "Primal forest: \n\n" << mPrimalForest << "\n";
        ost << "Control point set: \n\n";
        for(const auto& p : mCPts)
            std::cout << p << "\n";
        std::cout << "\n- end geometry object -\n";
    }
    
    void Geometry::rescale(const double sf)
    {
        for(uint ipoint = 0; ipoint < controlPtN(); ++ipoint) {
            auto& pt = point(ipoint);
            for(uint i = 0; i < 3; ++i)
                pt[i] *= sf;
        }
    }
    
    void Geometry::normalise()
    {
        const double max = std::numeric_limits<double>::max();
        const double min = std::numeric_limits<double>::lowest();
        Point3D minpt(max, max, max);
        Point3D maxpt(min, min, min);
        
        for(unsigned ipoint = 0; ipoint < controlPtN(); ++ipoint)
        {
            const auto& cpt = point(ipoint);
            for(unsigned i = 0; i < 3; ++i)
            {
                if(cpt[i] < minpt[i])
                    minpt[i] = cpt[i];
                if(cpt[i] > maxpt[i])
                    maxpt[i] = cpt[i];
            }
        }
        const auto diff = maxpt - minpt;
        std::vector<double> diffvec = diff.asVec();
        auto result = std::max_element(diffvec.begin(), diffvec.end());
        const double scalefactor = 1.0 / (*result);
        rescale(scalefactor);
        
    }
    
    std::istream& operator>>(std::istream& ist, Geometry& g)
    {
        g.load(ist);
        return ist;
    }
    
    /// Overload outut operator
    std::ostream& operator<<(std::ostream& ost, const Geometry& g)
    {
        g.print(ost);
        return ost;
    }
}
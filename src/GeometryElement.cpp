#include "GeometryElement.h"
#include "Geometry.h"

#include <boost/icl/closed_interval.hpp>

namespace nurbs {
    
    /// Definitino of static mutex
//    std::mutex GeometryElement::sMutex;
    
    GeometryElement::GeometryElement(const Geometry& g,
                                     const uint sp,
                                     const DoublePairVec& knots)
    : mGeom(g),
      mSpaceI(sp),
      mMutex(std::make_shared<std::mutex>()),
      mSize(std::make_pair(false, 0.0))
    {
        assert(knots.size() >= 2);
        for(const auto& kp : knots)
            mKnotIntervals.emplace_back(boost::icl::construct<boost::icl::continuous_interval<double>>(kp.first,
                                                                                                       kp.second,
                                                                                                       boost::icl::interval_bounds::closed()));
    }
    
    /// Get jacobian determinant from parent to physical space
    double GeometryElement::jacDet(const double u, const double v) const
    {
        auto it = mJDetCache.find(std::make_pair(u,v));
        if(it != mJDetCache.end())
            return it->second;
        ParamCoord c = paramCoord(u,v);
        const double jdet =  geometry().jacDet(c.s, c.t, spaceI()) * jacDetParam(u,v);
        insertCachedJDet(std::make_pair(u, v), jdet);
        return jdet;
    }
    
    Point3D GeometryElement::eval(const double u, const double v) const
    {
        auto it = mPointCache.find(std::make_pair(u, v));
        if(it != mPointCache.end())
            return it->second;
        ParamCoord c = paramCoord(u,v);
        const auto pt = geometry().eval(c.s, c.t, spaceI());
        insertCachedPoint(std::make_pair(u, v), pt);
        return pt;
    }
    
    Point3D GeometryElement::evalVertex(const Vertex v) const
    {
        switch (v) {
            case Vertex::VERTEX0:
                return eval(-1.0, -1.0);
                break;
            case Vertex::VERTEX1:
                return eval(1.0,-1.0);
                break;
            case Vertex::VERTEX2:
                return eval(-1.0, 1.0);
                break;
            case Vertex::VERTEX3:
                return eval(1.0, 1.0);
                break;
            default:
                error("Bad vertex enum");
                break;
        }
    }
    
    /// Get jacobian
    DoubleVecVec GeometryElement::jacob(const double u, const double v) const
    {
        auto it = mJacobCache.find(std::make_pair(u, v));
        if(it != mJacobCache.end())
            return it->second;
        ParamCoord c = paramCoord(u,v);
        const auto jacob_param = geometry().jacob(c.s, c.t, spaceI());
        const auto jacob_parent = jacobParam(u, v);
        DoubleVecVec r{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }};
        for(uint i = 0; i < 2; ++i)
            for(uint j = 0; j < 3; ++j)
                for(uint k = 0; k < 2; ++k)
                    r[i][j] += jacob_parent[i][k] * jacob_param[k][j];
        insertCachedJacob(std::make_pair(u, v), r);
        return r;
    }
    
    /// Get the normal
    Point3D GeometryElement::normal(const double u, const double v) const
    {
        ParamCoord c = paramCoord(u,v);
        return geometry().normal(c.s, c.t, spaceI());
    }
    
    /// Get the tangent
    Point3D GeometryElement::tangent(const double u, const double v, const ParamDir d) const
    {
        ParamCoord c = paramCoord(u,v);
        return geometry().tangent(c.s, c.t, spaceI(), d);
    }
    
    UIntVec GeometryElement::geometryDegree() const
    { return geometry().degree(spaceI()); }
    
    std::pair<bool, GPt2D> GeometryElement::containsParamPt(const double s, const double t) const
    {
        if(boost::icl::contains(boostInterval(S),s) && boost::icl::contains(boostInterval(T), t))
            return std::make_pair(true, parentCoord(s, t));
        else
            return std::make_pair(false, GPt2D());
    }
    
    bool GeometryElement::contains(const GeometryElement& e) const
    {
        if(containsParamPt(e.lowerBound()).first && containsParamPt(e.upperBound()).first)
            return true;
        return false;
    }
    
    void GeometryElement::print(std::ostream& ost) const
    {
        ost << "Geometry: " << &mGeom << "\n";
        ost << "Space index: " << mSpaceI << "\n";
        ost << "Knot intervals: ";
        for(const auto& i : mKnotIntervals)
            ost << "(" << i.lower() << ", " << i.upper() << ")  ";
    
    }
    
    double dist(const GeometryElement& e1, const GeometryElement& e2)
    {
        return dist(e1.eval(0.0,0.0), e2.eval(0.0, 0.0));
    }
}
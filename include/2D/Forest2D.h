#ifndef FOREST2D_H
#define FOREST2D_H

#include <vector>
#include <string>
#include <utility>
#include <memory>
#include <map>
#include <cassert>
#include <istream>

#include "BSplineSpace2D.h"
#include "base.h"
#include "AnalysisElement2D.h"

namespace nurbs {
    
    /// A forest (2d) is basically a set of bspline spaces
    /// with connectivity information. No geometry information
    /// is contained in a forest. This is instead represnted by
    /// a Geometry2D instance.
    
    /// A forest does however have a reference to a geometry instance.
    
    /// Forward declaration
    class Geometry2D;
    
    class Forest2D {
        
    public:
        
        /// Default constructor
        Forest2D()
        :
        mGeom(nullptr) { clear(); }
        
        /// Copy constructor. Required due to use of unique_ptr.
        Forest2D(const Forest2D& f);
        
        /// Assignment operator
        Forest2D& operator=(const Forest2D& f);
        
        /// Move constructor. Required due to use of unique_ptr.
        Forest2D(Forest2D&& f) = default;
        
        /// Move assignment operator
        Forest2D& operator=(Forest2D&& f) = default;
        
        /// Construct with geometry object
        Forest2D(const Geometry2D* g);
        
        /// Construct with vector of bspline spaces and
        /// vector of connectivites.
        Forest2D(const std::vector<BSplineSpace2D>& space_vec,
                 const std::vector<UIntVec>& conn_vec);
        
        /// clear all member data
        void clear();
        
        /// Number of bspline spaces
        inline uint spaceN() const { return mSpaces.size(); }
        
        /// Geometry getter
        inline const Geometry2D* geometry() const { return mGeom; }
        
        /// Geometry setter
        inline void setGeometry(const Geometry2D* g) { mGeom = g; }
        
        /// Bspline space getter
        inline const BSplineSpace2D& space(const uint i) const
        { return mSpaces.at(i); }
        
        /// Connectivity vector getter
        inline const UIntVec& connectivityVec(const uint i) const { return mNodalConn.at(i); }
        
        /// Get total number of elements in this forest
        uint elemN() const
        {
            uint n = 0;
            for(uint s = 0; s < spaceN(); ++s)
                n += space(s).nonzeroKnotSpanN();
            return n;
        }
        
        /// Element getter
        const AnalysisElement2D* element(const uint i) const;
        
        /// Get global index give a space index and a local nodal index
        inline uint globalI(const uint sp, const uint i) const
        {
            assert(sp < spaceN() && i < space(sp).basisFuncN());
            return mNodalConn.at(sp)[i];
        }
        
        /// get parametric interval given a bspline space index and local element
        /// index
        std::pair<double, double> knotInterval(const int sp, const uint iel) const
        {
            assert(sp < spaceN() && iel < space(sp).nonzeroKnotSpanN());
            const BSplineSpace2D& bs = space(sp);
            return std::make_pair(bs.uniqueKnot(iel), bs.uniqueKnot(iel + 1));
        }
        
        /// global dof for this forest
        uint globalDofN() const
        {
            uint n = 0;
            for(const auto& s : spaceVec())
                n += s.basisFuncN();
            return n;
        }
        
        /// load from a given input stream
        void load(std::istream& ist);
        
        /// print to given output stream
        void print(std::ostream& ost) const;
        
    protected:

    private:
        
        /// Const getter for vector of spaces
        inline const std::vector<BSplineSpace2D>& spaceVec() const { return mSpaces; }
        
        /// Non-const getter for vector of spaces
        inline std::vector<BSplineSpace2D>& spaceVec() { return mSpaces; }
        
        /// Const space map getter
        inline const std::map<std::string, uint>& spaceMap() const { return mSpaceMap; }
        /// space map non-const getter
        inline std::map<std::string, uint>& spaceMap() { return mSpaceMap; }
        
        /// Const nodal vector getter
        inline const std::map<uint, std::vector<uint>>& nodalConnVec() const { return mNodalConn; }
        
        /// Non-const nodal vector getter
        inline std::map<uint, std::vector<uint>>& nodalConnVec()  { return mNodalConn; }
        
        /// Add a space to the forest
        void addSpace(const BSplineSpace2D& s)
        {
            auto i = spaceMap().find(s.name());
            if(i != spaceMap().end())
                error("Cannot add two spaces to a forest with identical names");
            spaceVec().push_back(s);
            spaceMap().insert(std::make_pair(s.name(), spaceVec().size() - 1));
        }
        
        /// Does a space exist with the given string id? If so, return its
        /// index
        std::pair<bool, uint> validSpace(const std::string& s)
        {
            auto i = spaceMap().find(s);
            if(i == spaceMap().end())
                return std::make_pair(false, INVALID_UINT);
            return std::make_pair(true, i->second);
        }
        
        void addConnectivityVec(const std::string& s, const std::vector<uint>& v)
        {
            auto p = validSpace(s);
            if(!p.first)
                error("Bad connectivity vector specified to Forest2D instance");
            nodalConnVec().insert(std::make_pair(p.second, v));
        }
        
        /// Non-owning reference to geometry instance
        const Geometry2D* mGeom;
        
        /// Set of B-spline spaces
        std::vector<BSplineSpace2D> mSpaces;
        
        /// B-spline space mapping from B-spline space string to
        /// global space index
        std::map<std::string, uint> mSpaceMap;
        
        /// Tree connectivity.
        std::map<uint, std::vector<uint>> mNodalConn;
        
        /// Mapping from global element index to space index and local element
        /// index
        mutable std::map<uint, std::pair<uint, uint>> mElemIndexMap;
        
        /// Map from map from global element index to element instance
        /// This is a container for the element instances.
        mutable std::map<uint, std::unique_ptr<AnalysisElement2D>> mElems;
        
        /// Overload input operator
        friend std::istream& operator>>(std::istream& ist, Forest2D& f);
        
        /// Overload input operator
        friend std::ostream& operator<<(std::ostream& ost, const Forest2D& f);
    };
    
}
#endif
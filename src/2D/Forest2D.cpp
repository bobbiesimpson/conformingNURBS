#include <utility>
#include <memory>

#include "Forest2D.h"
#include "Geometry2D.h"
#include "base.h"
#include "InputDataStructures.h"
#include "AnalysisElement2D.h"

namespace nurbs {
    
    Forest2D::Forest2D(const Geometry2D* g)
    :
    mGeom(g)
    {
        if(g) { // g can be a nullptr, so check g exists
            clear();
            const Forest2D& primal = g->primalForest();
            mNodalConn = primal.mNodalConn;
            mSpaceMap = primal.mSpaceMap;
            mSpaces = primal.mSpaces;
        }
    }
    
    void Forest2D::clear()
    {
        mElemIndexMap.clear();
        mElems.clear();
        mNodalConn.clear();
        mSpaceMap.clear();
        mSpaces.clear();
    }
    
    Forest2D::Forest2D(const Forest2D& f)
    {
        mGeom = f.mGeom;
        mElemIndexMap = f.mElemIndexMap;
        /// perform deep copy of element instances
        for(const auto& e : f.mElems)
            mElems.insert(std::make_pair(e.first, make_unique<AnalysisElement2D>(*e.second)));
        mNodalConn = f.mNodalConn;
        mSpaceMap = f.mSpaceMap;
        mSpaces = f.mSpaces;
    }
    
    Forest2D& Forest2D::operator=(const Forest2D& f)
    {
        if(this == &f)  // check for self-assignment
            return *this;
        clear();
        mGeom = f.mGeom;
        mElemIndexMap = f.mElemIndexMap;
        /// perform deep copy of element instances
        for(const auto& e : f.mElems)
            mElems.insert(std::make_pair(e.first, make_unique<AnalysisElement2D>(*e.second)));
        mNodalConn = f.mNodalConn;
        mSpaceMap = f.mSpaceMap;
        mSpaces = f.mSpaces;
        return *this;
    }
    
    Forest2D::Forest2D(const std::vector<BSplineSpace2D>& space_vec,
                       const std::vector<UIntVec>& conn_vec)
    :
    mGeom(nullptr),
    mSpaces(space_vec)
    {
        if(space_vec.size() != conn_vec.size())
            error("Bad input data to Forest2D constructor: number of spaces"
                  "and size of connectivity vector not equal");
        for(uint s = 0; s < spaceN(); ++s) {
            mNodalConn.insert(std::make_pair(s, conn_vec[s]));
            mSpaceMap.insert(std::make_pair("space" + std::to_string(s), s));
        }
    }
    
    void Forest2D::load(std::istream& ist)
    {
        while(true) {
            BSplineSpace2D bs;
            if(!(ist >> bs))
                break;
            addSpace(bs);
        }
        endOfLoop(ist, '}', "Bad end of bspline space input");
        while (true) {
            ConnVecInput cv;
            if(!(ist >> cv))
                break;
            addConnectivityVec(cv.name, cv.data);
        }
        endOfLoop(ist, '}', "Bad end of connectivity vector input");
    }
    
    void Forest2D::print(std::ostream& ost) const
    {
        ost << "Forest: " << spaceN() << " trunks\n";
        for(const auto& n : spaceMap()) {
            ost << space(n.second) << "\n";
        }
        ost << "Nodal connectivity:\n";
        for(const auto& c : nodalConnVec())
            ost << c.first << ": " << c.second << "\n";
    }
    
    const AnalysisElement2D* Forest2D::element(const uint i) const
    {
        auto e = mElems.find(i);
        if(e != mElems.end())
            return e->second.get();
        uint start_i = 0;
        for(uint sp = 0; sp < spaceN(); ++sp) {
            const uint el_n = space(sp).nonzeroKnotSpanN();
            if((i - start_i) > (el_n - 1)) {
                start_i += el_n;
                continue;
            }
            const uint local_i = i - start_i;
            auto r = mElems.insert(std::make_pair(i, make_unique<AnalysisElement2D>(this, sp, local_i)));
            if(!r.second)
                error("Failed attempt to create element");
            auto r_i = mElemIndexMap.insert(std::make_pair(i, std::make_pair(sp, local_i)));
            if(!r_i.second)
                error("Error inserting element index mapping");
            return mElems[i].get();
        }
        error("Failed to create element"); return nullptr;
    }
    
    std::istream& operator>>(std::istream& ist, Forest2D& f)
    {
        f.load(ist);
        return ist;
    }
    
    std::ostream& operator<<(std::ostream& ost, const Forest2D& f)
    {
        f.print(ost);
        return ost;
    }
}
#include "MultiForest.h"
#include "Forest.h"
#include "Geometry.h"
#include "HCElement.h"
#include "BezierVectorElement.h"
#include "NedelecVectorElement.h"

#include <cassert>
#include <algorithm>
#include <bitset>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace nurbs
{
    /// Copy constructor
    MultiForest::MultiForest(const MultiForest& f)
    {
        mpGeom = f.mpGeom;
        mSpaceS = f.mSpaceS;
        mSpaceT = f.mSpaceT;
        mConn = f.mConn;
        mSignedConn = f.mSignedConn;
        mSignConn = f.mSignConn;
        mGlobalDofN = f.mGlobalDofN;
        mGlobalNonDegenerateDofN = f.mGlobalNonDegenerateDofN;
        for(const auto& e : f.mElems)
            mElems.insert(std::make_pair(e.first, e.second->copy()));
        for(const auto& e : f.mBezierElems)
            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        mElemN = f.mElemN;
    }
    
    /// Copy assignment operator
    MultiForest& MultiForest::operator=(const MultiForest& f)
    {
        if(this == &f)
            return *this;
        clear();
        mpGeom = f.mpGeom;
        mSpaceS = f.mSpaceS;
        mSpaceT = f.mSpaceT;
        mConn = f.mConn;
        mSignedConn = f.mSignedConn;
        mSignConn = f.mSignConn;
        mGlobalDofN = f.mGlobalDofN;
        mGlobalNonDegenerateDofN = f.mGlobalNonDegenerateDofN;
        for(const auto& e : f.mElems)
            mElems.insert(std::make_pair(e.first, e.second->copy()));
        for(const auto& e : f.mBezierElems)
            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        mElemN = f.mElemN;
        return *this;
    }
    
    void MultiForest::clear()
    {
        mSpaceS.clear();
        mSpaceT.clear();
        mConn.clear();
        mSignedConn.clear();
        mSignConn.clear();
        mGlobalDofN = 0;
        mGlobalNonDegenerateDofN = 0;
        mElems.clear();
        mBezierElems.clear();
        mElemN = std::make_pair(false, -1);
    }
    
    void MultiForest::initSpaces(const std::vector<uint>& sreduce,
                                 const std::vector<uint>& treduce)
    {
        assert(sreduce.size() >= 2);
        assert(treduce.size() >= 2);
        const Forest& gforest = geometry()->primalForest();
        const std::vector<std::vector<uint>> r_vecs{sreduce, treduce};
        const uint ndim = r_vecs.size();
        
        // first degree reduce according to reduce vectors
        for(uint b_dir = 0; b_dir < ndim; ++b_dir) {
            for(uint ispace = 0; ispace < gforest.spaceN(); ++ispace) {
                BSplineSpace scopy = gforest.space(ispace);
                
                // HARD CODED DEGREE ELEVATION
//                scopy.degreeElevate(S);
//                scopy.degreeElevate(T);
                // -- END HARD CODE
                
                for(uint dir = 0; dir < ndim; ++dir) {
                    for(uint reduce = 0; reduce < r_vecs[b_dir][dir]; ++reduce)
                        scopy.degreeReduce(ParamDirType(dir));
                }
                insertSpace(scopy, ParamDirType(b_dir));
            }
        }
        initConnectivity();
        initNedelecConnectivity();
        
    }
    
    uint MultiForest::elemN() const
    {
        if(mElemN.first)
            return mElemN.second;
        uint n = 0;
        for(uint s = 0; s < spaceN(); ++s)
            n += space(s, S).nonzeroKnotSpanN();
        mElemN = std::make_pair(true, n);
        return n;
    }
    
    const VAnalysisElement* MultiForest::element(const uint i) const
    {
        initEl(i);
        return mElems.at(i).get();
    }
    
    const VAnalysisElement* MultiForest::bezierElement(const uint i) const
    {
        initBezierEl(i);
        return mBezierElems.at(i).get();
    }
    
    void MultiForest::prefine(const uint nrefine)
    {
        if(0 == nrefine)
            return;
        
        uint count = 0;
        while(count < nrefine)
        {
            for(auto& s_space : mSpaceS)
            {
                s_space.degreeElevate(ParamDir::S);
                s_space.degreeElevate(ParamDir::T);
            }
            for(auto& t_space : mSpaceT) {
                t_space.degreeElevate(ParamDir::S);
                t_space.degreeElevate(ParamDir::T);
            }
            ++count;
        }
        initConnectivity();
        initNedelecConnectivity();
    }
    
    /// Apply uniform h-refinment to the multiforest
    void MultiForest::hrefine(const uint nrefine)
    {
        if(0 == nrefine)
            return;
        for(auto& s_space : mSpaceS)
            s_space.hrefine(nrefine);
        for(auto& t_space : mSpaceT)
            t_space.hrefine(nrefine);
        initConnectivity();
        initNedelecConnectivity();
        
    }
    
    Point3D MultiForest::collocPt(const uint ispace,
                                  const ParamDir d,
                                  const uint i) const
    {
        const GPt2D param_pt = space(ispace, d).grevilleAbscissaPt(i);
        return geometry()->eval(param_pt.s, param_pt.t, ispace);
    }
    
    DoublePairVec MultiForest::knotIntervals(const uint sp, const uint iel) const
    {
        const ParamDir dir = S; // choose S space for computations
        assert(sp < spaceN() && iel < space(sp, dir).nonzeroKnotSpanN());
        const BSplineSpace& bspace = space(sp, dir);
        DoublePairVec p_vec;
        const uint nel_s = bspace.uniqueKnotN(S) - 1;
        const uint i_index = iel % nel_s;
        const uint j_index = iel / nel_s;
        p_vec.emplace_back(std::make_pair(bspace.uniqueKnot(i_index, S), bspace.uniqueKnot(i_index + 1, S)));
        p_vec.emplace_back(std::make_pair(bspace.uniqueKnot(j_index, T), bspace.uniqueKnot(j_index + 1, T)));
        return p_vec;
    }
    
    double MultiForest::h() const
    {
        double max_h = 0.0;
        double total_a = 0.0;
        for(uint ielem = 0; ielem < elemN(); ++ielem) {
            const auto el = element(ielem);
            double a = 0.0;
            for(IElemIntegrate igpt(el->equalIntegrationOrder()); !igpt.isDone(); ++igpt) {
                const GPt2D gpt = igpt.get();
                a += el->jacDet(gpt.s, gpt.t) * igpt.getWeight();
            }
            total_a += a;
            max_h = (a > max_h) ? a : max_h;
        }
        return max_h / total_a;
    }
    
    void MultiForest::initConnectivity()
    {
        mConn.clear();
        mSignedConn.clear();
        mSignConn.clear();
        mGlobalDofN = 0;
        mGlobalNonDegenerateDofN = 0;
        mLocalElemIMap.clear();
        mElems.clear();
        mBezierElems.clear();
        mElemN = std::make_pair(false, 0);
        
        const int UNASSIGNED = -2; // flag for unassigned dof
        
        // Now assign connectivity
        std::cout << "Constructing multiforest connectivity\n";
        const Forest& gforest = geometry()->primalForest();
        uint g_index = 0;
        
        // loop over spaces
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            const BSplineSpace& sspace = space(ispace, S);
            const BSplineSpace& tspace = space(ispace, T);
            const auto space_pair = std::make_pair(sspace, tspace);
            const uint nbasis = sspace.basisFuncN() + tspace.basisFuncN();
            std::vector<int> g_indices(nbasis, UNASSIGNED); // -9 signifies unassigned index
            std::vector<Sign> g_sign(nbasis, Sign::POSITIVE); // default positive
            
            // loop over edges
            for(uint e = 0; e < NEDGES; ++e)
            {
                const uint iedge = gforest.globalEdgeI(ispace, e);
                
                for(const auto& s : gforest.connectedSpacesOnEdge(iedge))
                {
                    if(ispace == s) continue;
                    auto it = mConn.find(s);
                    if(it != mConn.end())
                    {
                        // the connectivity is assigned
                        Edge e_ref, e_in;
                        if(!gforest.connectedEdges(ispace, s, iedge, e_in, e_ref))
                            throw std::runtime_error("Bad connectivity");
                        const ParamDir param_ref = paramDir(e_ref, continuityType()); // get parametric direction given normal/tang. continuity
                        auto g_ivec = globalBasisIVec(s, e_ref, param_ref);
                        auto s_pair = vectorBasisEdgeDir(e_ref, param_ref, e_in);
                        if(Sign::NEGATIVE == edgeOrientation(e_in, e_ref))
                            std::reverse(g_ivec.begin(), g_ivec.end());
                        auto l_ivec = localBasisIVec(e_in, s_pair.second, space_pair);
                        assert(l_ivec.size() == g_ivec.size());
                        
                        for(uint i = 0; i < l_ivec.size(); ++i)
                        {
                            const uint l_i = l_ivec[i];
                            if(g_indices[l_i] != UNASSIGNED) continue;
                            g_indices[l_i] = g_ivec[i];
                            g_sign[l_i] = s_pair.first;
                        }
                        break; // we've dealt with the neighbouring space. Leave loop.
                    }
                    
                }
            }
            // Now fill in the remaining unassigned connectivities.
            // The sign of these indicies will be default positive.
            std::vector<uint> unsigned_givec(nbasis);
            for(uint i = 0; i < nbasis; ++i)
            {
                const int g_i = g_indices[i];
                if(UNASSIGNED == g_i) unsigned_givec[i] = g_index++;
                else unsigned_givec[i] = static_cast<uint>(g_i);
            }
            // and add the connectivity vectors to the global map
            insertConn(ispace, unsigned_givec);
            insertSignConn(ispace, g_sign);
            
            // add a copy for the degenerate connectivity
            std::vector<int> signed_copy;
            for(const auto& val : unsigned_givec)
                signed_copy.push_back(static_cast<int>(val));
            insertSignedConn(ispace, signed_copy);
            
        }
        mGlobalDofN = g_index; // assign no. of global dof.
        
        // now check for degenerate edges and set any dof on such edges equal to -1
        // while also recalculating the global number of dof.
        
        const int DEGENERATE = -1; // flag for degenerate dof
        uint current_index = 0;
        std::map<int, int> old2newdof_map;
        
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            const auto& sspace = space(ispace, ParamDir::S);
            const auto& tspace = space(ispace, ParamDir::T);
            
            for(uint e = 0; e < NEDGES; ++e)
            {
                if(gforest.degenerateEdge(ispace, e))
                {
                    const ParamDir param_ref = paramDir(edgeType(e), continuityType()); // get parametric direction given normal/tang. continuity
                    
                    if(ParamDir::S == param_ref)
                    {
                        for(const auto& lindex : localBasisIVec(edgeType(e), sspace))
                            setSignedConn(ispace, lindex, DEGENERATE);
                    }
                    else
                    {
                        for(const auto& lindex : localBasisIVec(edgeType(e), tspace))
                            setSignedConn(ispace, lindex + sspace.basisFuncN(), DEGENERATE);
                    }
                }
            }
            // degenerate dofs are flagged, so now assign new dof indices
            for(auto& index : mSignedConn[ispace])
            {
                if(DEGENERATE == index)
                    continue;
                auto search = old2newdof_map.find(index);
                if(search != old2newdof_map.end())
                    index = search->second;
                else
                {
                    old2newdof_map[index] = current_index;
                    index = current_index;
                    ++current_index;
                }
            }
        }
        mGlobalNonDegenerateDofN = current_index;


        // now assign edge and vertex adjacent of elements (for Galerkin integration purposes)
//        for(uint ispace = 0; ispace < spaceN(); ++ispace)
//        {
//            std::map<Edge, UIntVec> ge_map; // global element indices connected at edges
//            std::map<Vertex, uint> gv_map; // global element indices connected at vertices
//            
//            // populate ge_map
//            for(uint e = 0; e < NEDGES; ++e)
//            {
//                const uint iedge = gforest.globalEdgeI(ispace, e);
//                for(const auto& ie_space : gforest.connectedSpacesOnEdge(iedge))
//                {
//                    if(ispace == ie_space) continue;
//                    //std::cout << ispace << "\t" << ie_space << "\t" << e << "\n";
//                    Edge e1, e2;
//                    if(!gforest.connectedEdges(ispace, ie_space, iedge, e1, e2))
//                        throw std::runtime_error("Bad space connectivity");
//                    for(const auto& iel : localElIndices(space(ie_space,S), e2, edgeOrientation(e1, e2)))
//                        ge_map[edgeType(e)].push_back(globalElI(ie_space, iel));
//                }
//            }
//            
//            // now populate gv_map
//            for(uint v = 0; v < NVERTICES; ++v)
//            {
//                const UIntVec vspace_ivec = gforest.vertexConnectedSpaces(ispace, vertexType(v));
//                if(vspace_ivec.size() > 1)
//                    throw std::runtime_error("cannot deal with extraordinary points with valency > 4");
//                for(const auto& iv_space : vspace_ivec)
//                    gv_map[vertexType(v)] = globalElI(iv_space, localElIndex(space(iv_space, S), vertexType(v)));
//            }
//            
//            // now set vertex and element adjacency for all elements
//            const uint nel_s = space(ispace, S).uniqueKnotN(S) - 1;
//            const uint nel_t = space(ispace, S).uniqueKnotN(T) - 1;
//            for(int i = 0; i < nel_t; ++i)
//            {
//                for(int j = 0; j < nel_s; ++j)
//                {
//                    auto curr_el = element(ispace, i, j);
//                    auto curr_bel = bezierElement(ispace, i, j);
//                    
//                    curr_el->setVertexConnectedEl(Vertex::VERTEX0, element(ispace, i-1, j-1));
//                    curr_el->setVertexConnectedEl(Vertex::VERTEX1, element(ispace, i-1, j+1));
//                    curr_el->setVertexConnectedEl(Vertex::VERTEX2, element(ispace, i+1, j-1));
//                    curr_el->setVertexConnectedEl(Vertex::VERTEX3, element(ispace, i+1, j+1));
//                    curr_el->setEdgeConnectedEl(Edge::EDGE0, element(ispace, i-1, j));
//                    curr_el->setEdgeConnectedEl(Edge::EDGE1, element(ispace, i+1, j));
//                    curr_el->setEdgeConnectedEl(Edge::EDGE2, element(ispace, i, j-1));
//                    curr_el->setEdgeConnectedEl(Edge::EDGE3, element(ispace, i, j+1));
//                    
//                    curr_bel->setVertexConnectedEl(Vertex::VERTEX0, bezierElement(ispace, i-1, j-1));
//                    curr_bel->setVertexConnectedEl(Vertex::VERTEX1, bezierElement(ispace, i-1, j+1));
//                    curr_bel->setVertexConnectedEl(Vertex::VERTEX2, bezierElement(ispace, i+1, j-1));
//                    curr_bel->setVertexConnectedEl(Vertex::VERTEX3, bezierElement(ispace, i+1, j+1));
//                    curr_bel->setEdgeConnectedEl(Edge::EDGE0, bezierElement(ispace, i-1, j));
//                    curr_bel->setEdgeConnectedEl(Edge::EDGE1, bezierElement(ispace, i+1, j));
//                    curr_bel->setEdgeConnectedEl(Edge::EDGE2, bezierElement(ispace, i, j-1));
//                    curr_bel->setEdgeConnectedEl(Edge::EDGE3, bezierElement(ispace, i, j+1));
//                    
//                    // we're on Edge0
//                    if(0 == i)
//                    {
//                        curr_el->setEdgeConnectedEl(Edge::EDGE0, element(ge_map[Edge::EDGE0][j]));
//                        curr_bel->setEdgeConnectedEl(Edge::EDGE0, bezierElement(ge_map[Edge::EDGE0][j]));
//                        if(0 == j)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX0);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX0, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX0, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX0, element(ge_map[Edge::EDGE0][j-1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX0, bezierElement(ge_map[Edge::EDGE0][j-1]));
//                        }
//                        if(nel_s -1 == j)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX1);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX1, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX1, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX1, element(ge_map[Edge::EDGE0][j+1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX1, bezierElement(ge_map[Edge::EDGE0][j+1]));
//                        }
//                    }
//                    
//                    // we're on Edge1
//                    if(nel_t - 1 == i)
//                    {
//                        curr_el->setEdgeConnectedEl(Edge::EDGE1, element(ge_map[Edge::EDGE1][j]));
//                        curr_bel->setEdgeConnectedEl(Edge::EDGE1, bezierElement(ge_map[Edge::EDGE1][j]));
//                        if(0 == j)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX2);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX2, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX2, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX2, element(ge_map[Edge::EDGE1][j-1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX2, bezierElement(ge_map[Edge::EDGE1][j-1]));
//                            
//                        }
//                        if(nel_s - 1 == j)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX3);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX3, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX3, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX3, element(ge_map[Edge::EDGE1][j+1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX3, bezierElement(ge_map[Edge::EDGE1][j+1]));
//                        }
//                    }
//                    if(0 == j)
//                    {
//                        curr_el->setEdgeConnectedEl(Edge::EDGE2, element(ge_map[Edge::EDGE2][i]));
//                        curr_bel->setEdgeConnectedEl(Edge::EDGE2, bezierElement(ge_map[Edge::EDGE2][i]));
//                        if(0 == i)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX0);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX0, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX0, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX0, element(ge_map[Edge::EDGE2][i-1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX0, bezierElement(ge_map[Edge::EDGE2][i-1]));
//                        }
//                        if(nel_t - 1 == i)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX2);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX2, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX2, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX2, element(ge_map[Edge::EDGE2][i+1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX2, bezierElement(ge_map[Edge::EDGE2][i+1]));
//                        }
//                    }
//                    if(nel_s - 1 == j)
//                    {
//                        curr_el->setEdgeConnectedEl(Edge::EDGE3, element(ge_map[Edge::EDGE3][i]));
//                        curr_bel->setEdgeConnectedEl(Edge::EDGE3, bezierElement(ge_map[Edge::EDGE3][i]));
//                        if(0 == i)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX1);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX1, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX1, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX1, element(ge_map[Edge::EDGE3][i-1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX1, bezierElement(ge_map[Edge::EDGE3][i-1]));
//                        }
//                        if(nel_t - 1 == i)
//                        {
//                            auto it = gv_map.find(Vertex::VERTEX3);
//                            if(it != gv_map.end())
//                            {
//                                curr_el->setVertexConnectedEl(Vertex::VERTEX3, element(it->second));
//                                curr_bel->setVertexConnectedEl(Vertex::VERTEX3, bezierElement(it->second));
//                            }
//                        }
//                        else
//                        {
//                            curr_el->setVertexConnectedEl(Vertex::VERTEX3, element(ge_map[Edge::EDGE3][i+1]));
//                            curr_bel->setVertexConnectedEl(Vertex::VERTEX3, bezierElement(ge_map[Edge::EDGE3][i+1]));
//                        }
//                    }

                    // now get the bezier element and assign the same
                    // vertex and edge connectivities
//                    auto p_bel = bezierElement(ispace,i,j);
//                    
//                    // first assign vertex connected elements
//                    for(uint ivertex = 0; ivertex < NVERTICES; ++ivertex)
//                    {
//                        
//                        p_bel->setVertexConnectedEl(vertexType(ivertex), curr_el->getVertexConnectedEl(vertexType(ivertex)));
//                    }
//                    
//                    // and now assign edege connected elements
//                    for(uint iedge = 0; iedge < NEDGES; ++iedge)
//                    {
//                        p_bel->setEdgeConnectedEl(edgeType(iedge), curr_el->getEdgeConnectedEl(edgeType(iedge)));
//                    }
//                }
//            }
//        }
    }
    
    void MultiForest::initNedelecConnectivity()
    {
        mNedelecElems.clear();
        
        // The general idea is to copy the Bspline spaces and convert them
        // to Bezier form. The connectivity construction is then exactly
        // the same as in the Bspline case
        const std::vector<BSplineSpace> sspace_copy = mSpaceS;
        const std::vector<BSplineSpace> tspace_copy = mSpaceT;
        
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            mSpaceS[ispace].convertToBezier();
            mSpaceT[ispace].convertToBezier();
        }
        initConnectivity();
        
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            const auto& sp = space(ispace, ParamDir::S);
            const auto& sign_vec = this->globalDirVec(ispace);
            
            for(uint iel = 0; iel < sp.nonzeroKnotSpanN(); ++iel)
            {
                const auto bezier_el = bezierElement(ispace, iel);
                std::vector<Sign> sign_conn;
                
                for(const auto& ilocal : bezier_el->localBasisFuncI())
                    sign_conn.push_back(sign_vec[ilocal]);
                
                mNedelecElems.push_back(make_unique<NedelecVectorElement>(this,
                                                                          ispace,
                                                                          iel,
                                                                          bezier_el->globalBasisFuncI(),
                                                                          bezier_el->signedGlobalBasisFuncI(),
                                                                          sign_conn));
                
                auto nedelec_el = mNedelecElems[globalElI(ispace, iel)].get();
                
                for(uint ielem = 0; ielem < geometry()->primalForest().space(ispace).nonzeroKnotSpanN(); ++ielem) {
                    const auto pel = geometry()->primalForest().bezierElement(ispace, ielem);
                    if(pel->contains(*nedelec_el))
                        nedelec_el->setParent(pel);
                }
            }
        }
        
        mGlobalNedelecDofN  = globalDofN();
        
        // and now revert back to original discretisation
        mSpaceS = sspace_copy;
        mSpaceT = tspace_copy;
        initConnectivity();
    }
    
    UIntVec MultiForest::globalBasisIVec(const uint ispace,
                                         const Edge e,
                                         const ParamDir d) const
    {
        const BSplineSpace& space = this->space(ispace, d);
        const uint ns = space.basisFuncN(S);
        const uint nt = space.basisFuncN(T);
        std::vector<uint> g_vec; // return vector
        auto l_vec = localBasisIVec(e, ns, nt);
        for(const auto& l : l_vec)
            g_vec.push_back(globalI(ispace, l, d));
        return g_vec;
    }
    
    uint MultiForest::globalBasisI(const uint ispace,
                                   const Vertex v,
                                   const ParamDir d) const
    {
        const BSplineSpace& space = this->space(ispace, d);
        const uint ns = space.basisFuncN(S);
        const uint nt = space.basisFuncN(T);
        return globalI(ispace, localBasisI(v, ns, nt), d);
    }
    
    uint MultiForest::globalElI(const uint ispace,
                                const uint iel) const
    {
        assert(ispace < spaceN());
        auto it = mLocalElemIMap.find(std::make_pair(ispace, iel));
        if(it != mLocalElemIMap.end())
            return it->second;
        else {
            uint el_count = 0;
            for(int s = 0; s < ispace; ++s)
                el_count += space(s, S).nonzeroKnotSpanN();
            assert(iel < space(ispace,S).nonzeroKnotSpanN());
            const uint gel_index = el_count + iel;
            mLocalElemIMap[std::make_pair(ispace, iel)] = gel_index;
            return gel_index;
        }
    }
    
    uint MultiForest::globalElI(const uint ispace,
                                const uint i,
                                const uint j) const
    {
        const uint ns = space(ispace,S).uniqueKnotN(S) -1;
        return globalElI(ispace, ns * i + j );
    }
    
    VAnalysisElement* MultiForest::element(const uint i)
    {
        initEl(i);
        return mElems[i].get();
    }
    
    VAnalysisElement* MultiForest::bezierElement(const uint i)
    {
        initBezierEl(i);
        return mBezierElems[i].get();
    }
    
    void MultiForest::initEl(const uint i) const
    {
        auto e = mElems.find(i);
        if(e != mElems.end())
            return;
        else {
            uint start_i = 0;
            for(uint s = 0; s < spaceN(); ++s) {
                const uint el_n = space(s, S).nonzeroKnotSpanN(); // use S space
                if((i - start_i) > (el_n - 1)) {
                    start_i += el_n;
                    continue;
                }
                const uint local_i = i - start_i;
                auto r = mElems.insert(std::make_pair(i, make_unique<HCElement>(this, s, local_i)));
                if(!r.second) {
                    std::cout << i << "\n";
                    error("Failed attempt to create element");
                }
                return;
            }
        }
    }
    
    void MultiForest::initBezierEl(const uint i) const
    {
        auto e = mBezierElems.find(i);
        if(e != mBezierElems.end())
            return;
        else {
            uint start_i = 0;
            for(uint s = 0; s < spaceN(); ++s) {
                const uint el_n = space(s, S).nonzeroKnotSpanN(); // use S space
                if((i - start_i) > (el_n - 1)) {
                    start_i += el_n;
                    continue;
                }
                const uint local_i = i - start_i;
                auto r = mBezierElems.insert(std::make_pair(i, make_unique<BezierVectorElement>(this, s, local_i)));
                if(!r.second) {
                    std::cout << i << "\n";
                    error("Failed attempt to create element");
                }
                
                auto el = mBezierElems[i].get();
                
                // finally, create a reference to the parent element in the primal forest
                
                // search for parent element
                for(uint ielem = 0; ielem < geometry()->primalForest().space(s).nonzeroKnotSpanN(); ++ielem) {
                    const auto pel = geometry()->primalForest().bezierElement(s, ielem);
                    if(pel->contains(*el))
                        el->setParent(pel);
                }
                assert(el->parent() != nullptr);
                return;
            }
        }
    }
    
    VAnalysisElement* MultiForest::element(const uint ispace,
                                           const int i,
                                           const int j)
    {
        const uint ns = space(ispace,S).uniqueKnotN(S) - 1;
        if(j < 0 || j >= ns)
            return nullptr;
        const uint nt = space(ispace,T).uniqueKnotN(T) - 1;
        if(i < 0 || i >= nt)
            return nullptr;
        return element(globalElI(ispace, static_cast<uint>(i), static_cast<uint>(j)));
    }
    
    VAnalysisElement* MultiForest::bezierElement(const uint ispace,
                                                 const int i,
                                                 const int j)
    {
        const uint ns = space(ispace,S).uniqueKnotN(S) - 1;
        if(j < 0 || j >= ns)
            return nullptr;
        const uint nt = space(ispace,T).uniqueKnotN(T) - 1;
        if(i < 0 || i >= nt)
            return nullptr;
        return bezierElement(globalElI(ispace, static_cast<uint>(i), static_cast<uint>(j)));
    }
    
    void MultiForest::print(std::ostream& ost) const
    {
        ost << "Printing MultiForest data....\n\n";
        printImpl(ost);
    }

    
    UIntVec localBasisIVec(const Edge e,
                           const uint ns,
                           const uint nt)
    {
        std::vector<uint> l_vec;
        switch(e) {
            case Edge::EDGE2:
                for(uint c = 0; c < nt; ++c)
                    l_vec.push_back(c * ns);
                break;
            case Edge::EDGE3:
                for(uint c = 0; c < nt; ++c)
                    l_vec.push_back(c * ns + (ns - 1));
                break;
            case Edge::EDGE0:
                for(uint c = 0; c < ns; ++c)
                    l_vec.push_back(c);
                break;
            case Edge::EDGE1:
                for(uint c = 0; c < ns; ++c)
                    l_vec.push_back(c + (nt - 1) * ns);
                break;
        }
        return l_vec;
    }
    
    uint localBasisI(const Vertex v,
                     const uint ns,
                     const uint nt)
    {
        switch(v) {
            case Vertex::VERTEX0: return 0;
            case Vertex::VERTEX1: return ns - 1 ;
            case Vertex::VERTEX2: return (nt - 1) * ns;
            case Vertex::VERTEX3: return ns * nt - 1;
        }
    }
    
    UIntVec localBasisIVec(const Edge e, const BSplineSpace& s)
    { return localBasisIVec(e, s.basisFuncN(S), s.basisFuncN(T)); }
    
    uint localBasisI(const Vertex v,
                     const BSplineSpace& s)
    { return localBasisI(v, s.basisFuncN(S), s.basisFuncN(T)); }
    
    std::pair<uint, uint> localBasisIPair(const Edge e,
                                          const BSplineSpace& s)
    {
        switch (e) {
            case Edge::EDGE0:
                return std::make_pair(localBasisI(Vertex::VERTEX0, s),
                                      localBasisI(Vertex::VERTEX1, s));
                break;
            case Edge::EDGE1:
                return std::make_pair(localBasisI(Vertex::VERTEX2, s),
                                      localBasisI(Vertex::VERTEX3, s));
                break;
            case Edge::EDGE2:
                return std::make_pair(localBasisI(Vertex::VERTEX0, s),
                                      localBasisI(Vertex::VERTEX2, s));
                break;
            case Edge::EDGE3:
                return std::make_pair(localBasisI(Vertex::VERTEX1, s),
                                      localBasisI(Vertex::VERTEX3, s));
                break;
        }
    }
    
    UIntVec localBasisIVec(const Edge e,
                           const ParamDir d,
                           const std::pair<BSplineSpace, BSplineSpace>& spacepair)
    {
        const auto& thisspace = (ParamDir::S == d) ? spacepair.first : spacepair.second;
        
        auto l_vec = localBasisIVec(e, thisspace);
        if(ParamDir::T == d)
        {
            const uint n_sbasis = spacepair.first.basisFuncN();
            std::for_each(l_vec.begin(), l_vec.end(), [&n_sbasis](uint& v){ v += n_sbasis; });
        }
        return l_vec;
    }
    
    uint localBasisI(const Vertex v,
                     const ParamDir d,
                     const std::pair<BSplineSpace, BSplineSpace>& spacepair)
    {
        const auto& thisspace = (ParamDir::S == d) ? spacepair.first : spacepair.second;
        
        uint l_i = localBasisI(v, thisspace);
        if(ParamDir::T == d)
            l_i += spacepair.first.basisFuncN();
        return l_i;
    }
    
    bool interiorLocalI(const BSplineSpace& space,
                        const uint lindex)
    {
        assert(lindex < space.basisFuncN());
        const uint ns = space.basisFuncN(S);
        const uint nt = space.basisFuncN(T);
        
        // Are we on edge 2 or 3?
        if((lindex % ns) == 0 || (lindex % ns) == ns-1)
            return false;
        // Are we on edge 0 or 1?
        const uint tindex = lindex / ns;
        if((tindex % nt) == 0 || (tindex % nt) == nt-1)
            return false;
        
        // Otherwise we must be at an interior point
        return true;
    }

	std::pair<Sign, ParamDir> vectorBasisEdgeDir(const Edge edge_ref,
												 const ParamDir d_ref,
												 const Edge e_in)
	{
		switch(edge_ref) {
			case Edge::EDGE0:
				switch(e_in) {
					case Edge::EDGE0:
						if(S == d_ref)
							return std::make_pair(Sign::NEGATIVE, S);
						else
                            return std::make_pair(Sign::NEGATIVE, T);
						break;
					case Edge::EDGE1:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
						break;
					case Edge::EDGE2:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
						break;
					case Edge::EDGE3:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
						break;
				}
				break;
			case Edge::EDGE1:
				switch(e_in) {
					case Edge::EDGE0:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
						break;
					case Edge::EDGE1:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, S);
                        else
                            return std::make_pair(Sign::NEGATIVE, T);
						break;
					case Edge::EDGE2:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
						break;
					case Edge::EDGE3:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
						break;
				}
				break;
			case Edge::EDGE2:
				switch(e_in) {
					case Edge::EDGE0:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
						break;
					case Edge::EDGE1:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
						break;
					case Edge::EDGE2:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, S);
                        else
                            return std::make_pair(Sign::NEGATIVE, T);
						break;
					case Edge::EDGE3:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
						break;
				}
				break;
			case Edge::EDGE3:
				switch(e_in) {
					case Edge::EDGE0:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
						break;
					case Edge::EDGE1:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
						break;
					case Edge::EDGE2:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
						break;
					case Edge::EDGE3:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, S);
                        else
                            return std::make_pair(Sign::NEGATIVE, T);
						break;
				}
				break;
		}
	}

	std::pair<Sign, ParamDir> vectorBasisVertexConn(const Vertex vert_ref,
													const ParamDir d_ref,
													const Vertex v_in)
	{
        switch(vert_ref) {
            case Vertex::VERTEX0:
                switch(v_in) {
                    case Vertex::VERTEX0:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, S);
                        else
                            return std::make_pair(Sign::NEGATIVE, T);
                        break;
                    case Vertex::VERTEX1:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
                        break;
                    case Vertex::VERTEX2:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
                        break;
                    case Vertex::VERTEX3:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
                        break;
                }
                break;
            case Vertex::VERTEX1:
                switch(v_in) {
                    case Vertex::VERTEX0:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
                        break;
                    case Vertex::VERTEX1:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, S);
                        else
                            return std::make_pair(Sign::NEGATIVE, T);
                        break;
                    case Vertex::VERTEX2:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
                        break;
                    case Vertex::VERTEX3:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
                        break;
                }
                break;
            case Vertex::VERTEX2:
                switch(v_in) {
                    case Vertex::VERTEX0:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
                        break;
                    case Vertex::VERTEX1:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
                        break;
                    case Vertex::VERTEX2:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, S);
                        else
                            return std::make_pair(Sign::NEGATIVE, T);
                        break;
                    case Vertex::VERTEX3:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
                        break;
                }
                break;
            case Vertex::VERTEX3:
                switch(v_in) {
                    case Vertex::VERTEX0:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, S);
                        else
                            return std::make_pair(Sign::POSITIVE, T);
                        break;
                    case Vertex::VERTEX1:
                        if(S == d_ref)
                            return std::make_pair(Sign::POSITIVE, T);
                        else
                            return std::make_pair(Sign::NEGATIVE, S);
                        break;
                    case Vertex::VERTEX2:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, T);
                        else
                            return std::make_pair(Sign::POSITIVE, S);
                        break;
                    case Vertex::VERTEX3:
                        if(S == d_ref)
                            return std::make_pair(Sign::NEGATIVE, S);
                        else
                            return std::make_pair(Sign::NEGATIVE, T);
                        break;
                }
                break;
        }
	}
    
    UIntVec localElIndices(const BSplineSpace& space,
                           const Edge e,
                           const Sign d)
    {
        const uint ns = space.uniqueKnotN(S) - 1;
        const uint nt = space.uniqueKnotN(T) - 1;
        const uint nel = space.nonzeroKnotSpanN();
        UIntVec rvec;
        switch (e) {
            case Edge::EDGE0:
                rvec.resize(ns);
                std::iota(rvec.begin(), rvec.end(), 0);
                break;
            case Edge::EDGE1:
                rvec.resize(ns);
                std::iota(rvec.begin(), rvec.end(), nel - ns);
                break;
            case Edge::EDGE2:
                rvec.resize(nt);
                for(uint j = 0; j < nt; ++j)
                    rvec[j] = j * ns;
                break;
            case Edge::EDGE3:
                rvec.resize(nt);
                for(uint j = 0; j < nt; ++j)
                    rvec[j] = (j+1) * ns - 1;
                break;
        }
        if(Sign::NEGATIVE == d)
            std::reverse(rvec.begin(), rvec.end());
        return rvec;
    }
    
    uint localElIndex(const BSplineSpace& space,
                      const Vertex v)
    {
        const uint ns = space.uniqueKnotN(S) - 1;
        const uint nt = space.uniqueKnotN(T) - 1;
        switch (v) {
            case Vertex::VERTEX0:
                return 0; break;
            case Vertex::VERTEX1:
                return ns - 1; break;
            case Vertex::VERTEX2:
                return ns * (nt - 1); break;
            case Vertex::VERTEX3:
                return ns * nt - 1;
                break;
        }
    }
    
    UIntVecVec subEdgeIVecs(const BSplineSpace& space,
                            const Edge e)
    {
        // march along edge generating vectors of local
        // basis function indices that are broken by lines
        // of C^0 continuity
        UIntVecVec C0_edgevec; // the vector we're populating
        const ParamDir dir = paramDir(e, ContinuityType::TANGENT);
        uint start_i = localBasisIVec(e, space).front();
        const uint degree = space.degree(dir);
        const uint unique_knot_n = space.uniqueKnotN(dir);
        uint knot_count = 0;
        for(uint k = 0; k < unique_knot_n; ++k) {
            const uint knot_repeats = space.knotRepeats(k, dir);
            knot_count += knot_repeats;
            
            if(knot_repeats >= degree && k != 0) {
                std::vector<uint> subedge_vec;
                if(unique_knot_n - 1 == k)
                    subedge_vec.resize(knot_count - degree - 1);
                else
                    subedge_vec.resize(knot_count - degree);
                if(Edge::EDGE0 == e || Edge::EDGE1 == e) {
                    std::iota(subedge_vec.begin(), subedge_vec.end(), start_i);
                    start_i += subedge_vec.size() - 1;
                }
                else {
                    uint index = start_i;
                    const uint nbasis_s = space.basisFuncN(S);
                    for(auto& i : subedge_vec) {
                        i = index;
                        index += nbasis_s;
                    }
                    start_i += (subedge_vec.size() - 1) * nbasis_s;
                }
                C0_edgevec.push_back(subedge_vec);
                knot_count = degree + 1; // simulate C^-1 continuity
            }
        }
        return C0_edgevec;
    }
    
    Sign edgeOrientation(Edge e1, Edge e2)
    {
        std::bitset<1> b1, b2; // value of 1 = clockwise, 0 anti-clockwise
        if(Edge::EDGE1 == e1 || Edge::EDGE2 == e1)
            b1.set(0,true);
        if(Edge::EDGE1 == e2 || Edge::EDGE2 == e2)
            b2.set(0,true);
        auto r = b1 ^ b2;
        if(r.test(0))
            return Sign::POSITIVE;
        else return Sign::NEGATIVE;
    }
    
    ParamDir paramDir(const Edge e, const ContinuityType c)
    {
        switch (c) {
            case ContinuityType::NORMAL:
                switch (e) {
                    case Edge::EDGE0:
                        return T;
                        break;
                    case Edge::EDGE1:
                        return T;
                        break;
                    case Edge::EDGE2:
                        return S;
                        break;
                    case Edge::EDGE3:
                        return S;
                        break;
                }
                break;
            case ContinuityType::TANGENT:
                switch (e) {
                    case Edge::EDGE0:
                        return S;
                        break;
                    case Edge::EDGE1:
                        return S;
                        break;
                    case Edge::EDGE2:
                        return T;
                        break;
                    case Edge::EDGE3:
                        return T;
                        break;
                }
                break;
        }

    }
    
    void MultiForest::printImpl(std::ostream& ost) const
    {
        for(uint ispace = 0; ispace < spaceN(); ++ispace)
        {
            ost << "Space index: " << ispace << "\n";
            ost << "vector basis connectivity: ";
            for(uint ibasis = 0; ibasis < basisFuncN(ispace); ++ibasis) {
                ost << globalDir(ispace, ibasis) << globalI(ispace, ibasis) << " ";
            }
            ost << "\n";
        }
    }
    
    GPt2D projectPt(const GPt2D p, const Edge e1, const Edge e2)
    {
        switch (e1) {
            case Edge::EDGE0:
                switch (e2) {
                    case Edge::EDGE0:
                        return GPt2D(-p.s, -1.0);
                        break;
                    case Edge::EDGE1:
                        return GPt2D(p.s, 1.0);
                        break;
                    case Edge::EDGE2:
                        return GPt2D(-1.0, p.s);
                        break;
                    case Edge::EDGE3:
                        return GPt2D(1.0, -p.s);
                        break;
                    default:
                        throw std::runtime_error("Bad edge type");
                        break;
                }
                break;
            case Edge::EDGE1:
                switch (e2) {
                    case Edge::EDGE0:
                        return GPt2D(p.s, -1.0);
                        break;
                    case Edge::EDGE1:
                        return GPt2D(-p.s, 1.0);
                        break;
                    case Edge::EDGE2:
                        return GPt2D(-1.0, -p.s);
                        break;
                    case Edge::EDGE3:
                        return GPt2D(1.0, p.s);
                        break;
                    default:
                        throw std::runtime_error("Bad edge type");
                        break;
                }
                break;
                
            case Edge::EDGE2:
                switch (e2) {
                    case Edge::EDGE0:
                        return GPt2D(p.t, -1.0);
                        break;
                    case Edge::EDGE1:
                        return GPt2D(-p.t, 1.0);
                        break;
                    case Edge::EDGE2:
                        return GPt2D(-1.0, -p.t);
                        break;
                    case Edge::EDGE3:
                        return GPt2D(1.0, p.t);
                        break;
                    default:
                        throw std::runtime_error("Bad edge type");
                        break;
                }
                break;
                
            case Edge::EDGE3:
                switch (e2) {
                    case Edge::EDGE0:
                        return GPt2D(-p.t, -1.0);
                        break;
                    case Edge::EDGE1:
                        return GPt2D(p.t, 1.0);
                        break;
                    case Edge::EDGE2:
                        return GPt2D(-1.0, p.t);
                        break;
                    case Edge::EDGE3:
                        return GPt2D(1.0, -p.t);
                        break;
                    default:
                        throw std::runtime_error("Bad edge type");
                        break;
                }
                break;
                
            default:
                throw std::runtime_error("Bad edge type");
                break;
        }
    }
    
    GPt2D paramPt(const Vertex v)
    {
        switch (v) {
            case Vertex::VERTEX0:
                return GPt2D(-1.0, -1.0);
                break;
            case Vertex::VERTEX1:
                return GPt2D(1.0, -1.0);
                break;
            case Vertex::VERTEX2:
                return GPt2D(-1.0, 1.0);
                break;
            case Vertex::VERTEX3:
                return GPt2D(1.0, 1.0);
                break;
            default:
                throw std::runtime_error("Bad vertex type");
                break;
        }
    }
}

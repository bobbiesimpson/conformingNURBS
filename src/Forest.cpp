#include <ios>
#include <string>
#include <utility>
#include <algorithm>

#include "Forest.h"
#include "BSplineSpace.h"
#include "Geometry.h"
#include "NodalElement.h"
#include "MultiForest.h"
#include "NodalElement.h"
#include "BezierNodalElement.h"

namespace nurbs
{
    
    Forest::Forest(const Geometry& g) : mGeom(&g)
    {
        const Forest& f = g.primalForest();
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mGlobalDofN = f.mGlobalDofN;
        initEdgeConn();
        initCollocConn();
    }
    
    Forest::Forest(const Forest& f)
    {
        clear();
        mGeom = f.mGeom;
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mEdgeConn = f.mEdgeConn;
        mEdgeSpaceMap = f.mEdgeSpaceMap;
        mCVertexSpaceMap = f.mCVertexSpaceMap;
        mSubEdgeMap = f.mSubEdgeMap;
        mSuperEdgeMap = f.mSuperEdgeMap;
        mGlobalCollocConn = f.mGlobalCollocConn;
        mLocalCollocConn = f.mLocalCollocConn;
        mElemIndexMap = f.mElemIndexMap;
        for(const auto& e : f.mElems)
            mElems.insert(std::make_pair(e.first, e.second->copy()));
        for(const auto& e : f.mBezierElems)
            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        mGlobalDofN = f.mGlobalDofN;
    }
    
    /// Assignment operator
    Forest& Forest::operator=(const Forest& f)
    {
        if(this == &f) // check for self assignment
            return *this;
        clear();
        mGeom = f.mGeom;
        mSpaces = f.mSpaces;
        mSpaceMap = f.mSpaceMap;
        mNodalConn = f.mNodalConn;
        mEdgeConn = f.mEdgeConn;
        mEdgeSpaceMap = f.mEdgeSpaceMap;
        mCVertexSpaceMap = f.mCVertexSpaceMap;
        mSubEdgeMap = f.mSubEdgeMap;
        mSuperEdgeMap = f.mSuperEdgeMap;
        mGlobalCollocConn = f.mGlobalCollocConn;
        mLocalCollocConn = f.mLocalCollocConn;
        mElemIndexMap = f.mElemIndexMap;
        for(const auto& e : f.mElems)
            mElems.insert(std::make_pair(e.first, e.second->copy()));
        for(const auto& e : f.mBezierElems)
            mBezierElems.insert(std::make_pair(e.first, e.second->copy()));
        mGlobalDofN = f.mGlobalDofN;
        return *this;
    }
    
    /// clear all data
    void Forest::clear()
    {
        mSpaces.clear();
        mSpaceMap.clear();
        mNodalConn.clear();
        mEdgeConn.clear();
        mEdgeSpaceMap.clear();
        mCVertexSpaceMap.clear();
        mSubEdgeMap.clear();
        mSuperEdgeMap.clear();
        mGlobalCollocConn.clear();
        mLocalCollocConn.clear();
        mElemIndexMap.clear();
        mElems.clear();
        mBezierElems.clear();
        mGlobalDofN = std::make_pair(false, 0);
    }
    
	void Forest::load(std::istream& ist)
	{
		while(true) {
			BSplineSpace s;
			if(!(ist >> s))
				break;
            std::cout << "found space\n";
			addSpace(s);
		}
		endOfLoop(ist, '}', "Could not read B-spline space\n");

		while(true) {
			ConnVecInput cv;
			if(!(ist >> cv)) 
				break;
			addConnectivityVec(cv.name, cv.data);
		}
		endOfLoop(ist, '}', "Could not read connectivity vectors\n");
        initEdgeConn();
        initCollocConn();
	}

	void Forest::print(std::ostream& ost) const
	{
		ost << "Forest: " << spaceN() << " trunks\n";
		for(const auto& n : mSpaceMap) {
			ost << mSpaces[n.second] << "\n"; 
		}
        ost << "Nodal connectivity:\n";
        for(const auto& c : mNodalConn)
            ost << c.first << ": " << c.second << "\n";
	}
    
    const NAnalysisElement* Forest::element(const uint i) const
    {
        auto e = mElems.find(i);
        if(e != mElems.end())
            return e->second.get();
        uint start_i = 0;
        for(uint s = 0; s < spaceN(); ++s) {
            const uint el_n = space(s).nonzeroKnotSpanN();
            if((i - start_i) > (el_n - 1)) {
                start_i += el_n;
                continue;
            }
            const uint local_i = i - start_i;
            auto r = mElems.insert(std::make_pair(i, make_unique<NodalElement>(this, s, local_i)));
            if(!r.second)
                error("Failed attempt to create element");
            auto r_i = mElemIndexMap.insert(std::make_pair(i, std::make_pair(s, local_i)));
            if(!r_i.second)
                error("Failure inserting element index mapping");
            auto el = mElems[i].get();
            
            // search for parent
            for(uint ielem = 0; ielem < geometry()->primalForest().space(s).nonzeroKnotSpanN(); ++ielem) {
                const auto pel = geometry()->primalForest().element(s, ielem);
                if(pel->contains(*el))
                    el->setParent(pel);
            }
            assert(el->parent() != nullptr);
            return el;
        }
        error("Failed to create element"); return nullptr;
    }
    
    const NAnalysisElement* Forest::bezierElement(const uint i) const
    {
        auto e = mBezierElems.find(i);
        if(e != mBezierElems.end())
            return e->second.get();
        uint start_i = 0;
        for(uint s = 0; s < spaceN(); ++s) {
            const uint el_n = space(s).nonzeroKnotSpanN();
            if((i - start_i) > (el_n - 1)) {
                start_i += el_n;
                continue;
            }
            const uint local_i = i - start_i;
            
            // create the element
            auto r = mBezierElems.insert(std::make_pair(i, make_unique<BezierNodalElement>(this,
                                                                                           s,
                                                                                           local_i)));
            if(!r.second)
                error("Failed attempt to create element");
            
            // create mapping from global element index to space and local element index
            mElemIndexMap.insert(std::make_pair(i, std::make_pair(s, local_i)));
//            if(!r_i.second)
//                std::cout << "Failure inserting bezier element index mapping. Carrying on...\n";
            
            auto el = mBezierElems[i].get();
            
            // finally, create a referece to the parent element in the primal forest
            
            // search for parent element
            for(uint ielem = 0; ielem < geometry()->primalForest().space(s).nonzeroKnotSpanN(); ++ielem) {
                const auto pel = geometry()->primalForest().bezierElement(s, ielem);
                if(pel->contains(*el))
                    el->setParent(pel);
            }
            assert(el->parent() != nullptr);
            return el;
        }
        error("Failed to create Bezier element"); return nullptr;
    }
    
    bool Forest::degenerateEdge(const uint ispace,
                                const uint iedge) const
    {
        // first get the set of nodal indices along the given edge
        const auto localvec  = localBasisIVec(edgeType(iedge), space(ispace));
        auto gvec = globalIVec(ispace, localvec);
        std::sort(gvec.begin(), gvec.end());
        auto last = std::unique(gvec.begin(), gvec.end());
        gvec.erase(last, gvec.end());
        
        // if there is only one unique nodal dof, it must be a degenerate edge
        if(gvec.size() == 1)
        {
//            std::cout << "Degenerate node: " << gvec[0] << "\n";
            return true;
        }
        else
            return false;
    }
    
    uint Forest::globalVertexI(const uint ispace, const Vertex v) const
    {
        const BSplineSpace& s = space(ispace);
        const uint nbasis_s = s.basisFuncN(S);
        const uint nbasis_t = s.basisFuncN(T);
        switch(v) {
           case Vertex::VERTEX0: return globalI(ispace,0,0);
           case Vertex::VERTEX1: return globalI(ispace, nbasis_s - 1,0);
           case Vertex::VERTEX2: return globalI(ispace,0, nbasis_t - 1);
           case Vertex::VERTEX3: return globalI(ispace, nbasis_s - 1,nbasis_t - 1);
        }
    }
    
    std::pair<uint, uint> Forest::globalVertexPairI(const uint ispace, const Edge e) const
    {
        switch(e) {
            case Edge::EDGE0:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX0),
                                      globalVertexI(ispace, Vertex::VERTEX1));
            case Edge::EDGE1:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX2),
                                      globalVertexI(ispace, Vertex::VERTEX3));
            case Edge::EDGE2:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX0),
                                      globalVertexI(ispace, Vertex::VERTEX2));
            case Edge::EDGE3:
                return std::make_pair(globalVertexI(ispace, Vertex::VERTEX1),
                                      globalVertexI(ispace, Vertex::VERTEX3));
        }
    }
    
    Point3D Forest::collocPt(const uint sp, const uint i) const
    {
        const GPt2D param_pt = space(sp).grevilleAbscissaPt(i);
        return geometry()->eval(param_pt.s, param_pt.t, sp);
    }
    
    std::pair<bool, uint> Forest::connectedCollocPtI(const uint icspace,
                                                     const uint icpt,
                                                     const uint ifspace) const
    {
        // if space indices are indentical, we simply return the given local
        // index
        if(icspace == ifspace)
            return std::make_pair(true, icpt);
        
        // If the given local point is an interior point (to icspace)
        // then it will not be connected to ifspace.
        if(interiorLocalI(space(icspace), icpt))
            return std::make_pair(false, INVALID_UINT);
        
        // Now search for connected space vertices
        Vertex cspace_v, fspace_v;
        bool v_connected = connectedVertices(icspace, ifspace, cspace_v, fspace_v);
        const uint iglobal = globalCollocI(icspace, icpt);
        if(v_connected) { // loop through vertices and return if global index matches. Carry on to edge search if not v. connected
            for(uint v = 0; v < NVERTICES; ++v) {
                const uint ilocal = localBasisI(vertexType(v), space(ifspace));
                if(iglobal == globalCollocI(ifspace, ilocal))
                    return std::make_pair(true, ilocal);
            }
        }
        
        // And now search for connected edges
        const auto conn_pair = multiConnectedEdges(icspace, ifspace);
        if(conn_pair.first) {
            for(const auto& edge_pair : conn_pair.second) { // loop over edge pairs that are connected
                for(const auto& ilocal : localBasisIVec(edge_pair.second, space(ifspace)))
                    if(iglobal == globalCollocI(ifspace, ilocal))
                        return std::make_pair(true, ilocal);
            }
        }
//        Edge cspace_e, fspace_e;
//        bool e_connected = connectedEdges(icspace, ifspace, cspace_e, fspace_e);
//        if(e_connected) {
//            
//            for(const auto& ilocal : localBasisIVec(fspace_e, space(ifspace)))
//                if(iglobal == globalCollocI(ifspace, ilocal))
//                    return std::make_pair(true, ilocal);
//        }
        
        // If we reach here, the spaces are not connected.
        return std::make_pair(false, INVALID_UINT);
    }
    
    std::pair<bool, std::map<Edge, Edge>> Forest::multiConnectedEdges(const uint ispace1,
                                                                      const uint ispace2) const
    {
        std::map<Edge, Edge> emap;
        bool connected_flag = false;
        for(uint e1 = 0; e1 < NEDGES; ++e1) {
            Edge edge1 = edgeType(e1);
            const uint g1 = globalEdgeI(ispace1, e1);
            for(uint e2 = 0; e2 < NEDGES; ++e2) {
                Edge edge2 = edgeType(e2);
                if(g1 == globalEdgeI(ispace2, e2)) {
                    connected_flag = true;
                    emap.insert(std::make_pair(edge1, edge2));
                }
            }
        }
        return std::make_pair(connected_flag, emap);
    }
    
    bool Forest::connectedEdges(const uint ispace1,
                                const uint ispace2,
                                Edge& edge1,
                                Edge& edge2) const
    {
        for(uint e1 = 0; e1 < NEDGES; ++e1) {
            edge1 = edgeType(e1);
            const uint g1 = globalEdgeI(ispace1, e1);
            for(uint e2 = 0; e2 < NEDGES; ++e2) {
                edge2 = edgeType(e2);
                if(g1 == globalEdgeI(ispace2, e2))
                    return true;
            }
        }
        return false;
    }
    
    bool Forest::connectedEdges(const uint ispace1,
                                const uint ispace2,
                                const uint iedge,
                                Edge& edge1,
                                Edge& edge2) const
    {
        bool e1_connected = false;
        for(uint e1 = 0; e1 < NEDGES; ++e1) {
            if(globalEdgeI(ispace1, e1) == iedge) {
                edge1 = edgeType(e1);
                e1_connected = true;
                break;
            }
        }
        bool e2_connected = false;
        for(uint e2 = 0; e2 < NEDGES; ++e2) {
            if(globalEdgeI(ispace2, e2) == iedge) {
                edge2 = edgeType(e2);
                e2_connected = true;
                break;
            }
        }
        return e1_connected && e2_connected;
    }
    
    bool Forest::connectedVertices(const uint ispace1,
                                   const uint ispace2,
                                   Vertex& vert1,
                                   Vertex& vert2) const
    {
        for(uint v1 = 0; v1 < NVERTICES; ++v1) {
            vert1 = vertexType(v1);
            const uint glb_v1 = globalVertexI(ispace1, v1);
            for(uint v2 = 0; v2 < NVERTICES; ++v2) {
                vert2 = vertexType(v2);
                if(glb_v1 == globalVertexI(ispace2, v2))
                    return true;
            }
        }
        return false;
    }
    
    Sign Forest::globalEdgeDir(const uint ispace, const Edge e) const
    {
        auto v_pair = globalVertexPairI(ispace, e);
        if(v_pair.second > v_pair.first)
            return Sign::POSITIVE; // increasing vertex index == positive edge
        else return Sign::NEGATIVE; // negative otherwise
    }
    
    void Forest::initNodalConn()
    {
        // Now construct the new connectivity using the edge connectivy that
        // is independent of refinement
        
        const Forest& primal = geometry()->primalForest();
        std::map<uint, std::vector<uint>> temp_conn;
        std::map<uint, std::vector<uint>> edge_map; // map from global edge index to global node indices
        std::map<uint, uint> vx_map; // map from geometry vertex index to new vertex index
        uint current_index = 0; // the current global node index
        
        for(uint ispace = 0; ispace < spaceN(); ++ispace) {
            const BSplineSpace& s = space(ispace);
            std::vector<int> gnode_vec(s.basisFuncN(), -1); // -1 means unassigned
            
            // Use vertex connectivity from primal space
            for(uint ivertex = 0; ivertex < NVERTICES; ++ivertex)  {
                const Vertex vertex = vertexType(ivertex);
                const uint geom_ivertex = primal.globalVertexI(ispace, vertex);
                auto find = vx_map.find(geom_ivertex);
                if(find != vx_map.end())
                    gnode_vec[localBasisI(vertex, s)] = find->second;
                else {
                    const uint new_index = current_index++;
                    vx_map[geom_ivertex] = new_index;
                    gnode_vec[localBasisI(vertex, s)] = new_index;
                }
            }
            
            // Now assign edges that have been assigned previously
            for(uint iedge = 0; iedge < NEDGES; ++iedge) {
                const Edge edge = edgeType(iedge);
                const uint global_iedge = globalEdgeI(ispace, iedge);
                auto find = edge_map.find(global_iedge);
                auto vpair = localBasisIPair(edge, s); // ordered local vertex indices on edge
                assert(vpair.first != vpair.second);
                const Sign sign = (gnode_vec[vpair.second] > gnode_vec[vpair.first]) ? Sign::POSITIVE : Sign::NEGATIVE;
                //const Sign sign = globalEdgeDir(ispace, edge);
                auto local_indices = localBasisIVec(edge, s);
                if(Sign::NEGATIVE == sign)
                    std::reverse(local_indices.begin(), local_indices.end());
                if(find != edge_map.end()) { // the collocation indices are assigned on this edge
                    const auto glb_indices = find->second;
                    assert(local_indices.size() == glb_indices.size());
                    for(uint i = 0; i < local_indices.size(); ++i) {
                        if(gnode_vec[local_indices[i]] != -1)
                            continue;
                        else
                            gnode_vec[local_indices[i]] = glb_indices[i];
                    }
                }
                else { // the edge node indices have not been assigned. Let's do it now.
                    std::vector<uint> gvec;
                    for(const auto& i : local_indices) {
                        if(gnode_vec[i] != -1)
                            gvec.push_back(gnode_vec[i]);   
                        else {
                            const uint gindex = current_index++;
                            gvec.push_back(gindex);
                            gnode_vec[i] = gindex;
                        }
                    }
//                    if(Sign::NEGATIVE == sign) // make sure to store indices in positive order
//                        std::reverse(gvec.begin(), gvec.end());
                    edge_map[global_iedge] = gvec;
                }
            }
            // now fill up unassigned indices
            std::vector<uint> unsigned_gvec;
            for(auto& i : gnode_vec) {
                if(i != -1) {
                    unsigned_gvec.push_back(static_cast<uint>(i));
                    continue;
                }
                i = current_index++;
                unsigned_gvec.push_back(i);
            }
            temp_conn[ispace] = unsigned_gvec;
        }
        mNodalConn = temp_conn;
    }
    
    void Forest::initEdgeConn()
    {
        // tidy up before initialising
        mEdgeConn.clear();
        mEdgeSpaceMap.clear();
        mCVertexSpaceMap.clear();
        mSubEdgeMap.clear();
        mSuperEdgeMap.clear();
        
        std::map<std::vector<uint>, uint> edge_map; // mapping from edge vertices to global edge index
        std::vector<Edge> e_types{Edge::EDGE0, Edge::EDGE1, Edge::EDGE2, Edge::EDGE3};
        
        uint edge_count = 0;
        for(uint s = 0; s < spaceN(); ++s)
        {
            std::vector<uint> edge_conn; // edge connectivity for this space
            std::vector<uint> c_vertices;
            for(const auto& e : e_types)
            {
                auto v_pair = globalVertexPairI(s, e);
                auto e_conn = globalIVec(s, localBasisIVec(e, space(s)));
                std::sort(e_conn.begin(), e_conn.end());

                c_vertices.push_back(v_pair.first);
                c_vertices.push_back(v_pair.second);
                
                if(v_pair.first > v_pair.second) // make sure ordering is +ve
                    v_pair = std::make_pair(v_pair.second, v_pair.first);
                
                auto search = edge_map.find(e_conn);
                
                if(search != edge_map.end())
                    edge_conn.push_back(search->second);
                else
                {
                    edge_map[e_conn] = edge_count;
                    edge_conn.push_back(edge_count);
                    ++edge_count;
                }
            }
            std::sort(c_vertices.begin(), c_vertices.end());
            auto last = std::unique(c_vertices.begin(), c_vertices.end());
            c_vertices.erase(last, c_vertices.end());
            for(const auto& c : c_vertices)
                mCVertexSpaceMap[c].push_back(s);
            for(const auto& e : edge_conn)
                mEdgeSpaceMap[e].push_back(s);
            mEdgeConn[s] = edge_conn;
        }
        
        // realised the the code below is redundant. We can simply use the
        // nodal connectivity matrix to obtain the connectivity of the
        // Greville abscissa
        
//        // Now construct the Greville abscissa connectivity.
//        std::map<uint, std::vector<uint>> cedge_map; // map from global edge index to global colloc pt indices along edge
//        uint current_index = 0; // the current global colloc. pt index
//        
//        for(uint ispace = 0; ispace < spaceN(); ++ispace) {
//            const BSplineSpace& s = space(ispace);
//            std::vector<int> colloc_ivec(s.grevilleAbscissaPtN(), -1); // -1 means unassigned
//            
//            // First assign edges that have been assigned previously
//            for(uint iedge = 0; iedge < NEDGES; ++iedge) {
//                const Edge edge = edgeType(iedge);
//                const uint global_iedge = globalEdgeI(ispace, iedge);
//                auto find = cedge_map.find(global_iedge);
//                auto local_indices = localBasisIVec(edge, s);
//                if(Sign::NEGATIVE == globalEdgeDir(ispace, edge)) // reverse local indices if -ve edge orientation
//                    std::reverse(local_indices.begin(), local_indices.end());
//                if(find != cedge_map.end()) { // the collocation indices are assigned on this edge
//                    const auto glb_indices = find->second;
//                    assert(local_indices.size() == glb_indices.size());
//                    for(uint i = 0; i < local_indices.size(); ++i)
//                        colloc_ivec[local_indices[i]] = glb_indices[i];
//                }
//                else { // the edge collocation indices have not been assigned. Let's do it now.
//                    std::vector<uint> gvec;
//                    for(const auto& i : local_indices) {
//                        const uint gindex = current_index++;
//                        gvec.push_back(gindex);
//                        colloc_ivec[i] = gindex;
//                    }
//                    cedge_map[global_iedge] = gvec;
//                }
//            }
//            // now fill up unassigned indices
//            for(auto& i : colloc_ivec) {
//                if(i != -1)
//                    continue;
//                i = current_index++;
//            }
//        }
//        mCollocPtN = current_index; // and finally set the number of collocation points
        
        // Now construct maps that detail super- and sub-edges.
        // These only arise in the case of T-junctions in a NURBS discretisations
        // where more than two spaces may be connected to an edge.
//        for(uint ispace = 0; ispace < spaceN(); ++ispace) {
//            for(uint iedge = 0; iedge < NEDGES; ++iedge) {
//                const uint global_iedge = globalEdgeI(ispace, iedge);
//                std::cout << "space: " << ispace << ", edge: " << iedge << " " << globalIVec(ispace, localBasisIVec(edgeType(iedge), space(ispace))) << "\n";
//                const auto space_vec = connectedSpacesOnEdge(global_iedge);
//                assert(space_vec.size() != 0);
//                if(space_vec.size() > 1) continue; // edges already assigned
//                for(const auto& ebasis_ivec : subEdgeIVecs(space(ispace), edgeType(iedge))) {
//                   // std::cout << ebasis_ivec << "\n";
//                    auto glb_edge_ivec = globalIVec(ispace, ebasis_ivec);
//                    std::sort(glb_edge_ivec.begin(), glb_edge_ivec.end());
//                   // std::cout << glb_edge_ivec << "\n";
//                    auto search = edge_map.find(glb_edge_ivec);
//                    if(search != edge_map.end()) {
//                        if(search->second == global_iedge) continue;
//                        mSubEdgeMap[global_iedge].push_back(search->second);
//                        mSuperEdgeMap[search->second] = global_iedge;
//                    }
//                }
//            }
//        }
    }
    
    void Forest::initCollocConn()
    {
        
        // Loop through spaces and determine which collocation points lie in each
        // element and assign the indices of these points to the relevant entry
        // in the connectivity array
        
        // The aim is to construct the following maps...
        std::map<uint, std::vector<uint>> gconn; // global connectivity
        std::map<uint, std::vector<uint>> lconn; // local connectivity
        
        for(uint ispace = 0; ispace < spaceN(); ++ispace) {
            const auto sp = space(ispace);
            
            std::vector<ParamDir> paramvec{ParamDir::S, ParamDir::T};
            std::map<ParamDir,std::vector<std::vector<uint>>> param_conn;   // colloc. connectivity along each
                                                                            //parametric direction for this space
            for(const auto& dir : paramvec) {
                const auto gvec = sp.grevilleAbscissa(dir);
                const auto uknots = sp.uniqueKnotVec(dir);
                const uint nel = uknots.size() - 1;
                
                auto& conn = param_conn[dir];
                conn.resize(nel);

                uint iel = 0; // current local element index
                uint igrev = 0; // current greville abscissa index
                double lower = uknots[iel];
                double upper = uknots[iel + 1];
                
                while(true) {
                    const double currentgpt = gvec[igrev];
                    if(approximatelyEqual(currentgpt, upper, 1.e-9)) {
                        conn[iel].push_back(igrev);
                        ++iel; // move to next element
                        if(iel > nel - 1)
                            break; // we're done
                        lower = uknots[iel];
                        upper = uknots[iel + 1];
                    }
                    else if(currentgpt < upper) {
                        conn[iel].push_back(igrev);
                        ++igrev; // move to next greville point
                    }
                    else { // point must lie within next element
                        ++iel;
                        if(iel > nel - 1)
                            error("Bad element data while constructing Forest collocation connectivity");
                        lower = uknots[iel];
                        upper = uknots[iel + 1];
                    }
                }

            }
            // Now determine the global and local connectivity arrays
            // using the connectivity along each parametric direction
            for(uint j = 0; j < param_conn[ParamDir::T].size(); ++j) {
                const auto& jconn = param_conn[ParamDir::T][j];
                for(uint i = 0; i < param_conn[ParamDir::S].size(); ++i) {
                    const auto& iconn = param_conn[ParamDir::S][i];
                    const uint ielem_g = globalElI(ispace, i, j);
                    for(const auto& jcpt : jconn) {
                        for(const auto& icpt : iconn) {
                            const uint space_index = sp.globalGrevillePtI(icpt, jcpt);
                            const uint global_index = globalCollocI(ispace, space_index);
                            gconn[ielem_g].push_back(global_index);
                            lconn[ielem_g].push_back(space_index);
                        }
                    }
                }
            }

        }
//        for(const auto& l : lconn)
//            std::cout << l.first << "\t" <<l.second << "\n";
//        for(const auto& l : gconn)
//            std::cout << l.first << "\t" <<l.second << "\n";
        
        // and finally assign global collocation connectivity
        mGlobalCollocConn = gconn;
        mLocalCollocConn = lconn;
    }
    
    std::istream& operator>>(std::istream& ist, Forest& f)
    {
        f.load(ist);
        return ist;
    }
    
    std::ostream& operator<<(std::ostream& ost, const Forest& f)
    {
        f.print(ost);
        return ost;
    }
	
}

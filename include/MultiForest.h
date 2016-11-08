#ifndef NURBS_MULTI_FOREST_H
#define NURBS_MULTI_FOREST_H

#include "Forest.h"
#include "base.h"
#include "AnalysisElement.h"

#include <tuple>
#include <utility>
#include <cassert>
#include <numeric>

namespace nurbs {
   
    /// A multiforest is the interface for div- and curl-conforming elements
    /// It encapsulates the relevant degree reduction algorithms for constructing
    /// the correct spaces.
    
    /// The multiforest class is abstract and must be subclassed accordingly.
    /// Current plans are only to subclass to create H(div) and H(curl) spaces
    /// for emag BEM analysis. i.e. implementation of the EFIE and MFIE.
    
    class MultiForest {
        
    public:
        
        /// Destructor
        virtual ~MultiForest() {}
        
        /// Geometry getter
        const Geometry* geometry() const { return mpGeom; }
        
        /// Const B-spline space accessor
        const BSplineSpace& space(const uint ispace, const ParamDir dir) const
        {
            assert(ispace < spaceN());
            if(S == dir)
                return mSpaceS.at(ispace);
            else
                return mSpaceT.at(ispace);
        }
        
        /// Total number of elements in forest
        uint elemN() const;
        
        /// Number of B-spline spaces. S- and T- spaces should
        /// be equal in size (if constructed correctly)
        inline uint spaceN() const
        {
            assert(mSpaceS.size() == mSpaceT.size());
            return mSpaceS.size();
        }
        
        /// Vector Element getter
        const VAnalysisElement* element(const uint i) const;
        
        /// Bezier element getter
        const VAnalysisElement* bezierElement(const uint i) const;
        
        /// Non-const element getter
        VAnalysisElement* element(const uint i);
        
        /// Non-const bezier element getter
        VAnalysisElement* bezierElement(const uint i);
        
        /// Get element given space index and local element index within this space
        const VAnalysisElement* element(const uint ispace,
                                        const uint ielem) const
        {
            return element(globalElI(ispace, ielem));
        }
        
        /// Get a bezierelement given space index and
        /// a local element index within this space
        const VAnalysisElement* bezierElement(const uint ispace,
                                              const uint ielem) const
        {
            return bezierElement(globalElI(ispace, ielem));
        }
        
        const VAnalysisElement* nedelecElement(const uint i) const
        {
            return mNedelecElems[i].get();
        }
        
        /// Get the number of elements on the given space
        /// The number of elements in each space is equal,
        /// therefore use S parametric space as default.
        const uint elemN(const uint ispace) const
        {
            return space(ispace, ParamDir::S).nonzeroKnotSpanN();
        }
        
        /// Apply prefinement
        void prefine(const uint nrefine);
        
        /// Apply uniform h-refinment to the multiforest
        void hrefine(const uint nrefine = 1);
        
        /// Apply graded h-refinement
//        void hrefineGraded(const uint nrefine,
//                           const double sratio,
//                           const double tratio);
        
        /// Get global basis function index
        uint globalI(const uint ispace, const uint ibasis, const ParamDir dir) const
        {
            if(ParamDir::T == dir)
                return mConn.at(ispace)[space(ispace, S).basisFuncN() + ibasis];
            else return mConn.at(ispace)[ibasis];
        }
        
        /// Get global non-degenerate global dof index for given space, basis and parameteric direction
        int signedGlobalI(const uint ispace,
                                 const uint ibasis,
                                 const ParamDir dir) const
        {
            if(ParamDir::T == dir)
                return mSignedConn.at(ispace)[space(ispace, S).basisFuncN() + ibasis];
            else return mSignedConn.at(ispace)[ibasis];
        }
        
        /// Get the global basis as above without specifing parametric
        /// direction. Assumed ordering of connectivity is S-direction indices
        /// the T-direction indices.
        uint globalI(const uint ispace,
                     const uint ibasis) const
        { return mConn.at(ispace)[ibasis]; }
        
        /// Get non-dengereate global index for given space and local basis index
        int signedGlobalI(const uint ispace,
                          const uint ibasis) const
        { return mSignedConn.at(ispace)[ibasis]; }
        
        
        /// Return the connectivity vector for this space
        const UIntVec& globalIVec(const uint ispace) const
        { return mConn.at(ispace); }
        
        /// Get the non-degenerate connectivity vector for this space
        const IntVec& signedGlobalIVec(const uint ispace) const
        { return mSignedConn.at(ispace); }

        /// Get the direction of the basis function in the given space and
        /// local basis function index
        Sign globalDir(const uint ispace, const uint ibasis) const
        { return mSignConn.at(ispace)[ibasis]; }
        
        /// Return the vector of basis function directions for this space
        const std::vector<Sign>& globalDirVec(const uint ispace) const
        { return mSignConn.at(ispace); }
        
        /// Get global basis index given local indices in each parametric direction
        uint globalI(const uint ispace, const uint i, const uint j, const ParamDir dir) const
        {
            const uint index = j * space(ispace, dir).basisFuncN(S) + i;
            return globalI(ispace, index, dir);
        }
        
        /// Get global basis index given local indices in each parametric direction
        int signedGlobalI(const uint ispace,
                          const uint i,
                          const uint j,
                          const ParamDir dir) const
        {
            const uint index = j * space(ispace, dir).basisFuncN(S) + i;
            return signedGlobalI(ispace, index, dir);
        }
        
        /// No. of basis functions for this space and parametric
        /// direction.
        uint basisFuncN(const uint ispace, const ParamDir d) const
        { return space(ispace, d).basisFuncN(); }
        
        /// Total number of non-zero basis functions in this space
        /// (sum of both directions)
        uint basisFuncN(const uint ispace) const
        {
            return space(ispace, S).basisFuncN()
            + space(ispace, T).basisFuncN();
        }
        
        /// Get the number of collocation points
        uint collocPtN() const
        {
            //return globalDofN(); // same as the number of global dof
            return nonDegenerateGlobalDofN();
        }
        
        /// Get the number of collocation point for this space index
        /// and parametric direction
        uint collocPtN(const uint ispace,
                       const ParamDir d) const
        {
            return basisFuncN(ispace,d);
        }
        
        /// Get the global collocation point indes for a given space
        /// and local index
        int globalCollocI(const uint ispace,
                           const ParamDir d,
                           const uint i) const
        {
            // for now, use the global basis function connectivity
            return signedGlobalI(ispace,i, d);
        }
        
        /// Get global element index given space and local element index
        uint globalElI(const uint ispace,
                       const uint iel) const;
        
        /// Get the global element index given a space index
        /// and local indices in the s- and t- directions
        uint globalElI(const uint ispace,
                       const uint i,
                       const uint j) const;
        
        /// Get the global collocation point for the given space and local
        /// collocation index
        Point3D collocPt(const uint ispace,
                         const ParamDir d,
                         const uint i) const;
        
        /// Get the knot interval pairs for this element. This is the same
        /// for both basis directions.
        DoublePairVec knotIntervals(const uint ispace, const uint iel) const;
        
        /// Number of global dof. Takes account of degenerate dof.
        uint globalDofN() const
        {
            return nonDegenerateGlobalDofN();
//            return mGlobalDofN;
        }
          
        uint nonDegenerateGlobalDofN() const { return mGlobalNonDegenerateDofN; }
        
        
        /// Number of global dofs for nedelec discretisation
        uint globalNedelecDofN() const { return mGlobalNedelecDofN; }
        
        /// Perform the appropriate Piola transform
        virtual DoubleVecVec transformBasis(const DoubleVecVec& basis,
                                            const DoubleVecVec& jacob,
                                            const double jdet) const = 0;
        
        /// Normalised element size
        double h() const;
        
        /// Print to output stream
        void print(std::ostream& ost) const;
        
    protected:
        
        /// Default constructor. Protected since this is an abstract class.
        MultiForest()
        : mpGeom(nullptr),
          mElemN(std::make_pair(false, 0)) {}
        
        /// Construct with a geometry object
        MultiForest(const Geometry& g)
        :
        mpGeom(&g),
        mElemN(std::make_pair(false, 0)) {}
        
        /// Copy constructor
        MultiForest(const MultiForest& f);
        
        /// Copy assignment operator
        MultiForest& operator=(const MultiForest& f);
        
        /// Move constructor
        MultiForest(MultiForest&& f) = default;
        
        /// Move assignment operator
        MultiForest& operator=(MultiForest&& f) = default;
        
        /// Clear data
        void clear();
        
        /// Perform appropriate degree reduction if necessary
        /// and assign connectivity.
        /// Input vectors specify degree reduction. E.g. H(div) {{0,1},{1,0}}
        void initSpaces(const std::vector<uint>& sreduce,
                        const std::vector<uint>& treduce);
        

    protected:
        
        /// Add a space for a specified parametric direction
        void insertSpace(const BSplineSpace& s, const ParamDir dir)
        {
            switch(dir) {
                case(S):
                    insertSSpace(s);
                    break;
                case(T):
                    insertTSpace(s);
                    break;
            }
        }
        
        /// Insert global connectivity vector
        /// comprising both s- and t- node connectivities
        void insertConn(const uint isp,
                        const std::vector<uint>& c_vec)
        { mConn[isp] = c_vec;}
        
        /// Insert vector of degenerate dof
        void insertSignedConn(const uint isp,
                              const std::vector<int>& c_vec)
        { mSignedConn[isp] = c_vec; }
        
        void setSignedConn(const uint isp,
                           const uint lindex,
                           const double val)
        {
            assert(lindex < mSignedConn[isp].size());
            mSignedConn[isp].at(lindex) = val;
        }
        
        /// Insert global sign connectivity vector
        /// comprising both s- and t- sign connectivities
        void insertSignConn(const uint isp,
                            const std::vector<Sign>& s_vec)
        { mSignConn[isp] = s_vec;}
        
    private:
        
        /// Assign element edge and vertex adjacencies
        void initConnectivity();
        
        /// Generate connectivity of nedelec elements
        void initNedelecConnectivity();
        
        /// Continuity type getter (must be implemented)
        virtual ContinuityType continuityType() const = 0;
        
        /// Add space (trunk) to S space vector
        void insertSSpace(const BSplineSpace& s)
        { mSpaceS.push_back(s); }
        
        /// Add space (trunk) to T space vector
        void insertTSpace(const BSplineSpace& s)
        { mSpaceT.push_back(s); }
        
        /// Get global basis function along edge in specified param. dir.
        /// Indices are ordered positively (according to +ve local edges)
        UIntVec globalBasisIVec(const uint ispace,
                                const Edge e,
                                const ParamDir d) const;
        
        /// Get the global basis function given a space index
        /// local vertex, and parametric direction.
        uint globalBasisI(const uint ispace,
                          const Vertex v,
                          const ParamDir d) const;
        

        
        /// Implementation of element getter (caches elements)
        void initEl(const uint ielem) const;
        
        /// Caches bezier element data
        void initBezierEl(const uint ielem) const;
        
        /// Get an element given a space index and local tensor product
        /// indices. Note: this returns a nullptr if either of the
        /// local indicies i or j lies outside valid bounds
        VAnalysisElement* element(const uint ispace,
                                        const int i,
                                        const int j);
        
        /// Get a bezier element with a given space and tensor product
        /// element indices
        VAnalysisElement* bezierElement(const uint ispace,
                                        const int i,
                                        const int j);
        
        /// Subclasses can implement this if required.
        virtual void printImpl(std::ostream& ost) const;
        
        /// Reference to the geometry instance
        const Geometry* mpGeom;
        
        /// B-spline spaces, s-direction
        std::vector<BSplineSpace> mSpaceS;
        
        /// B-spline spaces, t-direction
        std::vector<BSplineSpace> mSpaceT;
        
        /// Connectivity of spaces
        std::map<uint, std::vector<uint>> mConn;
        
        /// Non-degenerate connectivity
        std::map<uint, std::vector<int>> mSignedConn;
        
        /// Orientation of vector basis in each bspline space
        std::map<uint, std::vector<Sign>> mSignConn;
        
        /// Global dof number. Calculated during connectivity construction.
        uint mGlobalDofN;
        
        /// Global dof number of non-dengerate dof.
        uint mGlobalNonDegenerateDofN;
        
        /// Global dof number for nedelec discretisation. Degenerate discretisation not yet implemented.
        uint mGlobalNedelecDofN;
        
        /// Mapping from space index and local element index to global element index
        mutable std::map<std::pair<uint, uint>, uint> mLocalElemIMap;
        
        /// Map containing vector element instances.
        mutable std::map<uint, std::unique_ptr<VAnalysisElement>> mElems;
        
        /// Map containing bezier vector element instances.
        mutable std::map<uint, std::unique_ptr<VAnalysisElement>> mBezierElems;
        
        /// Vector of nedelec elements. Connectivity information is stored within the element data structure.
        mutable std::vector<std::unique_ptr<VAnalysisElement>> mNedelecElems;
        
        /// Number of elements. Cached for efficiency.
        mutable std::pair<bool, uint> mElemN;
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const MultiForest& f)
        { f.print(ost); return ost;}
        
        /// Allow HCElement to access private members
        friend class HCElement;
        
    };
    
    /// Get the local basis func. indices over a prescribed edge of a
    /// tensor product b-spline given number of basis functions in S
    /// and T directions.
    UIntVec localBasisIVec(const Edge e,
                           const uint ns,
                           const uint nt);
    
    /// Get the local basis function for given a local vertex
    /// given number of basis functions in S and T directions
    uint localBasisI(const Vertex v,
                     const uint ns,
                     const uint nt);
    
    /// Get local basis function indices on an edge of a Bspline space
    UIntVec localBasisIVec(const Edge e,
                           const BSplineSpace& s);
    
    /// Wrapper function for localBasisI() which passes ref. to space
    uint localBasisI(const Vertex v,
                     const BSplineSpace& s);
    
    /// Get the ordered pair of local basis function indices at vertices
    /// for a given edge.
    std::pair<uint, uint> localBasisIPair(const Edge e,
                                          const BSplineSpace& s);
    
    /// Return local basis function indices given local edge,
    /// parametric direction and s-space.
    UIntVec localBasisIVec(const Edge e,
                           const ParamDir d,
                           const std::pair<BSplineSpace, BSplineSpace>& spacepair);
    
    /// Get the local basis function for a given vertex and parametric direction
    uint localBasisI(const Vertex v,
                     const ParamDir d,
                     const std::pair<BSplineSpace, BSplineSpace>& spacepair);
    
    /// Is the given local index an interior index? I.e. does not lie on the space boundary.
    bool interiorLocalI(const BSplineSpace& space,
                        const uint lindex);

	/// Given a reference space and an edge in this space, determine
	/// the orientation of the shared nodes in a space which shares this edge.
	/// The edge in the space in question is given by e_in.
	std::pair<Sign, ParamDir> vectorBasisEdgeDir(const Edge edge_ref,
												 const ParamDir d_ref,
												 const Edge e_in);

	/// Given a space with a given Vertex vert_ref and another space with vertex v_in,
	/// determine the shared local node and its orientation to the other space
	std::pair<Sign, ParamDir> vectorBasisVertexConn(const Vertex vert_ref,
													const ParamDir d_ref,
													const Vertex v_in);
    
    /// Return the local element indices along an edge. Direction specifies order.
    UIntVec localElIndices(const BSplineSpace& space,
                           const Edge e,
                           const Sign d);
    
    /// Helper function to get local element index given a vertex
    uint localElIndex(const BSplineSpace& space,
                      const Vertex v);
    
    /// Given a space and local edge, return vectors of local basis function indices that
    /// are broken by lines of C^0 continuity.
    UIntVecVec subEdgeIVecs(const BSplineSpace& space,
                            const Edge e);
    
    /// Given two edges, is the basis numbering along the edge +ve (same direction) or -ve (opposite direction)
    Sign edgeOrientation(Edge e1, Edge e2);
    
    /// Given a ContinuityType (normal or tangential) and edge
    /// return the appropriate parametric direction
    ParamDir paramDir(const Edge e, const ContinuityType c);
    
    /// Projects a point living in element 1 onto the adjacent edge
    /// in element 2.
    GPt2D projectPt(const GPt2D p, const Edge e1, const Edge e2);
    
    /// Get parametric coordinate of vertex
    GPt2D paramPt(const Vertex v);

    
}
#endif

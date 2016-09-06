#ifndef NURBS_NEDELEC_VECTOR_ELEMENT_H
#define NURBS_NEDELEC_VECTOR_ELEMENT_H

#include "AnalysisElement.h"
#include "BezierNodalElement.h"
#include "MultiForest.h"
#include "base.h"
#include <cassert>
#include <numeric>

namespace nurbs {
    
    /// Helper functions
    
    /// Get the set of Legendre bais functions at the given local coordinate
    /// xi \in [-1,1] and given degree p
    std::vector<double> LegendreBasis1D(const double xi,
                                        const uint p);
    
    std::vector<double> LegendreBasisDer1D(const double xi,
                                           const uint p);
    
    /// A representation of a Nedelec element that prescribes
    /// a vector basis of arbitrary degree
    
    class NedelecVectorElement : public VAnalysisElement
    {
        
    public:
        
        /// Construct with given multiforest, space index, element index and pointer to parent element
        NedelecVectorElement(const MultiForest* f,
                             const uint ispace,
                             const uint iel,
                             const UIntVec& unsignedconn,
                             const IntVec& signedconn,
                             const std::vector<Sign>& orientationconn,
                             const BezierNodalElement* pel = nullptr)
        :
        VAnalysisElement(*f->geometry(),
                         ispace,
                         f->knotIntervals(ispace, iel)),
        mMultiForest(f),
        mSpaceI(ispace),
        mSpaceElemI(iel),
        mpParent(pel),
        mUnsignedConn(unsignedconn),
        mSignedConn(signedconn),
        mDirectionConn(orientationconn) {}
        
        /// Copy function
        std::unique_ptr<AnalysisElement> copy() const override
        {
            return make_unique<NedelecVectorElement>(*this);
        }
        
        /// Number of parametric components
        uint componentN() const override { return 2; }
        
        /// Number of non-zero basis functions over this element.
        uint basisFuncN() const override { return localBasisFuncI().size(); }
        
        /// Evaluate physical coordinates using Bezier extraction
        virtual Point3D eval(const double u, const double v) const override;
        
        /// Evaluate tangent vector using Bezier extraction
        virtual Point3D tangent(const double u,
                                const double v,
                                const ParamDir dir) const override;
        
        /// Override normal evaluation
        virtual Point3D normal(const double u, const double v) const override
        {
            return cross(tangent(u,v,S), tangent(u,v,T)).normalise();
        }
        
        /// Override jacobian determinant evaluation
        virtual double jacDet(const double u, const double v) const override
        {
            return cross(tangent(u,v,S), tangent(u,v,T)).length() * jacDetParam(u,v);
        }
        
        /// Jacobian determinant with given tangent vectors
        virtual double jacDet(const double u,
                              const double v,
                              const Point3D& t1,
                              const Point3D& t2) const override
        {
            return cross(t1,t2).length() * jacDetParam(u,v);
        }
        
        /// Override jacobian evaluation
        virtual DoubleVecVec jacob(const double u, const double v) const override
        {
            DoubleVecVec jacob_param;
            jacob_param.push_back(tangent(u,v,S).asVec());
            jacob_param.push_back(tangent(u,v,T).asVec());
            const auto jacob_parent = jacobParam(u, v);
            DoubleVecVec r{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }};
            for(uint i = 0; i < 2; ++i)
                for(uint j = 0; j < 3; ++j)
                    for(uint k = 0; k < 2; ++k)
                        r[i][j] += jacob_parent[i][k] * jacob_param[k][j];
            return r;
        }
        
        /// Local basis function indices. Simply {0,1,...,n_b - 1}
        virtual UIntVec localBasisFuncI() const override
        {
            const uint nbasis = degree(ParamDir::S, 0) * degree(ParamDir::S, 1) + degree(ParamDir::T, 0) * degree(ParamDir::T, 1);
            UIntVec rvec(nbasis);
            std::iota(rvec.begin(), rvec.end(), 0);
            return rvec;
        }
        
        /// Simply call the connectivity stored in the multiforest
        virtual UIntVec globalBasisFuncI() const override
        {
            return mUnsignedConn;
        }
     
        /// Simply call the connectivity stored in the multiforest
        virtual IntVec signedGlobalBasisFuncI() const override
        {
            return mSignedConn;
        }
        
        /// Return the global basis function set
        virtual DoubleVecVec basis(const double u, const double v) const override
        {
            DoubleVecVec jacob_param;
            const auto& t1 = tangent(u,v,S);
            const auto& t2 = tangent(u,v,T);
            jacob_param.push_back(t1.asVec());
            jacob_param.push_back(t2.asVec());
            const auto jdet = cross(t1, t2).length();
            return multiForest()->transformBasis(localBasis(u,v), jacob_param, jdet);
        }
        
        /// Same as above, but with given tangent vectorss
        DoubleVecVec basis(const double u,
                           const double v,
                           const Point3D& t1,
                           const Point3D& t2) const override
        {
            DoubleVecVec jacob_param;
            jacob_param.push_back(t1.asVec());
            jacob_param.push_back(t2.asVec());
            const auto jdet = cross(t1, t2).length();
            return multiForest()->transformBasis(localBasis(u,v), jacob_param, jdet);
        }
        
        /// Return the set of local Legendre basis functions that
        /// define the vector basis WITHOUT applying the Piola transform
        virtual DoubleVecVec localBasis(const double u,
                                        const double v) const override
        {
            DoubleVecVec rvec;
            
            // loop over each vector component
            for(uint icomp = 0; icomp < componentN(); ++icomp)
            {
                const auto sbasis = LegendreBasis1D(u, degree(ParamDir::S, icomp));
                const auto tbasis = LegendreBasis1D(v, degree(ParamDir::T, icomp));
                
                const ParamDir compdir = ParamDirType(icomp);
                for(const auto& tval : tbasis)
                    for(const auto& sval : sbasis)
                        rvec.push_back(compdir == ParamDir::S ? DoubleVec{sval * tval, 0.0} : DoubleVec{0.0, sval * tval});
            }
            
            // and now apply global orientation
            assert(rvec.size() == mDirectionConn.size());
            
            for(uint ibasis = 0; ibasis < rvec.size(); ++ibasis)
                for(auto& b : rvec[ibasis])
                    b *= asDouble(mDirectionConn[ibasis]);
            
            return rvec;
        }
        
        /// Get the relevant basis function derivatives
        /// with the sign applied
        virtual DoubleVecVec localBasisDers(const double u,
                                            const double v,
                                            const DerivType dtype) const override
        {
            DoubleVecVec rvec;
            
            // loop over each vector component
            for(uint icomp = 0; icomp < componentN(); ++icomp)
            {
                
                DoubleVec sbasis, tbasis;
                if(DerivType::DS == dtype)
                {
                    sbasis = LegendreBasisDer1D(u, degree(ParamDir::S, icomp));
                    tbasis = LegendreBasis1D(v, degree(ParamDir::T, icomp));
                }
                else
                {
                    sbasis = LegendreBasis1D(u, degree(ParamDir::S, icomp));
                    tbasis = LegendreBasisDer1D(v, degree(ParamDir::T, icomp));
                }
                
                const ParamDir compdir = ParamDirType(icomp);
                for(const auto& tval : tbasis)
                    for(const auto& sval : sbasis)
                        rvec.push_back(compdir == ParamDir::S ? DoubleVec{sval * tval, 0.0} : DoubleVec{0.0, sval * tval});
                        }
            
            // and now apply global orientation
            assert(rvec.size() == mDirectionConn.size());
            
            for(uint ibasis = 0; ibasis < rvec.size(); ++ibasis)
                for(auto& b : rvec[ibasis])
                    b *= asDouble(mDirectionConn[ibasis]);
                    
                    return rvec;
        }
        
        /// Basis function degrees
        virtual uint degree(const ParamDir dir, const uint comp) const override
        {
            return space(ParamDirType(comp)).degree(dir);
        }
        
        /// Basis function degrees
        virtual UIntVec degree(const uint comp) const override
        {
            return space(ParamDirType(comp)).degree();
        }
        
        /// Get the number of connected collocation points on this element
        virtual uint collocPtN() const override
        {
            // TODO!
            return 0;
        }
        
        /// Get the global connected collocation indices
        virtual UIntVec globalCollocConn() const override
        {
            // TODO!
            return {};
        }
        
        /// Get the global collocation index given a local index
        virtual uint globalCollocI(const uint icpt) const override
        {
            // TODO!
            return 0;
        }
        
        /// Get the parent coordinate of a given collocation point index
        virtual GPt2D collocParentCoord(const uint icpt) const override
        {
            // TODO!
            return GPt2D();
        }
        
        /// Print out info for this element
        virtual void print(std::ostream& ost) const override
        {
            //            ost << "Nedelec vector element" << "\n";
            GeometryElement::print(ost);
        }
        
        /// Multiforest getter
        const MultiForest* multiForest() const { return mMultiForest; }
        
        
    private:
        
        /// Local element index getter
        uint spaceElementI() const { return mSpaceElemI; }
        
        /// Space index accessor
        uint spaceI() const { return mSpaceI; }
        
        /// Space getter for given parametric direction
        const BSplineSpace& space(const ParamDir d) const
        {
            return multiForest()->space(spaceI(), d);
        }
        
        /// Pointer to multiforest (where this element is stored)
        const MultiForest* mMultiForest;
        
        /// Bspline space index for this element
        const uint mSpaceI;
        
        /// Local 'element' index in the bspline space
        const uint mSpaceElemI;
        
        /// Pointer to parent element where geometry info is accessed
        const BezierNodalElement* mpParent;
        
        /// Unsigned connectivity (i.e. degenerate points not taken into account)n
        const UIntVec mUnsignedConn;
        
        /// Signed element connectivity (degenerate points taken into account)
        const IntVec mSignedConn;
        
        /// Orientation of basis functions (in global mesh)
        const std::vector<Sign> mDirectionConn;
        
    };

    
}

#endif
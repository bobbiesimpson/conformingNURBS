#ifndef NURBS_BEZIER_ELEMENT_H
#define NURBS_BEZIER_ELEMENT_H

#include "AnalysisElement.h"
#include "BezierNodalElement.h"
#include "MultiForest.h"

namespace nurbs {
    
    
    ///
    /// A representation of a vector bezier element
    /// used to represent vector fields primarily
    /// in electromagnetic analysis.
    ///
    class BezierVectorElement : public VAnalysisElement {
        
    public:
        
        BezierVectorElement(const MultiForest* f,
                            const uint ispace,
                            const uint iel,
                            const BezierNodalElement* pel = nullptr)
        : VAnalysisElement(*f->geometry(),
                           ispace,
                           f->knotIntervals(ispace, iel)),
        mMultiForest(f),
        mSpaceI(ispace),
        mElemI(iel),
        mpParentEl(pel) {}
        
        /// Virtual copy function
        std::unique_ptr<AnalysisElement> copy() const override
        {
            return make_unique<BezierVectorElement>(*this);
        }
        
        /// Scalar basis. Only one component.
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
        { return cross(t1,t2).length() * jacDetParam(u,v); }
        
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
        
        /// Jacobian with given tangent vectors (for efficiency)
        DoubleVecVec jacob(const double u,
                           const double v,
                           const Point3D& t1,
                           const Point3D& t2) const override
        {
            DoubleVecVec jacob_param;
            jacob_param.push_back(t1.asVec());
            jacob_param.push_back(t2.asVec());
            const auto jacob_parent = jacobParam(u, v);
            DoubleVecVec r{ { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }};
            for(uint i = 0; i < 2; ++i)
                for(uint j = 0; j < 3; ++j)
                    for(uint k = 0; k < 2; ++k)
                        r[i][j] += jacob_parent[i][k] * jacob_param[k][j];
            return r;
        }
        
        /// Return local basis function indices that are non-zero over the
        /// relevant Bspline space
        virtual UIntVec localBasisFuncI() const override
        {
            if(!mLocalBasisIVec.empty())
                return mLocalBasisIVec;
            UIntVec lvec;
            const ParamCoord c = paramCoord(lowerBound(S), lowerBound(T));
            auto s_indices = space(S).globalBasisFuncI(c.s, c.t);
            for(const auto& i : s_indices)
                lvec.push_back(i);
                
                // Now add basis indices in T-direction.
                const uint nbasis_s = multiForest()->basisFuncN(spaceI(), S);
                auto t_indices = space(T).globalBasisFuncI(c.s, c.t);
                for(const auto& i : t_indices)
                    lvec.push_back(i + nbasis_s);
                    mLocalBasisIVec = lvec;
            return mLocalBasisIVec;
        }
        
        /// Return the global (to the forest) basis function indices
        virtual UIntVec globalBasisFuncI() const override
        {
            if(!mGlobalBasisIVec.empty())
                return mGlobalBasisIVec;
            UIntVec temp;
            for(const auto& i : localBasisFuncI())
                temp.push_back(multiForest()->globalI(spaceI(), i));
                mGlobalBasisIVec = temp;
            return mGlobalBasisIVec;
        }
        
        /// Return the basis function values with appropriate piola.
        /// Since this computes the jacobian from scratch, this is not
        /// as efficient as the alternative version of this function
        /// which takes tangent vectors as inputs.
        virtual DoubleVecVec basis(const double u, const double v) const override
        {
            return multiForest()->transformBasis(localBasis(u,v), jacob(u, v), jacDet(u, v));
        }
        
        /// More efficient basis function evaluation with given tangent vectors
        DoubleVecVec basis(const double u,
                           const double v,
                           const Point3D& t1,
                           const Point3D& t2) const override
        {
            return multiForest()->transformBasis(localBasis(u,v), jacob(u,v,t1,t2), jacDet(u,v,t1,t2));
        }
        
        /// Get the untransformed (i.e. without Piola) vector basis functions
        /// using Bezier extraction
        virtual DoubleVecVec localBasis(const double u,
                                const double v) const override
        {
            // The return vector of (vector-valued) basis functions
            DoubleVecVec rvec;
            
            for(uint icomp = 0; icomp < componentN(); ++icomp)
            {
                const ParamDir compdir = ParamDirType(icomp);
                auto indices = space(compdir).localIndices(localElementI());
                
                // Using a reference makes a huge difference to speed
                const auto& op_u = space(compdir).extractionOperator(indices.first, S);
                const auto& op_v = space(compdir).extractionOperator(indices.second, T);
                //
//                
//                std::cout << op_u << "\n";
//                std::cout << op_v << "\n";
//                
                const auto b_u = nurbshelper::bernsteinPolynomial(u, degree(S, icomp));
                const auto b_v = nurbshelper::bernsteinPolynomial(v, degree(T, icomp));
                
                // now apply extraction operators to bernstein basis
                const uint n_u = op_u.size();
                std::vector<double> basis_u(n_u, 0.0);
                for(uint i = 0; i < n_u; ++i)
                    for(uint j = 0; j < op_u[0].size(); ++j)
                        basis_u[i] += op_u[i][j] * b_u[j];
                
                const uint n_v = op_v.size();
                std::vector<double> basis_v(n_v, 0.0);
                for(uint i = 0; i < n_v; ++i)
                    for(uint j = 0; j < op_v[0].size(); ++j)
                        basis_v[i] += op_v[i][j] * b_v[j];
                
                for(uint j = 0; j < basis_v.size(); ++j)
                    for(uint i = 0; i < basis_u.size(); ++i)
                    {
                        const double b = basis_u[i] * basis_v[j];
                        rvec.emplace_back(compdir == ParamDir::S ? DoubleVec{b, 0.0} : DoubleVec{0.0, b});
                    }
            }
            
            // now apply sign
            const auto& sign_vec = multiForest()->globalDirVec(spaceI());
            const auto l_ivec = localBasisFuncI();
            assert(l_ivec.size() == rvec.size());
            for(uint i = 0; i < basisFuncN(); ++i) {
                for(auto& b : rvec[i])
                    b *= asDouble(sign_vec[l_ivec[i]]);
            }
            return rvec;
        }
        
        /// Get the local basis derivatives either in S or T direction
        virtual DoubleVecVec localBasisDers(const double u,
                                            const double v,
                                            const DerivType dtype) const override
        {
            DoubleVecVec rvec;
            
            for(uint icomp = 0; icomp < componentN(); ++icomp)
            {
                const ParamDir compdir = ParamDirType(icomp);
                auto indices = space(compdir).localIndices(localElementI());
                
                // Using a reference makes a huge difference to speed
                const auto& op_u = space(compdir).extractionOperator(indices.first, S);
                const auto& op_v = space(compdir).extractionOperator(indices.second, T);
                
                // Get Bertnein basis in each parametric direction taking account
                // of derivatives.
                DoubleVec bu;
                DoubleVec bv;
                if(DerivType::DS == dtype)
                {
                    bu = nurbshelper::bernsteinPolynomialDeriv(u, degree(S, icomp));
                    bv = nurbshelper::bernsteinPolynomial(v, degree(T, icomp));
                }
                else
                {
                    bu = nurbshelper::bernsteinPolynomial(u, degree(S, icomp));
                    bv = nurbshelper::bernsteinPolynomialDeriv(v, degree(T, icomp));
                }
                
                const uint n_u = op_u.size();
                std::vector<double> basis_u(n_u, 0.0);
                for(uint i = 0; i < n_u; ++i)
                    for(uint j = 0; j < op_u[0].size(); ++j)
                        basis_u[i] += op_u[i][j] * bu[j];
                
                const uint n_v = op_v.size();
                std::vector<double> basis_v(n_v, 0.0);
                for(uint i = 0; i < n_v; ++i)
                    for(uint j = 0; j < op_v[0].size(); ++j)
                        basis_v[i] += op_v[i][j] * bv[j];
                
                const auto param_j = jacobParam(u, v);
//                const double jterm =
//                (DerivType::DS == dtype) ? 1.0/param_j[0][0] : 1.0 / param_j[1][1];
                
                for(uint j = 0; j < basis_v.size(); ++j)
                    for(uint i = 0; i < basis_u.size(); ++i)
                    {
                        rvec.push_back(ParamDir::S == compdir ? DoubleVec{basis_u[i] * basis_v[j], 0.0} : DoubleVec{0.0, basis_u[i] * basis_v[j]});
                    }
            }
            
            // now apply sign
            const auto& sign_vec = multiForest()->globalDirVec(spaceI());
            const auto& l_ivec = localBasisFuncI();
            assert(l_ivec.size() == rvec.size());
            for(uint i = 0; i < basisFuncN(); ++i)
            {
                for(auto& b : rvec[i])
                    b *= asDouble(sign_vec[l_ivec[i]]);
            }
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
            ost << "Bezier Vector element\n";
            GeometryElement::print(ost);
        }
        
        /// Multiforest getter
        const MultiForest* multiForest() const { return mMultiForest; }
        
        
    private:
        
        /// Local element index getter
        uint localElementI() const { return mElemI; }
        
        /// Space index accessor
        uint spaceI() const { return mSpaceI; }
        
        /// Space getter
        const BSplineSpace& space(const ParamDir d) const
        {
            return multiForest()->space(spaceI(), d);
        }
        
        /// Reference to multiforest
        const MultiForest* mMultiForest;
        
        /// The space index in the multiforest
        const uint mSpaceI;
        
        /// Element index local to spaces
        const uint mElemI;
        
        /// Reference to parent (geometry) bezier element
        const BezierNodalElement* mpParentEl;
        
        /// Cached local basis functions which are non-zero over this element
        mutable UIntVec mLocalBasisIVec;
        
        /// Cached global basis function indices which are non-zero over this element
        mutable UIntVec mGlobalBasisIVec;
        
        /// Cached of global basis function directions over this element
        mutable std::vector<Sign> mGlobalSignVec;
        
        
    };
    
}
#endif
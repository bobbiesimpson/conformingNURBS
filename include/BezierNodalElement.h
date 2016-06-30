#ifndef BEZIER_NODAL_ELEMENT_H
#define BEZIER_NODAL_ELEMENT_H

#include "AnalysisElement.h"
#include "BSplineSpace.h"
#include "Forest.h"
#include "base.h"
#include "NURBSCommon.h"

namespace nurbs {
    
    ///
    /// A reprenestation of a Bezier element that
    /// provides fast evaluation of (rational) B-spline
    /// basis functions and derivatives for analysis.
    ///
    
    class BezierNodalElement : public NAnalysisElement {
        
    public:
        
        /// Constructor
        BezierNodalElement(const Forest* f,
                           const uint ispace,
                           const uint iel,
                           const BezierNodalElement* pel = nullptr)
        : NAnalysisElement(*f->geometry(),
                           ispace,
                           f->knotIntervals(ispace,iel)),
        mForest(f),
        mSpace(&f->space(ispace)),
        mElemI(iel),
        mpParentEl(pel) {}
        
        /// Virtual copy function
        std::unique_ptr<AnalysisElement> copy() const override
        {
            return make_unique<BezierNodalElement>(*this);
        }
        
        /// Scalar basis. Only one component.
        uint componentN() const override { return 1; }
        
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
        
        /// Return local basis function indices that are non-zero over the
        /// relevant Bspline space
        UIntVec localBasisFuncI() const override
        {
            if(!mLocalBasisFuncI.empty())
                return mLocalBasisFuncI;
            ParamCoord c = paramCoord(lowerBound(S), lowerBound(T));
            mLocalBasisFuncI = space()->globalBasisFuncI(c.s, c.t);
            return mLocalBasisFuncI;
        }
        
        /// Return the global (to the forest) basis function indices
        UIntVec globalBasisFuncI() const override
        {
            if(!mGlobalBasisFuncI.empty())
                return mGlobalBasisFuncI;
            for(const auto& i : localBasisFuncI())
                mGlobalBasisFuncI.emplace_back(forest()->globalI(spaceI(),i));
                return mGlobalBasisFuncI;
        }
        
        /// REturn the basis function values
        DoubleVec basis(const double u, const double v) const override
        {
            auto indices = space()->localIndices(localElementI());
            
            // Using a reference makes a huge difference to speed
            const auto& op_u = space()->extractionOperator(indices.first, S);
            const auto& op_v = space()->extractionOperator(indices.second, T);
//
            const auto b_u = nurbshelper::bernsteinPolynomial(u, degree(S));
            const auto b_v = nurbshelper::bernsteinPolynomial(v, degree(T));
            
//            std::vector<double> rvec;
//            for(uint j = 0; j < b_v.size(); ++j)
//                for(uint i = 0; i < b_u.size(); ++i)
//                    rvec.push_back(b_u[i] * b_v[j]);
//            return rvec;
            
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
            
            std::vector<double> final;
            for(uint j = 0; j < basis_v.size(); ++j)
                for(uint i = 0; i < basis_u.size(); ++i)
                    final.push_back(basis_u[i] * basis_v[j]);
            return final;
        }
        
        /// Get the local basis derivatives either in S or T direction
        DoubleVec localBasisDers(const double u,
                                 const double v,
                                 const DerivType dtype) const override
        {
            auto indices = space()->localIndices(localElementI());
            
            // Using a reference makes a huge difference to speed
            const auto& op_u = space()->extractionOperator(indices.first, S);
            const auto& op_v = space()->extractionOperator(indices.second, T);
            
            // Get Bertnein basis in each parametric direction taking account
            // of derivatives.
            DoubleVec bu;
            DoubleVec bv;
            if(DerivType::DS == dtype) {
                bu = nurbshelper::bernsteinPolynomialDeriv(u, degree(S));
                bv = nurbshelper::bernsteinPolynomial(v, degree(T));
            }
            else {
                bu = nurbshelper::bernsteinPolynomial(u, degree(S));
                bv = nurbshelper::bernsteinPolynomialDeriv(v, degree(T));
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
            
            std::vector<double> final;
            const auto param_j = jacobParam(u, v);
            const double jterm = (DerivType::DS == dtype) ? 1.0/param_j[0][0] : 1.0 / param_j[1][1];
            
            for(uint j = 0; j < basis_v.size(); ++j)
                for(uint i = 0; i < basis_u.size(); ++i)
                    final.push_back(basis_u[i] * basis_v[j] * jterm);
                    
            return final;
        }
        
        /// Basis function degrees
        uint degree(const ParamDir dir, const uint comp = 0) const override
        {
            return space()->degree(dir);
        }
        
         /// Basis function degrees
        UIntVec degree(const uint comp = 0) const override
        {
            return space()->degree();
        }
        
        /// Get the number of connected collocation points on this element
        uint collocPtN() const override
        {
            return forest()->connectedCollocPtN(spaceI(), localElementI());
        }
        
        /// Get the global connected collocation indices
        UIntVec globalCollocConn() const override
        {
            return forest()->connectedGlobalCollocI(spaceI(), localElementI());
        }
        
        /// Get the global collocation index given a local index
        uint globalCollocI(const uint icpt) const override
        {
            return forest()->connectedGlobalCollocI(spaceI(), localElementI(), icpt);
        }
        
        /// Get the parent coordinate of a given collocation point index
        GPt2D collocParentCoord(const uint icpt) const override
        {
            const uint ilocal = forest()->connectedLocalCollocI(spaceI(), localElementI(), icpt);
            return parentCoord(space()->grevilleAbscissaPt(ilocal));
        }
        
        /// TODO
        void print(std::ostream& ost) const override
        {
            GeometryElement::print(ost);
            //            ost << "Forest: " << forest() << "\n";
            //            ost << "Space: " << space() << "\n";
            //            ost << "Local basis function indices: " << localBasisFuncI() << "\n";
            //            ost << "Global basis function indicies" << globalBasisFuncI() << "\n";
        }
        
        /// Local element index getter
        uint localElementI() const { return mElemI; }
        
        /// Forest accessor
        const Forest* forest() const {return mForest;}
        
        /// Accessor for space
        const BSplineSpace* space() const{ return mSpace; }
        
    
    protected:
        
    private:
        
        /// Reference to forest
        const Forest* mForest;
        
        /// Reference to space this element belongs to
        const BSplineSpace* mSpace;

        /// Element index (local to space)
        const uint mElemI;
        
        /// Reference to 'parent' element which exists
        /// in the primal forest
        const BezierNodalElement* mpParentEl;
        
        /// Cache the knot span indices for evaluating basis functions
        mutable UIntVec mLocalBasisFuncI;
        
        /// A vector containing the global basis function indices. (i.e. the global forest basis index)
        mutable UIntVec mGlobalBasisFuncI;
    };
    
}
#endif
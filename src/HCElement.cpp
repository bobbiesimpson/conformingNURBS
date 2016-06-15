#include <stdexcept>

#include "HCElement.h"

namespace nurbs {
    
    UIntVec HCElement::localBasisFuncI() const
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
    
    UIntVec HCElement::globalBasisFuncI() const
    {
        if(!mGlobalBasisIVec.empty())
            return mGlobalBasisIVec;
        UIntVec temp;
        for(const auto& i : localBasisFuncI())
            temp.push_back(multiForest()->globalI(spaceI(), i));
        mGlobalBasisIVec = temp;
        return mGlobalBasisIVec;
    }
    
    DoubleVecVec HCElement::basis(const double u, const double v) const
    {
        // Apply relevant Piola transform.
        return multiForest()->transformBasis(localBasis(u,v), jacob(u, v), jacDet(u, v));
    }
    
    DoubleVecVec HCElement::localBasis(const double u, const double v) const
    {
        ParamCoord p = paramCoord(u,v);
        DoubleVecVec rvec;
        //std::cout << space(S).basis(p.s, p.t, span(S)) << "\n";
        for(const auto& b : space(S).basis(p.s, p.t, span(S)))
            rvec.push_back({b, 0.0});
        for(const auto& b : space(T).basis(p.s, p.t, span(T)))
            rvec.push_back({0.0, b});
        // Now apply sign
        const auto& sign_vec = multiForest()->globalDirVec(spaceI());
        const auto l_ivec = localBasisFuncI();
        assert(l_ivec.size() == rvec.size());
        for(uint i = 0; i < basisFuncN(); ++i) {
            for(auto& b : rvec[i])
                b *= asDouble(sign_vec[l_ivec[i]]);
        }
        return rvec;
    }
    
    DoubleVecVec HCElement::localBasisDers(const double u,
                                           const double v,
                                           const DerivType dtype) const
    {
        const ParamCoord p = paramCoord(u,v);
        const DoubleVecVec jp = jacobParam(u,v); // parametric jacob.
        DoubleVecVec rvec; // vector to populate with derivatives
        
        if(DS == dtype) {
            for(const auto& bd : space(S).basisDers(p.s,p.t, dtype))
                rvec.push_back({bd * jp[0][0], 0.0});
            for(const auto& bd : space(T).basisDers(p.s,p.t,dtype))
                rvec.push_back({0.0, bd * jp[0][0]});
        }
        else if(DT == dtype) {
            for(const auto& bd : space(S).basisDers(p.s,p.t, dtype))
                rvec.push_back({bd * jp[1][1], 0.0});
            for(const auto& bd : space(T).basisDers(p.s,p.t,dtype))
                rvec.push_back({0.0, bd * jp[1][1]});
        }
        else
                throw std::runtime_error("Bad derivative specified.");
        // Now apply sign
        const auto& sign_vec = multiForest()->globalDirVec(spaceI());
        const auto l_ivec = localBasisFuncI();
        assert(l_ivec.size() == rvec.size());
        for(uint i = 0; i < basisFuncN(); ++i) {
            for(auto& b : rvec[i])
                b *= asDouble(sign_vec[l_ivec[i]]);
        }
        return rvec;
    }
    
    UIntVec HCElement::span(const ParamDir d) const
    {
        auto it = mSpan.find(d);
        if(it != mSpan.end())
            return it->second;
        mSpan[d] = space(d).span(lowerBound(S), lowerBound(T));
        return mSpan[d];
    }
    
}
#include "ISubElem.h"

namespace nurbs
{
    namespace elem {
        
        ISubElem::ISubElem(const uint ns,
                           const uint nt )
        :
        mNSubEls(ns * nt),
        mCurrentIndex(0)
        {
            // initialise sub element ranges
            const double s_inc = 2.0 / ns;
            const double t_inc = 2.0 / nt;
            for( uint i = 0; i < ns; ++i )
            {
                double s_l = -1.0 + (double) i * s_inc;
                double s_r = -1.0 + (double) ( i + 1 ) * s_inc;
                for( uint j = 0; j < nt; ++j )
                {
                    double t_l = -1.0 + (double) j * t_inc;
                    double t_r = -1.0 + (double) ( j + 1 ) * t_inc;
                    DoubleVec r;
                    r.push_back( s_l );
                    r.push_back( s_r );
                    r.push_back( t_l );
                    r.push_back( t_r );
                    mRanges.push_back( r );
                }
            }
        }
        
        ISubElem::ISubElem(const nurbs::DoubleVec& sknots,
                           const nurbs::DoubleVec& tknots)
        :
        mNSubEls((sknots.size() + 1) * (tknots.size() + 1)),
        mCurrentIndex(0)
        {
            const uint ns = sknots.size() + 1;
            const uint nt = tknots.size() + 1;
            
            // initialise sub element ranges
            for(uint i = 0; i < ns; ++i)
            {
                const double s_l = (0 ==i) ? -1.0 : sknots.at(i-1);
                const double s_r = (ns - 1 == i) ? 1.0 : sknots.at(i);
                
                if(s_r < s_l)
                    throw std::runtime_error("Bad knot value given for subelement integration domains.");
                if(std::abs(s_r) > 1.0 || std::abs(s_l) > 1.0)
                    throw std::runtime_error("Bad knot value given for subelement integration domains: value outside interval [-1.0, 1.0].");
                   
                
                for( uint j = 0; j < nt; ++j )
                {
                    const double t_l = (0 == j) ? -1.0 : tknots.at(j-1);
                    const double t_r = (nt - 1 == j) ? 1.0 : tknots.at(j);
                    
                    if(t_r < t_l)
                        throw std::runtime_error("Bad knot value given for subelement integration domains.");
                    if(std::abs(t_r) > 1.0 || std::abs(t_l) > 1.0)
                        throw std::runtime_error("Bad knot value given for subelement integration domains: value outside interval [-1.0, 1.0].");
                    
                    DoubleVec r;
                    r.push_back( s_l );
                    r.push_back( s_r );
                    r.push_back( t_l );
                    r.push_back( t_r );
                    mRanges.push_back( r );
                }
            }
        }
        
//        ISubElem::ISubElem(const double x,
//                           const nurbs::ParamDir d)
//        :
//        mNSubEls(2),
//        mCurrentIndex(0),
//        mTellesTransformApplied(std::make_pair(true, d))
//        {
//            if(ParamDir::S == d)
//            {
//                DoubleVec subcell1{-1.0, x, -1.0, 1.0};
//                DoubleVec subcell2{x, 1.0, -1.0, 1.0};
//                mRanges.push_back(subcell1);
//                mRanges.push_back(subcell2);
////                mTellesVec.push_back(std::make_pair(ParamDir::S, Sign::POSITIVE));
////                mTellesVec.push_back(std::make_pair(ParamDir::S, Sign::NEGATIVE));
//                
//            }
//            else
//            {
//                DoubleVec subcell1{-1.0, 1.0, -1.0, x};
//                DoubleVec subcell2{-1.0, 1.0, x, 1.0};
//                mRanges.push_back(subcell1);
//                mRanges.push_back(subcell2);
//            }
////            mTellesVec.push_back(std::make_pair(ParamDir::S, Sign::POSITIVE));
////            mTellesVec.push_back(std::make_pair(ParamDir::S, Sign::NEGATIVE));
//        }
        
        DoubleVec ISubElem::getRange() const
        {
            return mRanges[mCurrentIndex];
        }
        
        double ISubElem::jacob() const
        {
            DoubleVec current_range = getRange();
            return ( current_range[1] - current_range[0] ) / 2.0 * ( current_range[3] - current_range[2] ) / 2.0;
        }
        
//        double ISubElem::tellesJacob(const nurbs::GPt2D& gpt) const
//        {
//            // NEEDS TESTING!!!!
//            if(!tellesTransformApplied())
//                return jacob();
//            
//            assert(currentIndex() < 2);
//            
//            const double coord = gpt.get(tellesTransformDirection());
//            const double sign = (0 == currentIndex()) ? 1.0 : -1.0;
//            
//            return (1.0 - coord * sign) * jacob();
//            
//        }
        
        
        GPt2D ISubElem::get(const GPt2D& pt) const
        {
            auto pt_copy = pt;
            
//            if(tellesTransformApplied())
//            {
//                // TODO
//            }
            
            DoubleVec current_range = getRange();
            const double s =  (current_range[1] - current_range[0] ) / 2.0 * pt_copy.get(0) + (current_range[0] + current_range[1] ) / 2.0;
            const double t =  (current_range[3] - current_range[2] ) / 2.0 * pt_copy.get(1) + ( current_range[2] + current_range[3] ) / 2.0;
            return GPt2D(s,t);
        }
    }
}

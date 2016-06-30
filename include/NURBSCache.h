#ifndef NURBS_CACHE_H
#define NURBS_CACHE_H

#include "base.h"

#include <map>
#include <tuple>
#include <utility>
#include <mutex>


namespace nurbs {
    
    namespace nurbshelper {
        
        /// A singleton class responsible for caching basis functions,
        /// derivatives, spans etc. Bezier extraction may be a better
        /// future option
        
        class NURBSCache {
        
        public:
            
            /// Get the one and only instance
            static NURBSCache& Instance()
            {
                static NURBSCache theInstance;
                return theInstance;
            }
            
            /// Check is span has been cached and populate with relevant value
            std::pair<bool, uint> span(double s,
                                       const DoubleVec& knotvec,
                                       const uint p) const
            {
//                auto it = mSpanMap.find(std::make_tuple(s,knotvec,p));
//                if(it != mSpanMap.end())
//                    return std::make_pair(true, it->second);
                return std::make_pair(false, 0);
            }
            
            /// Cache the given span
            bool cacheSpan(double s,
                           const DoubleVec& knotvec,
                           const uint p,
                           const uint span)
            {
//                std::lock_guard<std::mutex> lock(mMutex);
//                auto insert = mSpanMap.insert(std::make_pair(std::make_tuple(s,knotvec,p), span));
//                return insert.second;
                return true;
            }
            
            /// Check if bspline basis has been cached and populate with relevant value
            std::pair<bool, DoubleVec> basis(const double s,
                                        const uint span,
                                        const DoubleVec& knotvec,
                                        const uint p) const
            {
//                auto it = mBasisMap.find(std::make_tuple(s, span, knotvec, p));
//                if(it != mBasisMap.end())
//                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVec{});
            }
            
            /// Check if bspline basis has been cached and populate with relevant value
            std::pair<bool, DoubleVec> bernsteinBasis(const double s,
                                                      const uint p) const
            {
                auto it = mBernsteinBasisMap.find(std::make_tuple(s,p));
                if(it != mBernsteinBasisMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVec{});
            }
            
            /// Cache the given span
            bool cacheBernsteinBasis(double s,
                                     const uint p,
                                     const DoubleVec& basis)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mBernsteinBasisMap.insert(std::make_pair(std::make_tuple(s,p), basis));
                return insert.second;
            }
            
            /// Check if bspline basis has been cached and populate with relevant value
            std::pair<bool, DoubleVec> bernsteinBasisDeriv(const double s,
                                                           const uint p) const
            {
                auto it = mBernsteinBasisDerivMap.find(std::make_tuple(s,p));
                if(it != mBernsteinBasisDerivMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVec{});
            }
            
            /// Cache the given span
            bool cacheBernsteinBasisDeriv(double s,
                                          const uint p,
                                          const DoubleVec& basis)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mBernsteinBasisDerivMap.insert(std::make_pair(std::make_tuple(s,p), basis));
                return insert.second;
            }
            
            /// Cache the given span
            bool cacheBasis(double s,
                            const uint span,
                            const DoubleVec& knotvec,
                            const uint p,
                            const DoubleVec& basis)
            {
//                std::lock_guard<std::mutex> lock(mMutex);
//                auto insert = mBasisMap.insert(std::make_pair(std::make_tuple(s,span, knotvec, p), basis));
//                return insert.second;
                return true;
            }
            
            /// Check if bspline basis has been cached and populate with relevant value
            std::pair<bool, DoubleVecVec> basisDer(const double s,
                                                   const uint span,
                                                   const DoubleVec& knotvec,
                                                   const int p,
                                                   const DerivOrder order) const
            {
//                auto it = mBasisDerMap.find(std::make_tuple(s, span, knotvec, p, order));
//                if(it != mBasisDerMap.end())
//                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVecVec{});
            }
            
            /// Cache the given span
            bool cacheBasisDer(double s,
                               const uint span,
                               const DoubleVec& knotvec,
                               const uint p,
                               const DerivOrder order,
                               const DoubleVecVec& basisder)
            {
//                std::lock_guard<std::mutex> lock(mMutex);
//                auto insert = mBasisDerMap.insert(std::make_pair(std::make_tuple(s,span,knotvec,p,order), basisder));
//                return insert.second;
                return true;
            }

            
        private:
            
            /// Disable constructor
            explicit NURBSCache() = default;
            
            /// Disable destructor
            ~NURBSCache() = default;
            
            /// Disable copy constructor
            NURBSCache(const NURBSCache& n)  = default;
            
            /// Disable copy assignment
            NURBSCache& operator=(const NURBSCache& n) = default;
            
            /// Disable move constructor
            NURBSCache(NURBSCache&& n) = default;
            
            /// Disable move assignment
            NURBSCache& operator=(NURBSCache&& n) = default;
            
            /// The mutex for the singleton
            mutable std::mutex mMutex;
            
            /// Span cache
            std::map<std::tuple<double, DoubleVec, uint>, uint> mSpanMap;
            
            /// Basis cache
            std::map<std::tuple<double, uint, DoubleVec, uint>, DoubleVec> mBasisMap;
            
            /// Basis derivatives cache
            std::map<std::tuple<double, uint, DoubleVec, uint, DerivOrder>, DoubleVecVec> mBasisDerMap;
            
            /// Bernstein basis cache
            std::map<std::tuple<double, uint>, DoubleVec> mBernsteinBasisMap;
            
            /// Bernstein derivative basis cache
            std::map<std::tuple<double, uint>, DoubleVec> mBernsteinBasisDerivMap;
            
        };
        
    }
}
#endif
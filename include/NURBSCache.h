#ifndef NURBS_CACHE_H
#define NURBS_CACHE_H

#include "base.h"
#include "Point3D.h"

#include <map>
#include <tuple>
#include <utility>
#include <mutex>
#include <unordered_map>


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
            
            bool cacheJacobDet(const uint ielem,
                                     const double u,
                                     const double v,
                                     const double jdet)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mJacobDetMap.insert(std::make_pair(std::make_tuple(ielem, u, v), jdet));
                return insert.second;
            }
            
            std::pair<bool, double> jacobDet(const uint ielem,
                                             const double u,
                                             const double v)
            {
                auto it = mJacobDetMap.find(std::make_tuple(ielem, u, v));
                if(it != mJacobDetMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, 0.0);
            }
            
            bool cacheJacob(const uint ielem,
                            const double u,
                            const double v,
                            const DoubleVecVec& jacob)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mJacobMap.insert(std::make_pair(std::make_tuple(ielem, u, v), jacob));
                return insert.second;
            }
            
            std::pair<bool, DoubleVecVec> jacob(const uint ielem,
                                                const double u,
                                                const double v)
            {
                auto it = mJacobMap.find(std::make_tuple(ielem, u, v));
                if(it != mJacobMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVecVec{});
            }
            
            bool cachePhysicalCoord(const uint ielem,
                                    const double u,
                                    const double v,
                                    const Point3D& p)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mPhysicalCoordMap.insert(std::make_pair(std::make_tuple(ielem, u, v), p));
                return insert.second;
            }
            
            std::pair<bool, Point3D> physicalCoord(const uint ielem,
                                                   const double u,
                                                   const double v)
            {
                auto it = mPhysicalCoordMap.find(std::make_tuple(ielem, u, v));
                if(it != mPhysicalCoordMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, Point3D());
            }
            
            bool cacheTangentDS(const uint ielem,
                                    const double u,
                                    const double v,
                                    const Point3D& p)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mTangentDSMap.insert(std::make_pair(std::make_tuple(ielem, u, v), p));
                return insert.second;
            }
            
            std::pair<bool, Point3D> tangentDS(const uint ielem,
                                                   const double u,
                                                   const double v)
            {
                auto it = mTangentDSMap.find(std::make_tuple(ielem, u, v));
                if(it != mTangentDSMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, Point3D());
            }
            
            bool cacheTangentDT(const uint ielem,
                                const double u,
                                const double v,
                                const Point3D& p)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mTangentDTMap.insert(std::make_pair(std::make_tuple(ielem, u, v), p));
                return insert.second;
            }
            
            std::pair<bool, Point3D> tangentDT(const uint ielem,
                                               const double u,
                                               const double v)
            {
                auto it = mTangentDTMap.find(std::make_tuple(ielem, u, v));
                if(it != mTangentDTMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, Point3D());
            }
            
            
            bool cacheVectorBasis(const uint ielem,
                                    const double u,
                                    const double v,
                                    const DoubleVecVec& basis)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mVectorBasisMap.insert(std::make_pair(std::make_tuple(ielem, u, v), basis));
                return insert.second;
            }
            
            std::pair<bool, DoubleVecVec> vectorBasis(const uint ielem,
                                                      const double u,
                                                      const double v)
            {
                auto it = mVectorBasisMap.find(std::make_tuple(ielem, u, v));
                if(it != mVectorBasisMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVecVec{});
            }
            
            bool cacheLocalVectorBasis(const uint ielem,
                                       const double u,
                                       const double v,
                                       const DoubleVecVec& basis)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mLocalVectorBasisMap.insert(std::make_pair(std::make_tuple(ielem, u, v), basis));
                return insert.second;
            }
            
            std::pair<bool, DoubleVecVec> localVectorBasis(const uint ielem,
                                                           const double u,
                                                           const double v)
            {
                auto it = mLocalVectorBasisMap.find(std::make_tuple(ielem, u, v));
                if(it != mLocalVectorBasisMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVecVec{});
            }
            
            bool cacheVectorBasisDer(const uint ielem,
                                     const double u,
                                     const double v,
                                     const DerivType deriv,
                                     const DoubleVecVec& basisder)
            {
                std::lock_guard<std::mutex> lock(mMutex);
                auto insert = mVectorBasisDerMap.insert(std::make_pair(std::make_tuple(ielem, u, v, deriv), basisder));
                return insert.second;
            }
            
            std::pair<bool, DoubleVecVec> vectorBasisDer(const uint ielem,
                                                         const double u,
                                                         const double v,
                                                         const DerivType deriv)
            {
                auto it = mVectorBasisDerMap.find(std::make_tuple(ielem, u, v, deriv));
                if(it != mVectorBasisDerMap.end())
                    return std::make_pair(true, it->second);
                return std::make_pair(false, DoubleVecVec{});
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
            
            std::map<std::tuple<uint, double, double>, std::tuple<Point3D, Point3D>> mTangentVectorMap;
            
            std::map<std::tuple<uint, double, double>, double> mJacobDetMap;
            
            std::map<std::tuple<uint, double, double>, Point3D> mPhysicalCoordMap;
            
            std::map<std::tuple<uint, double, double>, DoubleVecVec> mVectorBasisMap;
            
            std::map<std::tuple<uint, double, double>, DoubleVecVec> mLocalVectorBasisMap;
            
            std::map<std::tuple<uint, double, double>, DoubleVecVec> mJacobMap;
            
            std::map<std::tuple<uint, double, double, uint>, DoubleVecVec> mVectorBasisDerMap;
            
            std::map<std::tuple<uint, double, double>, Point3D> mTangentDSMap;
            
            std::map<std::tuple<uint, double, double>, Point3D> mTangentDTMap;
            
        };
        
    }
}
#endif
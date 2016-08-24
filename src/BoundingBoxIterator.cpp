#include <utility>

#include "BoundingBoxIterator.h"
#include "Forest.h"
#include "MultiForest.h"
#include "Geometry.h"
#include "Point3D.h"
#include "NodalElement.h"

namespace nurbs {
    
    BoundingBoxIterator::BoundingBoxIterator(const MultiForest& f)
    :
    mCurrentIndex(0),
    mBoundingBoxN(f.globalDofN())
    {
        std::vector<Point3D> pdata;
        pdata.resize(f.collocPtN());
        std::vector<bool> cached(f.collocPtN(), false);
        for(uint ispace = 0; ispace < f.spaceN(); ++ispace)
        {
            for(uint dir = 0; dir < 2; ++dir)
            {
                const ParamDir pdir = ParamDirType(dir);
                for(uint icpt = 0; icpt < f.collocPtN(ispace, pdir); ++icpt)
                {
                    const int gindex = f.globalCollocI(ispace, pdir, icpt);
                    
                    if(-1 == gindex) // degenerate point
                        continue;
                    
                    if(gindex > pdata.size())
                        error("Bad collocation data in BoundingBoxIterator");
                    
                    if(cached[gindex])
                        continue;
                    
                    pdata[gindex] = f.collocPt(ispace, pdir, icpt);
                    cached[gindex] = true;
                }
            }
        }
        mPtData = pdata;
        
        // Now compute bounding boxes of basis function spans
        const double minval = std::numeric_limits<double>::lowest();
        const double maxval = std::numeric_limits<double>::max();
        
        Point3D minpt(minval, minval,minval);
        Point3D maxpt(maxval, maxval, maxval);
        
        std::vector<std::pair<Point3D, Point3D>> bbdata(f.globalDofN(),
                                                        std::make_pair(maxpt, minpt));
        
        for(uint ispace = 0; ispace < f.spaceN(); ++ispace)
        {
            for(uint ielem = 0; ielem < f.elemN(ispace); ++ielem)
            {
                const auto elem = f.bezierElement(ispace, ielem);
                const Point3D upper = elem->approxUpperBound();
                const Point3D lower = elem->approxLowerBound();
                const auto globalBasisVec = elem->signedGlobalBasisFuncI();
                for(const auto& gindex : globalBasisVec)
                {
                    if(-1 == gindex) // degenerate point
                        continue;
                    
                    if(gindex > f.globalDofN())
                        error("Bad basis function index in bounding box iterator");
                    
                    auto& currentmin = bbdata[gindex].first;
                    auto& currentmax = bbdata[gindex].second;
                    currentmin = min(currentmin, lower);
                    currentmax = max(currentmax, upper);
                }
            }
            
        }
        mBBData = bbdata;
    }
    
    BoundingBoxIterator::BoundingBoxIterator(const Forest& f)
    :
    mCurrentIndex(0),
    mBoundingBoxN(f.globalDofN())
    {
        // First compute and store collocation data
        
        std::vector<Point3D> pdata;
        pdata.resize(f.collocPtN());
        std::vector<bool> cached(f.collocPtN(), false);
        for(uint ispace = 0; ispace < f.spaceN(); ++ispace)
        {
            for(uint icpt = 0; icpt < f.collocPtN(ispace); ++icpt)
            {
                const uint gindex = f.globalCollocI(ispace, icpt);
                if(gindex > pdata.size())
                    error("Bad collocation data in BoundingBoxIterator");
                if(cached[gindex])
                    continue;
                pdata[gindex] = f.collocPt(ispace, icpt);
                cached[gindex] = true;
            }
        }
        
        mPtData.resize(f.collocPtN());
        for(uint i = 0; i < pdata.size(); ++i)
            mPtData[i] = pdata[i];
        
        // Now compute bounding boxes of basis function spans
        const double minval = std::numeric_limits<double>::lowest();
        const double maxval = std::numeric_limits<double>::max();
        
        Point3D minpt(minval, minval,minval);
        Point3D maxpt(maxval, maxval, maxval);
        
        std::vector<std::pair<Point3D, Point3D>> bbdata(f.globalDofN(),
                                                        std::make_pair(maxpt, minpt));
        
        for(uint ispace = 0; ispace < f.spaceN(); ++ispace) {
            for(uint ielem = 0; ielem < f.elemN(ispace); ++ielem) {
                const auto elem = f.element(ispace, ielem);
                const Point3D upper = elem->approxUpperBound();
                const Point3D lower = elem->approxLowerBound();
                const auto globalBasisVec = elem->globalBasisFuncI();
                for(const auto& gindex : globalBasisVec) {
                    if(gindex > f.globalDofN())
                        error("Bad basis function index in bounding box iterator");
                    auto& currentmin = bbdata[gindex].first;
                    auto& currentmax = bbdata[gindex].second;
                    currentmin = min(currentmin, lower);
                    currentmax = max(currentmax, upper);
                }
            }
        }
        mBBData = bbdata;
    }
    
    BoundingBoxIterator& BoundingBoxIterator::operator++()
    {
        ++mCurrentIndex;
        return *this;
    }
    
    bool BoundingBoxIterator::isDone() const
    {
        return mCurrentIndex >= boundingBoxN();
    }
    
    Point3D BoundingBoxIterator::currentPt() const
    {
        return mPtData.at(currentIndex());
    }
    
    std::pair<Point3D, Point3D> BoundingBoxIterator::currentBoundingBox() const
    {
        return mBBData.at(currentIndex());
    }
    
    Point3D BoundingBoxIterator::currentLowerBound() const
    {
        return currentBoundingBox().first;
    }
    
    /// Get current upper bound
    Point3D BoundingBoxIterator::currentUpperBound() const
    {
        return currentBoundingBox().second;
    }
    
    
    
}
#ifndef NURBS_OUTPUT_VTK_H
#define NURBS_OUTPUT_VTK_H

#include <string>
#include <vector>


#include "base.h"


namespace nurbs {
    
    /// Forward declarations
    class Forest;
    class Point3D;
    
    /// A class responsible for creating output files for the
    /// 'mesh' objects defined by Geometry, Forest and Multiforest
    
    /// Responsible for outputing both geometry and solution fields
    
    class OutputVTK {
        
    public:
        
        /// Default constructor
        OutputVTK() : OutputVTK("unnamed_output") {}
        
        /// Construct with filename
        OutputVTK(const std::string& f,
                  const uint nsample = DEFAULT_NGRID_PTS)
        :
        mFilename(f),
        mSamplePtN(nsample) {}
        
        /// output the geometry of a forest to VTK
        void outputGeometry(const Forest& f) const;
        
        /// Output bounding boxes (primarily for checking bounding boxes
        /// for defining support of basis functions for H-matrices
        void outputBoundingBoxSet(const std::vector<std::pair<Point3D, Point3D>>& bdata) const;
        
        void outputComplexAnalysisField(const Forest& f,
                                        const std::string& fieldname,
                                        const std::vector<std::complex<double>>& soln) const;
        
        
        /// Sample point number setter
        void setSamplePtN(const uint n)
        {
            if(n < 1)
                error("Must specifiy non-zero sample points for output");
            mSamplePtN = n;
        }
        
        /// Filename getter
        const std::string& filename() const
        {
            return mFilename;
        }
        
        /// Sample point number getter
        uint samplePtN() const { return mSamplePtN; }
        
    private:
        
        /// Filename we write to
        const std::string mFilename;
        
        /// Number of sample points
        uint mSamplePtN;
    };
    

    
    
}
#endif
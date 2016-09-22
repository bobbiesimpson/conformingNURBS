#ifndef NURBS_OUTPUT_VTK_H
#define NURBS_OUTPUT_VTK_H

#include <string>
#include <vector>


#include "base.h"


namespace nurbs {
    
    /// Forward declarations
    class Forest;
    class MultiForest;
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
        
        /// Write a complex nodal field to a vtu file
        void outputComplexNodalField(const Forest& f,
                                     const std::string& fieldname,
                                     const std::vector<std::complex<double>>& soln) const;
        
        /// Output a complex vector field with nedelec elements
        void outputComplexVectorFieldNedelec(const MultiForest& f,
                                                        const std::string& fieldname,
                                             const std::vector<std::complex<double>>& soln) const;
        
        /// Write a complex vector field to a vtu file
        void outputComplexVectorField(const MultiForest& f,
                                      const std::string& fieldname,
                                      const std::vector<std::complex<double>>& soln) const;
                                      
        
        /// Write the Mie series solution of the surface current given a wavenumber k
        /// It is assumed the wave is travelling in the x-direction and polarised in the
        /// z-direction (i.e solution in Harrington book)
        void outputAnalyticalMieComplexVectorField(const MultiForest& f,
                                                   const std::string& fieldname,
                                                   const double k) const;
    
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
        
        /// Filename setter
        void setFilename(const std::string& newfile)
        {
            mFilename = newfile;
        }
        
        /// Sample point number getter
        uint samplePtN() const { return mSamplePtN; }

    private:
        
        /// Filename we write to
        std::string mFilename;
        
        /// Number of sample points
        uint mSamplePtN;
    };
    
    /// analytical mie surface current expression
    std::vector<std::complex<double> > mieSurfaceCurrent(const double k,
                                                         const double theta,
                                                         const double phi);
    
}
#endif
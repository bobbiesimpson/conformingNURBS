#include "OutputVTK.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"
#include <vtkHexahedron.h>

#include "Forest.h"
#include "IElem.h"

namespace nurbs {

    void OutputVTK::outputGeometry(const Forest& f) const
    {
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1; // number of cells in each parametric direction
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        
        uint sample_offset = 0;
        
        for(uint i = 0; i < f.elemN(); ++i) {
            const auto e = f.element(i);

            uint count = 0;
            for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                const ParamPt samplept = isamplept.getCurrentPt();
                const Point3D phys_coord = e->eval(samplept.s, samplept.t);
                points->InsertPoint(sample_offset + count, phys_coord.data());
                ++count;
            }
            for( uint t = 0; t < ncell; ++t ) {
                for( uint s = 0; s < ncell; ++s ) {
                    vtkSmartPointer< vtkCell > cell = vtkQuad::New();
                    cell->GetPointIds()->SetId(0, sample_offset + t * nsample + s );
                    cell->GetPointIds()->SetId(1, sample_offset + t * nsample + s + 1 );
                    cell->GetPointIds()->SetId(2, sample_offset + ( t + 1 ) * nsample + s + 1 );
                    cell->GetPointIds()->SetId(3, sample_offset + ( t + 1 ) * nsample + s );
                    grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds() );
                }
            }
            sample_offset += nsample * nsample;
        }
        grid->SetPoints(points);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_geometry.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        if(!writer->Write())
            error( "Cannot write vtk file" );
    }
    
    
    
    void OutputVTK::outputBoundingBoxSet(const std::vector<std::pair<Point3D, Point3D>>& bdata) const
    {
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        
        uint counter = 0;
        
        // Now loop over all the bounding boxes and add the vertices to the grid
        for(const auto& box : bdata) {
            const auto& v0 = box.first;
            const auto& v6 = box.second;
            const Point3D v1(v6[0],v0[1],v0[2]);
            const Point3D v2(v6[0], v6[1], v0[2]);
            const Point3D v3(v0[0], v6[1], v0[2]);
            const Point3D v4(v0[0], v0[1], v6[2]);
            const Point3D v5(v6[0], v0[1], v6[2]);
            const Point3D v7(v0[0], v6[1], v6[2]);
            const std::vector<Point3D> pvec{v0,v1,v2,v3,v4,v5,v6,v7};
            
            vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
            for(auto i = 0; i < pvec.size(); ++i) {
                const auto& p = pvec[i];
                points->InsertNextPoint(p.data());
                hex->GetPointIds()->SetId(i, counter++);
            }
            grid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
            
        }
        grid->SetPoints(points);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_bbdata.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        if(!writer->Write())
            error( "Cannot write vtk file for bounding box data" );
    }
    
        /// specialisation for complex types
        void OutputVTK::outputComplexAnalysisField(const Forest& f,
                                                   const std::vector<std::complex<double>>& soln) const
        {
            const uint nsample = samplePtN();   // number of sample points in each parametric direction
            const uint ncell = nsample - 1; // number of cells in each parametric direction
            vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
            vtkSmartPointer<vtkPoints> points = vtkPoints::New();
            vtkSmartPointer<vtkDoubleArray> solndata = vtkDoubleArray::New();
            solndata->SetNumberOfComponents(3);
            solndata->SetName("acoustic pressure");
            solndata->SetComponentName(0, "real");
            solndata->SetComponentName(1, "imag");
            solndata->SetComponentName(2, "abs");
    
            uint sample_offset = 0;
    
            for(uint i = 0; i < f.elemN(); ++i) {
                const auto e = f.element(i);
    
                const auto gbasisivec = e->globalBasisFuncI();
                uint count = 0;
                for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept) {
                    const ParamPt samplept = isamplept.getCurrentPt();
                    const Point3D phys_coord = e->eval(samplept.s, samplept.t);
                    points->InsertPoint(sample_offset + count, phys_coord.data());
    
                    // now interpolate solution
                    const auto basisvec = e->basis(samplept.s, samplept.t);
    
                    std::complex<double> val;
                    for(uint ibasis = 0; ibasis < basisvec.size(); ++ibasis)
                        val += soln[gbasisivec[ibasis]] * basisvec[ibasis];
                    solndata->InsertComponent(sample_offset + count, 0, val.real());
                    solndata->InsertComponent(sample_offset + count, 1, val.imag());
                    solndata->InsertComponent(sample_offset + count, 2, std::abs(val));
                    ++count;
                }
                for( uint t = 0; t < ncell; ++t ) {
                    for( uint s = 0; s < ncell; ++s ) {
                        vtkSmartPointer< vtkCell > cell = vtkQuad::New();
                        cell->GetPointIds()->SetId(0, sample_offset + t * nsample + s );
                        cell->GetPointIds()->SetId(1, sample_offset + t * nsample + s + 1 );
                        cell->GetPointIds()->SetId(2, sample_offset + ( t + 1 ) * nsample + s + 1 );
                        cell->GetPointIds()->SetId(3, sample_offset + ( t + 1 ) * nsample + s );
                        grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds() );
                    }
                }
                sample_offset += nsample * nsample;
            }
            grid->SetPoints(points);
            grid->GetPointData()->AddArray(solndata);
            vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
            const std::string fname = filename() + "_complex_soln.vtu";
            writer->SetFileName(fname.c_str());
            writer->SetInputData(grid);
            if(!writer->Write())
                error( "Cannot write vtk file" );
        }
}
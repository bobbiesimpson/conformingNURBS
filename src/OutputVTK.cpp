#include "OutputVTK.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"
#include <vtkHexahedron.h>
#include <complex>
#include "Forest.h"
#include "MultiForest.h"
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
            const auto e = f.bezierElement(i);
            
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
    void OutputVTK::outputComplexNodalField(const Forest& f,
                                            const std::string& fieldname,
                                            const std::vector<std::complex<double>>& soln) const
    {
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1; // number of cells in each parametric direction
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        vtkSmartPointer<vtkDoubleArray> solndata = vtkDoubleArray::New();
        solndata->SetNumberOfComponents(3);
        solndata->SetName(fieldname.c_str());
        solndata->SetComponentName(0, "real");
        solndata->SetComponentName(1, "imag");
        solndata->SetComponentName(2, "abs");
        
        uint sample_offset = 0;
        
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto e = f.bezierElement(i);
            
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
    
    void OutputVTK::outputComplexVectorField(const MultiForest& f,
                                             const std::string& fieldname,
                                             const std::vector<std::complex<double>>& soln) const
    {
        // number of sample points and cells in each parametric direction
        const uint nsample = samplePtN();
        const uint ncell = nsample - 1;
        
        assert(f.globalDofN() == soln.size());
        
        // create the vtk grid, points array and solution array
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        vtkSmartPointer<vtkDoubleArray> vtk_realsoln = vtkDoubleArray::New();
        vtkSmartPointer<vtkDoubleArray> vtk_imagsoln = vtkDoubleArray::New();
        vtkSmartPointer<vtkDoubleArray> vtk_abssoln = vtkDoubleArray::New();
        
        vtk_realsoln->SetNumberOfComponents(3);
        std::string name = fieldname + "_real";
        vtk_realsoln->SetName(name.c_str());
        
        vtk_imagsoln->SetNumberOfComponents(3);
        name = fieldname + "_imag";
        vtk_imagsoln->SetName(name.c_str());
        
        vtk_abssoln->SetNumberOfComponents(1);
        name = fieldname + "_abs";
        vtk_abssoln->SetName(name.c_str());
        
        // now loop over elements and sample solution
        uint sample_offset = 0;
        
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto el = f.bezierElement(i);
            const auto gbasisivec = el->globalBasisFuncI();
            
            uint count = 0;
            
            for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                const ParamPt samplept = isamplept.getCurrentPt();
            
                const Point3D phys_coord = el->eval(samplept.s, samplept.t);
                
                points->InsertPoint(sample_offset + count, phys_coord.data());
                const auto basis = el->basis(samplept.s, samplept.t);
                
                std::vector<std::complex<double>> val(3);
                for(size_t ibasis = 0; ibasis < basis.size(); ++ibasis)
                    for(unsigned i = 0; i < 3; ++i)
                        val[i] += soln[gbasisivec[ibasis]] * basis[ibasis][i];
                
                // now put this complex vector into the vtk arrays
                double absval = 0.0;
                for(unsigned i = 0; i < 3; ++i)
                {
                    const double re = val[i].real();
                    const double im = val[i].imag();
                    vtk_realsoln->InsertComponent(sample_offset + count, i, re);
                    vtk_imagsoln->InsertComponent(sample_offset + count, i, im);
                    absval += re * re + im * im;
                }
                
                // and finally insert absolute value
                vtk_abssoln->vtkDataArray::InsertComponent(sample_offset + count, 0, std::sqrt(absval));
                
                ++count;
            }
            
            // create the cell connectivity
            for( uint t = 0; t < ncell; ++t )
            {
                for( uint s = 0; s < ncell; ++s )
                {
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
        
        // and finally add the points and solutions to the grid and write!
        grid->SetPoints(points);
        grid->GetPointData()->AddArray(vtk_realsoln);
        grid->GetPointData()->AddArray(vtk_imagsoln);
        grid->GetPointData()->AddArray(vtk_abssoln);
        
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_complexvector.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        
        if(!writer->Write())
            error( "Cannot write vtk file" );
    }

}
#include "OutputVTK.h"
#include "Forest.h"
#include "MultiForest.h"
#include "IElem.h"
#include "algebra.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"
#include "vtkHexahedron.h"

#include <complex>
#include <thread>
#include <atomic>
#include <vector>

#include <boost/math/special_functions.hpp>
#include <gsl/gsl_sf_legendre.h>

namespace nurbs {
    
    void OutputVTK::outputGeometry(const Forest& f) const
    {
        const uint nsample = samplePtN();   // number of sample points in each parametric direction
        const uint ncell = nsample - 1; // number of cells in each parametric direction
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        
        uint sample_offset = 0;
        
        std::vector<uint> evec{8,17};
//        for(uint i = 0; i < 18; ++i)
//            evec.push_back(i);
        
//        for(uint i = 0; i < f.elemN(); ++i)
        for(const auto& i : evec)
        {
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
        for(size_t i = 0; i < bdata.size(); ++i) {
            const auto& box = bdata[i];
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
                
                // temp code
//                const Point3D p(0.0, 1.0, 0.0);
//                const Point3D k(5.0, 0.0, 0.0);
//                
//                std::vector<std::complex<double>> result{p[0], p[1], p[2]};
//                const auto wave = std::exp(-std::complex<double>(0.0, dot(k, phys_coord)));
//                for(auto& r : result)
//                    r *= wave;
//                
//                const auto& exact_val = result;
//                
//                const auto& t1 = e->tangent(samplept.s, samplept.t, nurbs::ParamDir::S);
//                const auto& t2 = e->tangent(samplept.s, samplept.t, nurbs::ParamDir::T);
//                const double jpiola = nurbs::cross(t1, t2).length();
//                
//                DoubleVecVec j;
//                j.push_back(t1.asVec());
//                j.push_back(t2.asVec());
//                j.push_back(
//                            {
//                                1.0/jpiola * (j[0][1] * j[1][2] - j[0][2] * j[1][1]),
//                                1.0/jpiola * (j[0][2] * j[1][0] - j[0][0] * j[1][2]),
//                                1.0/jpiola * (j[0][0] * j[1][1] - j[0][1] * j[1][0])
//                            });
//                auto jinv = inv3x3Mat(j);
//                for(auto& row : jinv)
//                    row.erase(row.begin() + 2);
//                
//                
//                std::vector<std::complex<double>> fexact(2);
//                std::vector<std::complex<double>> f_h(2);
//                for(size_t i = 0; i < 2; ++i)
//                    for(size_t j = 0; j < 3; ++j)
//                        fexact[i] += jpiola * jinv[j][i] * exact_val[j];
//                
//                const auto basisvec = e->localBasis(samplept.s, samplept.t);
//                
//                std::complex<double> val;
//                for(uint ibasis = 0; ibasis < basisvec.size(); ++ibasis)
//                    val += soln[gbasisivec[ibasis]] * basisvec[ibasis];
                
                // end temp code
                
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
    
    void OutputVTK::outputComplexVectorFieldNedelec(const MultiForest& f,
                                                    const std::string& fieldname,
                                                    const std::vector<std::complex<double>>& soln) const
    {
        // number of sample points and cells in each parametric direction
        const uint nsample = samplePtN();
        const uint ncell = nsample - 1;
        
        assert(f.globalNedelecDofN() == soln.size());
        
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
        const double degenerate_shift = 1.0e-6; // tolerance to shift sample points away from degenerate edges
        

        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto el = f.nedelecElement(i);
            //            const auto parent_el = el->parent();
            const auto gbasisivec = el->signedGlobalBasisFuncI();
            //            if(el->degenerate())
            //                continue;
            
            uint count = 0;
            
            for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                ParamPt samplept = isamplept.getCurrentPt();
                if(el->degenerate())
                {
                    samplept.s *= (1.0 - degenerate_shift);
                    samplept.t *= (1.0 - degenerate_shift);
                }
                
                const Point3D phys_coord = el->eval(samplept.s, samplept.t);
                
                points->InsertPoint(sample_offset + count, phys_coord.data());
                const auto basis = el->basis(samplept.s, samplept.t);
                
                std::vector<std::complex<double>> val(3, 0.0);
                for(size_t ibasis = 0; ibasis < basis.size(); ++ibasis)
                {
                    if(-1 == gbasisivec[ibasis]) // degenerate point
                        continue;
                    
                    for(unsigned i = 0; i < 3; ++i)
                        val[i] += soln[gbasisivec[ibasis]] * basis[ibasis][i];
                }
                
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
        const double degenerate_shift = 1.0e-6; // tolerance to shift sample points away from degenerate edges
        
//        std::vector<uint> evec{1};
        
//        for(const auto& i : evec)
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto el = f.bezierElement(i);
//            const auto parent_el = el->parent();
            const auto gbasisivec = el->signedGlobalBasisFuncI();
//            if(el->degenerate())
//                continue;
            
            uint count = 0;
            
            for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                ParamPt samplept = isamplept.getCurrentPt();
                if(el->degenerate())
                {
                    samplept.s *= (1.0 - degenerate_shift);
                    samplept.t *= (1.0 - degenerate_shift);
                }
            
                const Point3D phys_coord = el->eval(samplept.s, samplept.t);
                
                points->InsertPoint(sample_offset + count, phys_coord.data());
                const auto basis = el->basis(samplept.s, samplept.t);
                
                std::vector<std::complex<double>> val(3, 0.0);
                for(size_t ibasis = 0; ibasis < basis.size(); ++ibasis)
                {
                    if(-1 == gbasisivec[ibasis]) // degenerate point
                        continue;
                    
                    for(unsigned i = 0; i < 3; ++i)
                        val[i] += soln[gbasisivec[ibasis]] * basis[ibasis][i];
                }
                
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
    
    void OutputVTK::outputAnalyticalMieComplexVectorField(const MultiForest& f,
                                                          const std::string& fieldname,
                                                          const double k) const
    {
        // number of sample points and cells in each parametric direction
        const uint nsample = samplePtN();
        const uint ncell = nsample - 1;
        
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
        const double degenerate_shift = 1.0e-6; // tolerance to shift sample points away from degenerate edges
        
        for(uint i = 0; i < f.elemN(); ++i)
        {
            const auto el = f.bezierElement(i);
            //            const auto parent_el = el->parent();
            const auto gbasisivec = el->signedGlobalBasisFuncI();
            //            if(el->degenerate())
            //                continue;
            
            uint count = 0;
            
            for(ISamplePt isamplept(nsample); !isamplept.isDone(); ++isamplept)
            {
                ParamPt samplept = isamplept.getCurrentPt();
                if(el->degenerate())
                {
                    samplept.s *= (1.0 - degenerate_shift);
                    samplept.t *= (1.0 - degenerate_shift);
                }
                
                const Point3D phys_coord = el->eval(samplept.s, samplept.t);
                
                // BEWARE: HARDCODED ANALYTICAL SOLUTION!!
                const auto& p = phys_coord;
                double r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                double theta = acos(p[0] / r);
                double phi = atan2(p[2], p[1]);
                
                const auto j_vec = mieSurfaceCurrent(k, theta, phi);
                
                points->InsertPoint(sample_offset + count, phys_coord.data());
                
                // now put this complex vector into the vtk arrays
                double absval = 0.0;
                for(unsigned i = 0; i < 3; ++i)
                {
                    const double re = j_vec[i].real();
                    const double im = j_vec[i].imag();
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
    
    void OutputVTK::ouputQuadratureData(const std::string& fieldname,
                                        const std::vector<double>& rawdata,
                                        const std::vector<nurbs::GPt2D>& pts,
                                        const uint ns,
                                        const uint nt) const
    {
        // assume that the raw data is stored as j * ns + i
        const uint ncell_s = ns - 1;
        const uint ncell_t = nt - 1;
        
        // create the vtk grid, points array and solution array
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkUnstructuredGrid::New();
        vtkSmartPointer<vtkPoints> points = vtkPoints::New();
        vtkSmartPointer<vtkDoubleArray> vtk_soln = vtkDoubleArray::New();
        
        vtk_soln->SetNumberOfComponents(1);
        vtk_soln->SetName(fieldname.c_str());

        // now loop over elements and sample solution
        uint currentpt_index = 0;
        
        for(uint j = 0; j < nt; ++j)
        {
            for(uint i = 0; i < ns; ++i)
            {
                const uint index = j * ns + i;
                nurbs::Point3D point(pts[index].s, pts[index].t, 0.0);
                points->InsertPoint(currentpt_index, point.data());
                vtk_soln->vtkDataArray::InsertComponent(currentpt_index, 0, rawdata[currentpt_index]);
                ++currentpt_index;
            }
        }
        
        // and construct cells
        // create the cell connectivity
        for(uint t = 0; t < ncell_t; ++t)
        {
            for(uint s = 0; s < ncell_s; ++s)
            {
                vtkSmartPointer<vtkCell> cell = vtkQuad::New();
                cell->GetPointIds()->SetId(0, t * ns + s);
                cell->GetPointIds()->SetId(1, t * ns + s + 1 );
                cell->GetPointIds()->SetId(2, ( t + 1 ) * ns + s + 1);
                cell->GetPointIds()->SetId(3, ( t + 1 ) * ns + s );
                grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
            }
        }
        
        // and finally add the points and solutions to the grid and write!
        grid->SetPoints(points);
        grid->GetPointData()->AddArray(vtk_soln);
        
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkXMLUnstructuredGridWriter::New();
        const std::string fname = filename() + "_integrand.vtu";
        writer->SetFileName(fname.c_str());
        writer->SetInputData(grid);
        
        if(!writer->Write())
            error( "Cannot write vtk file" );
        
    }
    
    double OutputVTK::computeRCS(const MultiForest& f,
                                 const nurbs::Point3D& sample,
                                 const double k,
                                 const nurbs::Point3D& rhat,
                                 const double mu,
                                 const double omega,
                                 const std::vector<std::complex<double>>& soln) const
    {

        std::vector<std::complex<double>> result(3);
        
        const int nthreads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        auto boundsvec = bounds(nthreads, f.elemN());
        
        for(int i = 0; i < nthreads; ++i)
        {
            threads.push_back(std::thread(&OutputVTK::rcsThreadWorkerFunction,
                                          this,
                                          std::ref(f),
                                          sample,
                                          k,
                                          rhat,
                                          mu,
                                          omega,
                                          soln,
                                          boundsvec[i],
                                          boundsvec[i+1],
                                          std::ref(result)));
            
        }
        for(auto &t : threads)
            t.join();
        
        const double r = sample.length();
        ComplexDouble j(0.0,1.0);
        ComplexDouble coeff1 = - (j * omega * mu) / (4.0 * PI);
        ComplexDouble coeff2 = std::exp(-j*k*r) / r;
        
        std::vector<ComplexDouble> Er(3,0);
        double rcs = 0.0;
        for(unsigned i = 0; i < 3; ++i)
        {
            Er[i] = coeff1 * coeff2 * result[i];
            rcs += real(Er[i] *real(Er[i]) + imag(Er[i]) * imag(Er[i]));
        }
        rcs = 4.0 * PI * r * r * rcs;
        return rcs;
        

    }
    
    void OutputVTK::rcsThreadWorkerFunction(const MultiForest& f,
                                            const nurbs::Point3D& sample,
                                            const double k,
                                            const nurbs::Point3D& rhat,
                                            const double mu,
                                            const double omega,
                                            const std::vector<std::complex<double>>& soln,
                                            const uint start,
                                            const uint end,
                                            std::vector<std::complex<double>>& result) const
    {

        const ComplexDouble j(0.0,1.0);

        std::vector<ComplexDouble> integral(3,0);
        for(uint ielem = start; ielem < end; ++ielem)
        {
            const auto el = f.bezierElement(ielem);
            const auto gbasisivec = el->signedGlobalBasisFuncI();
            
            // Loop over gauss points
            for(nurbs::IElemIntegrate igpt(el->equalIntegrationOrder()); !igpt.isDone(); ++igpt)
            {
                const auto gpt = igpt.get();
                const auto w = igpt.getWeight();
                const auto jdet = el->jacDet(gpt);
                const auto basis = el->basis(gpt.s, gpt.t);
                const auto gpphys_coord = el->eval(gpt.s, gpt.t);
                
                std::vector<ComplexDouble> Jgpt(3,0);
                
                for(size_t ibasis = 0; ibasis < basis.size(); ++ibasis)
                {
                    if(-1 == gbasisivec[ibasis]) // degenerate point
                        continue;
                    
                    for(unsigned i = 0; i < 3; ++i)
                        Jgpt[i] += soln[gbasisivec[ibasis]] * basis[ibasis][i];
                }
                
                double rgprhat = 0.0;
                for(size_t i = 0; i < 3; ++i)
                    rgprhat += gpphys_coord[i]*rhat[i];
                
                ComplexDouble ejkrr = std::exp(j*k*rgprhat);
                
                nurbs::Point3D jreal(Jgpt[0].real(), Jgpt[1].real(), Jgpt[2].real());
                nurbs::Point3D jimag(Jgpt[0].imag(), Jgpt[1].imag(), Jgpt[2].imag());
                
                const auto cross_real = nurbs::cross(rhat, nurbs::cross(jreal, rhat));
                const auto cross_imag = nurbs::cross(rhat, nurbs::cross(jimag, rhat));
                
                for(size_t i = 0; i < 3; ++i)
                    Jgpt[i] = std::complex<double>(cross_real[i], cross_imag[i]);
                
                for(size_t i = 0; i < 3; ++i)
                    integral[i] += Jgpt[i] * ejkrr * w * jdet;
            }
        }

        std::lock_guard<std::mutex> block(mutex());
        for(int i = 0; i < 3; ++i)
            result[i] += integral[i];
    }
    
    std::vector<std::complex<double> > mieSurfaceCurrent( const double k,
                                                         const double theta,
                                                         const double phi)
    {
        const double PI = atan( 1.0 ) * 4.0;
        
        double error = 100.0;
        const double tol = 1.0e-10;
        std::complex<double> jtheta, jphi;
        double curr_mag = 0.0;
        double prev_mag;
        unsigned n = 1;
        const double ct = std::cos(theta);
        const double st = std::sin(theta);
        
        while(error > tol) {
            const auto i = std::complex<double>(0.0, 1.0);
            const auto an = 1.0 / k * std::pow(i, -n) * (2.0 * n + 1.0) / (n * (n + 1.0));
            const auto h = boost::math::sph_hankel_2(n, k);
            const auto Hbar = k * h;
            const auto h_d = n / k * boost::math::sph_hankel_2(n,k) - boost::math::sph_hankel_2(n + 1, k);
            const auto Hbar_d = k * h_d + h;
            
            
            if(std::abs(theta) < 1.0e-7) {
                const double tri_n = 0.5 * (n * n + n);
                jtheta += std::cos(phi) * an * ( tri_n / Hbar_d
                                                + (i * -tri_n) / Hbar );
                jphi += std::sin(phi) * an * (-tri_n /  Hbar_d
                                              - tri_n / (i * Hbar) );
            }
            else if(std::abs(theta - PI) < 1.0e-7) {
                const double tri_n = 0.5 * (n * n + n) * std::pow(-1.0, n);
                jtheta += std::cos(phi) * an * ( tri_n / Hbar_d
                                                + (i * tri_n) / Hbar );
                jphi += std::sin(phi) * an * (tri_n /  Hbar_d
                                              - tri_n / (i * Hbar) );
            }
            else {
                double l_d[n+1];
                double l[n+1];
                gsl_sf_legendre_Plm_deriv_array(n, 1, ct, l, l_d);
                const auto P = l[n-1];
                const auto Pd = l_d[n-1];
                
                
                jtheta += std::cos(phi) * an * ( (st * Pd) / Hbar_d
                                                + (i * P) / (st * Hbar) );
                jphi += std::sin(phi) * an * (P / (st * Hbar_d)
                                              - (st * Pd) / (i * Hbar) );
            }
            prev_mag = curr_mag;
            curr_mag = std::sqrt( std::abs(jtheta) * std::abs(jtheta)
                                 + std::abs(jphi) * std::abs(jphi) );
            error = std::abs(prev_mag - curr_mag);
            ++n;
        }
        const double cp = std::cos(phi);
        const double sp = std::sin(phi);
        
        //    return {ct * cp * jtheta - sp * jphi,
        //        ct * sp * jtheta + cp * jphi,
        //        -st * jtheta};
        
        return {-st * jtheta,
            ct * cp * jtheta - sp * jphi,
            ct * sp * jtheta + cp * jphi
        };
    }
}
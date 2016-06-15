//
//  Geometry2D.cpp
//  nurbslib
//
//  Created by Robert Simpson on 05/12/2014.
//
//
#include <iostream>

#include "Geometry2D.h"
#include "NURBSCommon.h"
#include "BSplineSpace2D.h"
#include "Forest2D.h"
#include "base.h"

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkSmartPointer.h>

namespace nurbs {
    
    double Geometry2D::jacDet(const double s, const uint sp) const
    {
        return tangent(s, sp).length();
    }
    
    /// Get normal at given parametric coord and space index
    Point3D Geometry2D::normal(const double s, const uint sp) const
    {
        Point3D t = tangent(s, sp);
        return Point3D{t[1], -t[0], t[2]}.asNormal();
    }
    
    /// Get tangent at given parametric coord and space index
    Point3D Geometry2D::tangent(const double s, const uint sp) const
    {
        const BSplineSpace2D& bs = primalForest().space(sp);
        UIntVec indices = bs.localBasisFuncI(s);
        const uint span = bs.span(s);
        DoubleVec basis = bs.basis(s, span);
        DoubleVec ders = bs.basisDers(s, span);
        Point3D a, ader;
        double w = 0.0; double wder = 0.0;
        for(uint i = 0; i < bs.degree() + 1; ++i) {
            const Point4D& p = controlPt(indices[i]);
            a += p.asUnweighted() * basis[i];
            ader += p.asUnweighted() * ders[i];
            w += p.getWeight() * basis[i];
            wder += p.getWeight() * ders[i];
        }
        return nurbshelper::getNonRationalDeriv({a, ader}, {w, wder});
    }
    
    /// Interpolate geometry at given parametric coord and space index
    Point3D Geometry2D::eval(const double s, const uint sp) const
    {
        const BSplineSpace2D& bs = primalForest().space(sp);
        UIntVec indices = bs.localBasisFuncI(s);
        DoubleVec basis = bs.basis(s);
        assert(basis.size() == indices.size());
        Point4D p;
        for(uint i = 0; i < basis.size(); ++i)
            p += controlPt(indices[i]) * basis[i];
        return p.asCartesian();
    }
    
    /// Write to vtu file.
    void Geometry2D::writeVTPOutput(const std::string& file,
                                    const uint nsample) const
    {
        vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer< vtkCellArray > cells = vtkSmartPointer< vtkCellArray >::New();
        uint ncount = 0;
        for(uint s = 0; s < spaceN(); ++s) {
            vtkSmartPointer<vtkPolyLine> polyline = vtkSmartPointer<vtkPolyLine>::New();
            polyline->GetPointIds()->SetNumberOfIds(nsample);
            const BSplineSpace2D& bspace = primalForest().space(s);
            for(uint c = 0; c < nsample; ++c) {
                double s_param = (bspace.upperParamLimit() - bspace.lowerParamLimit()) / (nsample - 1) * c + bspace.lowerParamLimit();
                Point3D p = eval(s_param, s);
                points->InsertNextPoint(p.data());
                polyline->GetPointIds()->SetId(c, ncount);
                ++ncount;
            }
            cells->InsertNextCell(polyline);
        }
        vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
        data->SetPoints(points);
        data->SetLines(cells);
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(file.c_str());
        writer->SetInputData(data);
        writer->Write();
    }
    
    /// Load from a file stream
    bool Geometry2D::load(std::istream& ist)
    {
        char ch;
        if(ist >> ch && ch != '{') {
            ist.unget();
            ist.clear(std::ios_base::failbit);
            return false;
        }
        clear();
        Forest2D f;
        if(!(ist >> f))
            error("Cannot load forest into geometry data structure");
        mPrimalForest = f;
        mPrimalForest.setGeometry(this);
        
        if(ist >> ch && ch != '{') { // We might be reading something else
            ist.unget();
            ist.clear(std::ios_base::failbit);
            return false;
        }
        std::vector<Point4D> cp_vec;
        Point4D p;
        while(ist >> p)
            cp_vec.push_back(p);
        endOfLoop(ist, '}', "Bad control point set");
        setCPtVec(cp_vec);
        return true;
    }
    
    /// Print to a file stream
    void Geometry2D::print(std::ostream& ost) const
    {
        ost << "- Geometry object -\n\n";
        ost << "Primal forest: \n\n" << primalForest() << "\n";
        ost << "Control point set: \n\n";
        for(const auto& p : cPtVec())
            std::cout << p << "\n";
        std::cout << "\n- end geometry object -\n";
    }

    std::istream& operator>>(std::istream& ist, Geometry2D& g)
    {
        g.load(ist);
        return ist;
    }
    
    /// Output operator
    std::ostream& operator<<(std::ostream& ost, const Geometry2D& g)
    {
        g.print(ost); return ost;
    }
    
    
}

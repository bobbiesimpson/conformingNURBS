//
//  IGalerkinIntegrate.cpp
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#include "IGalerkinIntegrate.h"

#include <fstream>
#include <stdexcept>

namespace nurbs {
    
    GPt4D GPt4D::rotate(const DoubleVecVec& smat,
                        const DoubleVecVec& fmat) const
    {
        GPt4D r_pt;
        for (uint i = 0; i < 2; ++i)
            for(uint j = 0; j < 2; ++j)
                r_pt[i] += smat[i][j] * data[j];
        for (uint i = 0; i < 2; ++i)
            for(uint j = 0; j < 2; ++j)
                r_pt[i + 2] += fmat[i][j] * data[j + 2];
        return r_pt;
    }
    
    void IGalerkinIntegrate::incrementImpl()
    {
        mDidUpdateSrcI = false; // assume false to start with
        elem::IElemIntegrate& f_quad = quadrature(FIELD);
        ++f_quad;
        if(f_quad.isDone()) {
            f_quad.restart();
            elem::IElemIntegrate& s_quad = quadrature(SOURCE);
            ++s_quad; mDidUpdateSrcI = true; // we updated the source index
            if(s_quad.isDone()) {
                s_quad.restart();
                ++mSubCellI;
            }
        }
    }
    
    GPt4D IGalerkinIntegrate::getCurrentUnitIntervalPt() const
    {
        const GPt2D s = getBasePt(SOURCE);
        const GPt2D f = getBasePt(FIELD);
        const Range ps_interval = make_interval(PARENT_SPACE_INTERVAL);
        const Range ups_interval = make_interval(UNIT_PARENT_SPACE_INTERVAL);
        return convert_interval(GPt4D(s,f), ps_interval, ups_interval);
    }
    
    GPt4D IGalerkinIntegrate::getCurrentBiUnitIntervalPt() const
    {
        return GPt4D(getBasePt(SOURCE), getBasePt(FIELD));
    }
    
    void IGalerkinIntegrate::outputDatFile(const std::string& sfile,
                                           const std::string& ffile)
    {
        std::ofstream ofs(sfile);
        if(!ofs)
            throw std::runtime_error("Cannot open file for writing quadrature points");
        std::ofstream off(ffile);
        if(!off)
            throw std::runtime_error("Cannot open file for writing quadrature points");
        for(restart(); !isDone(); ++(*this)) {
            const auto gpt = get();
            ofs << gpt[0] << "\t" << gpt[1] << "\t" << getWeight() << "\n";
            off << gpt[2] << "\t" << gpt[3] << "\t" << getWeight() << "\n";
        }
    }
    void IGalerkinIntegrate::print(std::ostream& ost) const
    {
        const elem::IElemIntegrate& s = quadrature(SOURCE);
        ost << "Source quadrature index: " << s.currentIndex() + 1
        << " / " << s.pointN() << "\n";
        const elem::IElemIntegrate& f = quadrature(FIELD);
        ost << "Field quadrature index: " << f.currentIndex() + 1
        << " / " << f.pointN() << "\n";
        ost << "Sub element index: " << subCellI() + 1 << " / "
        << subCellN();
        printImpl(ost);
    }
    
    std::ostream& operator<<(std::ostream& ost, const IGalerkinIntegrate& g)
    {
        g.print(ost);
        return ost;
    }
    
    Range make_interval(IntervalType type)
    {
        switch(type) {
            case PARENT_SPACE_INTERVAL:
                return Range(-1.0, 1.0);
                break;
            case UNIT_PARENT_SPACE_INTERVAL:
                return Range(0.0, 1.0);
                break;
            default:
                std::cerr << "unrecognised interval type: aborting\n";
                exit(EXIT_FAILURE);
        }
    }
    
    GPt4D convert_interval(const GPt4D& p,
                           const Range& from_range,
                           const Range& to_range)
    {
        std::vector<double> v = p.data; // take a copy of the coordinates
        const double c1 = from_range.second - from_range.first;
        const double c2 = to_range.second - to_range.first;
        const double j = c2 / c1;
        const double c3 = (to_range.second + to_range.first) / 2.0;
        const double c4 = (from_range.second + from_range.first) / 2.0 * j;
        std::for_each(v.begin(), v.end(), [=](double& x)
                      {
                          x = j * x + (c3 - c4);
                      });
        return GPt4D(v);
    }
    
    GPt4D convert_interval_Tri(const GPt4D& p)
    {
        const double xsi = p[0];
        const double eta1 = p[1] / p[0];
        const double eta2 = p[2];
        const double eta3 = p[3] / p[2];
        double s1 = 2 * xsi - 1;
        double t1 = 2 * eta1 - 1;
        double s2 = 2 * eta2 - 1;
        double t2 = 2 * eta3 - 1;
        
        return GPt4D {s1,t1,s2,t2};
    }
    
    std::ostream& operator<<(std::ostream& ost, const GPt4D& p)
    {
        ost << "(" << p.data.at(0) << ", " << p.data.at(1) << ", " << p.data.at(2)
        << ", " << p.data.at(3) << ")";
        return ost;
    }
}
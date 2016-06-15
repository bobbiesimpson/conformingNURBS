//
//  BSplineSpace2D.cpp
//  nurbslib
//
//  Created by Robert Simpson on 05/12/2014.
//
//

#include "BSplineSpace2D.h"
#include "InputDataStructures.h"

namespace nurbs {
    
    void BSplineSpace2D::load(std::istream& ist)
    {
        clear();
        char ch;
        if(ist >> ch && ch != '{') {
            ist.unget();
            ist.clear(std::ios_base::failbit);
            return;
        }
        std::string s;
        uint d;
        InputVec<double> k_vec;
        ist >> s >> d >> k_vec;
        if(!ist )
            error("Cannot read b-spline space");
        setName(s); setDegree(d); setKnotVec(k_vec.data);
        init(); // recalculate unique knots + intervals
    }
    
    void BSplineSpace2D::printData(std::ostream& ost) const
    {
        ost << "--------- B-SPLINE SPACE DATA ----------\n\n";
        ost << "Degree: " << degree() << "\n";
        ost << "Knot vector: ";
        printVector(knotVec(), ost);
        ost << "\n";
        ost << "Unique knot vector: ";
        printVector(uniqueKnotVec(), ost);
        ost << "\n---------\n";
        ost << "\n--------------------------------------\n";
    }
}

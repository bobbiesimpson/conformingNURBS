//
//  IVertexQuadrature.cpp
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#include "IVertexQuadrature.h"

namespace nurbs
{
	GPt4D IVertexQuadrature::getImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
		const Range unit_interval = make_interval(UNIT_PARENT_SPACE_INTERVAL);
		const Range parent_interval = make_interval(PARENT_SPACE_INTERVAL);
        
		const double xsi = p[0];
		const double eta1 = p[1];
		const double eta2 = p[2];
		const double eta3 = p[3];
        
        const double aux1 = xsi * eta1;
        const double aux2 = xsi * eta2;
        const double aux3 = xsi * eta3;
        
		switch(subCellI())
        {
			case 0:
				return convert_interval(GPt4D(xsi, aux1, aux2, aux3),
										unit_interval,
										parent_interval).rotate(sourceRotMat(), fieldRotMat());
				break;
			case 1:
				return convert_interval(GPt4D(aux1, xsi, aux2, aux3),
										unit_interval,
										parent_interval).rotate(sourceRotMat(), fieldRotMat());
				break;
			case 2:
				return convert_interval(GPt4D(aux1, aux2, xsi, aux3),
										unit_interval,
										parent_interval).rotate(sourceRotMat(), fieldRotMat());
				break;
			case 3:
				return convert_interval(GPt4D(aux1, aux2, aux3, xsi),
										unit_interval,
										parent_interval).rotate(sourceRotMat(), fieldRotMat());
				break;
			default:
				std::cerr << "Reached bad subelement. Aborting.\n";
				std::exit(EXIT_FAILURE);
		}
		
	}
    
	double IVertexQuadrature::getWeightImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
		return getBaseWeight(SOURCE) * getBaseWeight(FIELD) *
        p[0] * p[0] * p[0];
	}	
}

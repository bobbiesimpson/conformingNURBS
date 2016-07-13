//
//  IEqualQuadrature.cpp
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#include "IEqualQuadrature.h"

namespace nurbs {
    
    GPt4D IEqualQuadrature::getImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
		const Range unit_interval = make_interval(UNIT_PARENT_SPACE_INTERVAL);
		const Range parent_interval = make_interval(PARENT_SPACE_INTERVAL);
		const double xsi = p[0];
		const double eta1 = p[1];
		const double eta2 = p[2];
		const double eta3 = p[3];
        
        // some additional defitiions
        const double aux1 = ( 1. - xsi ) * eta3;
        const double aux2 = ( 1. - xsi * eta1 ) * eta2;
        const double aux3 = xsi + aux1;
        const double aux4 = xsi * eta1 + aux2;
        
		switch(subCellI()) {
			case 0:
				return convert_interval(GPt4D(aux1,
                                              aux2,
                                              aux3,
                                              aux4),
                                        unit_interval,
                                        parent_interval);
				break;
			case 1:
				return convert_interval(GPt4D(aux2,
                                              aux1,
                                              aux4,
                                              aux3),
										unit_interval,
										parent_interval);
				break;
			case 2:
                return convert_interval(GPt4D(aux1,
                                              aux4,
                                              aux3,
                                              aux2),
                                        unit_interval,
                                        parent_interval);
				break;
			case 3:
                return convert_interval(GPt4D(aux2,
                                              aux3,
                                              aux4,
                                              aux1),
                                        unit_interval,
                                        parent_interval);
				break;
			case 4:
                return convert_interval(GPt4D(aux3,
                                              aux2,
                                              aux1,
                                              aux4),
                                        unit_interval,
                                        parent_interval);
				break;
			case 5:
                return convert_interval(GPt4D(aux4,
                                              aux1,
                                              aux2,
                                              aux3),
                                        unit_interval,
                                        parent_interval);
				break;
			case 6:
                return convert_interval(GPt4D(aux3,
                                              aux4,
                                              aux1,
                                              aux2),
                                        unit_interval,
                                        parent_interval);
				break;
			case 7:
                return convert_interval(GPt4D(aux4,
                                              aux3,
                                              aux2,
                                              aux1),
                                        unit_interval,
                                        parent_interval);
				break;
			default:
				std::cerr << "Reached bad subelement. Aborting.\n";
				exit(EXIT_FAILURE);
		}
	}
    
	double IEqualQuadrature::getWeightImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
		return getBaseWeight(SOURCE) * getBaseWeight(FIELD)
        * p[0] * (1.0 - p[0]) * (1.0 - p[0] * p[1]);
	}
}
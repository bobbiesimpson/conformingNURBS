//
//  IEdgeQuadrature.cpp
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#include "IEdgeQuadrature.h"

namespace nurbs
{
	GPt4D IEdgeQuadrature::getImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
		const Range unit_interval = make_interval(UNIT_PARENT_SPACE_INTERVAL);
		const Range parent_interval = make_interval(PARENT_SPACE_INTERVAL);
        
		const double xsi = p[0];
		const double eta1 = p[1];
		const double eta2 = p[2];
		const double eta3 = p[3];
        
		switch(subCellI())
        {
            case 0:
                return convert_interval(GPt4D(( 1. - xsi ) * eta3 + xsi,
                                              xsi * eta2,
                                              ( 1. - xsi ) * eta3,
                                              xsi * eta1),
                                        unit_interval,
                                        parent_interval).rotate(sourceRotMat(), fieldRotMat());
                break;
            case 1:
                return convert_interval(GPt4D((1.0 - xsi) * eta3,
                                                xsi * eta2,
                                                xsi + (1.0 - xsi) * eta3,
                                                xsi * eta1),
                                        unit_interval,
                                        parent_interval).rotate(sourceRotMat(), fieldRotMat());
                break;
            case 2:
                return convert_interval(GPt4D((1.0 - xsi * eta1) * eta3 + xsi * eta1,
                                               xsi * eta2,
                                              (1.0 - xsi * eta1) * eta3,
                                               xsi),
                                        unit_interval,
                                        parent_interval).rotate(sourceRotMat(), fieldRotMat());
                break;
            case 3:
                return convert_interval(GPt4D((1.0 - xsi * eta1) * eta3 + xsi * eta1,
                                              xsi,
                                              (1.0 - xsi * eta1) * eta3,
                                              xsi * eta2),
                                        unit_interval,
                                        parent_interval).rotate(sourceRotMat(), fieldRotMat());
                break;
            case 4:
                return convert_interval(GPt4D((1.0 - xsi * eta1) * eta3,
                                              xsi * eta2,
                                              (1.0 - xsi * eta1) * eta3 + xsi * eta1,
                                              xsi),
                                        unit_interval,
                                        parent_interval).rotate(sourceRotMat(), fieldRotMat());
                break;
            case 5:
                return convert_interval(GPt4D((1.0 - xsi * eta1) * eta3,
                                              xsi,
                                              (1.0 - xsi * eta1) * eta3 + xsi * eta1,
                                              xsi * eta2),
                                        unit_interval,
                                        parent_interval).rotate(sourceRotMat(), fieldRotMat());
				break;
			default:
				std::cerr << "Reached bad subelement. Aborting.\n";
				std::exit(EXIT_FAILURE);
		}
		
	}
    
	double IEdgeQuadrature::getWeightImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
		const double s = p[0];
		const double s1 = p[1];
		const double ws = getBaseWeight(SOURCE);
		const double wf = getBaseWeight(FIELD);
		if(subCellI() < 2)
			return ws * wf * s * s * (1.0 - s);
		else if(subCellI() < 6 && subCellI() > 1)
			return ws * wf * s * s * (1.0 - s * s1);
		else {
			std::cerr << "Reached bad subelement. Aborting.\n";
			std::exit(EXIT_FAILURE);
		}
	}
	
}

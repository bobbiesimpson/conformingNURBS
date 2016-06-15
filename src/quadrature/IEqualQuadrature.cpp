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
		const double s = p[0];
		const double s1 = p[1];
		const double s2 = p[2];
		const double s3 = p[3];
		switch(subCellI()) {
			case 0:
				return convert_interval(GPt4D((1.0 - s) * s3,
												(1.0 - s * s1) * s2,
												s + (1.0 - s) * s3,
												s * s1 + (1.0 - s * s1) * s2),
										unit_interval,
										parent_interval);
				break;
			case 1:
				return convert_interval(GPt4D((1.0 - s * s1) * s2,
												(1.0 - s) * s3,
												s * s1 + (1.0 - s * s1) * s2,
												s + (1.0 - s) * s3),
										unit_interval,
										parent_interval);
				break;
			case 2:
				return convert_interval(GPt4D((1.0 - s) * s3,
												s * s1 + (1.0 - s * s1) * s2,
												s + (1.0 - s) * s3,
												(1.0 - s * s1) * s2),
										unit_interval,
										parent_interval);
				break;
			case 3:
				return convert_interval(GPt4D((1.0 - s * s1) * s2,
												s + (1.0 - s) * s3,
												s * s1 + (1.0 - s * s1) * s2,
												(1.0 - s) * s3),
										unit_interval,
										parent_interval);
				break;
			case 4:
				return convert_interval(GPt4D(s + (1.0 - s) * s3,
												(1.0 - s * s1) * s2,
												(1.0 - s) * s3,
												s * s1 + (1.0 - s * s1) * s2),
				 						unit_interval,
				 						parent_interval);
				break;
			case 5:
				return convert_interval(GPt4D(s * s1 + (1.0 - s * s1) * s2,
												(1.0 - s) * s3,
												(1.0 - s * s1) * s2,
												s + (1.0 - s) * s3),
										unit_interval,
										parent_interval);
				break;
			case 6:
				return convert_interval(GPt4D(s + (1.0 - s) * s3,
												s * s1 + (1.0 - s * s1) * s2,
												(1.0 - s) * s3,
												(1.0 - s * s1) * s2),
										unit_interval,
										parent_interval);
				break;
			case 7:
				return convert_interval(GPt4D(s * s1 + (1.0 - s * s1) * s2,
												s + (1.0 - s) * s3,
												(1.0 - s * s1) * s2,
												(1.0 - s) * s3),
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
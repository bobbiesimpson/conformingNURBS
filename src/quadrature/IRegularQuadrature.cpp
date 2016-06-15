//
//  IRegularQuadrature.cpp
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#include "IRegularQuadrature.h"

namespace nurbs
{
	GPt4D IRegularQuadrature::getImpl() const
	{
        return getCurrentBiUnitIntervalPt();
	}
    
	double IRegularQuadrature::getWeightImpl() const
	{
		return getBaseWeight(SOURCE) * getBaseWeight(FIELD);
	}
}
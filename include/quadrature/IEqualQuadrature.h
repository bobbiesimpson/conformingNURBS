//
//  IEqualQuadrature.h
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#ifndef NURBS_IEQUAL_QUADRATURE_H
#define NURBS_IEQUAL_QUADRATURE_H

#include <iostream>

#include "IGalerkinIntegrate.h"

namespace nurbs {
    
	class IEqualQuadrature : public IGalerkinIntegrate {
		
    public:
        
		/// Constructor
		IEqualQuadrature(const UIntVec& sourceorders,
                         const UIntVec& fieldorders)
        : IGalerkinIntegrate(sourceorders, fieldorders) {}
		
    private:
        
		/// For the current set of quadrature indices and sub element,
		/// return the corresponding 4-d quadrature point
		GPt4D getImpl() const;
        
		/// For the current 4-d quadrature point return the quadrature weight
		/// including any necessary jacobian determinants.
		double getWeightImpl() const;
		
		/// 8 subcells for Sauter and Schwab transformation
		uint subCellN() const { return 8; };
        
		/// Print function implementation
		void printImpl(std::ostream& ost) const {}
	};
}

#endif

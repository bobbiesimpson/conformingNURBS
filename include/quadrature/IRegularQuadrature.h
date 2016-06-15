//
//  IRegularQuadrature.h
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#ifndef NURBS__IREGULAR_QUADRATURE_H
#define NURBS__IREGULAR_QUADRATURE_H

#include <iostream>

#include "IGalerkinIntegrate.h"

namespace nurbs {
    
    class IRegularQuadrature : public IGalerkinIntegrate {
        
    public:
        
		IRegularQuadrature(const UIntVec& sourceorders,
                           const UIntVec& fieldorders)
        : IGalerkinIntegrate(sourceorders, fieldorders) {}
        
    private:
        
		/// For the current set of quadrature indices and sub element,
		/// return the corresponding 4-d quadrature point
		GPt4D getImpl() const;
        
		/// For the current 4-d quadrature point return the quadrature weight
		/// including any necessary jacobian determinants.
		double getWeightImpl() const;
		
		/// 6 subcells for this transformation (see p.318 of Sauter and Schwab)
		uint subCellN() const { return 1; };
        
		/// Print function implementation (Currently does nothing)
		void printImpl(std::ostream& ost) const {}
		
	};
	
}
#endif /* defined(__nurbslib__IRegularQuadrature__) */

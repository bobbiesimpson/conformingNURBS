//
//  IEqualQuadratureTri.h
//  nurbslib
//
//  Created by Zhaowei Liu on 10/11/2016.
//
//

#ifndef NURBS_IEQUAL_TRI_QUADRATURE_H
#define NURBS_IEQUAL_TRI_QUADRATURE_H

#include <iostream>

#include "IGalerkinIntegrate.h"

namespace nurbs {
    
	class IEqualQuadratureTri : public IGalerkinIntegrate {
		
    public:
        
		/// Constructor
		IEqualQuadratureTri(const UIntVec& sourceorders,
                         const UIntVec& fieldorders,
                            const Edge dedge)
        :
        mdedge(dedge),
        IGalerkinIntegrate(sourceorders, fieldorders) {}
		
    private:
        
		/// For the current set of quadrature indices and sub element,
		/// return the corresponding 4-d quadrature point
		GPt4D getImpl() const;
        
		/// For the current 4-d quadrature point return the quadrature weight
		/// including any necessary jacobian determinants.
        
		double getWeightImpl() const;
        		
		/// 6 subcells for Sauter and Schwab transformation
		uint subCellN() const { return 6; };
        
		/// Print function implementation
		void printImpl(std::ostream& ost) const {}
        
        Edge mdedge;
	};
}

#endif

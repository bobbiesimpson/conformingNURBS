//
//  IEdgeQuadrature.h
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#ifndef NURBS_IEDGE_QUADRATURE_H
#define NURBS_IEDGE_QUADRATURE_H

#include <iostream>

#include "algebra.h"
#include "base.h"
#include "IGalerkinIntegrate.h"

namespace nurbs {
    
    /// The class used for generating 4D quadrature points and weights
    /// in the case of two elements that have an adjacent edge.
    /// If no edge is specified, then it is assumed the adjacent edge lies on
    /// (s, 0) x (s, 0) which corresponds to Edge1 in the enumeration specified
    /// in base.h
    
    class IEdgeQuadrature : public IGalerkinIntegrate {
        
    public:
        
		IEdgeQuadrature(const UIntVec& sourceorders,
                        const UIntVec& fieldorders,
                        const Edge esrc,
                        const Edge efield)
        :
        IGalerkinIntegrate(sourceorders, fieldorders)
        {
            switch (esrc) {
                case Edge::EDGE0:
                    mSrcRotMatrix = rotMatrix(0.0);
                    break;
                case Edge::EDGE1:
                    mSrcRotMatrix = rotMatrix(PI);
                    break;
                case Edge::EDGE2:
                    mSrcRotMatrix = rotMatrix(-PI * 0.5);
                    break;
                case Edge::EDGE3:
                    mSrcRotMatrix = rotMatrix(PI * 0.5);
                    break;
            }
            
            switch (efield) {
                case Edge::EDGE0:
                    mFieldRotMatrix = rotMatrix(0.0);
                    break;
                case Edge::EDGE1:
                    mFieldRotMatrix = rotMatrix(PI);
                    break;
                case Edge::EDGE2:
                    mFieldRotMatrix = rotMatrix(-PI * 0.5);
                    break;
                case Edge::EDGE3:
                    mFieldRotMatrix = rotMatrix(PI * 0.5);
                    break;
            }
        }
        
    private:
        
		/// For the current set of quadrature indices and sub element,
		/// return the corresponding 4-d quadrature point
		GPt4D getImpl() const;
        
		/// For the current 4-d quadrature point return the quadrature weight
		/// including any necessary jacobian determinants.
		double getWeightImpl() const;
		
		/// 6 subcells for this transformation (see p.318 of Sauter and Schwab)
		uint subCellN() const { return 6; };
        
		/// Print function implementation (Currently does nothing)
		void printImpl(std::ostream& ost) const {}
        
        /// Source rotation matrix getter
        const DoubleVecVec& sourceRotMat() const { return mSrcRotMatrix; }
        
        /// Field rotation matrix getter
        const DoubleVecVec& fieldRotMat() const { return mFieldRotMatrix; }
        
        /// The rotation matrix used for handling different source edge
        /// singularities
        DoubleVecVec mSrcRotMatrix;
        
        /// The rotation matrix used for handling different field edge
        /// singularities
        DoubleVecVec mFieldRotMatrix;
		
	};
	
}
#endif

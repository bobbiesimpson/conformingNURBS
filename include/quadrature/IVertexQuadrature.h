//
//  IVertexQuadrature.h
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#ifndef NURBS_IVERTEX_QUADRATURE_H
#define NURBS_IVERTEX_QUADRATURE_H

#include <iostream>

#include "algebra.h"
#include "base.h"
#include "IGalerkinIntegrate.h"

namespace nurbs {
    
    /// The interface for computing 4D quadrature points and weights
    /// for Galerkin BEM quadrature in the case that the source
    /// and field elements have an adjacent edge. If no vertices are specified
    /// during constructing then the default is to assume the singular vertex is
    /// located at (-1,-1) in both the source and field elements.
    
    class IVertexQuadrature : public IGalerkinIntegrate {
        
    public:
        
        /// Constructor
		IVertexQuadrature(const UIntVec& sourceorders,
                          const UIntVec& fieldorders,
                          const Vertex svertex,
                          const Vertex fvertex)
        : IGalerkinIntegrate(sourceorders, fieldorders)
        {
            switch (svertex) {
                case Vertex::VERTEX0:
                    mSrcRotMatrix = rotMatrix(0.0);
                    break;
                case Vertex::VERTEX1:
                    mSrcRotMatrix = rotMatrix(PI * 0.5);
                    break;
                case Vertex::VERTEX2:
                    mSrcRotMatrix = rotMatrix(-PI * 0.5);
                    break;
                case Vertex::VERTEX3:
                    mSrcRotMatrix = rotMatrix(PI);
                    break;
            }
            
            switch (fvertex) {
                case Vertex::VERTEX0:
                    mFieldRotMatrix = rotMatrix(0.0);
                    break;
                case Vertex::VERTEX1:
                    mFieldRotMatrix = rotMatrix(PI * 0.5);
                    break;
                case Vertex::VERTEX2:
                    mFieldRotMatrix = rotMatrix(-PI * 0.5);
                    break;
                case Vertex::VERTEX3:
                    mFieldRotMatrix = rotMatrix(PI);
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
		uint subCellN() const { return 4; };
        
		/// Print function implementation (Currently does nothing)
		void printImpl(std::ostream& ost) const {}
        
        /// Source rotation matrix getter
        const DoubleVecVec& sourceRotMat() const { return mSrcRotMatrix; }
        
        /// Field rotation matrix getter
        const DoubleVecVec& fieldRotMat() const { return mFieldRotMatrix; }
        
        /// Rotation matrix for source element
        DoubleVecVec mSrcRotMatrix;
        
        /// Rotation matrix for field element
        DoubleVecVec mFieldRotMatrix;
		
	};
	
}
#endif

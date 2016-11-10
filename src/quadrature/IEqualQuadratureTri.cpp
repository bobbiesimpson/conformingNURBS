//
//  IEqualQuadrature.cpp
//  nurbslib
//
//  Created by Robert Simpson on 18/07/2014.
//
//

#include "IEqualQuadratureTri.h"

namespace nurbs {
    
    GPt4D IEqualQuadratureTri::getImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
        double xsi,eta1,eta2,eta3;
        if (mdedge == Edge::EDGE2) {
             xsi = p[0];
             eta1 = p[1];
             eta2 = p[2];
             eta3 = p[3];
        }
        else if (mdedge == Edge::EDGE3){
             xsi = 1.0 - p[0];
             eta1 = 1.0 - p[1];
             eta2 = 1.0 - p[2];
             eta3 = 1.0 - p[3];
        }
        else{
            throw std::runtime_error("Degenerated edge is not edge 2 or 3!");
        }
        // some additional defitiions
        const double aux1 = xsi;
        const double aux2 = xsi * (1 - eta1 + eta1 * eta2);
        const double aux3 = xsi * (1 - eta1 * eta2 * eta3);
        const double aux4 = xsi * (1 - eta1);
        const double aux5 = xsi * eta1 * (1 - eta2 + eta2 * eta3);
        const double aux6 = xsi * (1 - eta1 * eta2);
        const double aux7 = xsi * eta1 * (1 - eta2);
        const double aux8 = xsi * eta1 * (1 - eta2 * eta3);
        if (mdedge == Edge::EDGE2) {
            switch(subCellI()) {
                case 0:
                    return convert_interval_Tri_Edge2(GPt4D(aux1,
                                                      aux2,
                                                      aux3,
                                                      aux4));
                    break;
                case 1:
                    return convert_interval_Tri_Edge2(GPt4D(aux3,
                                                      aux4,
                                                      aux1,
                                                      aux2));
                    break;
                case 2:
                    return convert_interval_Tri_Edge2(GPt4D(aux1,
                                                      aux5,
                                                      aux6,
                                                      aux7));
                    break;
                case 3:
                    return convert_interval_Tri_Edge2(GPt4D(aux6,
                                                      aux7,
                                                      aux1,
                                                      aux5));
                    break;
                case 4:
                    return convert_interval_Tri_Edge2(GPt4D(aux3,
                                                      aux8,
                                                      aux1,
                                                      aux7));
                    break;
                case 5:
                    return convert_interval_Tri_Edge2(GPt4D(aux1,
                                                      aux7,
                                                      aux3,
                                                      aux8));
                    break;
                    
                default:
                    std::cerr << "Reached bad subelement. Aborting.\n";
                    exit(EXIT_FAILURE);
            }
        }
        else if (mdedge == Edge::EDGE3){
            switch(subCellI()) {
                case 0:
                    return convert_interval_Tri_Edge3(GPt4D(aux1,
                                                      aux2,
                                                      aux3,
                                                      aux4));
                    break;
                case 1:
                    return convert_interval_Tri_Edge3(GPt4D(aux3,
                                                      aux4,
                                                      aux1,
                                                      aux2));
                    break;
                case 2:
                    return convert_interval_Tri_Edge3(GPt4D(aux1,
                                                      aux5,
                                                      aux6,
                                                      aux7));
                    break;
                case 3:
                    return convert_interval_Tri_Edge3(GPt4D(aux6,
                                                      aux7,
                                                      aux1,
                                                      aux5));
                    break;
                case 4:
                    return convert_interval_Tri_Edge3(GPt4D(aux3,
                                                      aux8,
                                                      aux1,
                                                      aux7));
                    break;
                case 5:
                    return convert_interval_Tri_Edge3(GPt4D(aux1,
                                                      aux7,
                                                      aux3,
                                                      aux8));
                    break;
                    
                default:
                    std::cerr << "Reached bad subelement. Aborting.\n";
                    exit(EXIT_FAILURE);
            }

        }
        else{
            throw std::runtime_error("Degenerated edge is not edge 2 or 3!");
        }
    }
    
	double IEqualQuadratureTri::getWeightImpl() const
	{
		GPt4D p = getCurrentUnitIntervalPt();
        double xsi,eta1,eta2,eta3;
        if (mdedge == Edge::EDGE2) {
            xsi = p[0];
            eta1 = p[1];
            eta2 = p[2];
            eta3 = p[3];
        }
        else if (mdedge == Edge::EDGE3){
            xsi = 1.0 - p[0];
            eta1 = 1.0 - p[1];
            eta2 = 1.0 - p[2];
            eta3 = 1.0 - p[3];
        }
        else{
            throw std::runtime_error("Degenerated edge is not edge 2 or 3!");
        }
        // some additional defitiions
        const double aux1 = xsi;
        const double aux3 = xsi * (1 - eta1 * eta2 * eta3);
        const double aux6 = xsi * (1 - eta1 * eta2);
        double w = getBaseWeight(SOURCE) * getBaseWeight(FIELD) * xsi * xsi* xsi  * eta1 * eta1* eta2;
        switch(subCellI()) {
            case 0:
                return w / aux1 / aux3;
                break;
            case 1:
                return w / aux1 / aux3;
                break;
            case 2:
                return w / aux1 / aux6;
                break;
            case 3:
                return w / aux1 / aux6;
                break;
            case 4:
                return w / aux1 / aux3;
                break;
            case 5:
                return w / aux1 / aux3;
                break;
                
            default:
                std::cerr << "Reached bad subelement. Aborting.\n";
                exit(EXIT_FAILURE);
        }
    }
}
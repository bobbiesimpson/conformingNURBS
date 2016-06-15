//
//  InputDataStructures.h
//  nurbslib
//
//  Created by Robert Simpson on 05/12/2014.
//
//

#ifndef INPUT_DATA_STRUCTURES_H
#define INPUT_DATA_STRUCTURES_H

#include <iostream>

namespace nurbs {
    
    /// A simple structure for inputting connectivity vectors for a b-spline space
    template<typename T>
    struct InputVec {
        
        std::vector<T> data;
        
        friend std::istream& operator>>(std::istream& ist, InputVec<T>& v)
        {
            char ch;
            if(ist >> ch && ch != '{') {
                ist.unget();
                ist.clear( std::ios_base::failbit );
                return ist;
            }
            while(true) {
                T i;
                if(!(ist >> i >> ch))
                    break;
                //std::cout << i << " " << ch << "\n";
                if(ch == ',' || ch == '}') {
                    v.data.push_back(i);
                    if(ch == '}') {
                        ist.unget();
                        ist.clear(std::ios::failbit);
                        break;
                    }
                }
                else {
                    error("Bad vector end");
                }
            }
            endOfLoop(ist, '}', "Could not read vector");
            return ist;
        }
        
        friend std::ostream& operator<<(std::ostream& ost, const InputVec& v)
        {
            ost << "{";
            for(uint i = 0; i < v.data.size(); ++i) {
                ost << v.data[i];
                if(i < v.data.size() - 1)
                    ost <<  " ";
            }
            ost << "}";
            return ost;
        }
    };
    
    /// A simple struct used for holding input data when reading in connectivity info.
    struct ConnVecInput {
        
        std::string name;
        
        std::vector<uint> data;
        
        friend std::istream& operator>>(std::istream& ist, ConnVecInput& c)
        {
            char ch;
            if(ist >> ch && ch != '{') {
                ist.unget();
                ist.clear( std::ios_base::failbit );
                return ist;				
            }
            std::string s;
            if(!(ist >> s))
                error("Bad connectivity string");
            c.name = s;
            InputVec<uint> v;
            if(!(ist >> v))
                error("Bad connectivity vector reading");
            c.data = v.data;
            return ist;			
        }
        
        friend std::ostream& operator<<(std::ostream& ost, const ConnVecInput& c)
        {
            ost << c.name << " ";
            ost << "{";
            for(uint i = 0; i < c.data.size(); ++i) {
                ost << c.data[i];
                if(i < c.data.size() - 1)
                    ost <<  " ";
            }
            ost << "}";
            return ost;
        }	
    };
}
#endif

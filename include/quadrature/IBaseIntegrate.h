#ifndef IBASE_INTEGRATE_H
#define IBASE_INTEGRATE_H

#include "IElemIntegrate.h"

namespace nurbs {
    
    namespace elem {
        
        /// A base class that defines the interface for an integrator
        /// over a surface manifold element
        
        class IBaseIntegrate {
            
        public:
            
            /// Indicates that iterator is finished
            virtual bool isDone() const = 0;
            
            /// Get the two-dimensional guass point
            virtual GPt2D get() const = 0;
            
            /// Get specified component
            double get(uint icomp) const { return get().get(icomp); }
            
            /// Quadrature weight getter
            virtual double getWeight() const = 0;
            
            /// Restart the iterator
            virtual void restart() = 0;
            
            /// Prefix increment
            IBaseIntegrate& operator++()
            {
                incrementImpl();
                return *this;
            }
            
        protected:
            
            /// Protected constructor
            IBaseIntegrate() = default;
            
        private:
            
            /// Increment implementation function that must be implemented
            virtual void incrementImpl() = 0;
            
        };
    }
}


#endif
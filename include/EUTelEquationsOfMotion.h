/*
 * File:   EUTelMagneticFieldFinder.h
 *
 * Created on July 2, 2013, 12:53 PM
 */

#ifndef EUTELEQUATIONSOFMOTION_H
#define EUTELEQUATIONSOFMOTION_H

// EUTELESCOPE
#include "EUTelUtilityRungeKutta.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif

namespace eom {
        /** 
         * @class Implementation of particles differential
         * equation of motion
         */
        class EOMODE : public ODE {
          private:
            DISALLOW_COPY_AND_ASSIGN( EOMODE )
          public:
            explicit EOMODE( int neq ) : 
            ODE(neq),
            _h() {}
            
            virtual ~EOMODE() {};
            
            virtual TVectorD evalRHS( const TVectorD& point ) {
                
                streamlog_out( DEBUG2 ) << "EOMODE::evalRHS()" << std::endl;
                
                streamlog_out( DEBUG0 ) << "Input vector" << std::endl;
                streamlog_message( DEBUG0, point.Print();, std::endl; );
                
                TVectorD result( this->getNEquations() );
                
                const double mm = 1000.;
                const double k = 0.299792458/mm;
                
                const double tx = point[ 2 ];
                const double ty = point[ 3 ];
                const double q  = point[ 4 ];
                
                TVector2 a = A( tx, ty );
                const double dxdz  = tx;
                const double dydz  = ty;
                const double dtxdz = q * k * a.X();
                const double dtydz = q * k * a.Y();
                
                result[ 0 ] = dxdz;
                result[ 1 ] = dydz;
                result[ 2 ] = dtxdz;
                result[ 3 ] = dtydz;
                result[ 4 ] = 0;
                
                streamlog_out( DEBUG0 ) << "Result vector" << std::endl;
                streamlog_message( DEBUG0, result.Print();, std::endl; );
                
                streamlog_out( DEBUG2 ) << "-----------------------------EOMODE::evalRHS()------------------------------" << std::endl;
                
                return result;
            }
            
            void setBField( const TVector3& h ) {
                _h = h;
            }
            
          private:
              /**
               * Calculation of A vector necessary for rhs of particle's eom
               * @param tx particle's tx parameter
               * @param ty particle's ty parameter
               * 
               * @return 2d vector A
               */
            TVector2 A( double tx, double ty ) const {
               const double Bx = _h.X();
               const double By = _h.Y();
               const double Bz = _h.Z();
                
               const double sqrtFactor = sqrt( 1. + tx*tx + ty*ty );
               const double Ax = sqrtFactor * (  ty * ( tx * Bx + Bz ) - ( 1. + tx*tx ) * By );
               const double Ay = sqrtFactor * ( -tx * ( ty * By + Bz ) + ( 1. + ty*ty ) * Bx );
               
               return TVector2( Ax, Ay );
            }
            
            /** Magnetic field vector */
            TVector3 _h;
        };
}

#endif //EUTELEQUATIONSOFMOTION_H

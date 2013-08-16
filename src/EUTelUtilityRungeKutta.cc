/* 
 * File:   EUTelUtilityRungeKutta.cc
 * 
 * Created on August 15, 2013, 1:51 PM
 */

#include "EUTelUtilityRungeKutta.h"

/**
 * Runge-Kutta ODE solver constructor
 */
EUTelUtilityRungeKutta::EUTelUtilityRungeKutta( ) :
_safetyFactor(0.95),
_nIterations(0),
_absTolerance(1E-6),
_relTolerance(1E-6),
_butcherTableau(0),
_ode(0),
_step(0.)
{
}

/**
 * Runge-Kutta ODE solver constructor
 * 
 * @param atol absolute tolerance
 * @param rtol relative tolerance
 * @param safety safety factor
 */
EUTelUtilityRungeKutta::EUTelUtilityRungeKutta( double atol = 1E-6, double rtol = 1E-6, double safety = 0.95 ) : 
_safetyFactor(safety),
_nIterations(0),
_absTolerance(atol),
_relTolerance(rtol),
_butcherTableau(0),
_ode(0),
_step(0.)
{
}

EUTelUtilityRungeKutta::~EUTelUtilityRungeKutta( ) {
    delete _butcherTableau;
    delete _ode;
}

void EUTelUtilityRungeKutta::setInitValue( const TVectorD& initValue ) {
    this->_ode->setInitValue( initValue );
}

TVectorD EUTelUtilityRungeKutta::getInitValue( ) const {
    return this->_ode->getInitValue();
}

int EUTelUtilityRungeKutta::getNEquations( ) const {
    return this->_ode->getNEquations();
}

double EUTelUtilityRungeKutta::getNIterations( ) const {
    return _nIterations;
}

void EUTelUtilityRungeKutta::setSafetyFactor( double safetyFactor ) {
    this->_safetyFactor = safetyFactor;
}

double EUTelUtilityRungeKutta::getSafetyFactor( ) const {
    return _safetyFactor;
}

void EUTelUtilityRungeKutta::setRelTolerance( double relTolerance ) {
    this->_relTolerance = relTolerance;
}

double EUTelUtilityRungeKutta::getRelTolerance( ) const {
    return _relTolerance;
}

void EUTelUtilityRungeKutta::setAbsTolerance( double absTolerance ) {
    this->_absTolerance = absTolerance;
}

double EUTelUtilityRungeKutta::getAbsTolerance( ) const {
    return _absTolerance;
}

void EUTelUtilityRungeKutta::setStep( double step ) {
    this->_step = step;
}

double EUTelUtilityRungeKutta::getStep( ) const {
    return _step;
}

void EUTelUtilityRungeKutta::setRhs( ODE* ode ) {
    this->_ode = ode;
}

void EUTelUtilityRungeKutta::setButcherTableau( ButcherTableau* bt ) {
    this->_butcherTableau = bt;
}

ODE* EUTelUtilityRungeKutta::getRhs( ) const {
    return _ode;
}

/**
 *  Integrate equation of motion
 * 
 * @param h desired step
 * 
 * @return Solution Y(h)
 */
TVectorD EUTelUtilityRungeKutta::integrate( double h ) const {
    
    streamlog_out( DEBUG2 ) << "EUTelUtilityRungeKutta::integrate()" << std::endl;
    
    streamlog_out( DEBUG0 ) << "Step size: " << h << std::endl;
    
    unsigned int nComponents = _ode->getNEquations();
    
    streamlog_out( DEBUG0 ) << "N equations: " << nComponents << std::endl;
    
    TVectorD result( nComponents );
    
    if ( _butcherTableau->isEmbedded() ) {
        
        unsigned int nStages = _butcherTableau->getNStages();

        TVectorD *pm = new TVectorD[ nStages ];
        TVectorD *km = new TVectorD[ nStages ];
        for ( unsigned int m = 0; m < nStages; ++m ) { 
            pm[m].ResizeTo( 0, nComponents-1 );
            km[m].ResizeTo( 0, nComponents-1 ); 
        }

        TVectorD temp( nComponents );
        // ODE integration
        {
            // Calculation of pm and km
            for ( int m = 0; m < nStages; ++m ) {
                pm[ m ] = this->_ode->getInitValue();
                temp.Zero();
                for ( int n = 0; n <= m - 1 ; ++n ) {
                    temp = km[ n ];
                    const double bmn = _butcherTableau->_rungeKutta[ m ][ n ];
                    streamlog_out( DEBUG0 ) << "b[" << m << "][" << n << "]= " << bmn << std::endl;
                    temp *= bmn;
                    pm[ m ] += temp;
                }
                
                streamlog_out( DEBUG0 ) << "pm[ " << m << "]" << std::endl;
                streamlog_message( DEBUG0, pm[ m ].Print();, std::endl; );
                
                km[ m ] = _ode->evalRHS( pm[ m ] ); km[ m ] *= h;
                streamlog_out( DEBUG0 ) << "km[ " << m << "]" << std::endl;
                streamlog_message( DEBUG0, km[ m ].Print();, std::endl; );
            }
            
            for ( unsigned int m = 0; m < nStages; ++m ) {
                temp.Zero();
                const double cHO = _butcherTableau->_weightsHigherOrder[ m ];
                temp = km[ m ]; temp *= cHO;
                result += temp;
            }
            
            result += pm[0];
            
            streamlog_out( DEBUG0 ) << "Solution:" << std::endl;
            streamlog_message( DEBUG0, result.Print();, std::endl; );
        }
        
        // Estimate error
        TVectorD delta( nComponents );
        delta.Zero();
        {
            for ( unsigned int m = 0; m < nStages; ++m ) {
                temp.Zero();
                const double cHO = _butcherTableau->_weightsHigherOrder[ m ];
                const double cLO = _butcherTableau->_weightsLowerOrder[ m ];
                streamlog_out( DEBUG0 ) << "Weights HO/LO:\t" << cHO << "\t" << cLO << std::endl;
                temp = km[ m ]; temp *= (  cHO - cLO );
                delta += temp;
            }
            streamlog_out( DEBUG0 ) << "Error estimate:" << std::endl;
            streamlog_message( DEBUG0, delta.Print();, std::endl; );
        }
        
        delete[] pm;
        delete[] km;
        
    } else {
        streamlog_out( WARNING1 ) << "Not embedded Runge-Kutta is not supported!" << std::endl;
    }
    
    streamlog_out( DEBUG2 ) << "--------------------------------EUTelUtilityRungeKutta::integrate()-------------------------------" << std::endl;
    
    return result;
}



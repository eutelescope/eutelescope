/* 
 * File:   EUTelUtilityRungeKutta.h
 * Author: diont
 *
 * Created on August 15, 2013, 1:50 PM
 */

#ifndef EUTELUTILITYRUNGEKUTTA_H
#define	EUTELUTILITYRUNGEKUTTA_H

#include "EUTELESCOPE.h"
#include "EUTelUtility.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVectorD.h"
#include "TMatrixD.h"
#endif

#include <vector>

/**
 * @Class ButcherTableau
 * Class encapsulates Butcher tableau data structure 
 */
class ButcherTableau {
private:
    DISALLOW_COPY_AND_ASSIGN(ButcherTableau) // prevent users from making (default) copies of processors   
public:
    ButcherTableau() : 
    _isEmbedded(false),
    _nStages(0),
    _rungeKutta(),
    _nodes(),
    _weightsHigherOrder(),
    _weightsLowerOrder()
    {}
    
    virtual ~ButcherTableau() {};

    virtual void init( ) {};
    
    virtual void init( const TMatrixD& rk, const TVectorD& nodes, const TVectorD& weightsHO ) {
        this->_isEmbedded = false;
        this->_rungeKutta = rk;
        this->_nodes = nodes;
        this->_weightsHigherOrder = weightsHO;
    }
    
    virtual void init( const TMatrixD& rk, const TVectorD& nodes, const TVectorD& weightsHO, const TVectorD& weightsLO ) {
        this->_isEmbedded = true;
        this->_rungeKutta = rk;
        this->_nodes = nodes;
        this->_weightsHigherOrder = weightsHO;
        this->_weightsLowerOrder = weightsLO;
    }

    virtual void setEmbedded( bool _isEmbedded ) {
        this->_isEmbedded = _isEmbedded;
    }
    
    virtual bool isEmbedded( ) const {
        return this->_isEmbedded;
    }

    virtual void setNStages( int nStages ) {
        this->_nStages = nStages;
    }

    virtual int getNStages() const {
        return _nStages;
    }
    
private:
    /** Embedded Runge-Kutta flag*/
    bool _isEmbedded;
    
    /** Number of stages for solution evaluation */
    int _nStages;
    
public:
    /** Runge-Kutta weight matrix */
    TMatrixD _rungeKutta;
    
    /** Vector of nodes */
    TVectorD _nodes;
    
    /** Vector of weights */
    TVectorD _weightsHigherOrder;
    
    /** Vector of weight for lower order embedded Runge-Kutta */
    TVectorD _weightsLowerOrder;
            
};

/**
 * @class ButcherTableauCashKarp
 * Cash-Karp Butcher tableau for 5th order embedded Runge-Kutta integration
 */
class ButcherTableauCashKarp : public ButcherTableau {
private:
    DISALLOW_COPY_AND_ASSIGN(ButcherTableauCashKarp) // prevent users from making (default) copies of processors  
    
public:
    
    ButcherTableauCashKarp() : ButcherTableau() { init(); }
    
    virtual ~ButcherTableauCashKarp() {};
    
    virtual void init() {
        setEmbedded(true);
        
        const int nSolverStages = 6;
        
        _rungeKutta.ResizeTo( nSolverStages, nSolverStages );
        _nodes.ResizeTo( nSolverStages );
        _weightsHigherOrder.ResizeTo( nSolverStages );
        _weightsLowerOrder.ResizeTo( nSolverStages );
        
        setNStages( nSolverStages );
        
        // Fill Runge-Kutta matrix
        TMatrixD rk( nSolverStages, nSolverStages );
        rk.Zero();
        rk[1][0] = 1./5.;
        rk[2][0] = 3./40.;        rk[2][1] = 9./40.;
        rk[3][0] = 3./10.;        rk[3][1] = -9./10.;    rk[3][2] = 6./5.;
        rk[4][0] = -11./54.;      rk[4][1] = 5./2.;      rk[4][2] = -70./27.;     rk[4][3] = 35./27.;
        rk[5][0] = 1631./55296.;  rk[5][1] = 175./512.;  rk[5][2] = 575./13824.;  rk[5][3] = 44275./110592.; rk[5][4] = 253./4096.;
        
        // Fill nodes
        TVectorD nodes( nSolverStages );
        nodes[0] = 0.; nodes[1] = 1./5.; nodes[2] = 3./10.; nodes[3] = 3./5.; nodes[4] = 1.; nodes[5] = 7./8.;
        
        // Fill 5th order weights
        TVectorD weights5( nSolverStages );
        weights5[0] = 37./378.; weights5[1] = 0.; weights5[2] = 250./621.; weights5[3] = 125./594.; weights5[4] = 0.; weights5[5] = 512./1771.; 
        
        // Fill 4th order weights
        TVectorD weights4( nSolverStages );
        weights4[0] = 2825./27648.; weights4[1] = 0.; weights4[2] = 18575./48384.; weights4[3] = 13525./55296.; weights4[4] = 277./14336.; weights4[5] = 1./4.; 
        
        // Initialise tableau
        ButcherTableau::init( rk, nodes, weights5, weights4 );
    }
};

/**
 * @class ButcherTableauDormandPrince
 * Dormand-Prince Butcher tableau for 5th order embedded Runge-Kutta integration
 */
class ButcherTableauDormandPrince : public ButcherTableau {
private:
    DISALLOW_COPY_AND_ASSIGN(ButcherTableauDormandPrince) // prevent users from making (default) copies of processors  
    
public:
    
    ButcherTableauDormandPrince() : ButcherTableau() { init(); }
    
    virtual ~ButcherTableauDormandPrince() {};
    
    virtual void init() {
      
        setEmbedded(true);
        
        const int nSolverStages = 7;

        _rungeKutta.ResizeTo( nSolverStages, nSolverStages );
        _nodes.ResizeTo( nSolverStages );
        _weightsHigherOrder.ResizeTo( nSolverStages );
        _weightsLowerOrder.ResizeTo( nSolverStages );
        
        setNStages( nSolverStages );
        
        // Fill Runge-Kutta matrix
        TMatrixD rk( nSolverStages, nSolverStages );
        rk.Zero();
        rk[1][0] = 1./5.;
        rk[2][0] = 3./40.;       rk[2][1] = 9./40.;
        rk[3][0] = 44./45.;      rk[3][1] = -56./15.;      rk[3][2] = 32./9.;
        rk[4][0] = 19372./6561.; rk[4][1] = -25360./2187.; rk[4][2] = 64448./6561.; rk[4][3] = -212./729.;
        rk[5][0] = 9017./3168.;  rk[5][1] = -355./33.;     rk[5][2] = 46732./5247.; rk[5][3] = 49./176.;    rk[5][4] = -5103./18656.;
        rk[6][0] = 35./384.;     rk[6][1] = 0.;            rk[6][2] = 500./1113.;   rk[6][3] = 125./192.;   rk[6][4] = -2187./6784.;   rk[6][5] = 11./84.;
        
        // Fill nodes
        TVectorD nodes( nSolverStages );
        nodes[0] = 0.; nodes[1] = 1./5.; nodes[2] = 3./10.; nodes[3] = 4./5.; nodes[4] = 8./9.; nodes[5] = 1.; nodes[6] = 1.;
        
        // Fill 5th order weights
        TVectorD weights5( nSolverStages );
        weights5[0] = 35./384.; weights5[1] = 0.; weights5[2] = 500./1113.; weights5[3] = 125./192.; weights5[4] = -2187./6784.; weights5[5] = 11./84.; weights5[6] = 0.;
        
        // Fill 4th order weights
        TVectorD weights4( nSolverStages );
        weights4[0] = 5179./57600.; weights4[1] = 0.; weights4[2] = 7571./16695.; weights4[3] = 393./640.; weights4[4] = -92097./339200.; weights4[5] = 187./2100.; weights4[6] = 1./40.;
        
        // Initialise tableau
        ButcherTableau::init( rk, nodes, weights5, weights4 );
    }
};

/**
 * @class Ordinary differential equation
 * 
 * Class encapsulates initial value problem for a system of ordinary 
 * differential equations specified by explicit equation
 *  (*)     Y' = f(Y)
 * In this equation Y and f are vectors
 */
class ODE {
private:
    DISALLOW_COPY_AND_ASSIGN(ODE) // prevent users from making (default) copies of processors
public:
    ODE() : 
    _nEquations(0),
    _initValue(0)
    {}
    
    explicit ODE( int neq ) : 
    _nEquations( neq ),
    _initValue( neq )
    {}
 
    virtual ~ODE() {};
    
    /**
     * Evaluates right hand side of the differential equation (*)
     * This function has to be implemented in users instance of ODE
     * 
     * @param argument of function f 
     * @return right hand side evaluated at point
     */
    virtual TVectorD evalRHS( const TVectorD& ) = 0;
    
    void setInitValue( const TVectorD& initVec ) {
        streamlog_out( DEBUG0 ) << "ODE::setInitValue()" << std::endl;
        int initValNComponents = initVec.GetNrows();
        if ( initValNComponents !=  this->_nEquations ) {
            streamlog_out( WARNING1 ) << "Dimension of supplied vector of initial values (" << initValNComponents 
                                      << ") doesn't correspond to the number of equations (" << this->_nEquations << ")" << std::endl;
        } else {
            this->_initValue = initVec;
        }
        streamlog_out( DEBUG0 ) << "--------------------------------ODE::setInitValue()---------------------------" << std::endl;
    }

    TVectorD getInitValue( ) const {
        return _initValue;
    }

    int getNEquations( ) const {
        return _nEquations;
    }
    
private:
    /** Number of equations */
    int _nEquations;
    
    /** Initial value */
    TVectorD _initValue;
};

/**
 * @class  EUTelUtilityRungeKutta
 * Interface for Runge-Kutta type ODE solvers
 */
class EUTelUtilityRungeKutta {
private:
        DISALLOW_COPY_AND_ASSIGN(EUTelUtilityRungeKutta) // prevent users from making (default) copies of processors
public:
    EUTelUtilityRungeKutta();
    
    EUTelUtilityRungeKutta( double, double, double );
    
    virtual ~EUTelUtilityRungeKutta();

    int getNEquations() const;
    
    double getNIterations() const;
    
    void setSafetyFactor( double );
    
    double getSafetyFactor() const;
    
    void setRelTolerance( double );
    
    double getRelTolerance() const;
    
    void setAbsTolerance( double );
    
    double getAbsTolerance() const;
    
    void setButcherTableau( ButcherTableau* );

    void setRhs( ODE* );
    
    ODE* getRhs() const;
    
    void setInitValue( const TVectorD& );
    
    TVectorD getInitValue() const;
    
    void setStep(double _step);
    
    double getStep() const;
    
    TVectorD integrate( double ) const;
    
private:

    /** Adaptive step size safety factor */
    double _safetyFactor;
    
    /** Number of iterations performed to get the solution */
    double _nIterations;
    
    /** Absolute tolerance */
    double _absTolerance;
    
    /** Relative tolerance */
    double _relTolerance;

    ButcherTableau* _butcherTableau;    

    /** Right hand side of the ODE */
    ODE* _ode;
        
    /** Step size */
    double _step;
    
};

#endif	/* EUTELUTILITYRUNGEKUTTA_H */


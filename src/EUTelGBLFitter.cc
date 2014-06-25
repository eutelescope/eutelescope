/* 
 * File:   EUTelGBLFitter.cc
 * Contact: denys.lontkovskyi@desy.de
 * 
 * Created on January 25, 2013, 2:53 PM
 */

#ifdef USE_GBL

// its own header 
#include "EUTelGBLFitter.h"

// eutelescope includes ".h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelUtilityRungeKutta.h"
#include "EUTELESCOPE.h"

// marlin util includes
#include "mille/Mille.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// GBL
#include "include/GblTrajectory.h"
#include "include/GblPoint.h"
#include "include/GblData.h"
#include "include/BorderedBandMatrix.h"
#include "include/MilleBinary.h"
#include "include/VMatrix.h"


// system includes <>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <iterator>
#include <algorithm>

namespace eutelescope {

    EUTelGBLFitter::EUTelGBLFitter() : EUTelTrackFitter("GBLTrackFitter"),
    _trackCandidates(),
    _gblTrackCandidates(),
    _parPropJac(5, 5),
    _beamQ(-1),
    _eBeam(4.),
    _hitId2GblPointLabel(),
    _hitId2GblPointLabelMille(),
		_PointToState(),
		_binaryname(),
    _alignmentMode( Utility::noAlignment ), // default - no Alignment 
    _mEstimatorType(),
    _mille(0),
    _parameterIdXShiftsMap(),
    _parameterIdYShiftsMap(),
    _parameterIdZShiftsMap(),
    _parameterIdXRotationsMap(),
    _parameterIdYRotationsMap(),
    _parameterIdZRotationsMap(),
    _excludeFromFit(),
    _chi2cut(1000.),
    _eomIntegrator( new EUTelUtilityRungeKutta() ),
    _eomODE( 0 )
    {
       streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter default constructor FitTrackVec: " << getFitTrackVec() << std::endl;
                // Initialise ODE integrators for eom and jacobian
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
       // reset parent members: 
       streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter constructor over" << std::endl;
    }

    EUTelGBLFitter::EUTelGBLFitter(std::string name) : EUTelTrackFitter(name),
    _trackCandidates(),
    _gblTrackCandidates(),
    _parPropJac(5, 5),
    _beamQ(-1),
    _eBeam(-1.),
    _hitId2GblPointLabel(),
    _hitId2GblPointLabelMille(),
		_PointToState(),
		_binaryname(),
    _alignmentMode( Utility::noAlignment ),
    _mEstimatorType(),
    _mille(0),
    _parameterIdXShiftsMap(),
    _parameterIdYShiftsMap(),
    _parameterIdZShiftsMap(),
    _parameterIdXRotationsMap(),
    _parameterIdYRotationsMap(),
    _parameterIdZRotationsMap(),
    _excludeFromFit(),
    _chi2cut(1000.),
    _eomIntegrator( new EUTelUtilityRungeKutta() ),
    _eomODE(  0 )
    {
        streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter " << name << " constructor FitTrackVec: " << getFitTrackVec() << std::endl;
               // Initialise ODE integrators for eom and jacobian
                {
                    _eomODE = new eom::EOMODE(5);
                    _eomIntegrator->setRhs( _eomODE );
                    _eomIntegrator->setButcherTableau( new ButcherTableauDormandPrince );
                }
        streamlog_out( MESSAGE4) << " EUTelGBLFitter::EUTelGBLFitter constructor over" << std::endl;
   }

    EUTelGBLFitter::~EUTelGBLFitter() {
    }
    
    void EUTelGBLFitter::setParamterIdPlaneVec( const std::vector<int>& vector)
    {
      _parameterIdPlaneVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdXResolutionVec( const std::vector<float>& vector)
    {
      _parameterIdXResolutionVec = vector;
    }
 
    void EUTelGBLFitter::setParamterIdYResolutionVec( const std::vector<float>& vector)
    {
      _parameterIdYResolutionVec = vector;
    }
       

    void EUTelGBLFitter::setParamterIdXRotationsMap( const std::map<int, int>& map ) {
        _parameterIdXRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdYRotationsMap( const std::map<int, int>& map ) {
        _parameterIdYRotationsMap = map;
    }
    
    void EUTelGBLFitter::setParamterIdZRotationsMap( const std::map<int, int>& map ) {
        _parameterIdZRotationsMap = map;
    }

    void EUTelGBLFitter::setParamterIdZShiftsMap( const std::map<int, int>& map ) {
        _parameterIdZShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdYShiftsMap( const std::map<int, int>& map ) {
        _parameterIdYShiftsMap = map;
    }

    void EUTelGBLFitter::setParamterIdXShiftsMap( const std::map<int, int>& map ) {
        _parameterIdXShiftsMap = map;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdXRotationsMap() const {
        return _parameterIdXRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdYRotationsMap() const {
        return _parameterIdYRotationsMap;
    }
    
    const std::map<int, int>& EUTelGBLFitter::getParamterIdZRotationsMap() const {
        return _parameterIdZRotationsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdZShiftsMap() const {
        return _parameterIdZShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdYShiftsMap() const {
        return _parameterIdYShiftsMap;
    }

    const std::map<int, int>& EUTelGBLFitter::getParamterIdXShiftsMap() const {
        return _parameterIdXShiftsMap;
    }

    void EUTelGBLFitter::setMEstimatorType( const std::string& mEstimatorType ) {
        std::string mEstimatorTypeLowerCase = mEstimatorType;
        std::transform( mEstimatorType.begin(), mEstimatorType.end(), mEstimatorTypeLowerCase.begin(), ::tolower);
        
        if ( mEstimatorType.size() != 1 ) {
            streamlog_out( WARNING1 ) << "More than one character supplied as M-estimator option" << std::endl;
            streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
            return;
        }
        
        if ( mEstimatorType.compare("t") == 0 ||
             mEstimatorType.compare("h") == 0 ||
             mEstimatorType.compare("c") == 0   ) this->_mEstimatorType = _mEstimatorType;
        else {
            streamlog_out( WARNING1 ) << "M-estimator option " << mEstimatorType << " was not recognized" << std::endl;
            streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
        }
    }

    std::string EUTelGBLFitter::getMEstimatorType( ) const {
        return _mEstimatorType;
    }

    const std::map<long, int>& EUTelGBLFitter::getHitId2GblPointLabel( ) const {
        return _hitId2GblPointLabel;
    }

    void EUTelGBLFitter::setExcludeFromFitPlanes( const std::vector<int>& excludedPlanes ) {
        this->_excludeFromFit = excludedPlanes;
    }

    std::vector<int> EUTelGBLFitter::getExcludeFromFitPlanes( ) const {
        return _excludeFromFit;
    }

    TMatrixD EUTelGBLFitter::propagatePar(double ds) {
        /* for GBL:
           Jacobian for straight line track
           track = q/p, x', y', x, y
                    0,   1,  2, 3, 4
         */
        _parPropJac.UnitMatrix();
        _parPropJac[3][1] = ds; // x = x0 + xp * ds
        _parPropJac[4][2] = ds; // y = y0 + yp * ds
        return _parPropJac;
    }
     
    /** Propagation jacobian in inhomogeneous magnetic field
     * 
     * 
     * @param ds        Z - Z0 propagation distance
     * 
     * @return          transport matrix
     */
    
    TMatrixD EUTelGBLFitter::PropagatePar( double ds, double invP, double tx0, double ty0, double x0, double y0, double z0 ) {
        // for GBL:
        //   Jacobian for helical track
        //   track = q/p, x'  y'  x  y
        //            0,   1,  2, 3, 4
        //
        
        streamlog_out( DEBUG2 ) << "EUTelGBLFitter::PropagatePar()" << std::endl;
	// The formulas below are derived from equations of motion of the particle in
        // magnetic field under assumption |dz| small. Must be valid for |dz| < 20-30 cm

        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
        const double Bx         = B.at( vectorGlobal ).x();
        const double By         = B.at( vectorGlobal ).y();
        const double Bz         = B.at( vectorGlobal ).z();
	const double mm = 1000.;
	const double k = 0.299792458/mm;

        const double sqrtFactor = sqrt( 1. + tx0*tx0 + ty0*ty0 );

	const double Ax = sqrtFactor * (  ty0 * ( tx0 * Bx + Bz ) - ( 1. + tx0*tx0 ) * By );
	const double Ay = sqrtFactor * ( -tx0 * ( ty0 * By + Bz ) + ( 1. + ty0*ty0 ) * Bx );

	// Partial derivatives
	const double dAxdty0 = ty0 * Ax / (sqrtFactor*sqrtFactor) + sqrtFactor*( tx0*Bx + Bz );
	const double dAydtx0 = tx0 * Ay / (sqrtFactor*sqrtFactor) + sqrtFactor*( -ty0*By - Bz );

	const double dxdtx0 = ds;
	const double dxdty0 = 0.5 * invP * k * ds*ds * dAxdty0;

	const double dydtx0 = 0.5 * invP * k * ds*ds * dAydtx0;
	const double dydty0 = ds;

	const double dtxdty0 = invP * k * ds * dAxdty0;
	const double dtydtx0 = invP * k * ds * dAydtx0;

	const double dxdinvP0 = 0.5 * k * ds*ds * Ax;
	const double dydinvP0 = 0.5 * k * ds*ds * Ay;

	const double dtxdinvP0 = k * ds * Ax;
	const double dtydinvP0 = k * ds * Ay;

	// Fill-in matrix elements
	_parPropJac.UnitMatrix();
	_parPropJac[1][0] = dtxdinvP0;  _parPropJac[1][2] = dtxdty0;	
	_parPropJac[2][0] = dtydinvP0;  _parPropJac[2][1] = dtydtx0;	
	_parPropJac[3][0] = dxdinvP0;   _parPropJac[3][1] = dxdtx0;	_parPropJac[3][2] = dxdty0;	
	_parPropJac[4][0] = dydinvP0;   _parPropJac[4][1] = dydtx0;	_parPropJac[4][2] = dydty0;	
        
        if ( streamlog_level(DEBUG0) ){
             streamlog_out( DEBUG0 ) << "Propagation jacobian: " << std::endl;
            _parPropJac.Print();
        }
	
        if( streamlog_level(DEBUG4) ) 
        {
             streamlog_out(DEBUG4) << "step:" << ds ;
             streamlog_out(DEBUG4) << " to point x: " << std::setw(8) << x0 ;
             streamlog_out(DEBUG4) << " y: " << std::setw(8) << y0 ;
             streamlog_out(DEBUG4) << " z: " << std::setw(8) << z0 << std::endl;
             streamlog_out(DEBUG4) << " B=0:                  " ;
             streamlog_out(DEBUG4) << "                              and genB>0 : " << std::endl;
             for(int i =0;i<5;i++)
             {
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][0] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][1] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][2] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][3] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << propagatePar(ds)[i][4] ;
             streamlog_out(DEBUG4) << "                 " ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][0] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][1] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][2] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][3] ;
             streamlog_out(DEBUG4) << " " << std::setw(8) << _parPropJac[i][4] <<std::endl;
             }
        }

        streamlog_out( DEBUG2 ) << "-----------------------------EUTelGBLFitter::PropagatePar()-------------------------------" << std::endl;   

        return _parPropJac;
    }

    /**
     * Get extrapolated position of the track in global coordinate system
     * Formulas are also valid for vanishing magnetic field.
     * Calculation is based on numerical integration of equations of motion
     * 
     * @param ts track state
     * @param dz propagation distance along z
     * @return vector of track parameters in the global coordinate system
     */
    TVectorD EUTelGBLFitter::getXYZfromDzNum( double invP, double tx, double ty, double x0, double y0, double z0, double dz ) const {
        streamlog_out(DEBUG2) << "EUTelGBLFitter::getXYZfromDzNum()" << std::endl;
        
        // Get magnetic field vector
        gear::Vector3D vectorGlobal( x0, y0, z0 );        // assuming uniform magnetic field running along X direction
	const gear::BField&   B = geo::gGeometry().getMagneticFiled();
        const double bx         = B.at( vectorGlobal ).x();
        const double by         = B.at( vectorGlobal ).y();
        const double bz         = B.at( vectorGlobal ).z();

        TVector3 hVec(bx,by,bz);
        
        // Setup the equation
        TVectorD trackStateVec(5);
	trackStateVec[0] = x0;
	trackStateVec[1] = y0;
	trackStateVec[2] = tx;
	trackStateVec[3] = ty;
	trackStateVec[4] = invP;
        static_cast< eom::EOMODE* >(_eomODE)->setInitValue( trackStateVec );
        static_cast< eom::EOMODE* >(_eomODE)->setBField( hVec );
        
        // Integrate
        TVectorD result = _eomIntegrator->integrate( dz );
        
        streamlog_out(DEBUG0) << "Result of the integration:" << std::endl;
        streamlog_message( DEBUG0, result.Print();, std::endl; );
        
        streamlog_out(DEBUG2) << "---------------------------------EUTelGBLFitter::getXYZfromDzNum()------------------------------------" << std::endl;
        
        streamlog_out(DEBUG4) << "getXYZfromDzNum   dz:" << dz ;
        streamlog_out(DEBUG4) << " to point x0: " << std::setw(8) << x0 ;
        streamlog_out(DEBUG4) << " y0: " << std::setw(8) << y0 ;
        streamlog_out(DEBUG4) << " z0: " << std::setw(8) << z0 ;
        streamlog_out(DEBUG4) << " tx: " << std::setw(8) << tx ;
        streamlog_out(DEBUG4) << " ty: " << std::setw(8) << ty << std::endl;
       
        return result;
    }
 
    void EUTelGBLFitter::SetTrackCandidates( vector<const IMPL::TrackImpl*>& trackCandidatesVec) {

        this->_trackCandidatesVec = trackCandidatesVec	;
        return;
    }
       
    void EUTelGBLFitter::SetTrackCandidates( const EVENT::TrackVec& trackCandidates) {

        this->_trackCandidates = trackCandidates;
        return;
    }


    void EUTelGBLFitter::Clear() {

        EUTelTrackFitter::Clear();

        std::map< int, gbl::GblTrajectory* >::iterator it;
        for (it = _gblTrackCandidates.begin(); it != _gblTrackCandidates.end(); ++it) delete it->second;
        _gblTrackCandidates.clear();
        
        _hitId2GblPointLabel.clear();
        _hitId2GblPointLabelMille.clear();
    }


    // @TODO iplane, xPred, yPred must not be here. consider refactoring

    /** Add alignment derivative necessary for MILLIPEDE
     * 
     * @param point
     * @param alDer matrix of global parameter (alignment constants) derivatives. This relates this particular point to this particular sensor orientation. This completely depends on the state variables in LOCAL frame. 
     * @param globalLabels vector of alignment parameters ids. These are used to associate each alignment constant (shift x, y z and rotations) on each plane with a unique number. This is done so that if a particular alignment parameter i.e particular number is changed. The MILLEPEDE knows that all points hit measurements must change by the amount that is calculated by there particular jacobain. Note that for each point this change will be different due to different incidence angles etc.
     * @param iPlane plane id
     * @param xPred predicted by hit x-position (first approximation)
     * @param yPred predicted by hit y-position (first approximation)
     * @param xSlope predicted x-slope
     * @param ySlope predicted y-slope
     */
    void EUTelGBLFitter::addGlobalParametersGBL(gbl::GblPoint& point, TMatrixD& alDer, std::vector<int>& globalLabels, int iPlane,
            const double* predpos, double xSlope, double ySlope) {

        streamlog_out(MESSAGE1) << " addGlobalParametersGBL  -------------------- BEGIN ------------------" << std::endl;

        alDer[0][0] = -1.0; // dx/dx
        alDer[0][1] =  0.0; // dx/dy
        alDer[1][0] =  0.0; // dy/dx
        alDer[1][1] = -1.0; // dy/dy
        globalLabels[0] = _parameterIdXShiftsMap[iPlane]; // dx
        globalLabels[1] = _parameterIdYShiftsMap[iPlane]; // dy


        if (_alignmentMode == Utility::XYShiftXYRot
                || _alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYShiftXZRotXYRot
                || _alignmentMode == Utility::XYShiftYZRotXYRot
                || _alignmentMode == Utility::XYShiftXZRotYZRotXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][2] = -predpos[1]; // dx/rot
            alDer[1][2] =  predpos[0]; // dy/rot
            globalLabels[2] = _parameterIdZRotationsMap[iPlane]; // rot z
        }


        if (_alignmentMode == Utility::XYZShiftXYRot
                || _alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][3] =   -xSlope; // dx/dz
            alDer[1][3] =   -ySlope; // dy/dz
            globalLabels[3] = _parameterIdZShiftsMap[iPlane]; // dz
        }

        if (_alignmentMode == Utility::XYZShiftXZRotYZRotXYRot) {
            alDer[0][4] =   predpos[0]*xSlope; // dx/rot y
            alDer[1][4] =   predpos[0]*ySlope; // dy/rot y
            globalLabels[4] = _parameterIdYRotationsMap[iPlane]; // drot y  - actually X?
            alDer[0][5] =  -predpos[1]*xSlope; // dx/rot x          
            alDer[1][5] =  -predpos[1]*ySlope; // dy/rot x         
            globalLabels[5] = _parameterIdXRotationsMap[iPlane]; // drot x  - actually Y?
        }


// partial alignment 
        if (_alignmentMode == Utility::XYShiftXZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
        }
        if (_alignmentMode == Utility::XYShiftYZRotXYRot) {
            alDer[0][3] = predpos[1]*xSlope; // dx/rot x
            alDer[1][3] = predpos[1]*ySlope; // dy/rot x
            globalLabels[3] = _parameterIdXRotationsMap[iPlane]; // drot x
        }
 
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _parameterIdXRotationsMap[iPlane]; // drot x
        }
        if (_alignmentMode == Utility::XYShiftXZRotYZRotXYRot) {
            alDer[0][3] = predpos[0]*xSlope; // dx/rot y
            alDer[1][3] = predpos[0]*ySlope; // dy/rot y
            globalLabels[3] = _parameterIdYRotationsMap[iPlane]; // drot y
            alDer[0][4] = predpos[1]*xSlope; // dx/rot x
            alDer[1][4] = predpos[1]*ySlope; // dy/rot x
            globalLabels[4] = _parameterIdXRotationsMap[iPlane]; // drot x
        }

        point.addGlobals(globalLabels, alDer);
 
        streamlog_out(MESSAGE1) << " addGlobalParametersGBL -------------------- END ------------------ " << std::endl;
   }

    void EUTelGBLFitter::pushBackPoint( std::vector< gbl::GblPoint >& pointListTrack, const gbl::GblPoint& pointTrack, int hitid ) {
        pointListTrack.push_back(pointTrack);
       
        streamlog_out(MESSAGE0) << endl << "pushBackPoint: " << hitid << std::endl;
        // store point's GBL label for future reference
        _hitId2GblPointLabel[ hitid ] = static_cast<int>(pointListTrack.size());
    }

    void EUTelGBLFitter::pushBackPointandState( std::vector< gbl::GblPoint >* pointListTrack, const gbl::GblPoint pointTrack, EUTelTrackStateImpl *state) {
        pointListTrack->push_back(pointTrack);
       
        streamlog_out(MESSAGE0) << endl << "pushBackPoint: " << pointListTrack->size() <<  std::endl;
        // store point's GBL label for future reference
				const gbl::GblPoint *pointer_to_gblpoint = &pointTrack;
        _PointToState[ pointer_to_gblpoint ] = state; 
    }


    void EUTelGBLFitter::pushBackPointMille( std::vector< gbl::GblPoint >& pointListMille, const gbl::GblPoint& pointMille, int hitid ) {
        pointListMille.push_back(pointMille);
        
        streamlog_out(MESSAGE0) << endl << "pushBackPointMille: " << hitid << std::endl;
        // store point's GBL label for future reference
        _hitId2GblPointLabelMille[ hitid ] = static_cast<int>(pointListMille.size());
    }




void EUTelGBLFitter::UpdateTrackFromGBLTrajectory (gbl::GblTrajectory* traj, std::vector< gbl::GblPoint >* pointList){

	int number_of_points = pointList->size();
	int pointNum = 0;
	for(pointNum=0; pointNum < number_of_points; ++pointNum){ //Must make sure that the label starts at ????????
		EUTelTrackStateImpl * state = 	_PointToState[ &(pointList->at(pointNum)) ]; //get the state associated with this point
		if(state != NULL){

			TVectorD corrections;
			TMatrixDSym correctionsCov;
      traj->getResults(pointNum, corrections, correctionsCov );

     	streamlog_out(MESSAGE0) << endl << "State before we have added corrections: " << std::endl;
			state->Print();
			state->setX( state->getX() + corrections[0]);
			state->setY( state->getY() + corrections[1]);
			state->setTx( state->getTx() + corrections[2]);
			state->setTy( state->getTy() + corrections[3]);
			state->setInvP( state->getInvP() + corrections[4]);
     	streamlog_out(MESSAGE0) << endl << "State after we have added corrections: " << std::endl;
			state->Print();

			
		}

	}//END OF LOOP OVER POINTS

}

    /**
     * merging too sources of information: a copy of a Track Candidate (*itTrkCand) gets updated with traj
     * 
     * @param LCIOtrack pointer to track object to be stored
     * @param chi2     Chi2 of the track fit
     * @param ndf      NDF of the track fit
     */
     void EUTelGBLFitter::prepareLCIOTrack( gbl::GblTrajectory* gblTraj, vector<const IMPL::TrackImpl*>::const_iterator& itTrkCand, 
                                          double chi2, int ndf ){ //, double omega, double d0, double z0, double phi, double tanlam ) {

        streamlog_out( DEBUG4 ) << " ----------------  EUTelGBLFitter::prepareLCIOTrack -- BEGIN ------------- " << std::endl;

// output  track

        IMPL::TrackImpl * LCIOtrack = new IMPL::TrackImpl( **itTrkCand); 
 
        unsigned int numData;
        TVectorD corrections(5);
        TMatrixDSym correctionsCov(5);
 
        TVectorD residual(2);
        TVectorD measErr(2);
        TVectorD residualErr(2);
        TVectorD downWeight(2);
 
        streamlog_out(MESSAGE1) << endl; 
 
        int nstates = LCIOtrack->getTrackStates().size();
        streamlog_out(MESSAGE1) << "states " << nstates << "    " << endl;

        for(int i=0;i < nstates; i++) 
            {
                const IMPL::TrackStateImpl* const_trkState = static_cast <const IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i) ) ;

                double fitPointLocal[] = {0.,0.,0.};
                fitPointLocal[0] = const_trkState->getReferencePoint()[0] ;
                fitPointLocal[1] = const_trkState->getReferencePoint()[1] ;
                fitPointLocal[2] = const_trkState->getReferencePoint()[2] ;

                double bd0        = const_trkState->getD0() ;
 	        double bphi       = const_trkState->getPhi() ;
                double bomega     = const_trkState->getOmega() ;
	        double btanlambda = const_trkState->getTanLambda() ;
	        double bz0        = const_trkState->getZ0() ;

                const int hitGblLabel = _hitId2GblPointLabel[ const_trkState->id() ];

                streamlog_out(MESSAGE1) << hitGblLabel << " corr: " << const_trkState->id() 
                                        << " [d0]" << std::setw(6) << bd0 << ":" 
                                        << " [phi]" << std::setw(6) << bphi << ":" 
                                        << " [ome]" << std::setw(6) << bomega << ":" 
                                        << " [tanl]" << std::setw(6) << btanlambda << ":" 
                                        << " [z0]" << std::setw(6) << bz0 << ":" << " ["<< setw(3) << const_trkState->getLocation() <<"]" 
                                        << " points:"  << setw(8) << fitPointLocal[0] << setw(8) << fitPointLocal[1] << setw(8) << fitPointLocal[2] << " "  ;   

               streamlog_out(MESSAGE1) << endl; 
            }
 
          streamlog_out(MESSAGE1) << endl; 
  
        const EVENT::TrackerHitVec& ihits = (*itTrkCand)->getTrackerHits();
        int itrk = (*itTrkCand)->id();
        int nhits =  ihits.size( ) ;
        int expec = _parameterIdPlaneVec.size();
        streamlog_out(MESSAGE1) <<  " track itrk:" <<  itrk  << " with " << nhits << " at least " << expec << "       " ;//std::endl;
        for (unsigned int i = 0; i< ihits.size(); i++ ) 
        { 
            EVENT::TrackerHit* ihit = ihits[i];
            int ic = ihit->id();
            streamlog_out(MESSAGE0) <<  ic << " ";
        }
        streamlog_out(MESSAGE1) << std::endl;
 
          for(int i=0;i < nstates; i++) 
            {
                const IMPL::TrackStateImpl* const_trkState = static_cast <const IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i) ) ;
                      IMPL::TrackStateImpl* trkState = static_cast <      IMPL::TrackStateImpl*> ( LCIOtrack->getTrackStates().at(i) ) ;

                // first get the GBL trajectory fit results:
                const int hitGblLabel = _hitId2GblPointLabel[ const_trkState->id() ];

                _hitId2GblPointLabel.insert( std::make_pair(   trkState->id() , hitGblLabel  ) );

                gblTraj->getResults( hitGblLabel, corrections, correctionsCov );

                float fitPointLocal[] = {0.,0.,0.};
                fitPointLocal[0] = const_trkState->getReferencePoint()[0] ;
                fitPointLocal[1] = const_trkState->getReferencePoint()[1] ;
                fitPointLocal[2] = const_trkState->getReferencePoint()[2] ;

                double bd0        = const_trkState->getD0() ;
 	        double bphi       = const_trkState->getPhi() ;
                double bomega     = const_trkState->getOmega() ;
	        double btanlambda = const_trkState->getTanLambda() ;
	        double bz0        = const_trkState->getZ0() ;

                // apply GBL fit results:
                trkState->setD0(        const_trkState->getD0() + corrections[0] ) ;
 	        trkState->setPhi(       const_trkState->getPhi() + corrections[1] ) ;
                trkState->setOmega(     const_trkState->getOmega() + corrections[2] ) ;
	        trkState->setTanLambda( const_trkState->getTanLambda() + corrections[3] ) ;
	        trkState->setZ0(        const_trkState->getZ0() + corrections[4] ) ;

                double ed0        = trkState->getD0() ;
 	        double ephi       = trkState->getPhi()  ;
                double eomega     = trkState->getOmega()  ;
	        double etanlambda = trkState->getTanLambda()  ;
	        double ez0        = trkState->getZ0()  ;



                streamlog_out(MESSAGE1) << hitGblLabel << " corr: " << trkState->id() 
                                        << " [d0]" << std::setw(6)  << ed0 <<  ":"
                                        << " [phi]" << std::setw(6)  << ephi<<  ":"
                                        << " [ome]" << std::setw(6)  << eomega<<  ":"
                                        << " [tanl]" << std::setw(6)  << etanlambda<<  ":"
                                        << " [z0]" << std::setw(6)  << ez0 <<  ":" << " ["<< setw(3) << trkState->getLocation() <<"]" ;      

                if( trkState->getLocation() >= 0 ) {
                  gblTraj->getMeasResults( hitGblLabel, numData, residual, measErr, residualErr, downWeight);
                  // correct original values to the fitted ones

                  bool foundAhit = false; 
                  for (unsigned int i = 0; i< ihits.size(); i++ ) 
                  { 
                      EVENT::TrackerHit* ihit = ihits[i];
                      const int planeID = Utility::getSensorIDfromHit( static_cast< IMPL::TrackerHitImpl* >(ihit) );
                      if( planeID == trkState->getLocation()) {
                         foundAhit = true;   
                         const float* opoint = Utility::HitCDoubleShiftCFloat(ihit->getPosition(), residual );

                         trkState->setReferencePoint( opoint );                        
                      } 
                  }
 

                  double PointLocal[] = {0.,0.,0.};
                  PointLocal[0] = trkState->getReferencePoint()[0] ;
                  PointLocal[1] = trkState->getReferencePoint()[1] ;
                  PointLocal[2] = trkState->getReferencePoint()[2] ;

                  streamlog_out(MESSAGE1) << " points:"  << setw(8) << PointLocal[0] << setw(8)  << PointLocal[1] << setw(8) << PointLocal[2] << " "  ;   
                  streamlog_out(MESSAGE1) << "  re[0]:"  << std::setw(7) << residual[0] << " re[1]:"  << std::setw(7) << residual[1] ;
                }

                streamlog_out(MESSAGE1) << endl;      

            }
            streamlog_out(MESSAGE1) << std::endl;

/*
        const EVENT::TrackerHitVec& chits = LCIOtrack->getTrackerHits();
        itrk = LCIOtrack->id();
        nhits =  chits.size( ) ;
        expec =  _parameterIdPlaneVec.size();
        streamlog_out(MESSAGE1) <<  " track itrk:" <<  itrk  << " with " << nhits << " at least " << expec  << "       ";

        for (int i = 0; i< chits.size(); i++ ) 
        { 
            EVENT::TrackerHit* ihit = chits[i];
            int ic = ihit->id();
            streamlog_out(MESSAGE0) <<  ic << " ";
        }
        streamlog_out(MESSAGE1) << std::endl;
 */
         // prepare track
        LCIOtrack->setChi2 ( chi2 );      // Chi2 of the fit (including penalties)
        LCIOtrack->setNdf  ( ndf );        // Number of planes fired (!)
        
        // add track to LCIO collection vector

        if( getFitTrackVec() != 0 ) getFitTrackVec()->push_back( LCIOtrack );
 
        streamlog_out( DEBUG4 ) << " ----------------  EUTelGBLFitter::prepareLCIOTrack -- END ------------- " << std::endl;

    }

    void EUTelGBLFitter::FitSingleTrackCandidate(EVENT::TrackVec::const_iterator& itTrkCand)
    {
            // relies on sane itTrkCand -> sanity checks ?

	    // Z axis points along beam direction.
	    double pt = ( 1./(*itTrkCand)->getOmega() ) * _beamQ;
	    double px = (*itTrkCand)->getTanLambda() * pt;
	    double py = pt * sin( (*itTrkCand)->getPhi() );
	    double pz = pt * cos( (*itTrkCand)->getPhi() );

	    double tx   = px / pz;
	    double ty   = py / pz;
	    double invP = _beamQ / sqrt( px*px + pt*pt );
 
            streamlog_out(DEBUG4) << "::FitSingleTrackCandidate FitTracks   ";
            streamlog_out(DEBUG4) << " omega: " << std::setw(8) << (*itTrkCand)->getOmega()  ;
            streamlog_out(DEBUG4) << " Q: " << std::setw(8) << _beamQ ;
            streamlog_out(DEBUG4) << " px: " << std::setw(8) << px ;
            streamlog_out(DEBUG4) << " py: " << std::setw(8) << py ;
            streamlog_out(DEBUG4) << " pz: " << std::setw(8) << pt << std::endl;
 
	    // Reference point is the last hit on a track
	    const float *refPoint = (*itTrkCand)->getReferencePoint();

            // remember last hit-point from the track candidate below
            double start[3] = { 0.,0.,-500.}; // 
 
            // Loop over hits on a track candidate
            const EVENT::TrackerHitVec& hits = (*itTrkCand)->getTrackerHits();
            EVENT::TrackerHitVec::const_reverse_iterator itHit;

      const map< int, int > sensorMap = geo::gGeometry().sensorZOrdertoIDs();
      int planeID     = sensorMap.at(0); // the first first plane in the array of the planes according to z direction. // assume not tilted plane. 

      const double thicknessSen       = geo::gGeometry().siPlaneZSize( planeID );
      start[2]    = geo::gGeometry().siPlaneZPosition(planeID) - thicknessSen ; //initial z position to the most-first plane

      double dir[3]   = {0.,0.,1.};  // as good as any other value along z axis.
      float dpoint[3] = {0.,0.,0.};  // initialise output point-vector (World frame)
      float lpoint[3] = {0.,0.,0.};  // initialise output point-vector (sensor frame)

      streamlog_out ( DEBUG0 ) << "::FitSingleTrackCandidate starting point: "  <<  start[0] << " "  << start[1] << " " << start[2] << " dir:"<< dir[0]<<" "<<dir[1]<<" "<<dir[2] << endl;
 
      map< int, int >::const_iterator iter = sensorMap.begin();
      while ( iter != sensorMap.end() ) {
        bool found = false;
        if( iter->first > -999 )
          {
            int getErrorID = geo::gGeometry( ).findNextPlaneEntrance(  start,  dir, iter->second, dpoint ) ;

            // get to local frame for sensor newSensorID:
//            geo::gGeometry().master2Local(dpoint, lpoint);


            if( getErrorID < 0 ) 
            {
              streamlog_out ( DEBUG4 ) << "no Entrance: " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] << " err:"<< getErrorID << endl;
            } else {     
              streamlog_out ( DEBUG0 ) << "identified NextPlane Entrance: " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] <<  " err:"<< getErrorID << endl;
              start[0] = dpoint[0];
              start[1] = dpoint[1];
              start[2] = dpoint[2];
              for ( itHit = hits.rbegin(); itHit != hits.rend(); ++itHit) {

                 const int planeID = Utility::getSensorIDfromHit( static_cast< IMPL::TrackerHitImpl* >(*itHit) );

                 if(planeID == getErrorID)
                 {
                   streamlog_out ( DEBUG4 ) << "match sensor found: extrapl@ " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] <<  " @"<< getErrorID << endl;
                  //
                  const double* hitPointLocal = (*itHit)->getPosition();
                  double hitPointGlobal[] = {0.,0.,0.};
                  geo::gGeometry().local2Master(planeID,hitPointLocal,hitPointGlobal);
   		  streamlog_out(DEBUG4) << "hitm= " << hitPointGlobal[0] << " " << hitPointGlobal[1] << " " << hitPointGlobal[2] << std::endl;
                  found = true;
                 }
              }
              if(!found) 
              {
                streamlog_out ( DEBUG4 ) << "identified NextPlane UNMATCHED " <<  dpoint[0] << " " <<  dpoint[1] << " " << dpoint[2] <<  " err:"<< getErrorID << endl;
              } 
            }
          }
        ++iter;
      } 
    }
 		//getHmatrix is local to global. So we need global to local or using gbl library terminology. local (our global_ to measurement frame (our local). So we take the inverse

// GBL Trajectory treatment ::  Fit and dump into LCIO
    void EUTelGBLFitter::PerformFitGBLTrajectory( gbl::GblTrajectory* gblTraj, vector<const IMPL::TrackImpl*>::const_iterator& itTrkCand  ) {
                
            streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::PerformFitGBLTrajectory -- starting " << endl;

            double loss = 0.;
            double chi2 = 0.;
            int ndf = 0;
            // perform GBL fit

                int ierr = 0;

		        if ( streamlog_level(MESSAGE0) ){
	        	  std::cout << "pre fit FitTrack - trajectory: " << std::endl;
//		 	  gblTraj->printTrajectory(1);
//			  gblTraj->printData();
//	 		  gblTraj->printPoints(1);
	        	  std::cout << "pre fit FitTrack - trajectory: end;" << std::endl;

		        }
	
                if ( !_mEstimatorType.empty( ) ) ierr = gblTraj->fit( chi2, ndf, loss, _mEstimatorType );
                else ierr = gblTraj->fit( chi2, ndf, loss );

                streamlog_out(MESSAGE0) << "ierr : "<< ierr << " and chi2: " << chi2 << std::endl;

                if ( ierr  )
                    {
		        if ( streamlog_level(MESSAGE0) ){
	        	  std::cout << "after fit FitTrack - trajectory: " << std::endl;
//		 	  gblTraj->printTrajectory(1);
//			  gblTraj->printData();
//	 		  gblTraj->printPoints(1);
	        	  std::cout << "after fit FitTrack - trajectory: " << std::endl;
  		        }
		    }
  
                    // for some reason (??) need to keep the same numbering for trajectories as for the Track Candidates
                    vector<const IMPL::TrackImpl*>::const_iterator begin = _trackCandidatesVec.begin();
                    _gblTrackCandidates.insert( std::make_pair( std::distance( begin, itTrkCand ), gblTraj ) );
                
                    // Write fit result
                    prepareLCIOTrack( gblTraj, itTrkCand, chi2, ndf);

                    streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::PerformFitGBLTrajectory -- finished " << endl;
           
    }


void EUTelGBLFitter::CreateTrajectoryandFit(std::vector< gbl::GblPoint >* pointList,  gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr){
	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::CreateTrajectoryandFit -- BEGIN " << endl;

	double loss = 0.;

	if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( *chi2, *ndf, loss, _mEstimatorType );
  else ierr = traj->fit( *chi2, *ndf, loss );

	if( ierr != 0 ){
		streamlog_out(MESSAGE0) << "Fit failed! Track error: "<< ierr << " and chi2: " << *chi2 << std::endl;
	}
	else{
  streamlog_out(MESSAGE0) << "Fit Successful! Track error; "<< ierr << " and chi2: " << *chi2 << std::endl;
	}
	

	streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::CreateTrajectoryandFit -- END " << endl;
}







// convert input TrackCandidates and TrackStates into a GBL Trajectory
void EUTelGBLFitter::FillInformationToGBLPointObject(EUTelTrackImpl* EUtrack, std::vector< gbl::GblPoint >* pointList){
	// sanity check. Mustn't happen in principle. That the number of hits is greater than the number of hits
  if ( EUtrack->getTrackerHits().size() > geo::gGeometry().nPlanes() ){
  	streamlog_out(ERROR) << "Sanity check. This should not happen in principle. Number of hits is greater then number of planes" << std::endl;
   	return;
  }
  //Create the jacobian
  TMatrixD jacPointToPoint(5, 5);
  jacPointToPoint.UnitMatrix();
 	////////////////////////////////////////////////////////////////////////////////////////////////// loop through all states.
  for(int i=0;i < EUtrack->getTrackStates().size(); i++){		
		/////////////////////////////////////////////////////////////////////////////////////////////BEGIN to create GBL point 
		gbl::GblPoint point(jacPointToPoint);
  		EUTelTrackStateImpl* state = EUtrack->getTrackStates().at(i); //get the state for this track. Static cast from EVENT::TrackState to derived class IMPL::TrackStateImpl.
  		streamlog_out(DEBUG3) << "This is the track state being used in creation of GBL points" << std::endl;
			state->Print();

		//Need to find hit that this state may be associated with. Note this is a problem for two reasons. Not all states have a hit. Furthermore we can not associate a hit with a state with the current LCIO format. This must be fixed
		EVENT::TrackerHit* hit = NULL; //Create the hit pointer
		FindHitIfThereIsOne(EUtrack, hit, state); //This will point the hit to the correct hit object associated with this state. If non exists then point it will remain pointed to NULL
		double fitPointLocal[] = {0.,0.,0.};
  	fitPointLocal [0] = state->getReferencePoint()[0] ;
  	fitPointLocal [1] = state->getReferencePoint()[1] ;
  	fitPointLocal [2] = state->getReferencePoint()[2] ;
		addSiPlaneScattererGBL(point, state->getLocation()); //This we still functions still assumes silicon is the thin scatterer. This can be easily changed when we have the correct gear file. However we will always assume that states will come with scattering information. To take into account material between states this will be dealt with latter. 
		double fitPointGlobal[3];
		geo::gGeometry().local2Master( state->getLocation(), fitPointLocal, fitPointGlobal);
		streamlog_out(DEBUG3) << "This is the global position of the track state. Should be the same x,y as above: " <<fitPointGlobal[0]<<","<<fitPointGlobal[1]<<","<<fitPointGlobal[2]<< std::endl;	
		if(hit != NULL){
			double cov[4] ;
			state->setTrackStateHitCov(cov); //This part should not be done in the way it has. MUST FIX! Hit cov should be part of hits own class. Must learn more about LCIO data format
			
			addMeasurementGBL(point, hit->getPosition(),  fitPointLocal, cov, state->getH()); 		
			pushBackPointandState(pointList, point, state);

		}else{
			pushBackPointandState(pointList, point, state);
		}

		////////////////////////////////////////////////////////////////////////////////START TO CREATE SCATTERS BETWEEN PLANES
		if(i != (EUtrack->getTrackStates().size()-1)){
  			EUTelTrackStateImpl* state_next = EUtrack->getTrackStates().at(i+1); //get the next tracks state to determine dz between the two states
  		streamlog_out(DEBUG3) << "This is the track state that is one state ahead" << std::endl;
			state->Print();
			double fitPointLocal_next[] = {0.,0.,0.};  //Need this since we dont save the z parameter as state variable
			fitPointLocal_next [0] = state_next->getReferencePoint()[0] ;
  			fitPointLocal_next [1] = state_next->getReferencePoint()[1] ;
	  		fitPointLocal_next [2] = state_next->getReferencePoint()[2] ;

			double fitPointGlobal_next[3];
			geo::gGeometry().local2Master( state_next->getLocation(), fitPointLocal_next, fitPointGlobal_next );
		streamlog_out(DEBUG3) << "This is the global position of the track state. Should be the same x,y as above: " <<fitPointGlobal_next[0]<<","<<fitPointGlobal_next[1]<<","<<fitPointGlobal_next[2]<< std::endl;	
			float rad = geo::gGeometry().findRadLengthIntegral( fitPointGlobal, fitPointGlobal_next, true ); //We need to skip the volumes that contain the hits since this has already been counted. Must check this functions as expect????
			streamlog_out(DEBUG3) << "This is the radiation length between the two points  " << fitPointGlobal[0]<<","<<fitPointGlobal[1]<<","<<fitPointGlobal[2]<<" and  " <<fitPointGlobal_next[0]<<","<<fitPointGlobal_next[1]<<","<<fitPointGlobal_next[2] <<"Radition length:  "<< rad <<std::endl;

			///////////////////////////////////////////////////////////////////////////////////////////////////////BEGIN THE FIRST SCATTERING PLANE
			//These distances are from the last state plane. There are where the next scatterer should be
			float distance1 = (fitPointGlobal_next[2] + fitPointGlobal[2])/2 - (fitPointGlobal_next[2] - fitPointGlobal[2])/sqrt(12); 
			streamlog_out(DEBUG3) << "This is the distance to the first scatterer: " << distance1 <<std::endl;	
			//Note the distance is used along the track and not from the scattering plane. How should this be dealt with?
			TMatrix jacobianScat1(5,5); jacobianScat1 = state->getPropagationJacobianF(distance1);
			gbl::GblPoint pointScat1(jacobianScat1);
			TVectorD scat(2);
			scat[0] = 0.0; //This should always be 0 right? If not then it should be given as a parameter
			scat[1] = 0.0; 

 			const double scatvariance  = Utility::getThetaRMSHighland(GetBeamEnergy(), rad/2);
			streamlog_out(DEBUG3) << "Scattering mean angle and variance. Average angle: " << scat[0] <<"," <<scat[1] << ".Variance " <<scatvariance<<std::endl;	
			TVectorD scatPrecSensor(2);
 			scatPrecSensor[0] = 1.0 / (scatvariance * scatvariance );
 			scatPrecSensor[1] = 1.0 / (scatvariance * scatvariance );

  		pointScat1.addScatterer(scat, scatPrecSensor);
			pushBackPointandState(pointList, pointScat1, NULL);

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////END THE FIRST SCATTERING PLANE
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////BEGIN THE SECOND SCATTERING PLANE
			float distance2 = (fitPointGlobal_next[2] + fitPointGlobal[2])/2 + (fitPointGlobal_next[2] - fitPointGlobal[2])/sqrt(12);

			TMatrix jacobianScat2(5,5); jacobianScat1 = state->getPropagationJacobianF(distance2);
			gbl::GblPoint pointScat2(jacobianScat1);


  		pointScat2.addScatterer(scat, scatPrecSensor);
			pushBackPointandState(pointList, pointScat2, NULL);
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////END OF SECOND SCATTERING PLANE
			jacPointToPoint = state->getPropagationJacobianF(fitPointGlobal_next[2] - fitPointGlobal[2]); 				

		}  
		/////////////////////////////////////////////////////////////////////////////////////////END OF CREATE SCATTERERS BETWEEN PLANES
	
	}//END OF LOOP THROUGH ALL PLANES
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
}


void EUTelGBLFitter::addSiPlaneScattererGBL(gbl::GblPoint& point, int iPlane) {

	streamlog_out(MESSAGE1) << " addSiPlaneScattererGBL ------------- BEGIN --------------  " << std::endl;


	TVectorD scatPrecSensor(2);
	TVectorD scat(2); 
	scat[0] = 0.0; scat[1]=0.0; //This should always be 0 right? If not then it should be given as a parameter
	const double radlenSi           = geo::gGeometry().siPlaneRadLength(iPlane);
	const double thicknessSi        = geo::gGeometry().siPlaneZSize(iPlane);

	const double X0Si  = thicknessSi / radlenSi; // Si 
               
	
	const double tetSi  = Utility::getThetaRMSHighland(GetBeamEnergy(), X0Si);

        scatPrecSensor[0] = 1.0 / (tetSi * tetSi );
        scatPrecSensor[1] = 1.0 / (tetSi * tetSi );

        point.addScatterer(scat, scatPrecSensor);

	streamlog_out(MESSAGE1) << " addSiPlaneScattererGBL  ------------- END ----------------- " << std::endl;
}



//This will add the measurement of the hit and predicted position. Using the covariant matrix of the hit. NOT! the residual.
void EUTelGBLFitter::addMeasurementGBL(gbl::GblPoint& point, const double *hitPos, const double *statePos, double hitCov[4], TMatrixD HMatrix){
     
	streamlog_out(MESSAGE1) << " addMeasurementsGBL ------------- BEGIN --------------- " << std::endl;

 	TVectorD meas;
	meas[0] = hitPos[0] - statePos[0];
        meas[1] = hitPos[1] - statePos[1];

	TMatrixDSym measPrec(2,2); //Precision matrix is symmetric. The vector option that was here was silly since there could be correlation between variance and x/y.
        measPrec[0][0] = 1. / hitCov[0];	// cov(x,x)
        measPrec[1][1] = 1. / hitCov[2];	// cov(y,y)
	measPrec[0][1] = 1. / hitCov[1];  //cov(x,y)

	streamlog_out(DEBUG4) << "Residuals and covariant matrix for the hit:" << std::endl;
        streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0][0] << measPrec[0][1] << std::endl;
        streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20) << measPrec[1][0] << measPrec[1][1] << std::endl;

        point.addMeasurement(HMatrix, meas, measPrec);

	streamlog_out(MESSAGE1) << " addMeasurementsGBL ------------- END ----------------- " << std::endl;
}

void EUTelGBLFitter::FindHitIfThereIsOne(EUTelTrackImpl* EUtrack, EVENT::TrackerHit* hit, EUTelTrackStateImpl* state){
	
	int state_location = state->getLocation(); //Get this states location

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////Loop through all the hits on this track and see of it any are on this states plane using location
	const EVENT::TrackerHitVec& HitOnTrack = EUtrack->getTrackerHits(); //point to these hits by reference
  EVENT::TrackerHitVec::const_iterator itrHit;
	for( itrHit = HitOnTrack.begin(); itrHit != HitOnTrack.end(); ++itrHit){
  	const int planeID = Utility::getSensorIDfromHit( static_cast< IMPL::TrackerHitImpl* >(*itrHit) ); //Get the sensor ID for this hit
			streamlog_out(DEBUG5) << "Hit was on plane " << planeID <<"This state location is: " << state_location << std::endl;
			if(planeID == state_location){
				streamlog_out(DEBUG5) << "Hit and Plane ID the same. Point the hit object to this hit " << std::endl;
				hit = *itrHit;
			}//END OF IF STATEMENT
	}//END OF HIT LOOP
		
}

    /** Add a measurement to GBL point
     * 
     * @param point
     * @param meas measuremet vector (residuals) to be calculated in this routine
     * @param measPrec residuals weights (1/unc^2) to be calculated in this routine
     * @param hitpos hit position
     * @param predpos predicted by hit x-position (approximation)
     * @param hitcov hit covariance matrix
     * @param proL2m projection matrix from track coordinate system onto measurement system
     */
    void EUTelGBLFitter::addMeasurementsGBL(gbl::GblPoint& point, TVectorD& meas, TVectorD& measPrec, const double* hitpos,
            const double* predpos, const EVENT::FloatVec& hitcov, TMatrixD& proL2m) {
     
        streamlog_out(DEBUG4) << " addMeasurementsGBL " << std::endl;
 
        meas[0] = hitpos[0] - predpos[0];
        meas[1] = hitpos[1] - predpos[1];
        measPrec[0] = 1. / hitcov[0];	// cov(x,x)
        measPrec[1] = 1. / hitcov[2];	// cov(y,y)

        streamlog_out(DEBUG4) << "Residuals:" << std::endl;
        streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] << std::endl;
        streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20) << measPrec[1] << std::endl;

        point.addMeasurement(proL2m, meas, measPrec);
    }

    void EUTelGBLFitter::TrackCandidatesToGBLTrajectory( vector<const IMPL::TrackImpl*>::const_iterator& itTrkCand) {


            // sanity check. Mustn't happen in principle.
            if ( (*itTrkCand)->getTrackerHits().size() > geo::gGeometry().nPlanes() )
            {
              streamlog_out(ERROR) << "Sanity check. This should not happen in principle. Number of hits is greater then number of planes" << std::endl;
              return;
            }



            // Z axis points along beam direction.
            double pt = ( 1./(*itTrkCand)->getOmega() ) * _beamQ;
            double px = (*itTrkCand)->getTanLambda() * pt;
   	    double py = pt * sin( (*itTrkCand)->getPhi() );
	    double pz = pt * cos( (*itTrkCand)->getPhi() );

	    double tx   = px / pz;
	    double ty   = py / pz;
//	    double invP = _beamQ / sqrt( px*px + pt*pt ); // leave this hear as an example
 
            double p = _eBeam; // beam momentum
            double invP = 1./p;

            //GBL trajectory construction
            // lsit of GBL points that goes eventually into the trajectory
            std::vector< gbl::GblPoint > pointList;

            // introduce new jacobian from point to point // not sure why we need it here ? 
            //
            TMatrixD jacPointToPoint(5, 5);
            jacPointToPoint.UnitMatrix();

            // get measurement points:
            const EVENT::TrackerHitVec& hits = (*itTrkCand)->getTrackerHits();
            EVENT::TrackerHitVec::const_reverse_iterator itHit;
            unsigned int imatch=0;
            streamlog_out( MESSAGE0 ) << "list requested planes: " ;
            for(unsigned int izPlane=0; izPlane<_parameterIdPlaneVec.size(); izPlane++) {
              int kPlaneID =  _parameterIdPlaneVec[izPlane];
              streamlog_out( MESSAGE0 ) << " [" << imatch << ":" << kPlaneID ;
              bool ifound = false;
              for ( itHit = hits.rbegin(); itHit != hits.rend(); ++itHit) {
              const int planeID = Utility::getSensorIDfromHit( static_cast< IMPL::TrackerHitImpl* >(*itHit) );
               if( kPlaneID == planeID )
                {  
                  imatch++;
                  streamlog_out( MESSAGE0 ) << " yes "   ;
                  ifound = true;
                  break;
                }
              }     
              if( !ifound)  streamlog_out( MESSAGE0 ) << " not "   ;      
              streamlog_out( MESSAGE0 ) << "] "   ;
            }
            streamlog_out( MESSAGE0 ) << std::endl;

            streamlog_out( MESSAGE1 ) << "from the list of requested planes and available hit on a track candidate: " << imatch << "/" <<_parameterIdPlaneVec.size() << " found" << std::endl;

            if( imatch != _parameterIdPlaneVec.size() ) 
            { 
              streamlog_out(MESSAGE1) << " Number of hits does not correspond to the number of planes selected for tracking" << std::endl;
              return;
            }



            TVectorD meas(2);
            TVectorD measPrec(2); // precision = 1/resolution^2

            TVectorD scat(2);
            scat.Zero();

            TVectorD scatPrecSensor(2);
            TVectorD scatPrec(2);
 

            // loop through scattering planes (including measurement planes) ::
            int nstates = (*itTrkCand)->getTrackStates().size();
            for(int i=0;i < nstates; i++) 
            {

                streamlog_out(MESSAGE1) << "state: at " << i << " of " << nstates ;
                IMPL::TrackStateImpl* trk = static_cast < IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i) ) ;
                int trkVolumeID =  trk->getLocation();

                

                if ( trkVolumeID < 0 )
                {
                   streamlog_out(DEBUG0) << "Sanity check. SensorID can not be negative. Negative TrackStates are kept for the beginning and end of the TrackStates on Pattern Recognition. skip this one.";
                   streamlog_out(MESSAGE1) << std::endl; // to be consistent with the MESSAGE1 level printouts in the for-loop
                   continue;
                } 

                streamlog_out(MESSAGE1) << "  [ " << trk->id() << " ] " ;
                streamlog_out(MESSAGE1) << " ref: ";
                streamlog_out(MESSAGE1) << " [0] " <<  trk->getReferencePoint()[0] ;
                streamlog_out(MESSAGE1) << " [1] " <<  trk->getReferencePoint()[1] ;
                streamlog_out(MESSAGE1) << " [2] " <<  trk->getReferencePoint()[2] ;
                streamlog_out(MESSAGE1) << std::endl;

                double fitPointLocal[] = {0.,0.,0.};
                fitPointLocal [0] = trk->getReferencePoint()[0] ;
                fitPointLocal [1] = trk->getReferencePoint()[1] ;
                fitPointLocal [2] = trk->getReferencePoint()[2] ;

                double fitPointGlobal[] = {0.,0.,0.};
                geo::gGeometry().local2Master( trkVolumeID, fitPointLocal, fitPointGlobal);

                // fill next GBL point:
                gbl::GblPoint point(jacPointToPoint);


                if(  trkVolumeID == -1 )
                {
                   // begining of the track - does not exist normally
                }
                else if(  trkVolumeID == -2 )
                {
                   // end of the track - can ignore
                   streamlog_out(MESSAGE1) << std::endl;
                }
                if(  trkVolumeID >= 0 )
                {


                 // look for the measurement point corresponding to the TrackState Location (fitted hit on the volume surface)
                 EVENT::TrackerHitVec::const_iterator itHit;
                 for ( itHit = hits.begin(); itHit != hits.end(); itHit++) {
                  const int planeID = Utility::getSensorIDfromHit( static_cast< IMPL::TrackerHitImpl* >(*itHit) );
                  if( trkVolumeID == planeID )
                  { 
                    EVENT::FloatVec hitcov(4);
                    hitcov[0]=0.01;
                    hitcov[1]=0.00;
                    hitcov[2]=0.01;
                    hitcov[3]=0.00;

                    // check Processor Parameters for plane resolution // Denys
                    if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 )
                    {
                      for(unsigned int izPlane=0;izPlane<_parameterIdPlaneVec.size();izPlane++)
                      {
                        if( _parameterIdPlaneVec[izPlane] == planeID )
                        {  
                          hitcov[0] =  _parameterIdXResolutionVec[izPlane];
                          hitcov[2] =  _parameterIdYResolutionVec[izPlane];

                          hitcov[0] *= hitcov[0]; // squared !
                          hitcov[2] *= hitcov[2]; // squared !
                          break;
                        } 
                      }
                    }
                    else
                    {  
                      hitcov = (*itHit)->getCovMatrix();            
                    }

 
                    const double* hitPointLocal = (*itHit)->getPosition();
                    double hitPointGlobal[] = {0.,0.,0.};
                    geo::gGeometry().local2Master(planeID,hitPointLocal,hitPointGlobal);
/*
                     streamlog_out(MESSAGE0) << " | " <<  setw(3)  << right 
                                               <<   planeID  << " " << fixed << setw(11) << setprecision (3) << noshowpos << setw(7) << " " <<
                                                                      hitPointGlobal[0] << " : " <<  fitPointGlobal[0] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointGlobal[1] << " : " <<  fitPointGlobal[1] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointGlobal[2] << " : " <<  fitPointGlobal[2] ;
*/
                     streamlog_out(MESSAGE1) << " | " <<  setw(3)  << right 
                                               <<   planeID  << " " << fixed << setw(11) << setprecision (3) << noshowpos << setw(7) << " " <<
                                                                      hitPointLocal[0] << " : " <<  fitPointLocal[0] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointLocal[1] << " : " <<  fitPointLocal[1] 
                                                             << " " <<  noshowpos << setw(7) << " " <<
                                                                      hitPointLocal[2] << " : " <<  fitPointLocal[2] << " glo " << setw(8) << fitPointGlobal[2] ;
                     streamlog_out(MESSAGE1) << std::endl;

                 
                     // Calculate projection matrix
                     // GBL language "Local" -> our language "Telescope"="global" system
                     // GLB language "measurement" -> our language "measurement"="detector" system
                     TMatrixD proL2m(2, 2);
			geo::gGeometry().CalculateProjMatrix(proL2m, hitPointGlobal);
                
                     bool excludeFromFit = false;
                     if ( std::find( _excludeFromFit.begin(), _excludeFromFit.end(), planeID ) != _excludeFromFit.end() ) excludeFromFit = true;
 
                     if ( !excludeFromFit )
                     {
                       	addMeasurementsGBL(point, meas, measPrec, hitPointLocal, fitPointLocal, hitcov, proL2m);
			if( _alignmentMode > 0 ) {

 			        TMatrixD alDer; // alignment derivatives
				alDer.ResizeTo(2, 6);
        			alDer.Zero();

			        std::vector<int> globalLabels;
			        globalLabels.resize(6);

			        TVectorD statevector(5);
                                statevector.Zero(); 

                                statevector[0] = trk->getD0       ();    // 0 - x   
                                statevector[1] = trk->getPhi      ();    // 1 - y
                                statevector[2] = trk->getOmega    ();     // 2 - tx
                                statevector[3] = trk->getZ0       ();    // 3 - ty
                                statevector[4] = trk->getTanLambda();    // 4 - invP
                             
				double trackDirGlobal[] = { statevector[2], statevector[3], 1.};
		            	double trackDirLocal[] = { 0., 0., 0. };
	    			geo::gGeometry().master2LocalVec( planeID, trackDirGlobal, trackDirLocal );


               			addGlobalParametersGBL( point, alDer, globalLabels, planeID, fitPointLocal, trackDirLocal[0], trackDirLocal[1] );
                        }
                     }

                  }                          
                 }

 		 if( _alignmentMode == 0 ) {
                   addSiPlaneScattererGBL(point, trk->getLocation() );
	 	 }
                 pushBackPoint( pointList, point, trk->id() );
 
// get jacobian to arrive to next point (if not at last point already):: 
                 if( i < nstates-1 ){
                   // get next point
                   IMPL::TrackStateImpl* trk = static_cast < IMPL::TrackStateImpl*> ( (*itTrkCand)->getTrackStates().at(i+1) ) ;
                   const float *fitPoint  = trk->getReferencePoint();
                   const int nextVolumeID = trk->getLocation();
                   double hitPointLocal[] = {static_cast<double>(fitPoint[0]), static_cast<double>(fitPoint[1]), static_cast<double>(fitPoint[2]) };
                   double hitPointGlobal[] = {0.,0.,0.};
                   geo::gGeometry().local2Master( nextVolumeID, hitPointLocal, hitPointGlobal);

                   double ds = hitPointGlobal[2] - fitPointGlobal[2]; // actually it should be a method which returns Jacobian matrix (from state1 to state2)
                   jacPointToPoint = propagatePar(ds);
                   streamlog_out(MESSAGE1) << " step[" << i << "] = " << ds << std::endl;
                 } else {
                   streamlog_out(MESSAGE1) << std::endl;
                 }

                }

           }

            // check magnetic field at (0,0,0)
            const gear::BField&   B = geo::gGeometry().getMagneticFiled();
            const double Bmag       = B.at( TVector3(0.,0.,0.) ).r2();
                
            gbl::GblTrajectory* traj = 0;
            if ( Bmag < 1.E-6 ) {
                traj = new gbl::GblTrajectory( pointList, false );
            } else {
                traj = new gbl::GblTrajectory( pointList, true );
            }

            if( _alignmentMode == 0 ) { 
              PerformFitGBLTrajectory( traj, itTrkCand  );
            } else {
              prepareMilleOut( traj ) ;
            }

    }

// convert input TrackCandidates and TrackStates into a GBL Trajectory
    void EUTelGBLFitter::TrackCandidatesToGBLTrajectories( ) {
  	//      Clear(); //  This should be done explictly outside the class. So it is not hidden away.

	////////////////////////////////////////////////////////////////////////////////////////////////////////
 
	//
        vector<const IMPL::TrackImpl*>::const_iterator itTrkCand;
        int trackcounter = 0;
        for ( itTrkCand = _trackCandidatesVec.begin(); itTrkCand != _trackCandidatesVec.end(); ++itTrkCand) {
            streamlog_out(MESSAGE) << " track candidate # " << trackcounter << " of " << _trackCandidatesVec.size() << std::endl; 
            const EVENT::TrackerHitVec& hits = (*itTrkCand)->getTrackerHits();
            for(unsigned int i=0;i< hits.size();i++){
                streamlog_out(MESSAGE) << "hit " << i << " of " << hits.size() << " at " << hits[i] << endl;
            }
            // take care of a track candidate
            TrackCandidatesToGBLTrajectory( itTrkCand  );         
            trackcounter++;
        }
 
    }

    bool EUTelGBLFitter:: PerformMille() {
       streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::PerformMille started " << endl;
       bool _mille_success = false;

       unsigned ntraj = _gblTrackCandidates.size();
       for(unsigned int i=0;i < ntraj; i++){
         streamlog_out ( MESSAGE1 ) << " reading Track Candidate : " << i << " of " << ntraj << " gblTraj:" << _gblTrackCandidates.at(i) << endl;
         prepareMilleOut( _gblTrackCandidates.at(i) );
       }
       streamlog_out ( MESSAGE1 ) << " EUTelGBLFitter::PerformMille -- finished " << endl;
       return _mille_success;
    }


 /* prepare Millepede II compatible binaries 
  * gblTraj  - input GBL trajectory
  * itTrkCan - input Track CAndidates vector (used) to build GBL trajectory
  */
    void EUTelGBLFitter::prepareMilleOut( gbl::GblTrajectory* gblTraj ) {
   // no GBL fit in here! assumes it has been done before and the results stored in TrackImpl object

 
        if ( streamlog_level(MESSAGE0) ){
          std::cout << "MilleOut - gblTrajectory: " << std::endl;
 	  gblTraj->printTrajectory(1);
	  gblTraj->printData();
 	  gblTraj->printPoints(1);
        }

        //gblTraj->milleOut( *_mille );
        

    }


// Jacobian from below looks like //   (-1   0  -y   -dx/dz   x*dx/dz   -y*dx/dz)( x local )
																  //   (0  -1 	x   -dy/dz   x*dy/dz   -y*dy/dz )( y local )
                                  //                                             ( anti clockwise rotation around z) 
                                  //                                             (moving the plane in the -z direction)?????????
																	//                                             (this is clockwise rotation in the y direction)
																	// 	                                           (this is clockwise rotations in the x direction)

                  		 

void EUTelGBLFitter::CreateAlignmentToMeasurementJacobian(std::vector< gbl::GblPoint >* pointList ){

	int number_of_points = pointList->size();
	int pointNum = 0;
	for(pointNum=0; pointNum < number_of_points; ++pointNum){ //Must make sure that the label starts at ????????
		EUTelTrackStateImpl * state = 	_PointToState[ &(pointList->at(pointNum)) ]; //get the state associated with this point
		
		if(state != NULL and state->getHit() != NULL){
		_MilleInterface->CreateAlignmentToMeasurementJacobian(state); //Get the jacobain for that state and sensor
		_MilleInterface->CreateGlobalLabels(state);  //Gets the coorect labels for that sensor
		TMatrixD * Jac = _MilleInterface->getAlignmentJacobian();
		std::vector<int> labels =  _MilleInterface->getGlobalParameters();
		

		(pointList->at(pointNum)).addGlobals(labels, *Jac);
		

		}//END OF IF STATEMENT IF THE THERE IS A STATE WITH HIT
	}


}

void EUTelGBLFitter::SetHitCovMatrixFromFitterGBL(EUTelTrackStateImpl *state){

	double hitcov[4];

	int izPlane = state->getLocation();
	if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 ){
		if( state->getHit() != NULL ){
  		hitcov[0] = _parameterIdXResolutionVec[izPlane];
    	hitcov[2] = _parameterIdYResolutionVec[izPlane];

   		hitcov[0] *= hitcov[0]; // squared !
   		hitcov[2] *= hitcov[2]; // squared !
  	}
	}
	else{
		hitcov[0]=0.1;
		hitcov[2]=0.1;
	}

state->setTrackStateHitCov(hitcov);

}

} // namespace eutelescope


#endif

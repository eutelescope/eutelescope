#include "EUTelTrackStateImpl.h"

namespace eutelescope {

  // the standard requires static const ints to be defined outside the class declaration
  // so we do this here :
  const int EUTelTrackStateImpl::AtOther     ; 
  const int EUTelTrackStateImpl::AtFirstHit  ; 							    
  const int EUTelTrackStateImpl::AtLastHit   ;							    
    
    EUTelTrackStateImpl::EUTelTrackStateImpl() : TrackStateImpl(),
        _location(0),
        _tx(0),
        _ty(0),
        _x(0),
        _y(0),
        _invp(0),
        _covMatrix(),
				_zparameter(0),
				_hit()   
    {

        for(int i=0 ; i < TRKSTATENCOVMATRIX ; i++ ) {
            _covMatrix.push_back( 0.0 ) ; 
        }

        for(int i=0 ; i < TRKSTATENREFSIZE ; i++ ) {
            _reference[i] = 0.0 ;
        }

    }

		//This is the constructor to change TrackStateImpl LCIO to EUTelTrackState
    EUTelTrackStateImpl::EUTelTrackStateImpl( const IMPL::TrackStateImpl& o ) : TrackStateImpl(),
        _location(o.getLocation()),
        _tx(o.getOmega()),
        _ty(o.getZ0()),
        _x(o.getD0()),
        _y(o.getPhi()),
        _invp(o.getTanLambda()),
        _covMatrix(o.getCovMatrix()),
				_zparameter(0),
				_beamQ(-1)	
    {
			//I dont think this is needed to add this. Also need to add the hit at state level!!!!
     //   for(int i=0 ; i < TRKSTATENCOVMATRIX ; i++ ) {
     //       _covMatrix.push_back( 0.0 ) ; 
     //   }
				//has to be called this way to change from const float to float. 
        for(int i=0 ; i < TRKSTATENREFSIZE ; i++ ) {
            _reference[i] = o.getReferencePoint()[i];
        }

    }

    EUTelTrackStateImpl::EUTelTrackStateImpl(int location, float tx, float ty, float x, float y, float invp, const float* covMatrix, const float* reference) :
        _location(0),
        _tx(tx),
        _ty(ty),
        _x(x),
        _y(y),
        _invp(invp),
        _covMatrix(),
				_zparameter(0)
    {

        setLocation( location );

        for(int i=0 ; i < TRKSTATENCOVMATRIX ; i++ ) {
            _covMatrix.push_back( covMatrix[i] ) ; 
        }

        setReferencePoint(reference);
    }


    EUTelTrackStateImpl::EUTelTrackStateImpl(int location, float tx, float ty, float x, float y, float invp, const EVENT::FloatVec& covMatrix, const float* reference) :
        _location(0),
	_tx(tx),
        _ty(ty),
        _x(x),
        _y(y),
        _invp(invp),
        _covMatrix(covMatrix),
				_zparameter(0)
    {

        setLocation( location );

        setReferencePoint(reference);
    }

    EUTelTrackStateImpl::EUTelTrackStateImpl( const EUTelTrackStateImpl &p ) :
//        AccessChecked(),
         TrackStateImpl(),
        _location(0),
        _tx(p.getTx()),
        _ty(p.getTy()),
        _x(p.getX()),
        _y(p.getY()),
        _invp(p.getInvP()),
        _covMatrix(p.getCovMatrix()),
				_zparameter(0)
    {

        setLocation( p.getLocation() );

        setReferencePoint( p.getReferencePoint() );
				setHit(p.getHit());
    }




    EUTelTrackStateImpl::~EUTelTrackStateImpl() { /* no-op */ } 

    int EUTelTrackStateImpl::getLocation() const { return _location ;}
    float EUTelTrackStateImpl::getTx() const { return _tx ;}
    float EUTelTrackStateImpl::getTy() const { return _ty ; }
    float EUTelTrackStateImpl::getX() const { return _x ; }
    float EUTelTrackStateImpl::getY() const { return _y ;}
    float EUTelTrackStateImpl::getInvP() const { return _invp ;}
    float EUTelTrackStateImpl::getZParameter() const { return _zparameter; }

		EVENT::TrackerHit* EUTelTrackStateImpl::getHit() const { return _hit; }

    const EVENT::FloatVec& EUTelTrackStateImpl::getCovMatrix() const { return _covMatrix ; }
    const float* EUTelTrackStateImpl::getReferencePoint() const { return _reference ; }

    void EUTelTrackStateImpl::Print() const {
      	streamlog_out(DEBUG0) << "location " << getLocation() << " Tx:"<<getTx() << " Ty:"<<getTy() << " X:"<<getX() << " Y:"<<getY() << " InvP:"<<getInvP() << std::endl;
				streamlog_out(DEBUG0) << "Reference hit: "  << _reference[0] << "," << _reference[1] << "," <<_reference[2]  << std::endl;
streamlog_out(DEBUG0) << "Zparameter " << _zparameter <<std::endl;
streamlog_out(DEBUG0) << "Hit pointer " << _hit << std::endl;
streamlog_out(DEBUG0) << "The covariance matrix of the state" << std::endl;
	        streamlog_message( DEBUG0, getTrackStateCov().Print();, std::endl; );

}

    void  EUTelTrackStateImpl::setLocation( int location ){
//        checkAccess("EUTelTrackStateImpl::setLocation") ;

        // if( location < 0 || location >= EUTelTrackState::Location::size() ){
        //     throw( Exception( " EUTelTrackStateImpl::setLocation called with an undefined Location" )) ;
        // }

        _location = location  ;
    } 
    void  EUTelTrackStateImpl::setTx( float tx ){
//        checkAccess("EUTelTrackStateImpl::setTx") ;
        _tx = tx  ;
    } 
    void  EUTelTrackStateImpl::setTy( float ty ){ 
//      checkAccess("EUTelTrackStateImpl::setTy") ;
        _ty = ty ; 
    } 
    void  EUTelTrackStateImpl::setX( float x ) { 
//      checkAccess("EUTelTrackStateImpl::setX") ;
        _x = x ;
    } 
    void  EUTelTrackStateImpl::setY( float y ){
//      checkAccess("EUTelTrackStateImpl::setY") ;
        _y = y ; 
    } 
    void  EUTelTrackStateImpl::setInvP( float invp ){
//        checkAccess("EUTelTrackStateImpl::setInvP") ;
        _invp = invp ; 
    } 
    void  EUTelTrackStateImpl::setCovMatrix( const float* cov ){ 
//        checkAccess("EUTelTrackStateImpl::setCovMatrix") ;
        for( int i=0 ; i<TRKSTATENCOVMATRIX ; i++ ){
            _covMatrix[i] = cov[i]  ; 
        }
    } 
    void  EUTelTrackStateImpl::setCovMatrix( const EVENT::FloatVec& cov ){ 
//        checkAccess("EUTelTrackStateImpl::setCovMatrix") ;
        for( int i=0 ; i<TRKSTATENCOVMATRIX ; i++ ){
            _covMatrix[i] = cov[i]  ; 
        }
    }


    void  EUTelTrackStateImpl::setReferencePoint( const float* rPnt ){ 
//        checkAccess("EUTelTrackStateImpl::setReferencePoint") ;
        for( int i=0 ; i<TRKSTATENREFSIZE ; i++ ){
            _reference[i] = rPnt[i]  ; 
        }
    }
	
		void EUTelTrackStateImpl::setZParameter(float z){
        _zparameter = z ;
    }

		void EUTelTrackStateImpl::setbeamQ(int beamQ){ //This really should be a static I think but then exactly how to use this I dont know. MUST FIX
			_beamQ=beamQ;
		}

    void EUTelTrackStateImpl::setHit( EVENT::TrackerHit* hit) {
        _hit = hit;
    }

		void EUTelTrackStateImpl::setTrackStateHitCov(double cov[4]){

			if(_hit != NULL){
				_covHitMatrix[0] = cov[0];
				_covHitMatrix[1] = cov[1];
				_covHitMatrix[2] = cov[2];
				_covHitMatrix[3] = cov[3];

			}
}

void EUTelTrackStateImpl::getTrackStateHitCov( double (&cov)[4]) {
				cov[0] = _covHitMatrix[0];
				cov[1] = _covHitMatrix[1];
				cov[2] = _covHitMatrix[2];
				cov[3] = _covHitMatrix[3];
}


//This function will output the momentum of the track in cartesian coordinates into a TVector structure
    /** Calculate track momentum from track parameters 
     * @param ts track state with specified tx,ty,x,y,invP
     * @return 3-vector of momentum
     */
TVector3 EUTelTrackStateImpl::getPfromCartesianParameters() const {
	streamlog_out(DEBUG2) << "EUTelTrackStateImpl::getPfromCartesianParameters()--------------------------BEGIN" << std::endl;
		//Ignore else where must fix:

	streamlog_out(DEBUG2) << "Input parameters: tx,ty, beamq,invp "<<_tx <<","<<_ty<<","<<_beamQ<<","<<_invp<<std::endl;
	const double p  =  1. / (_invp * -1);      //1. / (_invp * _beamQ);
  const double px = p*_tx / sqrt( 1. + _tx*_tx + _ty*_ty );
  const double py = p*_ty / sqrt( 1. + _tx*_tx + _ty*_ty );
  const double pz = p    / sqrt( 1. + _tx*_tx + _ty*_ty );

	streamlog_out(DEBUG2) << "Output parameters: px,py, pz "<<px <<","<<py<<","<<pz<<","<<std::endl;
        
  streamlog_out(DEBUG2) << "-------------------------------EUTelTrackStateImpl::getPfromCartesianParameters()-------------------------END" << std::endl;
        
  return TVector3(px,py,pz);
}
	//This function uses the class member findIntersectionWithCertainID. This find the point on the sensor specified by nextPlaneID and returns the global position in output.
	int EUTelTrackStateImpl::findIntersectionWithCertainID(int nextPlaneID, float* output ) const {
	streamlog_out(DEBUG5) << "EUTelTrackStateImpl::findIntersectionWithCertainID----------------------------BEGIN" << std::endl;
	TVector3 pVec = getPfromCartesianParameters();
	streamlog_out(DEBUG5) << "Momentum: " << pVec[0]<<","<<pVec[1]<<","<<pVec[2]<<","<< std::endl;
//	int sensorID = geo::gGeometry().findIntersectionWithCertainID( _x, _y, _zparameter, pVec[0],pVec[1],pVec[2], _beamQ, nextPlaneID, output ); //Input global positions and momentum in cartesian
//	streamlog_out(DEBUG5) << "SensorID" << sensorID << std::endl;
	streamlog_out(DEBUG5) << "EUTelTrackStateImpl::findIntersectionWithCertainID----------------------------END" << std::endl;
	return nextPlaneID;
}

TVector3 EUTelTrackStateImpl::getXYZfromArcLength( float s ) const {
//		streamlog_out(DEBUG2) << "EUTelTrackStateImpl::getXYZfromArcLength----------------------------BEGIN" << std::endl;

//	TVector3 pVec = getPfromCartesianParameters();
//	TVector3 pos = geo::gGeometry().getXYZfromArcLength( _x, _y, _zparameter, pVec[0], pVec[1], pVec[2], _beamQ, s);

//		streamlog_out(DEBUG2) << "EUTelTrackStateImpl::getXYZfromArcLength----------------------------END" << std::endl;	
//	return pos;

}	

//This function returns the H matrix of the state. This relates the state in global coordinates to local coordinates (The hit measurement). Note it assumes that z axis is a parameter so you only need the rotations of the plane you must be on. Note the H matrix from Geometry is Local->Global so we must take the inverse.
    /** Retrieve track state projection onto measurement space matrix
     * 
     * @param ts track state
     * @return 
     */
TMatrixD EUTelTrackStateImpl::getH() const {

	streamlog_out( DEBUG2 ) << "EUTelTrackStateImpl::getH()---------------------------------------BEGIN" << std::endl;
	
  	TMatrixD H(2,5);//2x5 matrix. 
  	H.Zero();
  	double trkPoint[] = { _x, _y, _zparameter };
  	const TGeoHMatrix* globalH = geo::gGeometry().getHMatrix( trkPoint );
        
  	if ( streamlog_level(DEBUG0) ) {
  		streamlog_out( DEBUG0 ) << "Local to global transformation matrix:" << std::endl;
		globalH->Print();
	}
        
	const TGeoHMatrix& globalHInv = globalH->Inverse();
	if ( streamlog_level(DEBUG0) ) {
		streamlog_out( DEBUG0 ) << "Global to local transformation matrix:" << std::endl;
	 	globalHInv.Print();
	}
        
	const double* rotation = globalHInv.GetRotationMatrix();

	// Fill necessary components
  	H[0][0] = rotation[0]; // x projection, xx
  	H[0][1] = rotation[1]; // y projection, xy
  	H[1][0] = rotation[3]; // x projection, yx
  	H[1][1] = rotation[4]; // y projection, yy

	if ( streamlog_level(DEBUG0) ) {
		streamlog_out( DEBUG0 ) << "Matrix H:" << std::endl;
		H.Print();
	}

 	streamlog_out( DEBUG2 ) << "EUTelTrackStateImpl::getH()---------------------------------------END" << std::endl;       
	return H;
}

    /** Convert track state to the vector object. Useful for matrix operations
     * 
     * @param ts track stare
     * @return vector of parameters
     */
TVectorD EUTelTrackStateImpl::getTrackStateVec() const {
	streamlog_out( DEBUG1 ) << "EUTelTrackStateImpl::getTrackStateVec()------------------------BEGIN" << std::endl;
        TVectorD x(5);
	  x[0] = _x;
	  x[1] = _y;
	  x[2] = _tx;
	  x[3] = _ty;
	  x[4] = _invp;
        
	if ( streamlog_level(DEBUG0) ){
		streamlog_out( DEBUG0 ) << "Track state:" << std::endl;
		Print();
	}
  	streamlog_out( DEBUG1 ) << "EUTelTrackStateImpl::getTrackStateVec()------------------------END" << std::endl;
 	return x;
}
    
/** Convert track state parameter covariances to the matrix object. Useful for matrix operations
* 
* @param ts track state
* @return covariance matrix
*/
TMatrixDSym EUTelTrackStateImpl::getTrackStateCov() const {

	streamlog_out( DEBUG1 ) << "EUTelTrackStateImpl::getTrackStateCov()----------------------------BEGIN" << std::endl;
	TMatrixDSym C(5);   
	const EVENT::FloatVec& trkCov = getCovMatrix();        
	C.Zero();
            
	C[0][0] = trkCov[0]; 
	C[1][0] = trkCov[1];  C[1][1] = trkCov[2]; 
	C[2][0] = trkCov[3];  C[2][1] = trkCov[4];  C[2][2] = trkCov[5]; 
	C[3][0] = trkCov[6];  C[3][1] = trkCov[7];  C[3][2] = trkCov[8];  C[3][3] = trkCov[9]; 
	C[4][0] = trkCov[10]; C[4][1] = trkCov[11]; C[4][2] = trkCov[12]; C[4][3] = trkCov[13]; C[4][4] = trkCov[14]; 
        
	if ( streamlog_level(DEBUG0) ){
  		streamlog_out( DEBUG0 ) << "Track state covariance matrix:" << std::endl;
	    	C.Print();
	}
        
	return C;
	streamlog_out( DEBUG1 ) << "EUTelTrackStateImpl::getTrackStateCov()----------------------------END" << std::endl;
}


TMatrix EUTelTrackStateImpl::getPropagationJacobianF( float dz ){

	TVector3 pVec = getPfromCartesianParameters();
	TMatrix jacobian(5,5);
	jacobian.Zero();
	streamlog_out( DEBUG1 ) << "These are the parameters being used in the calculation of the Jacobian. "<< "X= "<< _x << " Y= "<< _y << " Zparameter= "<< _zparameter << " Px= " << pVec[0] << " Py= " << pVec[1] << " Pz= " << pVec[2] << " Qbeam= " <<_beamQ  << std::endl;
	jacobian = geo::gGeometry().getPropagationJacobianF(  _x, _y, _zparameter, pVec[0],pVec[1],pVec[2], _beamQ, dz);

	return jacobian;
}

TVector3 EUTelTrackStateImpl::getIncidenceVectorInLocalFrame(){
	TVector3 pVec =	getPfromCartesianParameters();
	TVector3 pVecUnit = pVec;   //.Unit();
    streamlog_out(DEBUG2) << "Momentum in global coordinates  Px,Py,Pz= " << pVec[0]<<","<<pVec[1]<<","<<pVec[2] << std::endl;
	double globalVec[] = { pVecUnit[0],pVecUnit[1],pVecUnit[2] };
	double localVec[3];

	geo::gGeometry().master2LocalVec( _location ,globalVec, localVec );

	TVector3 pVecUnitLocal;
	pVecUnitLocal[0] = localVec[0]; 	pVecUnitLocal[1] = localVec[1]; 	pVecUnitLocal[2] = localVec[2]; 
  streamlog_out(DEBUG2) << "Momentum in local coordinates  Px,Py,Pz= " << pVecUnitLocal[0]<<","<<pVecUnitLocal[1]<<","<<pVecUnitLocal[2]<< std::endl;

	return pVecUnitLocal;
}
 			 

} // namespace eutelescope



#include "EUTelTrackStateImpl.h"

namespace eutelescope {

  // the standard requires static const ints to be defined outside the class declaration
  // so we do this here :
  const int EUTelTrackStateImpl::AtOther  ; 
  const int EUTelTrackStateImpl::AtFirstHit  ; 							    
  const int EUTelTrackStateImpl::AtLastHit  ;							    
    
    EUTelTrackStateImpl::EUTelTrackStateImpl() :
        _location(0),
        _tx(0),
        _ty(0),
        _x(0),
        _y(0),
        _invp(0),
        _covMatrix()
    {

        for(int i=0 ; i < TRKSTATENCOVMATRIX ; i++ ) {
            _covMatrix.push_back( 0.0 ) ; 
        }

        for(int i=0 ; i < TRKSTATENREFSIZE ; i++ ) {
            _reference[i] = 0.0 ;
        }

    }

    EUTelTrackStateImpl::EUTelTrackStateImpl(int location, float tx, float ty, float x, float y, float invp, const float* covMatrix, const float* reference) :
        _location(0),
        _tx(tx),
        _ty(ty),
        _x(x),
        _y(y),
        _invp(invp),
        _covMatrix()
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
        _covMatrix(covMatrix)
    {

        setLocation( location );

        setReferencePoint(reference);
    }

    EUTelTrackStateImpl::EUTelTrackStateImpl( const EUTelTrackStateImpl &p ) :
        AccessChecked(),
        _location(0),
        _tx(p.getTx()),
        _ty(p.getTy()),
        _x(p.getX()),
        _y(p.getY()),
        _invp(p.getInvP()),
        _covMatrix(p.getCovMatrix())
    {

        setLocation( p.getLocation() );

        setReferencePoint( p.getReferencePoint() );
    }




    EUTelTrackStateImpl::~EUTelTrackStateImpl() { /* no-op */ } 

    int EUTelTrackStateImpl::getLocation() const { return _location ;}
    float EUTelTrackStateImpl::getTx() const { return _tx ;}
    float EUTelTrackStateImpl::getTy() const { return _ty ; }
    float EUTelTrackStateImpl::getX() const { return _x ; }
    float EUTelTrackStateImpl::getY() const { return _y ;}
    float EUTelTrackStateImpl::getInvP() const { return _invp ;}

    const EVENT::FloatVec& EUTelTrackStateImpl::getCovMatrix() const { return _covMatrix ; }
    const float* EUTelTrackStateImpl::getReferencePoint() const { return _reference ; }


    void  EUTelTrackStateImpl::setLocation( int location ){
        checkAccess("EUTelTrackStateImpl::setLocation") ;

        // if( location < 0 || location >= EUTelTrackState::Location::size() ){
        //     throw( Exception( " EUTelTrackStateImpl::setLocation called with an undefined Location" )) ;
        // }

        _location = location  ;
    } 
    void  EUTelTrackStateImpl::setTx( float tx ){
        checkAccess("EUTelTrackStateImpl::setTx") ;
        _tx = tx  ;
    } 
    void  EUTelTrackStateImpl::setTy( float ty ){ 
        checkAccess("EUTelTrackStateImpl::setTy") ;
        _ty = ty ; 
    } 
    void  EUTelTrackStateImpl::setX( float x ) { 
        checkAccess("EUTelTrackStateImpl::setX") ;
        _x = x ;
    } 
    void  EUTelTrackStateImpl::setY( float y ){
        checkAccess("EUTelTrackStateImpl::setY") ;
        _y = y ; 
    } 
    void  EUTelTrackStateImpl::setInvP( float invp ){
        checkAccess("EUTelTrackStateImpl::setInvP") ;
        _invp = invp ; 
    } 
    void  EUTelTrackStateImpl::setCovMatrix( const float* cov ){ 
        checkAccess("EUTelTrackStateImpl::setCovMatrix") ;
        for( int i=0 ; i<TRKSTATENCOVMATRIX ; i++ ){
            _covMatrix[i] = cov[i]  ; 
        }
    } 
    void  EUTelTrackStateImpl::setCovMatrix( const EVENT::FloatVec& cov ){ 
        checkAccess("EUTelTrackStateImpl::setCovMatrix") ;
        for( int i=0 ; i<TRKSTATENCOVMATRIX ; i++ ){
            _covMatrix[i] = cov[i]  ; 
        }
    } 

    void  EUTelTrackStateImpl::setReferencePoint( const float* rPnt ){ 
        checkAccess("EUTelTrackStateImpl::setReferencePoint") ;
        for( int i=0 ; i<TRKSTATENREFSIZE ; i++ ){
            _reference[i] = rPnt[i]  ; 
        }
    } 

} // namespace eutelescope



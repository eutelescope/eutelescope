#include "EUTelTrackImpl.h"
#include "EUTelTrackStateImpl.h"

#include "lcio.h"

#include <cmath>
#include <sstream>

namespace eutelescope {
  
    EUTelTrackImpl::EUTelTrackImpl() :
        _chi2(0),
        _ndf(0),
        _dEdx(0),
        _dEdxError(0),
        _tracks(),
        _hits(),
        _trackStates() { 
        }

  // copy constructor
  EUTelTrackImpl::EUTelTrackImpl(const EUTelTrackImpl& o) : AccessChecked(),
  _chi2(o.getChi2()),
  _ndf(o.getNdf()),
  _dEdx(o.getdEdx()),
  _dEdxError(o.getdEdxError()),
  _tracks(o.getTracks()),
  _hits(o.getTrackerHits()),
  _trackStates(o.getTrackStates())
  { 

      *this = o ; // call operator =
    
  }

    const EUTelTrackImpl & EUTelTrackImpl::operator = ( const EUTelTrackImpl &o )
    {
    _chi2 = o._chi2 ;
    _ndf = o._ndf ;
    _dEdx = o._dEdx ;
    _dEdxError = o._dEdxError ;

    std::copy( o._hits.begin() ,  o._hits.end() , std::back_inserter( _hits ) ) ;

    std::copy( o._tracks.begin() ,  o._tracks.end() , std::back_inserter( _tracks ) ) ;

    for( unsigned int i=0; i< o._trackStates.size() ; i++ ){
      _trackStates.push_back( new EUTelTrackStateImpl(  *o._trackStates[i] ) ) ;
    }

    return *this ;
    }


    EUTelTrackImpl::~EUTelTrackImpl() { 

        // delete track states in here ?
        for( unsigned int i=0; i<_trackStates.size() ; i++ ){

	  delete _trackStates[i] ;
        }

    } 

    float EUTelTrackImpl::getTx() const {               return ( _trackStates.size()>0 ? _trackStates[0]->getTx()             : 0 ) ;  }
    float EUTelTrackImpl::getTy() const {               return ( _trackStates.size()>0 ? _trackStates[0]->getTy()            : 0 ) ;  }
    float EUTelTrackImpl::getX() const {                return ( _trackStates.size()>0 ? _trackStates[0]->getX()          : 0 ) ;  }
    float EUTelTrackImpl::getY() const {                return ( _trackStates.size()>0 ? _trackStates[0]->getY()             : 0 ) ;  }
    float EUTelTrackImpl::getInvP() const {             return ( _trackStates.size()>0 ? _trackStates[0]->getInvP()      : 0 ) ;  }
    const EVENT::FloatVec& EUTelTrackImpl::getCovMatrix() const {   
      static const EUTelTrackStateImpl dummy ;
      return ( _trackStates.size()>0 ? _trackStates[0]->getCovMatrix()      : dummy.getCovMatrix() ) ;  
    }
    const float* EUTelTrackImpl::getReferencePoint() const { 
      static const EUTelTrackStateImpl dummy ;
      return ( _trackStates.size()>0 ? _trackStates[0]->getReferencePoint() : dummy.getReferencePoint() ) ;  
    }


    float EUTelTrackImpl::getChi2() const { return _chi2 ;}
    int   EUTelTrackImpl::getNdf() const { return _ndf ;}
    float EUTelTrackImpl::getdEdx() const { return _dEdx ; }
    float EUTelTrackImpl::getdEdxError() const { return _dEdxError ; }

    const EVENT::TrackerHitVec & EUTelTrackImpl::getTrackerHits() const {
        return _hits ;
    }

    const EUTelTrackVec & EUTelTrackImpl::getTracks() const {
        return _tracks ;
    } 

    const EUTelTrackStateVec & EUTelTrackImpl::getTrackStates() const {
        return _trackStates ;
    } 

    const EUTelTrackStateImpl* EUTelTrackImpl::getClosestTrackState( float x, float y, float z ) const {
        EUTelTrackStateImpl* closest = _trackStates[0] ;
        const float * refP = _trackStates[0]->getReferencePoint() ;
        float shortest_distance_square = pow( ( x - refP[0] ) , 2 ) + pow( ( y - refP[1] ) , 2 ) + pow( ( z - refP[2] ) , 2 ) ;
        float current_distance_square = 0 ;

        for( unsigned int i=1 ; i < _trackStates.size() ; i++ ){
            refP = _trackStates[i]->getReferencePoint() ;
            current_distance_square = pow( ( x - refP[0] ) , 2 ) + pow( ( y - refP[1] ) , 2 ) + pow( ( z - refP[2] ) , 2 ) ;
            if( current_distance_square < shortest_distance_square ){
                closest = _trackStates[i] ;
                shortest_distance_square = current_distance_square ;
            }
        }
        return closest ;
    } 

    const EUTelTrackStateImpl* EUTelTrackImpl::getTrackState( int location ) const {

        for( unsigned int i=0 ; i < _trackStates.size() ; i++ ){
            if( _trackStates[i]->getLocation() == location ){
                return _trackStates[i] ;  
            }
        }
        return NULL ;
    } 


    void  EUTelTrackImpl::setTx( float tx ){

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setTx within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setTx( tx ) ;
    } 
    void  EUTelTrackImpl::setTy( float ty ){ 

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setTy within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setTy( ty ) ;
    } 
    void  EUTelTrackImpl::setX( float x ) { 

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setX within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setX( x ) ;
    } 
    void  EUTelTrackImpl::setY( float y ){

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setY within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setY( y ) ;
    } 
    void  EUTelTrackImpl::setInvP( float invp ){

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setInvP within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setInvP( invp ) ;
    } 

    void  EUTelTrackImpl::setCovMatrix( const float* cov ){ 

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setCovMatrix within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setCovMatrix( cov ) ;
    } 
    void  EUTelTrackImpl::setCovMatrix( const EVENT::FloatVec& cov ){ 

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setCovMatrix within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setCovMatrix( cov ) ;
    } 

    void  EUTelTrackImpl::setReferencePoint( const float* rPnt){ 

        if( _trackStates.size() == 0 ){
            // create a first TrackState for backwards compatibility
            EUTelTrackStateImpl* ts = new EUTelTrackStateImpl() ;
            _trackStates.push_back( ts ) ;
        }

        if( _trackStates.size() != 1 ){
            throw( lcio::Exception( " trying to use setReferencePoint within Track object containing more than one TrackState." )) ;
        }

        (static_cast< EUTelTrackStateImpl* >(_trackStates[0]) )->setReferencePoint( rPnt ) ;
    } 

    void  EUTelTrackImpl::setChi2( float chi2 ){ 
        checkAccess("EUTelTrackImpl::setChi2") ;
        _chi2 = chi2 ; 
    } 
    void  EUTelTrackImpl::setNdf( int ndf ){ 
        checkAccess("EUTelTrackImpl::setNdf") ;
        _ndf = ndf  ; 
    } 
    void  EUTelTrackImpl::setdEdx( float dEdx ){ 
        checkAccess("EUTelTrackImpl::setdEdx") ;
        _dEdx = dEdx ; 
    } 
    void  EUTelTrackImpl::setdEdxError( float dEdxError ){
        checkAccess("EUTelTrackImpl::setdEdxError") ;
        _dEdxError = dEdxError  ;
    }   

    void EUTelTrackImpl::addHit( EVENT::TrackerHit* hit) {
        _hits.push_back( hit ) ;
    }

    void  EUTelTrackImpl::addTrack( EUTelTrackImpl* trk ) {
        checkAccess("EUTelTrackImpl::addTrack") ;
        _tracks.push_back( trk ) ;
    }

    void  EUTelTrackImpl::addTrackState( EUTelTrackStateImpl* trkstate ) {
        checkAccess("EUTelTrackImpl::addTrackState") ;
        if( trkstate->getLocation() != EUTelTrackStateImpl::AtOther &&
            getTrackState( trkstate->getLocation() ) != NULL )
        {
            std::stringstream ss;
            ss << "another TrackState already exists with Location set to: " << trkstate->getLocation() ;
            throw( lcio::Exception( ss.str() )) ;
        }
        _trackStates.push_back( trkstate ) ;
    }

  EUTelTrackStateVec & EUTelTrackImpl::trackStates()  {
    checkAccess("EUTelTrackImpl::trackStates") ;
    return _trackStates ;
  } 

} // namespace eutelescope



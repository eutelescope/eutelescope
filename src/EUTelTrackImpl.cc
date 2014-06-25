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
  EUTelTrackImpl::EUTelTrackImpl(const EUTelTrackImpl& o) : EVENT::LCObject(), AccessChecked(),
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

	//This constructor will change IMPL::TrackImpl into EUTeltrack
	EUTelTrackImpl::EUTelTrackImpl(const IMPL::TrackImpl& o) : EVENT::LCObject(), AccessChecked(),
  _chi2(o.getChi2()),
  _ndf(o.getNdf()),
  _dEdx(o.getdEdx()),
  _dEdxError(o.getdEdxError())
  { 
		EVENT::TrackStateVec states =  o.getTrackStates();
		EVENT::TrackStateVec::const_iterator state;
		for(state = states.begin(); state != states.end(); ++state){
			EUTelTrackStateImpl* EUState = new EUTelTrackStateImpl(*(*state));
			EVENT::TrackerHitVec hits = 	o.getTrackerHits();
			EVENT::TrackerHitVec::const_iterator hit;
			for(hit = hits.begin();hit != hits.end(); ++hit){
				if((*state)->getLocation() == Utility::getSensorIDfromHit(*hit)){
					EUState->setHit(*hit);
					break;
				}else{
					EUState->setHit(NULL);
				}
			}
		addTrackState(EUState);
		}			
	    
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


    float EUTelTrackImpl::getTx(unsigned int i) const {               return ( _trackStates.size()>i ? _trackStates[i]->getTx()        : 0 ) ;  }
    float EUTelTrackImpl::getTy(unsigned int i) const {               return ( _trackStates.size()>i ? _trackStates[i]->getTy()        : 0 ) ;  }
    float EUTelTrackImpl::getX(unsigned int i) const {                return ( _trackStates.size()>i ? _trackStates[i]->getX()         : 0 ) ;  }
    float EUTelTrackImpl::getY(unsigned int i) const {                return ( _trackStates.size()>i ? _trackStates[i]->getY()         : 0 ) ;  }
    float EUTelTrackImpl::getInvP(unsigned int i) const {             return ( _trackStates.size()>i ? _trackStates[i]->getInvP()      : 0 ) ;  }


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

    const EUTelTrackStateImpl* EUTelTrackImpl::getFirstTrackState(  ) const {
        if( _trackStates.size() > 0 ) {
           int firstElement = 0; 
           return _trackStates[ firstElement ] ;  
        }
        return NULL ;
    } 

    const EUTelTrackStateImpl* EUTelTrackImpl::getTrackState( int location ) const {
        for( unsigned int i=0 ; i < _trackStates.size() ; i++ ){
            if( _trackStates[i]->getLocation() == location ){
                return _trackStates[i] ;  
            }
        }
        return NULL ;
    } 

    void  EUTelTrackImpl::setTypeBit( int  index, bool val){
        checkAccess("TrackImpl::setTypeBit") ;
        _type.set( index, val  )  ;
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

    void  EUTelTrackImpl::addTrack( EUTelTrackImpl* trk ) {
        checkAccess("EUTelTrackImpl::addTrack") ;
        _tracks.push_back( trk ) ;
    }

    void  EUTelTrackImpl::addTrackState( EUTelTrackStateImpl* trkstate ) {
        checkAccess("EUTelTrackImpl::addTrackState") ;
        Print();
//        std::cout << " addTrackState  " << trkstate << " of " << _trackStates.size() << std::endl;
/*        if( trkstate->getLocation() != EUTelTrackStateImpl::AtOther &&
            getTrackState( trkstate->getLocation() ) != NULL )
        {
            std::stringstream ss;
            ss << "another TrackState already exists with Location set to: " << trkstate->getLocation() ;
            throw( lcio::Exception( ss.str() )) ;
        }*/
        _trackStates.push_back( trkstate ) ;
    }

  EUTelTrackStateVec & EUTelTrackImpl::trackStates()  {
    checkAccess("EUTelTrackImpl::trackStates") ;
    return _trackStates ;
  }
	//We dont want to attach hits to a track in most cases. Since then you lose the information of what state it was part of. So you add hits to states and states to tracks. However you can access hit directy with this function.
 const EVENT::TrackerHitVec EUTelTrackImpl::getHitsOnTrack() const {
	EUTelTrackStateVec::const_iterator state;		
	EVENT::TrackerHitVec hits;
	for( state =  _trackStates.begin(); state != _trackStates.end(); ++state){
		if((*state)->getHit() == NULL){
			continue;
		}else{
		hits.push_back((*state)->getHit());
		}
	}
	return hits;
}

IMPL::TrackImpl* EUTelTrackImpl::CreateLCIOTrack(){
streamlog_out( DEBUG4 ) << " ---------------- EUTelTrackImpl::CreateLCIOTrack-- BEGIN ------------- " << std::endl;

	IMPL::TrackImpl* LCIOtrack = new IMPL::TrackImpl;
	//Loop over all state on the track and fill new track object
	EUTelTrackStateVec tracks = getTrackStates();
	streamlog_out( DEBUG0 ) << "The size of the state " << tracks.size() <<std::endl;  
	EUTelTrackStateVec::const_iterator trackstate;
	for(trackstate=tracks.begin(); trackstate !=tracks.end(); trackstate++){
		TVectorD statevector = (*trackstate)->getTrackStateVec();


		IMPL::TrackStateImpl* implstate     = static_cast <IMPL::TrackStateImpl*> (*trackstate ); //This is possible since EUTelTrack is derived from IMPL::TrackState

		////////Add our state variables into container. The covariant matrix is for our coordinate system
		implstate->setD0(statevector[0]); //x position global
		implstate->setPhi(statevector[1]); //y position global
		implstate->setOmega(statevector[2]); //tx position global
		implstate->setZ0(statevector[3]); //ty position global
		implstate->setTanLambda(statevector[4]); //invp position global		 

                streamlog_out(MESSAGE3) <<  "  " << (*trackstate ) -> id() << "  " << (*trackstate ) -> getLocation() << " " << (*trackstate )->getX() << " " << (*trackstate )->getY() << std::endl;
		LCIOtrack->addTrackState( implstate );
	}

	

   	// Assign hits to LCIO TRACK
	const EVENT::TrackerHitVec& trkcandhits = getHitsOnTrack();
	streamlog_out( DEBUG0 ) << "The size of hits " << trkcandhits.size() <<std::endl; 
    	EVENT::TrackerHitVec::const_iterator itrHit;
    	for ( itrHit = trkcandhits.begin(); itrHit != trkcandhits.end(); ++itrHit ){
    		LCIOtrack->addHit( *itrHit );
    	}
	
			return LCIOtrack;
streamlog_out( DEBUG4 ) << " ---------------- EUTelTrackImpl::CreateLCIOTrack-- END ------------- " << std::endl;
} 

	//int EUTelTrackImpl::getType(){
//		return  _type;
//	}


} // namespace eutelescope



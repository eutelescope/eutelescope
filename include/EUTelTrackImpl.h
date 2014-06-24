#ifndef IMPL_EUTELTRACKIMPL_H
#define IMPL_EUTELTRACKIMPL_H 1

#include "streamlog/streamlog.h"
#include <iostream>
#include <iomanip>      // std::setw

#include "IMPL/AccessChecked.h"
#include "EVENT/TrackerHit.h"

#include "LCIOSTLTypes.h"

#include <bitset>

#include "EUTelTrackStateImpl.h"
#include "EUTelUtility.h"

namespace eutelescope {

  class EUTelTrackImpl;
  /**Vector of (pointers to) Tracks.*/
  typedef std::vector<EUTelTrackImpl*> EUTelTrackVec;

//  class EUTelTrackStateImpl ;
  /**Vector of (pointers to) TrackStates.*/
  typedef std::vector<EUTelTrackStateImpl*> EUTelTrackStateVec ;
  
  using namespace std;

  class EUTelTrackImpl : public EVENT::LCObject, public IMPL::AccessChecked {
    
  public: 
    
    /** Default constructor, initializes values to 0.
     */
    EUTelTrackImpl() ;

		EUTelTrackImpl(const IMPL::TrackImpl& o);
    
    /** Copy constructor - creates shallow copy, i.e. all data members are copied but pointers to other LCObjects
     *  i.e. TrackerHits and Tracks are preserved.
     */
    EUTelTrackImpl( const EUTelTrackImpl& ) ; 
    const EUTelTrackImpl& operator = ( const EUTelTrackImpl& o ) ; // operator =
    

    /// Destructor.
    virtual ~EUTelTrackImpl() ; 

IMPL::TrackImpl* CreateLCIOTrack();

	
    
    virtual int id() const { return simpleUID() ; }
 
    virtual void Print(){
//      PrintTrackStates();   
      for(unsigned int i=0; i < this->getTrackStates().size() ; i++) {    
        streamlog_out(MESSAGE0) << "track " << id() << " state:" << i  << " id: " << setw(5) << (getTrackStates().at(i))->id() << " location: " ;
         streamlog_out(MESSAGE0) << setw(3) << (getTrackStates().at(i))->getLocation() << " at "<< getTrackStates().at(i) ;
         streamlog_out(MESSAGE0) << " at Tx: " << getTx(i) << " at Ty: " << getTy(i) << " getX() " << getX(i) << " getY() " << getY(i) << " "  ;
        const float*   point = (getTrackStates().at(i))->getReferencePoint();
        streamlog_out(MESSAGE0) << std::setw(7) << point[0] << " " << std::setw(7) << point[1] << " " << std::setw(7) << point[2] ;
        streamlog_out(MESSAGE0) << std::endl;
      }
        //IMPL::TrackImpl* track = static_cast< IMPL::TrackImpl*> (*itTrk);
        for ( size_t i = 0; i< _hits.size(); i++ ) 
        { 
            EVENT::TrackerHit* ihit = _hits[i];
            int ic = ihit->id();
            streamlog_out(DEBUG5) << ic << " Is the hit IDs " << std::endl;
        }
        streamlog_out(DEBUG5) << std::endl;

    }

    virtual float getTx(unsigned int i=0) const ;

    virtual float getTy(unsigned int i=0) const ;

    virtual float getX(unsigned int i=0) const ;

    virtual float getY(unsigned int i=0) const ;

    virtual float getInvP(unsigned int i=0) const ;

    /** Covariance matrix of the track parameters. Stored as lower triangle matrix where
     *  the order of parameters is:   d0, phi, omega, z0, tan(lambda).
     *  So we have cov(d0,d0), cov( phi, d0 ), cov( phi, phi), ...
     *  Information is stored in the first TrackState of this Track, @see TrackState.
     */
    virtual const EVENT::FloatVec & getCovMatrix() const ;

    /** Reference point of the track parameters.
     *  The default for the reference point is the point of closest approach.
     *  Information is stored in the first TrackState of this Track, @see TrackState.
     */
    virtual const float* getReferencePoint() const ;

    /** Chi^2 of the track fit.
     */
    virtual float getChi2() const ;

    /** Number of degrees of freedom  of the track fit.
     */
    virtual int getNdf() const ;

    /** dEdx of the track.
     */
    virtual float getdEdx() const;

    /** Error of dEdx.
     */
    virtual float getdEdxError() const;

    /** The tracks (as Track objects) that have been combined to this track.
     */
    virtual const EUTelTrackVec & getTracks() const ;

    /** Returns very first track state associated with this track. @see TrackState.
     */
    virtual const EUTelTrackStateImpl* getFirstTrackState(  ) const ;

    /** Returns track states associated with this track. @see TrackState.
     */
    virtual const EUTelTrackStateVec & getTrackStates() const ;

		virtual const EVENT::TrackerHitVec  getHitsOnTrack() const;
  
    virtual void PrintTrackStates(){
      streamlog_out(MESSAGE0) << " printing track states " << _trackStates.size() << std::endl;
//      for(unsigned int i=0; i<states.size(); i++){
//        EUTelTrackStateImpl* state = static_cast<EUTelTrackStateImpl*> (states.at(i));
//        streamlog_out(MESSAGE0) << state->getLocation() << " ";
//      } 
       for( unsigned int i=0 ; i < _trackStates.size() ; i++ ){
          streamlog_out(MESSAGE0) << i << " of " <<  _trackStates.size() << " location: " <<    _trackStates[i]->getLocation() << std::endl ;
        }
      streamlog_out(MESSAGE0) << std::endl;
  
    }

    /** Returns track state closest to the given point. @see TrackState.
     */
    virtual const EUTelTrackStateImpl* getClosestTrackState( float x, float y, float z ) const ;


    /** Returns track state for the given location - or NULL if not found. @see TrackState.
     *  location can be set to: AtIP, AtFirstHit, AtLastHit, AtCalorimeter, AtVertex, AtOther
     */
    virtual const EUTelTrackStateImpl* getTrackState( int location ) const ;


    /** Optionally ( check/set flag(LCIO::TRBIT_HITS)==1)  return the hits that have been used 
     *  to create this track.
     */
    virtual const EVENT::TrackerHitVec &  getTrackerHits() const ;

	//	int getType();
    

    // setters 
    virtual void  setTypeBit( int index , bool val=true) ;

    virtual void  setTx( float ) ;                          // stored in TrackState
    virtual void  setTy( float ) ;                          // stored in TrackState
    virtual void  setX( float ) ;                           // stored in TrackState
    virtual void  setY( float ) ;                           // stored in TrackState
    virtual void  setInvP( float ) ;                        // stored in TrackState
    virtual void  setCovMatrix( const float* ) ;            // stored in TrackState
    virtual void  setCovMatrix( const EVENT::FloatVec& ) ;  // stored in TrackState
    virtual void  setReferencePoint( const float* ) ;       // stored in TrackState

    virtual void  setChi2( float chi2 ) ;
    virtual void  setNdf( int ndf ) ;
    virtual void  setdEdx( float dEdx ) ;
    virtual void  setdEdxError( float dEdxError ) ;

    virtual void  addTrack( EUTelTrackImpl* trk ) ;
    virtual void  addTrackState( EUTelTrackStateImpl* trkstate ) ;

    // direct access to the track state vector 
    virtual  EUTelTrackStateVec & trackStates() ;


  protected:
      
      std::bitset<32> _type ;

    float _chi2 ;
    int   _ndf ; 
    float _dEdx ;
    float _dEdxError ;

    EUTelTrackVec _tracks ;
    EVENT::TrackerHitVec _hits ;

    EUTelTrackStateVec _trackStates ;



}; // class


} // namespace eutelescope
#endif /* ifndef EUTELTRACKIMLP_H */

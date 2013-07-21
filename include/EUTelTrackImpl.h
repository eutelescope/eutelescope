#ifndef IMPL_EUTELTRACKIMPL_H
#define IMPL_EUTELTRACKIMPL_H 1

#include "IMPL/AccessChecked.h"
#include "EVENT/TrackerHit.h"

#include "LCIOSTLTypes.h"

#include <map>



namespace eutelescope {

  class EUTelTrackImpl;
  /**Vector of (pointers to) Tracks.*/
  typedef std::vector<EUTelTrackImpl*> EUTelTrackVec;

  class EUTelTrackStateImpl ;
  /**Vector of (pointers to) TrackStates.*/
  typedef std::vector<EUTelTrackStateImpl*> EUTelTrackStateVec ;
  
  class EUTelTrackImpl : public IMPL::AccessChecked {
    
  public: 
    
    /** Default constructor, initializes values to 0.
     */
    EUTelTrackImpl() ;
    
    /** Copy constructor - creates shallow copy, i.e. all data members are copied but pointers to other LCObjects
     *  i.e. TrackerHits and Tracks are preserved.
     */
    EUTelTrackImpl( const EUTelTrackImpl& ) ; 
    const EUTelTrackImpl& operator = ( const EUTelTrackImpl& o ) ; // operator =
    

    /// Destructor.
    virtual ~EUTelTrackImpl() ; 
    
    virtual int id() const { return simpleUID() ; }


    virtual float getTx() const ;

    virtual float getTy() const ;

    virtual float getX() const ;

    virtual float getY() const ;

    virtual float getInvP() const ;

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


    /** Returns track states associated with this track. @see TrackState.
     */
    virtual const EUTelTrackStateVec & getTrackStates() const ;


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
    virtual const EVENT::TrackerHitVec & getTrackerHits() const ;
    

    // setters 

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
    virtual void  addHit( EVENT::TrackerHit* hit) ;

    // direct access to the track state vector 
    virtual  EUTelTrackStateVec & trackStates() ;


  protected:

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

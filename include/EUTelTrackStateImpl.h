#ifndef EUTELTRACKSTATEIMPL_H
#define EUTELTRACKSTATEIMPL_H 1

#include "IMPL/AccessChecked.h"

#include "LCIOSTLTypes.h"

#include <map>


#define TRKSTATENCOVMATRIX 15
#define TRKSTATENREFSIZE 3

namespace eutelescope {

  class EUTelTrackStateImpl : public IMPL::AccessChecked {

  public: 

    /** Default constructor, initializes values to 0.
     */
    EUTelTrackStateImpl() ;
    EUTelTrackStateImpl(int, float, float, float, float, float, const float*, const float* ) ;
    EUTelTrackStateImpl(int, float, float, float, float, float, const EVENT::FloatVec&, const float* ) ;
    /** Copy constructor which takes as an argument an EVENT::EUTelTrackState reference */
    EUTelTrackStateImpl(const EUTelTrackStateImpl &p );


    
    /// Destructor.
    virtual ~EUTelTrackStateImpl() ; 
    
    static const int AtOther = 0 ; // any location other than the ones defined below	     
    static const int AtFirstHit = 1 ; 							    
    static const int AtLastHit = 2 ;							        

    virtual int id() const { return simpleUID() ; }


    /** The location of the track state.
     *  Location can be set to: AtIP, AtFirstHit, AtLastHit, AtCalorimeter, AtVertex, AtOther
     */
    virtual int getLocation() const ;

    virtual float getTx() const ;

    virtual float getTy() const ;

    virtual float getX() const ;

    virtual float getY() const ;

    virtual float getInvP() const ;

    /** Covariance matrix of the track parameters. Stored as lower triangle matrix where
     * the order of parameters is:   x, y, tx, ty, q/p.
     * So we have cov(x,x), cov( y, x ), cov( y, y ), ...
     */
    virtual const EVENT::FloatVec & getCovMatrix() const ;

    /** Reference point of the track parameters, e.g. the origin at the IP, or the position
     *  of the first/last hits or the entry point into the calorimeter.
     */
    virtual const float* getReferencePoint() const ;
   

    // setters 
    virtual void  setLocation( int location ) ;
    virtual void  setTx( float ) ;
    virtual void  setTy( float ) ;
    virtual void  setX( float ) ;
    virtual void  setY( float ) ;
    virtual void  setInvP( float ) ;

    virtual void  setCovMatrix( const float* ) ;
    virtual void  setCovMatrix( const EVENT::FloatVec& ) ;

    virtual void  setReferencePoint( const float* ) ;


  protected:

    int _location ; // location defined by EUTelTrackStateLocationEnum
    float _tx ;
    float _ty ;
    float _x ;
    float _y ;
    float _invp ;

    EVENT::FloatVec _covMatrix ;
    float  _reference[TRKSTATENREFSIZE] ;

}; // class

} // namespace IMPL
#endif /* ifndef EUTELTRACKSTATEIMLP_H */


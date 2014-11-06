#ifndef EUTELTRACKSTATEIMPL_H
#define EUTELTRACKSTATEIMPL_H 1

#include "streamlog/streamlog.h"
#include <iostream>


#include "IMPL/AccessChecked.h"

#include "IMPL/TrackStateImpl.h"


#include "LCIOSTLTypes.h"

#include <map>

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#endif

#include "EUTelGeometryTelescopeGeoDescription.h"


#define TRKSTATENCOVMATRIX 15
#define TRKSTATENREFSIZE 3

namespace eutelescope {

  class EUTelTrackStateImpl : public IMPL::TrackStateImpl { //, IMPL::AccessChecked {

  public: 

    /** Default constructor, initializes values to 0.
     */
    EUTelTrackStateImpl() ;
		EUTelTrackStateImpl(const IMPL::TrackStateImpl& o );
    EUTelTrackStateImpl(int, float, float, float, float, float, const float*, const float*) ;
    EUTelTrackStateImpl(int, float, float, float, float, float, const EVENT::FloatVec&, const float* ) ;
    /** Copy constructor which takes as an argument an EVENT::EUTelTrackState reference */
    EUTelTrackStateImpl(const EUTelTrackStateImpl &p );


    
    /// Destructor.
    virtual ~EUTelTrackStateImpl() ; 
    
    static const int AtOther    =  0 ; // any location other than the ones defined below	     
    static const int AtFirstHit = -1 ; 							    
    static const int AtLastHit  = -2 ;							        

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

	virtual float getZParameter() const ;

    	TVectorD getTrackStateVec() const ;

	TMatrixDSym getTrackStateCov() const;

	virtual TVector3 getPfromCartesianParameters() const;

	virtual int findIntersectionWithCertainID(int , float*) const; 

	virtual TVector3 getXYZfromArcLength( float s ) const;

    	virtual TMatrixD getH() const;

	virtual TMatrix  getPropagationJacobianF( float dz );

		virtual EVENT::TrackerHit* getHit() const;

    /** Covariance matrix of the track parameters. Stored as lower triangle matrix where
     * the order of parameters is:   x, y, tx, ty, q/p.
     * So we have cov(x,x), cov( y, x ), cov( y, y ), ...
     */
    	virtual const EVENT::FloatVec & getCovMatrix() const ;

    /** Reference point of the track parameters, e.g. the origin at the IP, or the position
     *  of the first/last hits or the entry point into the calorimeter.
     */

    virtual const float* getReferencePoint() const ;

		virtual TVector3 getIncidenceVectorInLocalFrame();

		virtual void getTrackStateHitCov( double (&cov)[4]);

   
    	virtual void Print() const;

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

		virtual void setZParameter(float);

		virtual void setbeamQ(int);

		virtual void setHit( EVENT::TrackerHit* hit);

		virtual void setTrackStateHitCov(double cov[4]);



  protected:

		//Static since there is not point in using more memory than needed. Since the particle change will always be the same for a single run
		int _beamQ; //This seems a strange parameter to store here. However to state variable _invp is q/p. Therefore to determine p here you need Q. _invp should be renamed or changed? 

    int _location ; // location defined by EUTelTrackStateLocationEnum
    float _tx ;
    float _ty ;
    float _x ;
    float _y ;
    float _invp ;

    EVENT::FloatVec _covMatrix ;
    float  _reference[TRKSTATENREFSIZE] ;

		float _zparameter;
	
		EVENT::TrackerHit* _hit;

		double _covHitMatrix[4];

}; // class

} // namespace IMPL
#endif /* ifndef EUTELTRACKSTATEIMLP_H */


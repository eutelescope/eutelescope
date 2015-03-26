// Version: $Id$
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelUtility.h"

//#include "TrackerHitImpl2.h"
#include "IMPL/TrackerHitImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "mille/Mille.h"

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TMinuit.h>
#include <TSystem.h>
#include <TMath.h>
#include <TVector3.h>
class TMinuit;
#else
#error *** You need ROOT to compile this code.  *** 
#endif


namespace eutelescope {

 
  //Specify a Rectangular in a sensor
  class SensorRectangular {
	protected:
	// SensorID
	int sensor;
	// lowest pixel in X direction
	int A; 
	// lowest pixel in Y direction
	int B;
	// highest pixel in X direction
	int C;
	// highest pixel in Y direction
	int D;
	public:
	SensorRectangular(int s, int a, int b, int c, int d) : sensor(s), A(a), B(b), C(c), D(d) {};
	SensorRectangular() : sensor(0), A(0), B(0), C(0), D(0) {};
	int getSensor() const {return sensor; }
	//look if x and y are inside the foreseen rectangular
	bool isInside(int x, int y) const {  return (x >=A && x <=C && y >=B && y <=D); }
	void print() { streamlog_out(MESSAGE4) << "Sensor: " << sensor << ": (" << A << "|" << B << ") to (" << C << "|" << D << ")" << std::endl; }
	
  };
  
  class RectangularArray {
	protected:
		std::map<int,SensorRectangular > _rect;
	
	public:
	void addRectangular(SensorRectangular &s) { _rect[s.getSensor() ] = s;}
	
	bool isInside(int s, int x, int y) {
		std::map<int,SensorRectangular >::iterator it = _rect.find(s);
		if (it == _rect.end()) { // not in the map means no limit on this sensor -> always true
			return true;
		}
		SensorRectangular cSensor = _rect[s];
		return cSensor.isInside(x,y);
	}
		
  };


  //! Straight line fit processor
  /*!
   *
   */

  class EUTelMille : public marlin::Processor {

  public:
    class hit
    {
    public:
      hit()
      {
      }
      hit(double tx, double ty, double tz, double rx, double ry, double rz,int i)
      {
        x = tx;
        y = ty;
        z = tz;
        resolution_x = rx;
        resolution_y = ry;
        resolution_z = rz;
    
        planenumber = i;
      }
      double x;
      double y;
      double z;
      double resolution_x;
      double resolution_y;
      double resolution_z;
  
      int planenumber;
    };

   

    class trackfitter
    {
    public:
      ~trackfitter()
      {}
      trackfitter()
      {}
     //  trackfitter(std::vector<hit> &h)
//       {
//         hitsarray = h;
//       }
      trackfitter(hit *h, unsigned int num)
      {
        hitsarray = h;
        n = num;
      }

      double dot(const double *a, const double *b) const
      {
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
      }
      double fit (double *x)
      {
        double chi2 = 0.0;
   
        const unsigned int n = 3;

        const double b0 = x[0];
        const double b1 = x[1];
        const double b2 = 0.0;
    
        const double alpha = x[2];
        const double beta =  x[3];

        const double c0 = TMath::Sin(beta);
        const double c1 = -1.0*TMath::Cos(beta) * TMath::Sin(alpha);
        const double c2 = TMath::Cos(alpha) * TMath::Cos(beta);
    
        double c[n] = {c0, c1, c2}; 
        for(size_t i = 0; i < n;i++)
          {
            const double p0 = hitsarray[i].x;
            const double p1 = hitsarray[i].y;
            const double p2 = hitsarray[i].z;
        
            const double resol_x = hitsarray[i].resolution_x;
            const double resol_y = hitsarray[i].resolution_y;
            const double resol_z = hitsarray[i].resolution_z;
        
            const double pmb[n] = {p0-b0, p1-b1, p2-b2}; //p - b
        
            const double coeff = dot(c, pmb);
            const double t[n] = {
              b0 + c0 * coeff - p0,
              b1 + c1 * coeff - p1,
              b2 + c2 * coeff - p2
            }; 
        
            //sum of distances divided by resolution^2
            chi2 += t[0]*t[0] / pow(resol_x,2) 
              + t[1]*t[1] / pow(resol_y,2) 
              + t[2]*t[2] / pow(resol_z,2);
          }
 
        return chi2;
      }
    private:
      //std::vector<hit> hitsarray;
    hit *hitsarray;
    unsigned int n;
    };

    //! Variables for hit parameters
    class HitsInPlane {
    public:
      HitsInPlane(){
        measuredX = 0.0;
        measuredY = 0.0;
        measuredZ = 0.0;
      }
      HitsInPlane(double x, double y, double z)
      {
        measuredX = x;
        measuredY = y;
        measuredZ = z;
      }
      bool operator<(const HitsInPlane& b) const
      {
        return (measuredZ < b.measuredZ);
      }
      double measuredX;
      double measuredY;
      double measuredZ;
    };

    virtual void FitTrack(
                          unsigned int nPlanesFitter,
                          double xPosFitter[],
                          double yPosFitter[],
                          double zPosFitter[],
                          double xResFit[],
                          double yResFit[],
                          double chi2Fit[2],
                          double residXFit[],
                          double residYFit[],
                          double angleFit[2]
                          );


    //recursive method which searches for track candidates - with omits!
    virtual void findtracks2(
                            int missinghits,
                            std::vector<IntVec > &indexarray, //resulting vector of hit indizes
                            IntVec vec, //for internal use
                            std::vector<std::vector<EUTelMille::HitsInPlane> > &_hitsArray, //contains all hits for each plane
                            unsigned int i, //plane number
                            int y //hit index number
                            );


    //recursive method which searches for track candidates
    virtual void findtracks(
                            std::vector<IntVec > &indexarray, //resulting vector of hit indizes
                            IntVec vec, //for internal use
                            std::vector<std::vector<EUTelMille::HitsInPlane> > &_hitsArray, //contains all hits for each plane
                            int i, //plane number
                            int y //hit index number
                            );

    //! Returns a new instance of EUTelMille
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelMille.
     */
    virtual Processor * newProcessor() {
      return new EUTelMille;
    }

    virtual bool hitContainsHotPixels( TrackerHitImpl   * hit) ;


    //! Default constructor
    EUTelMille ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called for first event per run
    /*! Reads hotpixel information from hotPixelCollection into hotPixelMap
     * to be used in the sensor exclusion area logic 
     */
    virtual void  FillHotPixelMap(LCEvent *event);

    //! Called every event
    /*! This is called for each event in the file. Each element of the
     *  pulse collection is scanned and the center of the cluster is
     *  translated into the external frame of reference thanks to the
     *  GEAR geometry description.
     *
     *  The cluster center might be calculate using a standard linear
     *  charge center of gravity algortihm or applying a more
     *  sophisticated non linear eta function. This behaviour is
     *  regulated by the user from the steering file.
     *
     *  @throw UnknownDataTypeException if the cluster type is unknown
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);

    //! Track search method, fills class members _xPos and _trackResidX (Y,Z as well) 
    /** */
    void findMatchedHits( int &, Track* ) ;

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();


    //! Histogram booking
    /*! Some control histograms are filled during this procedure in
     *  order to be able to perform easy check on the quality of the
     *  output hits and also to understand if the frame of reference
     *  conversion has been properly done. Of course this method is
     *  effectively doing something only in the case MARLIN_USE_AIDA.
     */
    void bookHistos();

    TVector3 Line2Plane(int iplane, const TVector3& lpoint, const TVector3& lvector ); 

    virtual inline int getAllowedMissingHits(){return _allowedMissingHits;}
    virtual inline int getMimosa26ClusterChargeMin(){return _mimosa26ClusterChargeMin;}


  protected:


    //! Ordered sensor ID
    /*! Within the processor all the loops are done up to _nPlanes and
     *  according to their position along the Z axis (beam axis).
     *
     *  This vector is containing the sensorID sorted according to the
     *  same rule.
     */
    IntVec _orderedSensorID;
    IntVec _orderedSensorID_wo_excluded;


    //! reference HitCollection name 
    /*!
     */
    std::string      _referenceHitCollectionName;
    bool             _useReferenceHitCollection;
    LCCollectionVec* _referenceHitVec;    
 
    //! TrackerHit collection name
    /*! Input collection with hits.
     */
    StringVec _hitCollectionName;

    //! TRACK collection name
    /*! Output collection with fitted tracks.
     */
    std::string _trackCollectionName;

    //! Hot pixel collection name.
    /*! 
     * this collection is saved in a db file to be used at the clustering level
     */
    std::string _hotPixelCollectionName;

    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  first level key   sensor unique 
     *              value sensor map
     *  sensor map key    unique row number
     *             value  vector of column numbers.
     */
    
    std::map<std::string, bool > _hotPixelMap;

    //! Sensor ID vector
    IntVec _sensorIDVec;

    //! Sensor ID map (inverse sensorIDVec) 
    std::map< int, int > _sensorIDVecMap;
 
    // parameters

    float _distanceMax;
    FloatVec _distanceMaxVec;
    std::vector<unsigned int > _excludePlanes; //only for internal usage
    IntVec _excludePlanes_sensorIDs; //this is going to be
                                                //set by the user.
    IntVec _FixedPlanes; //only for internal usage
    IntVec _FixedPlanes_sensorIDs; //this is going to be
    //set by the user.
    
    StringVec _pedeSteerAddCmds; // allows user-added commands in the pede steering file


    int _maxTrackCandidates;
    int _maxTrackCandidatesTotal;

    std::string _binaryFilename;

    float _telescopeResolution;
    bool _onlySingleHitEvents;
    bool _onlySingleTrackEvents;
    int _alignMode;
    bool _useResidualCuts;

    FloatVec _residualsXMin;
    FloatVec _residualsYMin;
    FloatVec _residualsXMax;
    FloatVec _residualsYMax;

    FloatVec _resolutionX;
    FloatVec _resolutionY;
    FloatVec _resolutionZ;

    IntVec _FixParameter;


    int _generatePedeSteerfile;
    std::string _pedeSteerfileName;
    bool _runPede;
    int _usePedeUserStartValues;
    FloatVec _pedeUserStartValuesX;
    FloatVec _pedeUserStartValuesY;
    FloatVec _pedeUserStartValuesZ;
    
    FloatVec _pedeUserStartValuesAlpha;
    FloatVec _pedeUserStartValuesBeta;
    FloatVec _pedeUserStartValuesGamma;

    int _inputMode;
    int _allowedMissingHits;
    int _mimosa26ClusterChargeMin;

    float _testModeSensorResolution;
    float _testModeXTrackSlope;
    float _testModeYTrackSlope;

    FloatVec _testModeSensorZPositions;

    FloatVec _testModeSensorXShifts;
    FloatVec _testModeSensorYShifts;
    FloatVec _testModeSensorGamma;
    FloatVec _testModeSensorAlpha;
    FloatVec _testModeSensorBeta;

    std::string _alignmentConstantLCIOFile;
    std::string _alignmentConstantCollectionName;

    IntVec _useSensorRectangular;

  private:

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

    // Excluded planes
    int _nExcludePlanes;

    // Statistics
    int _nMilleDataPoints;
    int _nMilleTracks;

    // Mille
    Mille * _mille;

    //! Conversion ID map.
    /*! In the data file, each cluster is tagged with a detector ID
     *  identify the sensor it belongs to. In the geometry
     *  description, there are along with the sensors also "passive"
     *  layers and other stuff. Those are identify by a layerindex. So
     *  we need a conversion table to go from the detectorID to the
     *  layerindex.
     */
    std::map< int, int > _conversionIdMap;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;
    std::map<std::string, AIDA::IHistogram1D * > _aidaHistoMap1D;
    std::map<std::string, AIDA::IHistogram2D * > _aidaHistoMap2D;
    std::map<std::string, AIDA::IProfile1D * >   _aidaHistoMapProf1D;
 
    static std::string _numberTracksLocalname;

    static std::string _chi2XLocalname;
    static std::string _chi2YLocalname;

    static std::string _residualXLocalname;
    static std::string _residualYLocalname;
    static std::string _residualZLocalname;

    static std::string _residualXvsXLocalname;
    static std::string _residualXvsYLocalname;
    static std::string _residualYvsXLocalname;
    static std::string _residualYvsYLocalname;
    static std::string _residualZvsXLocalname;
    static std::string _residualZvsYLocalname;
 
#endif

    size_t _nPlanes;

    std::vector<DoubleVec > _xPos;
    std::vector<DoubleVec > _yPos;
    std::vector<DoubleVec > _zPos;

    std::vector<DoubleVec > _trackResidX;
    std::vector<DoubleVec > _trackResidY;
    std::vector<DoubleVec > _trackResidZ;


//    std::vector<DoubleVec > _xPosTrack;
//    std::vector<DoubleVec > _yPosTrack;
//    std::vector<DoubleVec > _zPosTrack;

    double * _xPosHere;
    double * _yPosHere;
    double * _zPosHere;
    double * _waferResidX;
    double * _waferResidY;
    double * _waferResidZ;
    double * _telescopeResolX;
    double * _telescopeResolY;
    double * _telescopeResolZ;
    double * _xFitPos;
    double * _yFitPos;

    DoubleVec _siPlaneZPosition;

    //! Fill histogram switch
    /*! Only for debug reason
     */
    bool _histogramSwitch;

    //! Limits the pixels on each sensor-plane to a sub-rectangular  
    RectangularArray _rect;

  };

  //! A global instance of the processor
  EUTelMille gEUTelMille;

}
#endif
#endif


/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCORRELATOR_H
#define EUTELCORRELATOR_H
#if defined(USE_GEAR)

// eutelescope includes ".h"

//ROOT includes
#include "TVector3.h"


// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
//#include <TrackerHitImpl2.h>
#include <IMPL/TrackerHitImpl.h>



// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#endif


// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {

  //! Hit and cluster correlator
  /*! This processor makes histograms that show the correlation
   *  between clusters of a detector and another.
   *  We study the correlation between centre of clusters' X and Y
   *  separately.
   *  At the end of our study we'll have (n*n - n) histograms for X
   *  and (n*n - n) for Y ( n is the number of our detectors ). We
   *  have this number of histograms because we don't create those
   *  histograms where there will be the correlation between a
   *  detector and himself.
   *  These histograms are put into different directories,
   *  one for X and one for Y, with cluster and hit correlation too.
   *  For each directory we have pairs of histograms that differ only
   *  in the order of DetectorID and these result with X and Y
   *  reversed, but otherwise they are the same.
   *  We do the same for hit correlation.
   *
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Cluster collection</b>: A collection with cluster
   *  produced by previous processors like EUTelClusteringProcessor or
   *  EUTelClusterFilter.
   *
   *
   *  @author Silvia Bonfanti, Uni. Insubria  <mailto:silviafisica@gmail.com>
   *  @author Loretta Negrini, Uni. Insubria  <mailto:loryneg@gmail.com>
   *  @version $Id$
   *
   */

  class EUTelCorrelator : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelCorrelator
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelCorrelator.
     */
    virtual Processor * newProcessor() {
      return new EUTelCorrelator;
    }

    //! Default constructor
    EUTelCorrelator ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! For each event, we loop over the input cluster collection and
     *  we fill in the correlation histograms using the x and y
     *  coordinate of the cluster center. Only clusters not belonging
     *  to the same sensors are used.
     *
     *  Histograms are booked in the first event.
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();


    //! Histogram booking
    /*! This method is used to book all the correlation histograms. It
     *  is called by processEvent when processing the first event.
     */
    void bookHistos();

    //! internal functtion: return the ID of a plane selected as a reference plane for correlation plots

    virtual int getFixedPlaneID(){return _fixedPlaneID;} 


  protected:

    //! set the plane you would like to use as a reference/strating point for correlation plots
    int _fixedPlaneID;

    //! vector of correlation band cuts in X (upper limit)
    std::vector< float  > _residualsXMax;
    //! vector of correlation band cuts in X (lower limit) 
    std::vector< float  > _residualsXMin;
    //! vector of correlation band cuts in Y (upper limit)       
    std::vector< float  > _residualsYMax;         
    //! vector of correlation band cuts in Y (lower limit) 
    std::vector< float  > _residualsYMin;

    int _minNumberOfCorrelatedHits;

    //! Input collection name.
    /*! This is the name of the output hit collection.
     */
    std::string _inputHitCollectionName;

    //! Output collection name.
    /*! This is the name of the output hit collection.
     */
    std::string _outputHitCollectionName;

    //! output collection for correlated 
    /*! 
     */ 
    LCCollectionVec* _outputCorrelatedHitCollectionVec;
 
    //! Input cluster charge cut
    /*!
     */
    int _clusterChargeMin;

    //! How many events are needed to get reasonable correlation plots 
    /*! (and Offset DB values) 
     *
     */
    int _events;

    //! Cluster collection list (EVENT::StringVec) 
    /*!
     */
    EVENT::StringVec  _clusterCollectionVec;

    std::vector<double> guessSensorOffset(int internalSensorID, int externalSensorID, std::vector<double> cluCenter );

  private:

    //! Initialization flag
    bool _isInitialize;

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

    //! Number of sensors
    int _noOfDetectors;

    //! First pixel along X
    /*! This is an associative map relating the sensorID to the first 
     *  pixel along X
     */
    std::map< int, int> _minX;

    //! Last pixel along X
    /*! This is an associative map relating the sensorID to the last 
     *  pixel along X
     */

    std::map< int, int> _maxX;

    //! First pixel along Y
    /*! This is an associative map relating the sensorID to the first 
     *  pixel along Y
     */

    std::map< int, int> _minY;

    //! Last pixel along Y
    /*! This is an associative map relating the sensorID to the last 
     *  pixel along Y
     */

    std::map< int, int> _maxY;

   //! First pixel along X
    /*! This is an associative map relating the sensorID to the first 
     *  pixel along X
     */
    std::map< float, float> _hitMinX;

    //! Last pixel along X
    /*! This is an associative map relating the sensorID to the last 
     *  pixel along X
     */

    std::map< float, float> _hitMaxX;

    //! First pixel along Y
    /*! This is an associative map relating the sensorID to the first 
     *  pixel along Y
     */

    std::map< float, float> _hitMinY;

    //! Last pixel along Y
    /*! This is an associative map relating the sensorID to the last 
     *  pixel along Y
     */

    std::map< float, float> _hitMaxY;

    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    //! vector of Rotation Matrix elements
    std::vector< std::map<int,double> > _siPlanesRotations;

    //! vector of Sensor Pitch X
    std::vector< double > _siPlanesPitchX;

    //! vector of Sensor Pitch Y
    std::vector< double > _siPlanesPitchY;

    //! vector of Sensor Offset X
    std::vector< double > _siPlanesOffsetX;

    //! vector of Sensor Offset Y
    std::vector< double > _siPlanesOffsetY;


    //! An array with the Z position of planes
    double * _siPlaneZPosition;

    //! Sensor ID map (inverse sensorIDVec) 
    std::map< int, int > _sensorIDVecMap;


    //! Sensor ID vector, 
    /*! it's position along Z axis
     */ 
    std::vector< int > _sensorIDVecZOrder;

    //! sensor ID to position along Z id
    /*!
     * 
     */
    std::map<int, int> _sensorIDtoZOrderMap;

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

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

    /** Histogram info file name */
    std::string _histoInfoFileName;


    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    //! Correlation histogram matrix
    /*! This is used to store the pointers of each histogram
     */
    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _clusterXCorrelationMatrix;
    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _clusterYCorrelationMatrix;

    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _clusterXCorrShiftMatrix;
    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _clusterYCorrShiftMatrix;
    std::map< unsigned int , AIDA::IHistogram1D*  > _clusterXCorrShiftProjection;
    std::map< unsigned int , AIDA::IHistogram1D*  > _clusterYCorrShiftProjection;

    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _hitXCorrelationMatrix;
    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _hitYCorrelationMatrix;
    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _hitXCorrShiftMatrix;
    std::map< unsigned int , std::map< unsigned int , AIDA::IHistogram2D* > > _hitYCorrShiftMatrix;
    std::map< unsigned int , AIDA::IHistogram1D*  > _hitXCorrShiftProjection;
    std::map< unsigned int , AIDA::IHistogram1D*  > _hitYCorrShiftProjection;


    //! Base name of the correlation histogram
    static std::string _clusterXCorrelationHistoName;
    static std::string _clusterYCorrelationHistoName;
    static std::string _clusterXCorrShiftHistoName;
    static std::string _clusterYCorrShiftHistoName;
    static std::string _clusterXCorrShiftProjectionHistoName;
    static std::string _clusterYCorrShiftProjectionHistoName;
    
    static std::string _hitXCorrelationHistoName;
    static std::string _hitYCorrelationHistoName;
    static std::string _hitXCorrShiftHistoName;
    static std::string _hitYCorrShiftHistoName;
    static std::string _hitXCorrShiftProjectionHistoName;
    static std::string _hitYCorrShiftProjectionHistoName;

#endif

    bool _hasClusterCollection;
    bool _hasHitCollection;

    std::vector<int> _sensorIDVec;
    std::map<int, int> _sensorIDtoZ;
  };

  //! A global instance of the processor
  EUTelCorrelator gEUTelCorrelator;


}

#endif // USE_GEAR
#endif

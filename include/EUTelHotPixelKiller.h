/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELHOTPIXELKILLER
#define EUTELHOTPIXELKILLER 1

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelSimpleSparsePixel.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/LCCollectionVec.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <map>


namespace eutelescope {

  //! Processor to mask hot pixels
  /*! This processor is used to keep hot matrix out from the analysis
   *  procedure. This processor is based on the idea that if a pixel
   *  is found to be a part of a cluster to often, probably it is a noisy
   *  pixel and should be removed from the game.
   *
   *  The input status collection has to be writable, so if it is read
   *  in from a file, it has to be copied locally (EUTelCopyPedestalProcessor).
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Status collection </b> the current status collection
   *
   *  <h4>Output collections</h4>
   *
   *  @param NoOfEventPerCycle The number of event to take before
   *  proceeding with the masking
   *  @param MaxAllowedFiringFreq This number [0,1] represents the
   *  maximum allowed firing frequency. Set it to a suitable value
   *  depending on the occupancy.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   *
   */

  class EUTelHotPixelKiller : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelHotPixelKiller
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelHotPixelKiller.
     */
    virtual Processor * newProcessor() {
      return new EUTelHotPixelKiller;
    }

    //! Default constructor
    EUTelHotPixelKiller ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and performs some asserts about
     *  the value of the provided parameters
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run LCRunHeader of the this current run
     *
     *  @throw InvalidParameterException if a paramter is wrongly set
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. If the current @c
     *  evt is flagged to be used for update, then the selected
     *  algorithm wrapper is called
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     *
     *  @throw InvalidParameterException if information in the cellID
     *  are inconsistence
     */
    virtual void processEvent (LCEvent * evt);

    //! Initialize geometry
    /*! Set the number of detectors in the setup and their boundaries.
     *
     *  @param evt The LCIO event
     */
    virtual void initializeGeometry( LCEvent * evt ) ;

    //! HotPixelFinder
    /*!
     */
    void HotPixelFinder(EUTelEventImpl *input);
    
    int getBuildHotPixelDatabase() const   { return _flagBuildHotPixelDatabase; }
    
    //! Check call back
    /*! This method is called every event just after the processEvent
     *  one. For the time being it is just calling the pixel
     *  monitoring protected method
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void check(LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished. Just printing a good bye message
     */
    virtual void end();


  protected:

    std::string _lcioWriteMode ; 


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! Histogram with the firing frequency 2D distribution
    static std::string _firing2DHistoName;

    //! Histogram with the firing cumulative 1D distribution
    static std::string _firing1DHistoName;

    //! book histogram method
    void bookAndFillHistos();


#endif

    //! Print the summary
    std::string printSummary() const ;



    //! Input collection name for ZS data
    /*! The input collection is the calibrated data one coming from
     *  the input data file. It is, usually, called
     *  "zsdata" and it is a collection of TrackerData
     */
    std::string _zsDataCollectionName;

    //! Noise collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _noiseCollectionName;

    //! Status collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _statusCollectionName;

    //! Hot pixel collection name.
    /*! 
     * this collection is saved in a db file to be used at the clustering level
     */
    std::string _hotPixelCollectionName;

 
    //! The excluded planes list
    /*! This is a list of sensor ids for planes that have to be
     *   excluded from the clustering.
     */
    std::vector<int> _ExcludedPlanes;


    //! Number of events for update cycle
    int _noOfEventPerCycle;

    //! Maximum allowed firing frequency
    float _maxAllowedFiringFreq;

    //! Map relating ancillary collection position and sensorID
    /*! The first element is the sensor ID, while the second is the
     *  position of such a sensorID in all the ancillary collections
     *  (noise, pedestal and status).
     */
    std::map< int, int > _ancillaryIndexMap;

 
    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  Key is the sequential (counter) id of a hit,
     *  Value - (sensor) unique pixel Id.
     */

    std::vector< std::map< int, int > > _hitIndexMapVec;
  
        
    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created (inverse
     *  mapping to _hitIndexMapVec). 
     *  Key is a (sensor) unique pixel Id 
     *  Value - sequential (counter) id of a hit.
     */
    
    std::vector< std::map< int, int > > _inverse_hitIndexMapVec;

    
    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  Key is a (sensor) unique pixel Id (to be addressed via
     *  _inverse_hitIndexMapVec)
     *  Value - EUTelSimpleSparsePixel pointer.
     */
    
    std::vector< std::map< int, EUTelSimpleSparsePixel* > > _pixelMapVec;



    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Number of detector planes in the run
    /*! It is the total number of elements in the status collection.
     */
    int _noOfDetectors;

    //! First pixel along X
    /*! An associative map between sensorID and the corresponding
     *  minimum X pixel.
     */
    std::map<int, int > _minX;

    //! Last pixel along X
    /*! An associative map between sensorID and the corresponding
     *  maximum X pixel.
     */
    std::map<int, int > _maxX;

    //! First pixel along Y
    /*! An associative map between sensorID and the corresponding
     *  minimum Y pixel.
     */
    std::map<int, int > _minY;

    //! Last pixel along Y
    /*! An associative map between sensorID and the corresponding
     *  maximum Y pixel.
     */
    std::map<int, int > _maxY;

    //! Sensor ID vector
    std::vector< int > _sensorIDVec;

    //! Total number of cycle
    int _totalNoOfCycle;

    //! Hot Pixel DB output file
    std::string _hotpixelDBFile;

  private:

    //! The current updating cycle number
    unsigned short _iCycle;

    //! A vector with the number of killed pixel
    std::vector< std::vector< unsigned short > > _killedPixelVec;

    //! A vector with the firing frequency value
    std::vector< std::vector< unsigned short > > _firingFreqVec;

    //! Simple data decoding and HotPixel database
    /*
     */
    int _flagBuildHotPixelDatabase; 
    
    //! write out the list of hot pixels
    /*!
     */
    void HotPixelDBWriter(LCEvent * event);


  };

  //! A global instance of the processor
  EUTelHotPixelKiller gEUTelHotPixelKiller;

}
#endif

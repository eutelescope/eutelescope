/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPROCESSORNOISYPIXELFINDER_H
#define EUTELPROCESSORNOISYPIXELFINDER_H

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelGenericSparsePixel.h"

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

//! Simple helper struct for a sensor pixel layout
/*  Contains the pixel count as well as the offset values
 *  this is necessay if fore example pixel indices range from
 *  50-99 instead of 0-50. The offset allows mapping to a
 *  range from 0 to (pixel count - 1) which is required to
 *  access the array in which hits are counted.
 */
struct sensor {
	int offX, offY;
	int sizeX, sizeY;
};

//! Processor to write out hot pixels 
/*! This processor is used to keep hot matrix out from the analysis
 *  procedure. It checks if pixels fired above a certain frequency
 *  and writes them into a TrackerData collection in an external 
 *  file if they did.
 *  You have to specify which plane IDs shall be processed (SensorIDVec)
 *  as well as which ones shall be excluded (ExcludedPlanes). This is
 *  necessary for histogramming and fast data processing.
 *
 *  @param NoOfEvents The amount of events to determine the firing frequency
 *
 *  @param SensorIDVec An integer vector containing the sensor IDs of the
 *  planes which are processed
 *
 *  @param MaxAllowedFiringFreq The firing frequency cut which shall be applied
 *
 *  @param HotPixelDBFile Name of the output file, currently appending to a file 
 *  does not work, use differnt files if multiple hot pixel collections are created
 *
 *  @param ExcludedPlanes Planes to be excluded from processing
 *
 *  @param HotPixelCollectionName The name of the collection in the output file
 */
class EUTelProcessorNoisyPixelFinder : public marlin::Processor {

public:
    //! Returns a new instance of EUTelProcessorNoisyPixelFinder
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelProcessorNoisyPixelFinder.
     */
    virtual Processor* newProcessor() {
      return new EUTelProcessorNoisyPixelFinder;
    }

    //! Default constructor
    EUTelProcessorNoisyPixelFinder();

    //! Default destructor
    ~EUTelProcessorNoisyPixelFinder();

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
    virtual void processEvent(LCEvent * evt);

    //! Initialize geometry
    /*! Set the number of detectors in the setup and their boundaries.
     *
     *  @param evt The LCIO event
     */
    void initializeHitMaps() ;

    //! HotPixelFinder
    void HotPixelFinder(EUTelEventImpl *input);
    
    //! Check call back
    /*! This method is called every event just after the processEvent
     *  one. For the time being it is just calling the pixel
     *  monitoring protected method
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void check( LCEvent* event );

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

    //! Input collection name for ZS data
    /*! The input collection is the calibrated data one coming from
     *  the input data file. It is, usually, called
     *  "zsdata" and it is a collection of TrackerData
     */
    std::string _zsDataCollectionName;

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
    int _noOfEvents;

    //! Maximum allowed firing frequency
    float _maxAllowedFiringFreq;
    
    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  Key is a (sensor) unique pixel Id (to be addressed via
     *  _inverse_hitIndexMapVec)
     */
    std::map<int, sensor> _sensorMap;

    //! Map holding the 2D-"array" which counts the hits
    /*! The key is the sensorID and the array is implemented as
     *  a vector of vectors. They are resized initially so there 
     *  is no real overhead with using std::vectors instead of
     * arrays.
     */
    std::map<int, std::vector<std::vector<int> >* > _hitVecMap;
    
    //! Map for storing the hot pixels in a std::vector as a value
    /*! The key is once again the sensorID.
     */
    std::map<int, std::vector<EUTelGenericSparsePixel> > _hotPixelMap;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Sensor ID vector
    /*! Passed as a argument via the steering file, here you
     *  specify for which sensors hot pixels should be determined
     */
    EVENT::IntVec _sensorIDVec;

    //! Hot Pixel DB output file
    std::string _hotpixelDBFile;

    //! write out the list of hot pixels
    void HotPixelDBWriter();

    //! Flag which will be set once we're done finding noisy pixels
    bool _finished;
};

//! A global instance of the processor
EUTelProcessorNoisyPixelFinder gEUTelProcessorNoisyPixelFinder;

}
#endif

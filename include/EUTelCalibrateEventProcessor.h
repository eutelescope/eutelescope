/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCALIBRATEEVENTPROCESSOR_H
#define EUTELCALIBRATEEVENTPROCESSOR_H 1

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// lcio includes <.h>

// system includes <>
#include <string>
#include <map>

namespace eutelescope {

  //! Calibration processor for the EUTelescope
  /*! This processor is used in a standard analysis sequence to apply
   *  the calibration (mainly pedestal value) as previously saved into
   *  a condition file or database.
   *
   *  The main things done by this processors are the following:
   *
   *  \li Cross check that the input data file and the condition
   *  calibration files are compatible. In other words, the number of
   *  detectors saved in the condition file has to be the same as the
   *  one in the input file. Also the number of pixels for each sensor
   *  has to be the same. In the case this cross check fails, the data
   *  analysis cannot continue and a
   *  eutelescope::IncompatibleDataSetException(std::string) is
   *  thrown.
   *
   *  \li In the case, the user would like to apply a common mode
   *  suppression procedure, this is done on a event and detector
   *  basis. If MARLIN_USE_AIDA is defined, a common mode distribution
   *  plot is also filled.
   *
   *  \li Each pixel raw signal is calibrated, in the meaning that the
   *  corresponding pedestal is subtracted and if the user switched it
   *  on, also the common mode is removed.
   *
   *  \li An output collection of TrackerData named "data" is created
   *  storing the calibrated information of each pixel. This has to be
   *  considered as the input collection for the cluster search
   *  processor
   *
   *  <h4>Input collections</h4>
   *  <br><b>RawDataCollection</b>. This is a collection of
   *  TrackerRawData containing all pixel signals as they are readout
   *  by the DAQ system. This collection is coming from the input
   *  slcio file.
   *
   *  <br><b>PedestalCollection</b>. This is a collection of
   *  TrackerData containing all pixel pedestal values as they are
   *  calculated by EUTelPedestalNoiseProcessor. This collection is
   *  coming in from the condition file.
   *
   *  <br><b>NoiseCollection</b>. This is a collection of TrackerData
   *  contaning all pixel noise values as they are calculated by
   *  EUTelPedestalNoiseProcessor. This collection is coming from
   *  the condition file and is used only if the user switched on the
   *  common mode calculation.
   *
   *  <br><b>StatusCollection</b>. This is a collection of
   *  TrackerRawData containing the status of each pixels. This
   *  collection is added to the event by the condition processor and
   *  is used only if common mode flag is active.
   *
   *  <h4>Output</h4>
   *
   *  <br><b>DataCollection</b>. This is a collection of TrackerData
   *  containing the calibrated signal of each pixel. This is the
   *  entry point for cluster search. The user can decide the name of
   *  this output collection via the steering parameter
   *  DataCollectionName
   *
   *  @param RawDataCollectionName Name of the input raw data collection
   *
   *  @param PedestalCollectionName Name of the input (condition)
   *  pedestal collection
   *
   *  @param NoiseCollectionName Name of the input (condition) noise
   *  collection
   *
   *  @param StatusCollectionName Name of the input (condition) status
   *  collection
   *
   *  @param DebugHistoFilling Flag to set the filling of debug
   *  histo on (true) or off (false)
   *
   *  @param PerformCommonMode Flag to set the common mode suppression on
   *  (true) or off (false)
   *
   *  @param HitRejectionCut Threshold in SNR to consider a pixel hit
   *  and exclude it from the common mode calculation.
   *
   *  @param MaxNoOfRejectedPixels Maximum allowed number of rejected
   *  pixels in a event. If this value is exceeded, than the event is
   *  skipped. If this parameter is set to -1 then the control is
   *  by-passed
   *
   *  @param DataCollectionName The name of the output calibrated data
   *  collection
   *
   *  @param HistoInfoFileName The name of the XML containing the
   *  histogram information file.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   *
   *
   */

  class EUTelCalibrateEventProcessor:public marlin::Processor {

  public:


    //! Returns a new instance of EUTelCalibrateEventProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelCalibrateEventProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelCalibrateEventProcessor;
    }

    //! Default constructor
    EUTelCalibrateEventProcessor ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and reset all needed data
     *  members. In the case the user set the _fillDebugHisto then
     *  she/he warned that the procedure is going to slow down
     *  considerably
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. From the run header the number of detector is
     *  retrieved.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. In the case this is
     *  the first event, then a cross-check is done on the
     *  compatibility of input raw data and pedestal collections. They
     *  have to contain the same number of detector planes and each of
     *  them should have exactly the same number of pixels. Of course
     *  this check is not exhaustive, but at least should avoid
     *  segmentation fault.
     *
     *  @throw IncompatibleDataSetException in the case the two
     *  collections are found to be incompatible
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. It can be used to fill check
     *  plots. For the time being there is nothing to check and do in
     *  this slot.
     *
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check (LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is finished. It
     *  prints only a goodbye message
     */
    virtual void end();

    //! Initialize the geometry information
    /*! This method is called to initialize the geometry information,
     *  namely the order in which the sensorID are organized into the
     *  input collection. In this function  the _ancillaryIndexMap
     *  is filled to get the correct pedestal frame object for the input data.
     *
     *  @param event The LCEvent to be used for geometry
     *  initialization.
     *
     *  @throw In case the @event does not contain all the needed
     *  information, a SkipEventException is thrown and the geometry
     *  will be initialize with the following event.
     */
    void initializeGeometry( LCEvent * event ) throw ( marlin::SkipEventException );

  protected:

    //! Input collection name.
    /*! For the time being we have just one collection that can be used as input
     */
    std::string _rawDataCollectionName;

    //! Pedestal collection name.
    /*! We have three different output collections that are going to be
     *  saved in the file within the first event possibly so that the
     *  next processor will have access to the pedestal, noise and
     *  status immediately.  We have a collection for the pedestal, one
     *  for the noise and another for the status.  The first two are
     *  self explaining, while few more words are needed for the status
     *  one.  This one is used, for example, to set a bad pixel mask or
     *  to flag pixel as already belonging to an already reconstructed
     *  cluster.
     */
    std::string _pedestalCollectionName;

    //! Noise collection name.
    /*! See EUTelCalibrateEventProcessor::_pedestalCollectionName for the detailed description
     */
    std::string _noiseCollectionName;

    //! Status collection name.
    /*! See EUTelCalibrateEventProcessor::_pedestalCollectionName for the detailed description
     */
    std::string _statusCollectionName;

    //! Calibrated data collection name.
    /*! The name of the output calibrated data collection. Usually
     *  simply "data"
     */
    std::string _calibratedDataCollectionName;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Filling of the debug histo
    /*! If the user wants to fill some low level detector debug
     *  histograms, she/he has to switch this flag on. An AIDA
     *  interface is also needed.
     */
    bool _fillDebugHisto;


    //! Perform common mode suppression
    /*! The user can perform a common mode suppression step by
     *  switching on this flag.  The common mode suppression is
     *  done on a detector base and it is rejecting possible hit
     *  pixels.
     */
    int _doCommonMode;

    /* Row wise common mode calculation */
    int _CommonModeRowWise;

    //! Hit rejection threshold
    /*! In the case the user wants to suppress the common mode, we
     *  need to remove from the calculation all pixels having a SNR
     *  such that they can be considered hit pixels.
     */
    float _hitRejectionCut;

    //! Maximum number of excluded pixels
    /*! This is the maximum allowed number of pixel exceeding the
     *  common mode cut.  If in a particular event there is a number
     *  of rejected pixels exceeding this value, then the event is
     *  rejected and not used for common mode suppression.  Usually
     *  this number if a not negligible fraction of the whole matrix
     *  (say 30%).  Having some many hit pixels is usually a good
     *  indication of a "crazy" event and it is better to skip it.
     *
     *  In the case the run is characterized by a very high flux, so
     *  this number has to be considerably high to avoid the rejection
     *  of all events. This cut can by-passed setting its value to -1.
     */
    int _maxNoOfRejectedPixels;

    //! Maximum number of excluded pixels per row
    /*! This is the maximum allowed number of pixel exceeding the
     *  common mode cut in a line. Used only when common mode
     *  algorithm is RowWise.
     */
    int _maxNoOfRejectedPixelPerRow;

    //! Maximum number of skipped row
    /*! This is the maximum allowed number of skipped rows per event
    //! in the common mode rejection. Used only when common mode
    //! algorithm is RowWise.
    */
    int _maxNoOfSkippedRow;

    //! The histogram information file
    /*! This string contain the name of the histogram information
     *  file. This is selected by the user in the steering file.
     *
     *  @see eutelescope::EUTelHistogramManager
     *  @see eutelescope::EUTelHistogramInfo
     */
    std::string _histoInfoFileName;

  private:

    //! First pixel along X
    /*! This array of int is used to store the number of the first
     *  pixel along the X direction
     */
    IntVec _minX;

    //! Last pixel along X
    /*! This array of int is used to store the number of the last
     *  pixel along the X direction
     */
    IntVec _maxX;

    //! First pixel along Y
    /*! This array of int is used to store the number of the first
     *  pixel along the Y direction
     */
    IntVec _minY;

    //! Last pixel along Y
    /*! This array of int is used to store the number of the last
     *  pixel along the Y direction
     */
    IntVec _maxY;

    //! Maximum number of consecutive missing events
    /*! This processor only applies to RAW data input collections, but
     *  not to break the generality, it will be active also in the
     *  case of ZS data analysis. In such a case after @a
     *  _maxNoOfConsecutiveMissing events, the warning message is not
     *  issued anymore.
     */
    static const unsigned short _maxNoOfConsecutiveMissing;

    //! Number of consecutive missing
    /*! This is a counter of events with missing input
     *  collection. This value is compared with
     *  _maxNoOfConsecutiveMissing to issue or not the warning
     *  message.
     */
    unsigned short _noOfConsecutiveMissing;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

      //! Name of the raw data histogram
      /*! This histogram contains the raw data value as they come in
      *  without any processing. It is useful only in the debug phase
      *  of the DAQ and for low level detector characterization.
      */
    static std::string _rawDataDistHistoName;


    //! Name of the pedestal corrected signal distribution
    /*! Very similar to the
     *  EUTelCalibrateEventProcessor::_rawDataDistHisto, the histogram
     *  named after this variable is filled with raw data corrected by
     *  the pedestal value. Again it is useful only in the debug phase
     *  of the DAQ and for low level detector characterization.
     */
    static std::string _dataDistHistoName;

    //! Common mode distribution histograms
    static std::string _commonModeDistHistoName;

    //! Skipped pixel histogram
    static std::string _skippedPixelDistHistoName;

    //! Skipped row histogram
    static std::string _skippedRowDistHistoName;

    //! Skipped pixel per row histogram
    static std::string _skippedPixelPerRowDistHistoName;

    //! AIDA histogram map
    /*! The histogram filling procedure may occur in many different
     *  places, while it is usually a good reason to keep the booking
     *  procedure in one place only. To recall an histogram pointer
     *  from a different place in the code they are stored within this
     *  histogram map. The second object of the pair has to a
     *  AIDA::IBaseHistogram since this is the base class of all
     *  different kind of histograms, profiles included. It has also
     *  to be a pointer since, this is a pure virtual class and we
     *  want to use <code>dynamic_cast</code> to convert them back to
     *  their original cast.
     */
    std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;
#endif

    //! Map relating ancillary collection position and sensorID
    /*! The first element is the sensor ID, while the second is the
     *  position of such a sensorID in all the ancillary collections
     *  (noise, pedestal and status).
     */
    std::map< int, int > _ancillaryIndexMap;

    //! Map relating the sensorID and pixels
    std::map< int, unsigned int > _noOfPixelMap;

    //! Map relating the sensorID and pixels per row
    std::map< int, unsigned int > _noOfPixelPerRowMap;

    //! Map relating the sensorID and the rows
    std::map< int, unsigned int > _noOfRowMap;

    //! Geometry ready switch
    /*! This boolean reveals if the geometry has been properly
     *  initialized or not.
     */
    bool _isGeometryReady;

  };

  //! A global instance of the processor
  EUTelCalibrateEventProcessor gEUTelCalibrateEventProcessor;

}
#endif

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELUPDATEPEDESTALNOISEPROCESSOR
#define EUTELUPDATEPEDESTALNOISEPROCESSOR 1

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>



namespace eutelescope {

  //! Update pedestal and noise values
  /*! This processor is used to keep the value of pedestal and noise
   *  update during the beam data analysis. The rationale for this
   *  processor is that, usually a specific "empty" run to calculate
   *  the pedestal and noise value is done at the beginning of the
   *  data taking and after that the acquisition may run for hours
   *  before another pedestal run. In order not to spoil the data
   *  quality because of thermal run away, or excessive degradation of
   *  pixel noise, this processor will take care to update the
   *  pedestal / noise (optionally status) matrix during the data
   *  analysis.
   *
   *  This processor can be called at any time after
   *  EUTelClusteringProcessor, because the information it needs are
   *  contained in the status matrix. All pixels neither BADPIXEL nor
   *  HIT are considered "empty" and then they can be used to update
   *  the initial value of pedestal and noise. Optionally one may also
   *  want to periodically change the status matrix (not yet
   *  implemented).
   *
   *  How often this update procedure is called is set by the user
   *  using the UpdateFrequency parameter. This corresponds to the
   *  number of events in between two updates. Set it to 1 update
   *  every single event.  Moreover the user can also select which
   *  algorithm to use for the update. Here comes a list of available
   *  algorithm:
   *
   *  \li <b>FixedWeight</b>: corresponding to
   *  EUTELESCOPE::FIXEDWEIGHT. In this specific algorithm, the user
   *  has to provide an integer number that is used as a weight; for
   *  simplicity say <code> W = 100 </code>. Given that:
   *  <code>P[i]</code> is the pedestal value after @c i interactions
   *  and that @c D is the current interaction pixel signal, then
   *  <code>P[i+1] = [(W - 1) / W] * P[i] + (1 / W) * D </code>. It
   *  looks like having a on-line average where the number of
   *  elements is kept constant to weight value. An analogous approach
   *  is used for the noise calculation.
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Raw data collection</b> the collection with the full raw data matrix.
   *
   *  <b>Pedestal collection</b> the current pedestal collection
   *
   *  <b>Noise collection </b> the current noise collection
   *
   *  <b>Status collection </b> the current status collection
   *
   *  <h4>Output collections</h4>
   *
   *  <b>Pedestal collection</b> the updated pedestal collection
   *
   *  <b>Noise collection </b> the updated noise collection
   *
   *  @param RawDataCollectionName name of the raw data collection
   *  @param PedestalCollectionName name of the pedestal collection
   *  @param NoiseCollectionName name of the noise collection
   *  @param StatusCollectionName name of the pixel status collection
   *  @param UpdateAlgorithm name of the algorithm to be used
   *  @param UpdateFrequency update frequency in events
   *  @param FixedWeightValue the value of the fixed weight
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */

  class EUTelUpdatePedestalNoiseProcessor : public marlin::Processor {

  public:


    //! Returns a new instance of EUTelUpdatePedestalNoiseProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelUpdatePedestalNoiseProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelUpdatePedestalNoiseProcessor;
    }

    //! Default constructor
    EUTelUpdatePedestalNoiseProcessor ();

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

    //! Fixed weight update algorithm
    /*! This method is called by the processEvent if the _updateAlgo
     *  has been set to EUTELESCOPE::FIXEDWEIGHT. It has been set as
     *  protected because it has to called from within the class
     *  itself and specifically from the processEvent(LCEvent*)
     *  method.
     *
     *  The algorithm can be described as follows:
     *
     *  \li All pixels that have a status equal to
     *  EUTELESCOPE::GOODPIXEL, i. e. they are neither hit nor bad,
     *  are then considered empty or in other words they can be used
     *  for pedestal calculation.
     *
     *  \li Their signal is added to previous pedestal value according
     *  to this formula: <code>P[i+1] = [(W - 1) / W] * P[i] + (1 / W)
     *  * D </code>, where W is the @c _fixedWeightValue and @c D is
     *  the pixel current value.
     *
     *  @param evt The current LCEvent event as passed by the
     *  processEvent
     */
    void fixedWeightUpdate(LCEvent * evt);

    //! Pixel monitoring
    /*! This method is used to collect some information about the
     *  pedestal and noise update. Updating pedestal values is of
     *  utmost importance when analysing long lasting run where a
     *  possible thermal run away may occurr. Anyway this procedure is
     *  not risk-free and bias can be induced by the update
     *  itself. Monitoring how the pedestal and noise values of a
     *  selected number of pixels are changing can be used as figure
     *  of merit to estimate the suitable value of updateFrequency and
     *  the consistency of the update.
     *
     *  This method has been set protected because it has to be called
     *  by the check(LCEvent *) call back. Moreover it is assumed to
     *  be algorithm independent, so that if other algorithm are
     *  implemented the monitoring task can be left untouched.
     *
     *  @param evt The LCEvent as passed by check(LCEvent*)
     *
     *  @throw InvalidParameterException if something bad happens with
     *  pixel index calculation.
     */
    void pixelMonitoring(LCEvent * evt);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! Fill in AIDA DataPointSet for pixel monitoring
    /*! If, and only if Marlin is using AIDA and an AIDAProcessor has
     *  been called before this processor (as it is expected to be)
     *  some dataPointSet are filled and saved into the AIDA output
     *  file.  Those dataPointSet will held the pedestal and noise
     *  evolution versus update iterations (read time).
     *
     *  This method is set protected since it is going to be called
     *  only by the end() call back.
     *
     *  //todo It is probably better to allow the user to save the
     *  monitoring info also if the Marlin is not using AIDA
     *
     *  @bug Because of a bug in RAIDA, this is not even compiling
     *  when using this AIDA implementation
     */
    void saveMonitoring();
#endif

    //! Raw data collection name
    /*! This is the input collection containing the current raw data
     */
    std::string _rawDataCollectionName;

    //! Pedestal collection name
    /*! Input pedestal collection name. Default value is pedestal.
     */
    std::string _pedestalCollectionName;

    //! Noise collection name
    /*! Input noise collection name. Default value is noise.
     */
    std::string _noiseCollectionName;

    //! Status collection name
    /*! Input status collection name. Default value is status.
     */
    std::string _statusCollectionName;

    //! Algorithm name
    /*! The name of the algorithm the user wants to apply. Suitable
     *  const string can be found within the EUTELESCOPE class.
     *
     *  @see EUTELESCOPE::FIXEDWEIGHT
     */
    std::string _updateAlgo;

    //! Monitored pixel list.
    /*! This vector of integer is used by the user in the steering
     *  file to provide a list of pixels to be monitored during the
     *  update process. Each pixel has to be specified by three
     *  integer numbers: the detectordID, the xCoord and the yCoord.
     *  If the size of _monitoredPixel is zero, then the full
     *  monitoring procedure will be switched off.
     */
    IntVec _monitoredPixel;

    //! Pedestal monitor
    /*! This is a vector of float vectors. There is a FloatVec for
     *  each monitored pixels and the size of this FloatVec is equal
     *  to the number of times the update pedestal algorithm has been
     *  called + 1.
     */
    std::vector<FloatVec > _monitoredPixelPedestal;

    //! Noise monitor
    /*! This is a vector of float vectors. There is a FloatVec for
     *  each monitored pixels and the size of this FloatVec is equal
     *  to the number of times the update pedestal algorithm has been
     *  called + 1.
     */
    std::vector<FloatVec > _monitoredPixelNoise;

    //! Update frequency
    /*! This is an integer number representing how often this
     *  algorithm has to be applied. It is set by the user via the
     *  processor parameter and has to be positive. An update
     *  frequency of 0 is meaningless and will be automatically re-set
     *  to 1, i.e. updated every single event.
     */
    int _updateFrequency;


    //! Fixed weight value
    /*! This integer number is used in the FixedWeight algorithm to
     *  update both pedestal and noise values.
     *
     *  @see fixedWeightUpdate(LCEvent*)
     */
    int _fixedWeight;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

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

  private:


  };

  //! A global instance of the processor
  EUTelUpdatePedestalNoiseProcessor gEUTelUpdatePedestalNoiseProcessor;

}
#endif

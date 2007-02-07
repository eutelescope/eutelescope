// -*- C++ -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPEDESTALNOISEPROCESSOR_H
#define EUTELPEDESTALNOISEPROCESSOR_H 1

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerRawDataImpl.h>

// system includes <>
#include <string>
#include <cmath>


namespace eutelescope
{

  /**  Pedestal and noise  processor for marlin.
   * 
   *  <h4>Input - Prerequisites</h4>
   *
   *  <h4>Output</h4> 
   * 
   *  @param RawDataCollectionName Name of the input data collection
   *  @param CalculationAlgorithm Name of the calculation algorithm used
   *  @param NoOfCMIteration Number of common suppression iterations
   *  @param HitRejectionCut Value of the SNR of a pixel to be considered a hit
   *  @param MaxNoOfRejectedPixels Maximum allowed number of rejected pixels per event
   *  @param BadPixelMaskCut Threshold for bad pixel identification
   *  @param FirstEvent First event to be used for pedestal calculation
   *  @param LastEvent Last event to be used for pedestal calculation
   *  @param PedestalCollectionName Name of the output pedestal collection
   *  @param NoiseCollectionName Name of the output noise collection
   *  @param StatusCollectionName Name of the output pixel status collection
   * 
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelPedestalNoiseProcessor.h,v 1.1.1.1 2007-02-07 10:53:12 bulgheroni Exp $ 
   */

  class EUTelPedestalNoiseProcessor:public marlin::Processor
  {

  public:

     
    //! Return a new instance of EUTelPedestalNoiseProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     */
    virtual Processor * newProcessor ()
    {
      return new EUTelPedestalNoiseProcessor;
    }

    //! Default constructor 
    EUTelPedestalNoiseProcessor ();

    /*! Called at the begin of the job before anything is read.  Use
     *  to initialize the processor, e.g.book histograms.
     */
    virtual void init ();

    /** Called for every run.
     */
    virtual void processRunHeader (LCRunHeader * run);

    /** Called for every event - the working horse.
     */
    virtual void processEvent (LCEvent * evt);


    //! Check event method
    /*! This method is called by the Marlin execution framework as
     * soon as the processEvent is over. It can be used to fill check
     * plots.
     */
    virtual void check (LCEvent * evt);

    //! Mask bad pixels
    /*! This method is called at the end of every loop and is used to
     *  identify too noisy pixels and remove them from the following
     *  analysis steps. Bad pixels are flagged as
     *  EUTELESCOPE::BADPIXEL in the
     *  EUTelPedestalNoiseProcessor::_status vector. <br>For the time
     *  being there is only one procedure implemented based on the
     *  noise distribution (<b>"NoiseDistribution"</b>). As a first
     *  thing both the <i>mean</i> and the <i>RMS</i> of the
     *  distribution are calculated. Then the region in the
     *  distribution with pixel having a noise exceeding
     *  _badPixelMaskCut * RMS + mean is considered bad. This
     *  algorithm has the advantage of being pretty robust and system
     *  independent since it does not use any system unit as ADC or
     *  ENC, but, at the same time, as the disadvantage of being
     *  difficult to tune. In fact, while the general user has a
     *  natural idea of what could be the noise of a bad pixel in
     *  her/his own system, it is not easy and clear which should be
     *  the _badPixelMaskCut to be used to obtain the same result.
     *
     *  @todo Ask other teams about their favorite masking algorithm.
     */
    void maskBadPixel();

    //! Calculation done in the first loop
    /*! This method is called within the
     *  EUTelPedestalNoiseProcessor::processEvent(LCEvent *) when the
     *  loop counter is 0. Within this method the first approximations
     *  of both pedestal and noise are calculated on a detector based.
     *  For the time being there is only algorithm ("MeanRMS")
     *  implemented based on the calculation of each pixel mean and
     *  RMS. The algorithm is chosen by the user from the steering
     *  file.
     *
     *  @param evt This is the current event passed by
     *  EUTelPedestalNoieProcessor::processEvent(LCEvent*)
     *
     *  @todo Ask other teams about their favorite calculation
     *  algorithm
     */
    void firstLoop(LCEvent * evt);

    //! Calculation done in all loops other than the first
    /*! This method is called by the
     *  EUTelPedestalNoiseProcessor::processEvent(LCEvent*) when the
     *  loop counter is different from 0 and it is repeated
     *  EUTelPedestalNoiseProcessor::_noOfCMIterations times.  The
     *  main difference from this method and
     *  EUTelPedestalNoiseProcessor::firstLoop(LCEvent*) is in the
     *  common mode suppression algorithm and accidental hit
     *  rejection. This two procedures cannot be applied into the
     *  first loop because they need at least a first approximation of
     *  the pedestal and the noise value of each pixel. For what the
     *  pedestal / noise algorithm is concerned, please have a look at
     *  the description given for the first loop.
     *
     *  @param evt This is the current event passed by the
     *  EUTelPedestalNoiseProcessor::processEvent(LCEvent*).
     * 
     *  @see EUTelPedestalNoiseProcessor::firstLoop(LCEvent*)
     */
    void otherLoop(LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is finished. The
     *  following operations are performed: \li Temporary pedestal and
     *  noise vectors are copied into the final ones.  \li Temporary
     *  vectors are cleared. \li The bad pixel masking procedure is
     *  invoked. \li The loop counter is incremented and crosschecked
     *  if other loop on events are needed or if the full processor is
     *  over. In the case the loop as to be re-executed, then the
     *  return value of the processor is set to false and a
     *  RewindDataFilesException is thrown. In the case, instead, the
     *  number of loops is the final one, then the return value is set
     *  to true and pedestal, noise and status are moved to a
     *  Tracker(Raw)Data object.
     *
     *  @bug RewindDataFilesException is properly thrown but it is not
     *  catched by the Marlin main program. In fact, this exception is
     *  implemented only during the event loop. But to throw
     *  RewindDataFilesException while the last event is processed you
     *  should know in advance how many events are in the file.
     *
     *  @todo For the time being, the three output vectors (_pedestal,
     *  _noise and _status) are not moved to a LCObject and then to
     *  collection.
     *
     *  @bug We can calculate the final values of pedestal and noise
     *  only after the last event has been processed. So the natural
     *  place to do it is within end(), but here we do not have any
     *  LCEvent to attach the collection we want to save.
     *
     */
    virtual void end();


  protected:

    //! Input collection name.
    /*! For the time being we have just one collection that can be used as input 
     */
    std::string _rawDataCollectionName;

    //! Pedestal collection name.
    /*! We have three different output collections that are going to be
     *  saved in the file within the first event possibly so that the
     *  next processor will have access to the pedestal, noise and
     *  status immeditalety.  We have a collection for the pedestal, one
     *  for the noise and another for the status.  The first two are
     *  self exaplaining, while few more words are needed for the status
     *  one.  This one is used, for example, to set a bad pixel mask or
     *  to flag pixel as already belonging to an already reconstructed
     *  cluster.
     */
    std::string _pedestalCollectionName;

    //! Noise collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _noiseCollectionName;

    //! Status collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _statusCollectionName;

    //! Pedestal calculation algorithm
    /*! People around the word are using several different algorithm
     *    to calculate the pedestal and noise values of pixel
     *    detector. To allow everyone to use her/his favorite method,
     *    this same processor can be used with different algorithm
     *    implementations. The user can select among the following
     *    methods already implemented: \li MeanRMS --> The pedestal is
     *    the mean pixel signal --> The noise is the sigma of the
     *    pixel signal distribution \todo Implement all other methods
     *    according to user wishes.
     */
    std::string _pedestalAlgo;

    //! Bad pixel masking algorithm
    /*! The selected algorithm for bad pixel masking. @see
     *  EUTelPedestalNoiseProcessor::maskBadPixel()
     */
    std::string _badPixelAlgo;

    //! Common mode suppression iterations
    /*! This is the number of times the common mode suppression
     *    algorithm has to be applied to the data. Usually one iteration
     *    is enough, but in case of a very noisy setup it can be
     *    repeated several times.  Set it to zero to skip common mode
     *    suppression at all.
     */
    int _noOfCMIterations;

    //! Event loop counter
    /*! This is a counter for the number of loops. The processor will
     * loop the first time (_iLoop == 0) for pedestal and noise
     * calculations and other _noOfCMIterations for common mode
     * suppression and next order correction to the pedestal and noise
     * values.
     */
    int _iLoop;

    //! Hit rejection cut
    /*! A.k.a. common mode cut. During the common mode suppression
     *  algorithm the average of all pixels signals is computed and
     *  then subtracted from their initial value. This is usually done
     *  on specific empty run where no signals from particles are
     *  expected. Whenever a hit occurs, it might harm the suppression
     *  algorithm. For that reason, all pixels with a SNR higher than
     *  this cut are excluded from the common mode calculation.  The
     *  SNR is calculated using the first approximation to the noise
     *  obtained in the first run where all pixels are used.
     */
    float _hitRejectionCut;

    //! Maximum number of excluded pixels
    /*! This is the maximum allowed number of pixel exceeding the
     *  common mode cut.  If in a particular event there is a number
     *  of rejected pixels exceeding this value, then the event is
     *  rejected and not used for common mode suppression.  Usually
     *  this number if a not negligigle fraction of the whole matrix
     *  (say 30%).  Having some many hit pixels is usually a good
     *  indication of a "crazy" event and it is better to skip it.
     */
    int _maxNoOfRejectedPixels;
    
    //! Counter for skipped event due to common mode
    /*! This counter is used to enumerate how many events in the run
     *  have been skipped because of a number of hit pixels exceeding
     *  _maxNoOfRejectedPixels. At the end of the procedure this
     *  number is shown to the user who can decide to keep or to drop
     *  the pedestal noise. As mentioned already, if properly set
     *  _maxNoOfRejectedPixels is a good indicator of an event in
     *  which something went badly wrong; consequently having a large
     *  number of skipped events is a good indicator that the run was
     *  somehow screw up! @see EUTelPedestalNoiseProcessor::_maxNoOfRejectedPixels.
     */
    int _noOfSkippedEvent;

    //! Bad pixel masking cut
    /*! This value is used to mask bad pixels.  After each loop on
     *  events the average (noiseAvg) and the variance (noiseSigma) of
     *  the noise distribution is calculated. To identify bad pixels
     *  the following criteria is adopted.  A pixel is masked as bad
     *  if:
     *    
     * \li its noise is == 0. This is usually a good indication a pixel
     * is dead!  \li its noise is exceeding noiseAvg + _badPixelMaskCut
     * * noiseSigma.
     *
     * \todo In principle also a cut to remove the left tail of the
     * distribution can be introduced, but probably not so important.
     */
    float _badPixelMaskCut;

    //! Number of detector planes in the run
    /*! This is the total number of detector saved into this input
     *  file <br>NOTE: Pedestal, noise and especially common mode
     *  suppression is * done on a detector base. It is retrieved from
     *  a runHeader * parameter.
     */
    int _noOfDetector;

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

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Temporary array to store the running average
    /*! This a vector of vectors of doubles. The outern layer represents
     * the detector, while the inner layer is used to temporary store
     * the running average of pixel output signal.  The outer layer is
     * properly sized in the runHeader(), when the number of detector
     * plane in the run is discovered.  The inner layer is done in the
     * firstEvent of each detector module.
     */
    std::vector < DoubleVec > _tempPede;

    //! Temporary array to store the running sigma
    /*! This a vector of vectors of doubles. The outern layer represents
     * the detector, while the inner layer is used to temporary store
     * the running sigma of pixel output signal.  The outer layer is
     * properly sized in the runHeader(), when the number of detector
     * plane in the run is discovered.  The inner layer is done in the
     * firstEvent of each detector module.
     */
    std::vector < DoubleVec > _tempNoise;

    //! Temporary array to store the number of entries
    /*! In order to have both the running average and sigma working,
     * another array to store how many "entries" a pixel has should be
     * recorded. This might be useless there is an entry for each event,
     * but this is not the case when applying "hit rejection" in common
     * mode suppression algorithm. The underlying phylosophy is the same
     * used for _tempPede and _tempNoise
     */
    std::vector < IntVec > _tempEntries;

    //! Array to store the intermediate/final pedestal value
    /*! At the end of the first loop on events, a first approximation
     *  of pedestal and noise values is available. This value is then
     *  stored into this array to allow further improvements in the
     *  following event loops
     */
    std::vector < DoubleVec > _pedestal;

    //! Array to store intermediate/final noise value
    /*! As for EUTelPedestalNoiseProcessor::_pedestal but now for the
     *  noise
     */
    std::vector < DoubleVec > _noise;

    //! Array to store the status 
    /*! This is the place where the status of pixel is saved
     */
    std::vector < IntVec > _status;
    
    //! First event for pedestal calculation
    /*! This is the first event to be used for pedestal
     *  calculation. If the input file is not a specific pedestal run
     *  and only some events are really empty, then the user has to
     *  select this range.
     */
    int _firstEvent;

    //! Last event for pedestal calculation
    /*! This is the last event to be used for pedestal calculation. If
     *  the imput file is not a specific pedestal run and only some
     *  events are really empty, then the user has to select this
     *  range. If the user set this value to -1 or to an event number
     *  in excess to the total number of events, then all events are
     *  used for pedestal calculation.
     */
    int _lastEvent;

    //! Boolean flag for pedestal calculation
    /*! This boolean flag is used to modify the
     *  EUTelPedestalNoiseProcessor::processEvent(LCEvent *)
     *  behaviour. This boolean is set to true during the first run
     *  over the selected events
     */
    bool _doPedestal;

  };
}
#endif

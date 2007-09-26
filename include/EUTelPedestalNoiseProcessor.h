// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
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

// AIDA includes <.h>
#ifdef MARLIN_USE_AIDA
#include <AIDA/IBaseHistogram.h>
#endif


// system includes <>
#include <string>
#include <cmath>



namespace eutelescope {

  //! Pedestal and noise  processor for Marlin.
  /*! This processor has the task to calculate the pedestal and noise
   *  value of each single pixel in the EUDET telescope. The input
   *  data, as they come from the DAQ, or converted from another
   *  format, are organized in a collection of TrackerRawData named
   *  rawdata. This collection has as many elements as the number of
   *  sensors in the telescope. This number may or may not correspond
   *  to the number of sensitive planes into the telescope geometry,
   *  because a plane can in principle be made by more than one
   *  sensor. This distinction between detectors and planes has been
   *  introduced because several algorithm like common mode
   *  suppression and bad pixel masking should be performed on a
   *  detector and not on a plane basis.  Moreover, in the case of a
   *  sensor like the MimoSTAR2 that is divided into two channels
   *  featuring different pixel characteristics, the sensor can be in
   *  principle be considered as two different detectors.  Each
   *  Tracker(Raw)Data has cellID0/1 set using the CellIDEncoder
   *  utility with EUTELESCOPE::DEFAULTMATRIXENCODING. This encoding
   *  offers the possibility to store the "sensorID" (5 bits), the
   *  "xMin", "xMax", "yMin", "yMax" with 12 bit each. The xMin and
   *  yMin are different from 0 only in the case there are more than
   *  one sensor per plane or the detector is actually a channel of a
   *  bigger detector.
   *  
   *  The user can choose which algorithm should be used for pedestal
   *  calculation. @see EUTelPedestalNoiseProcessor::_pedestalAlgo
   *  for a detailed description of the available methods.
   *
   *  The user can select the bad pixel masking algorithm, changing
   *  the value of _badPixelAlgo. @see
   *  EUTelPedestalNoiseProcessor::_badPixelAlgo
   *
   *  <h4>Input collection</h4>
   *  <b>Raw data</b> A collection of TrackerRawData 
   *
   *  <h4>Output file</h4> 
   *  The output of this processor will be saved into a separate file
   *  and can be reloaded using a condition processor. 
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
   *  @param OutputPedeFile Name of the output pedestal file
   *  @param ASCIIOutputSwitch To enable/disable the generation of
   *  ASCII output files
   *
   *  <h4>Typical steering file for pedestal production</h4> The
   *  following code can be used as a steering file for pedestal
   *  production.  
   *  
   *  @code 
   *  <marlin>
   *   <global> 
   *    <parameter name="LCIOInputFiles"> pedestal-run.slcio  </parameter> 
   *    <parameter name="GearXMLFile" value="gear-telescope.xml"/>
   *    <parameter name="MaxRecordNumber" value="5001"/> 
   *    <parameter name="SkipNEvents" value="0 "/> 
   *    <parameter name="SupressCheck" value="false"/> 
   *    <parameter name="Verbosity" value="MESSAGE"/>
   *   </global>
   * 
   *   <execute>
   *    <processor name="AIDAHistogrammingInterface"/>
   *    <processor name="PedestalAndNoiseCalculator"/>
   *   </execute>
   *
   *   <processor name="AIDAHistogrammingInterface" type="AIDAProcessor">
   *    <!--
   *      Processor that handles AIDA files. Creates on directory per
   *      processor. Processors only need to create and fill the
   *      histograms, clouds and tuples. Needs to be the first
   *      ActiveProcessor
   *    -->
   *    <!-- compression of output file 0: false >0: true (default) -->
   *    <parameter name="Compress" type="int" value="1"/>
   *    <!-- filename without extension-->
   *    <parameter name="FileName" type="string" value="pedestal-histo"/>
   *    <!-- type of output file xml (default) or root ( only OpenScientist)-->
   *    <parameter name="FileType" type="string" value="root"/>
   *   </processor>
   *  
   *   <processor name="PedestalAndNoiseCalculator" type="EUTelPedestalNoiseProcessor">
   *    <!--EUTelPedestalNoiseProcessor computes the pedestal and noise values of a pixel detector-->
   *    <!--Input raw data collection-->
   *    <parameter name="RawDataCollectionName" type="string" lcioInType="TrackerRawData"> rawdata </parameter>
   *    <!--Threshold for bad pixel identification-->
   *    <parameter name="BadPixelMaskCut" type="float" value="3.5"/>
   *    <!--Select the algorithm for bad pixel masking-->
   *    <parameter name="BadPixelMaskingAlgorithm" type="string" value="NoiseDistribution"/>
   *    <!--Select the algorithm for pede/noise calculation-->
   *    <parameter name="CalculationAlgorithm" type="string" value="MeanRMS"/>
   *    <!--First event for pedestal calculation-->
   *    <parameter name="FirstEvent" type="int" value="0"/>
   *    <!--Threshold for rejection of hit pixel (SNR units)-->
   *    <parameter name="HitRejectionCut" type="float" value="4"/>
   *    <!--Last event for pedestal calculation-->
   *    <parameter name="LastEvent" type="int" value="-1"/>
   *    <!--Maximum allowed number of rejected pixels per event-->
   *    <parameter name="MaxNoOfRejectedPixels" type="int" value="1000"/>
   *    <!--Number of common mode suppression iterations-->
   *    <parameter name="NoOfCMIteration" type="int" value="1"/>
   *    <!--Noise collection name-->
   *    <!--parameter name="NoiseCollectionName" type="string" value="noiseDB"/-->
   *    <!--The filename (w/o .slcio) to store the pedestal file-->
   *    <parameter name="OutputPedeFile" type="string" value="pedestal-db"/>
   *    <!--Pedestal collection name-->
   *    <!--parameter name="PedestalCollectionName" type="string" value="pedestalDB"/-->
   *    <!--Status collection name-->
   *    <!--parameter name="StatusCollectionName" type="string" value="statusDB"/-->
   *   </processor>
   *  </marlin>
   *  @endcode
   *
   *  Note that you don't need a LCIOOutputProcessor at the end since
   *  it is the EUTelPedestalNoiseProcessor itself taking care of
   *  saving the output pedestal file.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelPedestalNoiseProcessor.h,v 1.19 2007-09-26 15:15:52 bulgheroni Exp $ 
   *
   *  @todo For the time being the final pedestal/noise/status objects
   *  are stored into a LCIO and they will be successively accessed by
   *  a ConditionProcess; consider the possibility to store those
   *  objects also into a DB.
   *
   */

  class EUTelPedestalNoiseProcessor:public marlin::Processor
  {

  public:

     
    //! Returns a new instance of EUTelPedestalNoiseProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelPedestalNoiseProcess.
     */
    virtual Processor * newProcessor ()
    {
      return new EUTelPedestalNoiseProcessor;
    }

    //! Default constructor 
    EUTelPedestalNoiseProcessor ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and reset all needed data
     *  members. In principle this could also be the right place to
     *  book histograms, but since those are also grouped according to
     *  the detector numbers we need to have this parameter available.
     */
    virtual void init ();

    //! Called for every run.
    /*! At the beginning of every run the run header is read and
     *  processed by this method. As a first thing, the input run
     *  header is dynamically re-casted to a EUTelRunHeaderImpl and
     *  then important things like the number of detectors and the
     *  pixel detector boundaries are dumped from the file. After that
     *  the EUTelPedestalNoiseProcess::bookHistos() is called.
     *  
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! Since the behavior of the PedestalNoise processor is different
     *  if this is the first or one of the following loop, this method
     *  is just calling
     *  EUTelPedestalNoiseProcessor::firstLoop(LCEvent*) or
     *  EUTelPedestalNoiseProcessor::otherLoop(LCEvent*)
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
     *
     */
    virtual void check (LCEvent * evt);

    //! Mask bad pixels
    /*! This method is called at the end of every loop and is used to
     *  identify too noisy pixels and remove them from the following
     *  analysis steps. Bad pixels are flagged as
     *  EUTELESCOPE::BADPIXEL in the
     *  EUTelPedestalNoiseProcessor::_status vector. Here comes a list
     *  of all implemented algorithms:
     *
     *  \li <b>"NoiseDistribution"</b>. As a first thing both the
     *  <i>mean</i> and the <i>RMS</i> of the distribution are
     *  calculated. Then the region in the distribution with pixel
     *  having a noise exceeding <code>_badPixelMaskCut * RMS +
     *  mean</code> is considered bad. This algorithm has the
     *  advantage of being pretty robust and system independent since
     *  it does not use any system unit as ADC or ENC, but, at the
     *  same time, as the disadvantage of being difficult to tune. In
     *  fact, while the general user has a natural idea of what could
     *  be the noise of a bad pixel in her/his own system, it is not
     *  easy and clear which should be the _badPixelMaskCut to be used
     *  to obtain the same result.
     *  
     *  \li <b>"AbsoluteNoiseValue"</b>. This second algorithm, even
     *  if it is less general than the "NoiseDistribution" one it
     *  offers to the user a higher degree of sensibility. With this
     *  algorithm, the user has to provide as a _badPixelMaskCut the
     *  maximum acceptable value of noise in ADC units. At the
     *  pedestal calculation level, it is still impossible to provide
     *  any reasonable "calibration" able to provide a conversion
     *  factor from ADC units to equivalent noise charge.
     *
     *  @todo Ask other teams about their favorite masking algorithm.
     */
    void maskBadPixel();

    //! Calculation done in the first loop
    /*! This method is called within the
     *  EUTelPedestalNoiseProcessor::processEvent(LCEvent *) when the
     *  loop counter is 0. Within this method the first approximations
     *  of both pedestal and noise are calculated on a detector based.
     *  Here comes a list of all available processing algorithms:
     *
     *  \li <b>MeamRMS</b>. This method can be selected using the
     *  static const char EUTELESCOPE::MEANRMS, or typing MeanRMS in
     *  the steering file. This algorithm computes the pedestal of a
     *  pixel as its average value over the full range of selected
     *  events. The RMS of such signal distribution is used as the
     *  pixel noise.
     *
     *  \li <b>AIDAProfile</b>. This method can be selected using the
     *  static const char EUTELESCOPE::AIDAPROFILE, or typing
     *  AIDAProfile in the steering file. To use this algorithm, Marlin
     *  and consequently Eutelescope, have to be linked against an
     *  AIDA implementation. In the case MARLIN_USE_AIDA is undefined
     *  and the user selects this algorithm, then the program will
     *  automatically fall back to MeanRMS. This algorithm is based on
     *  the use of a AIDA::IProfile2D. Such an object can be booked
     *  having every pixel centered into a single bin and can be
     *  filled with all entries of the TrackerRawData. At the end of
     *  the event loop, each bin contains the pixel mean signal and
     *  each bin error is the noise.
     *
     *  @param evt This is the current event passed by
     *  EUTelPedestalNoieProcessor::processEvent(LCEvent*)
     *
     *  @todo A full debug of the AIDAProfile algorithm is required.
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
     *  @see EUTelPedestalNoiseProcessor::firstLoop(LCEvent*) for the
     *  description of available algorithms.
     */
    void otherLoop(LCEvent * evt);

    //! Additional loop for bad pixel masking
    /*! This method is called within the processEvent(LCEvent * evt)
     *  if the user wants to perform an additional (very fast) loop
     *  over the event looking for pixels firing too often. This has
     *  to be the last loop and both pedestal and noise are not
     *  changed, only the status is updated. 
     * 
     *  A histogram is filled with the firing frequency distribution
     *  for each detector.
     *
     *  @param evt The LCEvent containing all the input collections.
     */
    void additionalMaskingLoop(LCEvent * evt);

    //! Book histograms
    /*! This method is used to prepare the needed directory structure
     *  within the current ITree folder and books all required
     *  histograms. Histogram pointers are stored into
     *  EUTelPedestalNoiseProcess::_aidaHistoMap so that they can be
     *  recalled and filled from anywhere in the code.  Apart from the
     *  histograms listed in EUTelPedestalNoiseProcessor::fillHistos()
     *  there is also a common mode histo described here below:
     *
     *  \li commonModeHisto: 1D histogram to store the calculated
     *  common mode value for each event. This histogram is booked and
     *  filled only if the loop counter is greater-equal than 1,
     *  because for _iLoop == 0 there is no common mode suppression.
     *  This histo is not filled with the other because it needs to be
     *  updated every event.
     * 
     *  @see EUTelPedestalNoiseProcessor::fillHistos() for the todos
     */
    void bookHistos();

    //! Fill histograms
    /*! This method is used to fill in some control histograms. It is
     *  called from end() so it has available the final
     *  vectors. Results are grouped in order to have different loops
     *  and different detectors separated. Here comes a list of filled
     *  histograms:
     *
     *  \li pedestalDistribution: 1D histogram with all good pixel
     *  pedestal values.
     *
     *  \li noiseDistribution: 1D histogram with all good pixel noise
     *  values.
     *
     *  \li pedestalMap: 2D histogram binned as the detector number of
     *  pixels representing the geometrical distribution of pedestal
     *  values
     *
     *  \li noiseMap: as above but for the noise values.
     *
     *  \li statusMap: as above but for the pixel status. Remember that
     *  <code>GOODPIXEL == 0</code> and <code>BADPIXEL == 1</code>.
     *
     */
    void fillHistos();
    

    //! Called after data processing.
    /*! This method is called when the loop on events is finished. It
     *  is checking whether the calculation is properly finished or
     *  not. 
     *  A very common error occurs when the file finished without a
     *  EORE or when the MaxRecordNumber was set to low to loop over
     *  all the needed record. To check this is very easy because we
     *  just have to crosscheck if _iLoop is equal to noOfCMIterations.
     */
    virtual void end();

    //! True if this the current one is the last event in the files
    /*! This utility is used to know if the current events is the last
     *  one. This is currently used in this process to allow data
     *  rewind and multiple loops on events.
     *  
     *  @return true if this is the last event
     */
    inline virtual bool isLastEvent() const { return (_iEvt == _lastEvent); }

    //! Finishes up the current loop 
    /*! This method is called at the end of the event loop, just
     *  before falling into the end() callback. Here we can perform
     *  all the needed operations to finish up the current loop and,
     *  if needed, restart the loop once more. It the case this was
     *  the last loop, just save the output results in an external
     *  file. Here comes a more detailed description:
     *  
     *  \li In a way that depends on the calculation algorithm,
     *  pedestal and noise vectors are filled with the current loop
     *  calculation results. 
     *
     *  \li The bad pixel masking procedure is invoked. 
     *
     *  \li Check histograms are filled.
     *
     *  \li The loop counter is incremented and crosschecked if other
     *  loop on events are needed or if the full processor is over. In
     *  the case the loop as to be re-executed, then the return value
     *  of the processor is set to false and a
     *  RewindDataFilesException is thrown. In the case, instead, the
     *  number of loops is the final one, then the return value is set
     *  to true and pedestal, noise and status are moved to a
     *  Tracker(Raw)Data object.
     *
     *  @param fromMaskingLoop Boolean true when the finalizeProcessor
     *  is called from the additional masking loop.
     *
     *  @see EUTelPedestalNoiseProcessor::fillHistos() for a detailed
     *  description on histogram filling.
     *
     *  @throw RewindDataException to restart the event loop 
     *
     *  @throw StopProcessingException to stop the looping in the
     *  final number of looping is reached.
     *
     */
    virtual void finalizeProcessor(bool fromMaskingLoop = false);


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
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _noiseCollectionName;

    //! Status collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _statusCollectionName;

    //! Pedestal calculation algorithm
    /*! People around the word are using several different algorithms
     *  to calculate the pedestal and noise values of pixel
     *  detector. To allow everyone to use her/his favorite method,
     *  this same processor can be used with different algorithm
     *  implementations. The user can select among the following
     *  methods already implemented: 
     * 
     *  \li <b>MeanRMS</b>. This method is based on the on line
     *  calculation of the pedestal and noise value as the mean and
     *  the RMS of each pixel signal distribution. The actual values
     *  of pedestal and noise are calculated.
     *  
     *  \li <b>AIDAProfile</b>. This algorithm is the easiest one from
     *  the coding point of view since it relies on algorithm already
     *  coded into the used AIDA implementation. The idea behind is
     *  that all pixel output signals are used to fill in an AIDA
     *  profile. At the end of the event loop, each bin of the profile
     *  will be characterized by pixel mean value (pedestal) and
     *  spread (noise). The use of this algorithm requires that the
     *  MARLIN_USE_AIDA is set to 1 and that the AIDAProcessor is
     *  properly called before this processor. In the case one or both
     *  of these conditions are not fullfil the choice fall back on
     *  the very safe and always possible MeanRMS algorithm simply
     *  alerting the user of the change.  

     *  @bug All debug tests have been done using RAIDA as AIDA
     *  implementation and due to a bug in the
     *  AIDAIProfile2D::binRMS() method, the calculated noise is not
     *  correct. Moreover, this algorithm looks particularly slow. <br>
     *
     *  @todo Implement all other methods according to user wishes. In
     *  particular, we can imagine to implement the ROOTProfile as
     *  soon as the ROOTProcessor will be ready and fully tested. 
     */
    std::string _pedestalAlgo;

    //! Bad pixel masking algorithm
    /*! The selected algorithm for bad pixel masking. @see
     *  EUTelPedestalNoiseProcessor::maskBadPixel()
     */
    std::string _badPixelAlgo;

    //! Bad pixel masking cut
    /*! This value is used to mask bad pixels.  
     * 
     *  @see EUTelPedestalNoiseProcess::maskBadPixel() for a detailed
     *  description about the implemented algorithm and the use of
     *  _badPixelMaskCut in all cases.
     */
    float _badPixelMaskCut;

    //! Common mode suppression iterations
    /*! This is the number of times the common mode suppression
     *  algorithm has to be applied to the data. Usually one iteration
     *  is enough, but in case of a very noisy setup it can be
     *  repeated several times.  Set it to zero to skip common mode
     *  suppression at all.
     */
    int _noOfCMIterations;

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
     *  this number if a not negligible fraction of the whole matrix
     *  (say 30%).  Having some many hit pixels is usually a good
     *  indication of a "crazy" event and it is better to skip it.
     */
    int _maxNoOfRejectedPixels;
    
    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;
    
    //! First event for pedestal calculation
    /*! This is the first event to be used for pedestal
     *  calculation. If the input file is not a specific pedestal run
     *  and only some events are really empty, then the user has to
     *  select this range.
     */
    int _firstEvent;

    //! Last event for pedestal calculation
    /*! This is the last event to be used for pedestal calculation. If
     *  the input file is not a specific pedestal run and only some
     *  events are really empty, then the user has to select this
     *  range. If the user set this value to -1 or to an event number
     *  in excess to the total number of events, then all events are
     *  used for pedestal calculation.
     */
    int _lastEvent;

    //! Output pedestal/noise/status file name
    /*! At the end of this processor a LCIO file containing the
     *  results will be saved with this name
     */
    std::string _outputPedeFileName;

    //! Boolean flag for pedestal calculation
    /*! This boolean flag is used to modify the
     *  EUTelPedestalNoiseProcessor::processEvent(LCEvent *)
     *  behavior. This boolean is set to true during the first run
     *  over the selected events
     */
    bool _doPedestal;

    //! Boolean switch for ASCII output file
    /*! This switch is set in the steering file to enable the ASCII
     *  output file for pedestal, noise and status.
     *  These files are particularly useful when loading initial
     *  configuration for zero suppression operation mode.
     */
    bool _asciiOutputSwitch;

  private:

    //! Detector name
    /*! This string is used to copy the detector name from the run
     *  header to the event "header"
     */ 
    std::string _detectorName;

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

    //! Temporary array to store the running average
    /*! This a vector of vectors of doubles. The outern layer represents
     * the detector, while the inner layer is used to temporary store
     * the running average of pixel output signal.  The outer layer is
     * properly sized in the runHeader(), when the number of detector
     * plane in the run is discovered.  The inner layer is done in the
     * firstEvent of each detector module.
     */
    std::vector < FloatVec > _tempPede;

    //! Temporary array to store the running sigma
    /*! This a vector of vectors of doubles. The outern layer represents
     * the detector, while the inner layer is used to temporary store
     * the running sigma of pixel output signal.  The outer layer is
     * properly sized in the runHeader(), when the number of detector
     * plane in the run is discovered.  The inner layer is done in the
     * firstEvent of each detector module.
     */
    std::vector < FloatVec > _tempNoise;

    //! Temporary array to store the number of entries
    /*! In order to have both the running average and sigma working,
     * another array to store how many "entries" a pixel has should be
     * recorded. This might be useless there is an entry for each event,
     * but this is not the case when applying "hit rejection" in common
     * mode suppression algorithm. The underlying philosophy is the same
     * used for _tempPede and _tempNoise
     */
    std::vector < IntVec > _tempEntries;

    //! Array to store the intermediate/final pedestal value
    /*! At the end of the first loop on events, a first approximation
     *  of pedestal and noise values is available. This value is then
     *  stored into this array to allow further improvements in the
     *  following event loops
     */
    std::vector < FloatVec > _pedestal;

    //! Array to store intermediate/final noise value
    /*! As for EUTelPedestalNoiseProcessor::_pedestal but now for the
     *  noise
     */
    std::vector < FloatVec > _noise;

    //! Array to store the status 
    /*! This is the place where the status of pixel is saved
     */
    std::vector < ShortVec > _status;

    //! Additional bad pixel masking vector
    std::vector < ShortVec > _hitCounter;

    //! Event loop counter
    /*! This is a counter for the number of loops. The processor will
     * loop the first time (_iLoop == 0) for pedestal and noise
     * calculations and other _noOfCMIterations for common mode
     * suppression and next order correction to the pedestal and noise
     * values.
     */
    int _iLoop;

#ifdef MARLIN_USE_AIDA
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

    //! Name of the temporary AIDA 2D profile
    /*! The histogram pointed by this name is used in the case
     *  EUTELESCOPE::AIDAPROFILE pedestal calculation algorithm is
     *  selected. In fact, in such a case, a IProfile2D for each
     *  detector is booked in such a way that there each pixel is
     *  centered into a bin. In the end() call back the content of
     *  these profiles are moved back into the
     *  EUTelPedestalNoiseProcess::_pedestal and
     *  EUTelPedestalNoiseProcessor::_noise standard vectors. The
     *  profiles are cleared and ready to be re-used into an eventual
     *  following loop.
     */
    static std::string _tempProfile2DName;

    //! Pedestal distribution 1D histo base name
    /*! @see EUTelPedestalNoiseProcessor::_noiseHistoName;
     */
    static std::string _pedeDistHistoName;

    //! Noise distribution 1D histo base name
    /*! All histograms names defined therein should be actually
     *  considered as base names, in the meaning they are going to be
     *  padded with the loop and detector numbers. For example for the
     *  noise distribution histo on detector number 5 during loop 1 it
     *  is going to be <i>noiseHisto-d5-l1<i>
     */
    static std::string _noiseDistHistoName;

    //! Common mode distribution 1D histo base name
    /*! @see EUTelPedestalNoiseProcessor::_noiseHistoName;
     */
    static std::string _commonModeHistoName;

    //! Pedestal 2D map
    /*! @see EUTelPedestalNoiseProcessor::_noiseHistoName;
     */
    static std::string _pedeMapHistoName;

    //! Noise 2D map
    /*! @see EUTelPedestalNoiseProcessor::_noiseHistoName;
     */
    static std::string _noiseMapHistoName;

    //! Pixel status 2D map
    /*! @see EUTelPedestalNoiseProcessor::_noiseHistoName;
     */
    static std::string _statusMapHistoName;


    //! Firing frequency histogram
    /*! This histogram contains the cumulative distribution of firing
     *  frequency. It is filled only if the additional masking loop is
     *  done.
     */ 
    static std::string _fireFreqHistoName;


#endif 

    //! Histogram switch
    /*! Useful flag to switch on or off the histogramming
     *  operations. In case of problem with booking or filling this
     *  switch will be turned off even if marlin is using AIDA and the
     *  AIDAProcessor is available.
     *
     */
    bool _histogramSwitch;

    //! Additional bad masking loop
    bool _additionalMaskingLoop;

  };

  //! A global instance of the processor
  EUTelPedestalNoiseProcessor gEUTelPedestalNoiseProcessor;      

}
#endif

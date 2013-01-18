/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCALCULATEETAPROCESSOR_H
#define EUTELCALCULATEETAPROCESSOR_H 1

// eutelescope includes ".h"
#include "EUTelPseudo1DHistogram.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// lcio includes <.h>

// system includes <>
#include <string>
#include <map>
#include <set>


#undef MARLIN_USE_HISTOGRAM
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#define MARLIN_USE_HISTOGRAM
#endif
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#define MARLIN_USE_HISTOGRAM
#endif


namespace eutelescope {

  //! Eta function calculator for the EUTelescope
  /*! This processor is used to calculate the eta function of a set of
   *  detectors inside the telescope setup.
   *
   *  This function is used to give the value of the cluster center
   *  correcting the value obtained with the charge center of mass
   *  technique taking into account the effective shape of the
   *  collecting elements on the sensors.
   *
   *  The idea behind this non linear correction, is that in a normal
   *  experiment, the cluster center has to have a flat probability to
   *  occurr over the full pixel area.
   *
   *  As a general assumption we are considering the two geometric
   *  directions (x and y) as independent, so that the problem can be
   *  easily factorized along x and along y. Here comes a brief
   *  description of the steps for the eta function calculation:
   *
   *  \li For each cluster, the charge center of mass is calculated;
   *
   *  \li The (signed) distance between the seed pixel and the CoG is
   *  used to fill in a histogram (pseudo-histogram is used to avoid
   *  any dependency on a specific histogramming package). This
   *  histogram is binned from minum half pitch to plus half picth.
   *
   *  \li When all clusters have been scanned, we can obtain the eta
   *  function as the integral of the produced histogram
   *
   *  \li This eta function has to be normalized and shifted by half a
   *  pitch
   *
   *  <h4>Input collections </h4>
   *  <b>ClusterCollection</b>. A collection of clusters (TrackerPulse)
   *
   *  <h4>Output colelctions</h4>
   *  None
   *
   *  <h4>Output file</h4>
   *  <b>OutputEtaFileName</b>. This is the name of the output LCIO
   *  file containing the eta calculation results. This file can be
   *  reloaded into Marlin via a condition processor.
   *
   *  @param ClusterCollectionName The name of the input cluster
   *  collection.
   *
   *  @param EventNumber The number of events to use for the
   *  calculation
   *
   *  @param NumberOfBins  This vector contains the number
   *  of bins along x and y respectively to be used while binning the
   *  two eta function
   *
   *  @param ClusterQualitySelection This parameter is used to select
   *  the quality of clusters used for Eta calculation. Reccomended is
   *  0 == kGoodCluster. To apply other, more elaborated selection
   *  criteria see eutelescope::EUTelClusterFilter.
   *
   *  @param ClusterTypeSelection This string is
   *  used to choose how the charge center of gravity is
   *  calculated. Type "FULL" to use the full cluster, whatever number
   *  of pixels or shape it has; type "NxMPixel" to use only pixels of
   *  the cluster belonging to a NxM region around the seed; type
   *  NPixel to use only the N pixels of the cluster having the
   *  highest signal.
   *
   *  @param NxMPixelClusterSize This parameter is used only when
   *  ClusterTypeSelection=="NxMPixel" and represents the size along x
   *  and y of the cluster region of interest. In case this is bigger
   *  or equal to the full cluster, then "FULL" is
   *  used.
   *
   *  @param NPixelSize This parameter is used only when
   *  ClusterTypeSelection=="NPixel" and represents the number of
   *  pixels with the highest signal to be used in the CoG
   *  calculation. In case this number is exceeding the cluster
   *  multiplicity then "FULL" is used.
   *
   *  @param EtaXCollectionName This is the name of the collection of
   *  eta functions along the X axis. This collection is saved in the
   *  output file and not made available to the current event.
   *
   *  @param EtaYCollectionName This is the name of the collection of
   *  eta functions along the Y axis. This collection is saved in the
   *  output file and not made available to the current event.
   *
   *  @param OutputEtaFileName The name of the output file.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   *
   *
   */

  class EUTelCalculateEtaProcessor:public marlin::Processor {

  public:


    //! Returns a new instance of EUTelCalculateEtaProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelCalculateEtaProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelCalculateEtaProcessor;
    }

    //! Default constructor
    EUTelCalculateEtaProcessor ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters, resets event / run counters and
     *  clears the std vectors
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. From the run header the number of detector is
     *  retrieved and used to prepare the pseudo histogram vectors.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. Pseudohistograms
     *  are filled
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     *
     *  @throw UnknownDataTypeException if the cluster type stored in
     *  the TrackerPulse is unknown.
     */
    virtual void processEvent (LCEvent * evt);

    //! Finish the eta calculation
    /*! To calculate the eta function, a certain number of events has
     *  to accumulated, and only when this number
     *  (EUTelCalculateEtaProcessor::_nEvent) is reached the
     *  calculation can proceed. This method is making the real
     *  calculation and it is called within processEvent(LCEvent*)
     *  when the isLastEvent() is true.
     *
     *  @throw RewindDataFilesException when the calculation is over,
     *  so that the user can now re-loop on events and applying the
     *  calculated eta function
     */
    virtual void finalizeProcessor();


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
    /*! For the time being this is just printing a good bye message.
     */
    virtual void end();

    //! Identifier of the last event to be used.
    /*! This method can be profitably used to identify which is the
     *  last event to be used for the eta calculation.
     */
    inline bool isLastEvent() const { return ( _nEvent == _iEvt ) ; }

  protected:

    //! Input collection name.
    /*! For the time being we have just one collection that can be
     *  used as input
     */
    std::string _clusterCollectionName;

    //! Number of events for the calculation
    /*! This processor can be executed also on a fraction of the total
     *  number of events, say N events. After that N events have been
     *  processed, the eta funcion can be calculated, events are
     *  rewind to the beginning and this processor won't do anything
     *  more.
     *
     *  To allow the user to make conditional steering files, the
     *  processor return value "isEtaCalculationFinished" will be
     *  false until N events, and then it will become true.
     *
     *  A priori the user may don't know how many events are in the
     *  input files. So, if she/he wants to run the calculation over
     *  the full event range, this parameter has to set to -1. If
     *  _nEvent is exceeding the total number of events in the file or
     *  the MaxRecordParameter in the global section, then the
     *  calculation will stop just before the last available event.
     *
     */
    int _nEvent;

    //! Number of bins along x and y directions
    /*! These two integer numbers (stored into a IntVec) are playing a
     *  crucial role for Eta function calculation. They are used to
     *  divide the seed pixel into horizontal and vertical slices, and
     *  the number of times the charge center of mass is found in a
     *  certain slide is counted.
     *
     *  This number cannot be to low, because otherwise the Eta
     *  function will loose sensitivity and cannot be too high as well
     *  since we would like to have a statistically significant number
     *  of entries for each bin.
     *
     *  In the case of not squared pixels, the binning along the two
     *  directions can be different.
     */
    std::vector<int >  _noOfBin;

    //! Cluster quality selection
    /*! This parameter (int) is set by the user via the steering file. For
     *  the eta function calculation, it is a good attitude to use
     *  only complete and not merging clustering. This integer value
     *  has to be converted to ClusterQuality to be compared with the
     *  value stored into the cluster object.
     */
    int _clusterQuality;

    //! Cluster type selection
    /*! This parameter (string) is set by the user via the steering
     *  file. Independently of the cluster type, one may want to
     *  calculate the eta function not using all the pixels in the
     *  cluster; for example, one can decide to use just a smalled n x
     *  m sub-cluster, or only the first n pixels when the cluster is
     *  sorted according to pixel signal.
     *
     *  Here comes a summary of available option:<br> <b>Full</b>: all
     *  pixels in the cluster are used.<br> <b>NxMPixel</b>: only
     *  pixels belonging to a N x M sub-cluster centered around the
     *  seed.<br> <b>NPixel</b>: only the first N pixel with the
     *  highest signal.
     */
    std::string _clusterTypeSelection;

    //! Cluster max size along x and y
    /*! In the case the user selected the NxMCluster, this vector of
     *  integer is used to store the maximum size along x and y of the
     *  sub-cluster.
     */
    std::vector<int > _xyCluSize;

    //! Number of pixel in the cluster
    /*! In the case the user selected the NPixel cluster type
     *  selection, this parameter represents the number of pixels in
     *  the cluster.
     */
    int _nPixel;

    //! Eta X output collection name
    /*! This is the name of the output collection containing all the
     *  eta functions along x
     */
    std::string _etaXCollectionName;

    //! Eta Y output collection name
    /*! This is the name of the output collection containing all the
     *  eta functions along y
     */
    std::string _etaYCollectionName;

    //! The eta condition output file name
    /*! The eta collections are saved within a collection
     */
    std::string _outputEtaFileName;

  private:

    //! Boolean return value
    /*! This boolean is used as return value for conditional steering
     *  file. Eta function calculation requires to loop over a certain
     *  minimum number of events (_nEvent) before having a resonable
     *  result and to allow other processor to be executed. Until
     *  _nEvent, this boolean return value will be false, and will
     *  become true afterward.
     */
    bool _isEtaCalculationFinished;

    //! Pseudo histograms for CoG along x
    /*! This is an associative collection of pseudo histogram
     *  pointer. The key value is the sensor ID.
     */
    std::map< int, EUTelPseudo1DHistogram * > _cogHistogramX;


    //! Pseudo histograms for CoG along y
    /*! This is an associative collection of pseudo histogram
     *  pointer. The key value is the sensor ID.
     */
    std::map< int, EUTelPseudo1DHistogram * > _cogHistogramY;

    //! Pseudo histograms with the CoG integral along x
    /*! This is an associative collection of pseudo histogram
     *  pointer. The key value is the sensor ID.
     *  Those pseudo histograms
     *  will be filled during the end() with the integral function of
     *  _cogHistogramX
     */
    std::map< int,EUTelPseudo1DHistogram* > _integralHistoX;

    //! Pseudo histograms with the CoG integral along y
    /*! This is an associative collection of pseudo histogram
     *  pointer. The key value is the sensor ID.
     *  Those pseudo histograms
     *  will be filled during the end() with the integral function of
     *  _cogHistogramY
     */
    std::map<int, EUTelPseudo1DHistogram* > _integralHistoY;

    //! Number of detector planes in the run
    /*! This is the total number of detector saved into this input
     *  file
     */
    int _noOfDetector;

    //! The detector name
    /*! This is got from the input header and copied into the
     *  condition file header
     */
    std::string _detectorName;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! The left end of the pseudo histogram
    static const double _min;

    //! The right end of the pseudo histogram
    static const double _max;


#ifdef MARLIN_USE_HISTOGRAM
    //! Base name for CoG histogram along x
    static std::string _cogHistogramXName;

    //! Base name for CoG histogram along y
    static std::string _cogHistogramYName;

    //! Base name for CoG integral histogram along x
    static std::string _cogIntegralXName;

    //! Base name for CoG integral histogram along y
    static std::string _cogIntegralYName;

    //! Base name for the Eta histogram along x
    static std::string _etaHistoXName;

    //! Base name for the Eta histogram along x
    static std::string _etaHistoYName;


    //! Reject singple pixel cluster from the eta calculation
    int _rejectsingplepixelcluster;

    //! Base name for the CoG histo 2D
    /*! This 2D histogram is used to show where the charge center of
     *  gravity is found within the seed pixel. It is telling if there
     *  is a certain correlation between the x and y directions and
     *  how the CoG is biased
     */
    static std::string _cogHisto2DName;
#endif

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

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

    //! Already booked SensorID
    std::set< int > _alreadyBookedSensorID;

  };

  //! A global instance of the processor
  EUTelCalculateEtaProcessor gEUTelCalculateEtaProcessor;

}
#endif

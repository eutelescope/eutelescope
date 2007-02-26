// -*- C++ -*-
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

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 

// system includes <>
#include <string>
#include <map>

// forward declaration
class PseudoHistogram;

#undef MARLIN_USE_HISTOGRAM
#ifdef MARLIN_USE_AIDA
#define MARLIN_USE_HISTOGRAM
#endif
#ifdef MARLIN_USE_ROOT
#define MARLIN_USE_HISTOGRAM
#endif


#ifdef MARLIN_USE_AIDA
namespace AIDA {
  class IBaseHistogram;
}
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
   *  <h4>Input - Prerequisites</h4> <b>ClusterCollection</b>. A
   *  collection of clusters. The name of this collection has to be
   *  specified in the steering file.  <br><b>EventNumber</b>. This is
   *  the number of events will be used for the eta calculation. Set
   *  it to -1 to use all available
   *  events. <br><b>NumberOfBins</b>. This vector contains the number
   *  of bins along x and y respectively to be used while binning the
   *  two eta function
   *  projections. <br><b>ClusterQualitySelection</b>. This parameter
   *  is used to select the quality of clusters used for Eta
   *  calculation. Reccomended is 0 ==
   *  kGoodCluster. <br><b>ClusterTypeSelection</b>. This string is
   *  used to choose how the charge center of gravity is
   *  calculated. Type "FULL" to use the full cluster, whatever number
   *  of pixels or shape it has; type "NxMPixel" to use only pixels of
   *  the cluster belonging to a NxM region around the seed; type
   *  NPixel to use only the N pixels of the cluster having the
   *  highest signal. <br><b>NxMPixelClusterSize</b>. This parameter
   *  is used only when ClusterTypeSelection=="NxMPixel" and
   *  represents the size along x and y of the cluster region of
   *  interest. In case this is bigger or equal to the full cluster,
   *  then "FULL" is used. <br><b>NPixelSize</b>. This parameter is
   *  used only when ClusterTypeSelection=="NPixel" and represents the
   *  number of pixels with the highest signal to be used in the CoG
   *  calculation. In case this number is exceeding the cluster
   *  multiplicity then "FULL" is used.
   *
   *  <h4>Output</h4>
   *  
   *  @param ClusterCollectionName The name of the input cluster
   *  collection.
   *  @param EventNumber The number of events to use for the calculation
   *  @param NumberOfBins The number of bins (x and y) for the calculation
   *  @param ClusterQualitySelection The cluster quality for the calculation
   *  @param ClusterTypeSelection The way the CoG is calculated
   *  @param NxMPixelClusterSize Number of pixels when ClusterTypeSelection=="NxMPixel"
   *  @param NPixelSize Number of pixels when ClusterTypeSelection="NPixel"
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  
   *  @version $Id: EUTelCalculateEtaProcessor.h,v 1.1 2007-02-26 09:30:20 bulgheroni Exp $
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
    /*  This is called for each event in the file. Pseudohistograms
     *  are filled
     * 
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
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
    virtual void finishCalculation();


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
     *  _nEvent is exceeding the total number of events in the file,
     *  the calculation will stop when the last available event is
     *  used.
     *
     *  @bug _nEvent is never reached and the eta not calculated when
     *  for example the user from the steering file set a global
     *  MaxRecordNumber lesser than _nEvent. To fix it, there should
     *  be a way to have access to global parameter from within a
     *  processor. 
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
    /*! This is a vector of pointers to pseudo histogram objects one
     *  for each detector in the telescope. Those histograms are
     *  filled with signed distanced between the CoG and the seed
     *  coordinates
     */
    std::vector<PseudoHistogram* > _cogHistogramX;

    //! Pseudo histograms for CoG along y
    /*! This is a vector of pointers to pseudo histogram objects one
     *  for each detector in the telescope. Those histograms are
     *  filled with signed distanced between the CoG and the seed
     *  coordinates
     */
    std::vector<PseudoHistogram* > _cogHistogramY;

    //! Pseudo histograms with the CoG integral along x
    /*! This is a vector of pointers to pseudo histogram objects one
     *  for each detector in the telescope. Those pseudo histograms
     *  will be filled during the end() with the integral function of
     *  _cogHistogramX
     */
    std::vector<PseudoHistogram* > _integralHistoX;

    //! Pseudo histograms with the CoG integral along y
    /*! This is a vector of pointers to pseudo histogram objects one
     *  for each detector in the telescope. Those pseudo histograms
     *  will be filled during the end() with the integral function of
     *  _cogHistogramY
     */
    std::vector<PseudoHistogram* > _integralHistoY;

    //! Number of detector planes in the run
    /*! This is the total number of detector saved into this input
     *  file 
     */
    int _noOfDetector;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

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

    //! Base name for the CoG histo 2D
    /*! This 2D histogram is used to show where the charge center of
     *  gravity is found within the seed pixel. It is telling if there
     *  is a certain correlation between the x and y directions and
     *  how the CoG is biased
     */
    static std::string _cogHisto2DName;
#endif

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
#endif

  };

  //! A global instance of the processor
  EUTelCalculateEtaProcessor gEUTelCalculateEtaProcessor;      

}
#endif

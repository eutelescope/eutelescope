// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELHISTOGRAMMAKER_H
#define EUTELHISTOGRAMMAKER_H 1

// eutelescope includes ".h" 
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>

// AIDA includes <.h>
#ifdef MARLIN_USE_AIDA
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <vector>
#include <map>
#include <string>


namespace eutelescope {

  //! Histogram filler for EUTelescope
  /*! This processor has the task to fill many control histograms that
   *  can serve as benchmark to understand if the detectors are
   *  working or not. 
   *
   *  In the output histogram file, the user will find a folder for
   *  each of the detector in the telescope and a bunch of histograms
   *  described here below:
   *
   *  @li <b> Cluster Signal </b> (histo name = "clusterSignal"). This
   *  is a cluster spectrum using all the pixels in the cluster.
   *
   *  @li <b> Cluster Signal N x N </b> (histo name =
   *  "clusterSignal3x3" and "clusterSignal4x4"). Those two histograms
   *  are the cluster spectrum considering only the pixels contained
   *  into the 5x5 and 3x3 sub-matrix around the seed. Those pixels
   *  are not sorted!
   *
   *  @li <b> Cluster Signal N </b> (histo name = "clusterSignalN"
   *  with N an integer number). Those histograms are the cluster
   *  spectrum considering only the N most significant pixels.
   *
   *  @li <b> Seed pixel signal </b> (histo name = "seedSignal"). This
   *  is the signal spectrum of the highest pixel in the cluster. This
   *  is the place to look for the calibration peak.
   
   *  @li <b> Hit map</b> (histo name = "hitMap"). This is a 2D plot
   *  with each bin representing one pixel in the sensor. The number
   *  of entries for each bin corresponds to the number of times a
   *  seed pixel has been found on that pixel.
   *
   *  <h4>Input</h4>
   *  
   *  @param A TrackerPulse input collection
   *  @param A vector containing the cluster N x N spectra to be
   *  filled
   *  @param A vector containing the cluster N spectra to be filled.
   *  
   *  <h4>Output</h>
   *
   *  None
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version  $ 
   *
   */
  
  class EUTelHistogramMaker : public marlin::Processor {

  public:
    
    
    //! Returns a new instance of EUTelHistogramMaker
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelPedestalNoiseProcess.
     */
    virtual Processor * newProcessor () {
      return new EUTelHistogramMaker;
    }

    //! Default constructor 
    EUTelHistogramMaker ();

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
    /*  All histograms are filled with the current event information.
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Check event method
    /*! Nothing to check
     * 
     *  @param evt The LCEvent event as passed by the ProcessMgr
     *
     */
    virtual void check (LCEvent * evt);

    //! Book histograms
    /*! This method is used to prepare the needed directory structure
     *  within the current ITree folder and books all required
     *  histograms. Histogram pointers are stored into
     *  EUTelHistogramMaker::_aidaHistoMap so that they can be
     *  recalled and filled from anywhere in the code.  
     */
    void bookHistos();


    //! Called after data processing.
    /*! This method is called when the loop on events is finished. It
     *  prints out a good bye message and nothing more
     */
    virtual void end();

  protected:

    //! Input collection name.
    /*! This is the name of the TrackerPulse input collection.
     */
    std::string _pulseCollectionName;

  private:

    //! List of cluster spectra N
    /*! This vector contains a list of cluster spectra we want to fill
     *  in. 
     */ 
    std::vector<int > _clusterSpectraNVector;

    //! List of cluster spectra NxN
    /*! This vector contains a list of cluster spectra N x N we want
     *   to fill. For example, if it contains "3", then the cluster 3x3
     *   spectrum will be filled.
     */
    std::vector<int > _clusterSpectraNxNVector;

    //! The number of detectors
    /*! The number of sensors in the telescope. This is retrieve from
     *  the run header
     */ 
    int _noOfDetector;

    //! The first pixel on the x side
    /*! One value for each detector.
     */ 
    std::vector<int > _minX;
    
    //! The first pixel on the y side
    /*! One value for each detector.
     */ 
    std::vector<int > _minY;

    //! The last pixel on the x side
    /*! One value for each detector.
     */ 
    std::vector<int > _maxX;

    //! The last pixel on the y side
    /*! One value for each detector.
     */ 
    std::vector<int > _maxY;
    
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

    //! Cluster signal histogram base name.
    /*! This is the name of the cluster signal histogram. To this
     *  name, the detector number is added in order to make it
     *  unique. This is used also for the other cluster spectra.
     */
    static std::string _clusterSignalHistoName;

    //! Seed pixel signal histo name
    /*! This is the seed pixel spectrum histogram. To this name,
     *  the detector ID is added in order to make it unique.
     */ 
    static std::string _seedSignalHistoName;

    //! Hit map histogram name
    /*! This is the hit map in pixel coordinate histogram name.
     */ 
    static std::string _hitMapHistoName;

#endif 
    
    //! Event counter
    /*! This is used to count the processed events
     */ 
    int _iEvt;

  };

  //! A global instance of the processor
  EUTelHistogramMaker gEUTelHistogramMaker;      

}
#endif

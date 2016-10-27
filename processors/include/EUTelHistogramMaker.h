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
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
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
   *  <h4>Input collections </h4>
   *  A tracker pulse collection with clusters to be histogrammed.
   *
   *  @param PulseCollectionName The name of the TrackerPulse input
   *  collection.
   *
   *  @param HistoInfoFileName The name of the XML information file
   *  with the histogram booking information.
   *
   *  @param ClusterNxN A vector containing the cluster N x N spectra to be
   *  filled
   *
   *  @param ClusterN A vector containing the cluster N spectra to be filled.
   *
   *  <h4>Output collections </h4>
   *
   *  None
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
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

    //! Initialize the geometry
    /*! Get the geometry information, detector size and sensor ID from
     *  the noise collection. Since version v00-00-09, the user must
     *  provide valid noise and status collections. In case the real
     *  pedestal is missing, the user can generate a fake one by using
     *  the EUTelAutoPedestalNoiseProcessor.
     *
     *  @param event An LCIO event.
     */
    void initializeGeometry( LCEvent * event );

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

    //! Noise collection name.
    /*! This is the name of the TrackerData collection containing the
     *  noise information. The presence of this collection is not
     *  compulsory. If it is found in the event, then also the "noise"
     *  related histograms are filled, otherwise those are skipped.
     */
    std::vector<std::string> _noiseCollectionName;

    //! Status collection name.
    /*! This is the name of the TrackerRawData collection containing
     *  the pixel status information. As for the noise collection,
     *  also this one is not compulsory for the proper behaviour of
     *  this processor. Having the status collection will allow the
     *  filling the of noise related histograms.
     */
    std::vector<std::string> _statusCollectionName;

    //! The histogram information file
    /*! This string contain the name of the histogram information
     *  file. This is selected by the user in the steering file.
     *
     *  @see eutelescope::EUTelHistogramManager
     *  @see eutelescope::EUTelHistogramInfo
     */
    std::string _histoInfoFileName;

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

    //! The ancillary map
    /*! This is a map relating the sensorID and the position of such a
     *  sensorID in the noise / status collection.
     */
    std::map< int, int > _ancillaryMap;

    //! SensorID vector
    /*! This is a vector of sensorID
     */
    std::vector< int > _sensorIDVec;

    //! Geometry ready flag
    bool _isGeometryReady;

    //! The first pixel on the x side
    /*! One value for each detector.
     */
    std::map< int, int > _minX;

    //! The first pixel on the y side
    /*! One value for each detector.
     */
    std::map< int, int > _minY;

    //! The last pixel on the x side
    /*! One value for each detector.
     */
    std::map< int, int > _maxX;

    //! The last pixel on the y side
    /*! One value for each detector.
     */
    std::map< int, int > _maxY;

    //! Switch to turn on and off the noise histo filling
    /*! To fill noise and SNR related histograms, the presence of the
     *  noise and the status collections has to be verified first.
     */
    bool _noiseHistoSwitch;

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

    //! Seed pixel SNR name
    /*! This is the seed pixel SNR histogram name
     */
    static std::string _seedSNRHistoName;

    //! Cluster SNR histogram name
    /*! This is the name of the histogram containing the seed pixel
     *  SNR
     */
    static std::string _clusterSNRHistoName;

    //! Event multiplicity histogram name
    /*! There is one of this histogram for each plane in the telescope
     *  setup. For each event this histogram is filled with the number
     *  of clusters passing the thresholds found.
     */
    static std::string _eventMultiplicityHistoName;

    //! Cluster noise histogram name
    /*! This is the full cluster noise histogram name
     */
    static std::string _clusterNoiseHistoName;


    //! Number of hit pixel
    /*! This is a histogram showing the number of hit pixel inside
     *  a digital fixed frame cluster
     */
    static std::string _clusterNumberOfHitPixelName;
#endif

    //! Run counter
    /*! This is used to count the processed run header
     */
    int _iRun;

  };

  //! A global instance of the processor
  EUTelHistogramMaker gEUTelHistogramMaker;

}
#endif

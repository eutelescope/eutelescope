// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCLUSTERINGPROCESSOR_H
#define EUTELCLUSTERINGPROCESSOR_H 1

// eutelescope includes ".h" 
#include "EUTelExceptions.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// aida includes <.h>
#ifdef MARLIN_USE_AIDA
#include <AIDA/IBaseHistogram.h>
#endif

// lcio includes <.h> 
#include <IMPL/TrackerRawDataImpl.h>

// system includes <>
#include <string>
#include <map>
#include <cmath>
#include <vector>

namespace eutelescope {

  //! Clustering processor for the EUTelescope
  /*! This processor is used to search within the current data matrix
   *  (after pedestal subtraction and eventual common mode
   *  suppression) for clusters. In this contest a cluster is a group
   *  of neighboring pixels, fulfilling certain conditions and
   *  representing the signal produced in the detector by the passage
   *  of a particle. Hence, clustering is then considered the first
   *  step in a tracking procedure.
   *
   *  To obtain space points from clusters, the actual geometrical
   *  description of the sensors (pixel pitches, definition of the
   *  local frame of reference, relation to the global frame of
   *  reference are required. Moreover, a suitable algorithm to
   *  calculate the cluster center, in a linear (charge center of
   *  gravity) or non-linear (eta function) is needed.
   *
   *  There are different ways to build up a clusters: the user can
   *  choose which algorithm to use via the clusteringAlgo
   *  parameter. Different algorithms may requires different way to
   *  store, at least temporary, the cluster information in
   *  memory. For example, in the case of fixed frame NxM clustering
   *  the natural object to store a cluster is a EUTelFFClusterImpl
   *  object that is inheriting from TrackerData; other, algorithms
   *  can used different re-implementation of TrackerData. 
   *
   *  In order to write data on disk in the most general and
   *  consistent way, whatever object has been used during the
   *  clustering procedure to store cluster info, those are then moved
   *  to a TrackPulse before writing. The original information about
   *  the TrackerData is saved as a reference in the TrackerPulse for
   *  possible further use.
   *
   *  The TrackerPulse cell id encoding is very similar to the
   *  EUTELESCOPE::CLUSTERDEFAULTENCODING but instead of having the
   *  quality, it has another field named ClusterType used to identify
   *  the class used to store the cluster information.
   *
   *  @see processEvent(LCEvent*) for a detailed description of
   *  available algorithms.
   *
   *  <h4>Input collections</h4> 
   *
   *  <b>Data collection</b>: the input data TrackerData collection
   *  name. This collection is containing data pedestal subtracted
   *  and, eventually, common mode corrected
   *
   *  <b>Noise collections</b>: the noise TrackerData collection
   *  name as read from the Condition Processor. This object contains
   *  the noise information of all pixels in all telescope
   *  detectors. Those numbers are used to calculate seed and cluster
   *  signal to noise ratio.
   *  
   *  <b>Status collection</b>: the status TrackerData collection
   *  name as read from the Condition Processor. This object is used
   *  to exclude from the clustering procedure pixels defined as bad
   *  during the pedestal processor.
   *
   *  <h4>Output collections</h4>
   *
   *  <b>Pulse collection</b>: this is the TrackerPulse collection
   *  containing all clusters found in the event.
   *
   *  @param DataCollectionName The name of the input data collection.
   *
   *  @param NoiseCollectionName The name of the noise collection.
   *
   *  @param StatusCollectionName The name of the status collection.
   *
   *  @param PulseCollectionName The name of the output TrackerPulse collection.
   * 
   *  @param ClusteringAlgo  This is a string representing which algorithm
   *  should be used for the clustering procedure.
   *
   *  @param ClusterSizeX The maximum size of the cluster in pixel
   *  unit along the x direction. It is used only for
   *  EUTELESCOPE::FIXEDFRAME algorithm. It has to be an odd number
   *  since the seed pixel is bound to be the cluster center.
   *
   *  @param ClusterSizeY As the ClusterSizeX but along the y
   *  direction.
   *
   *  @param SeedPixelCut This is the SNR threshold used to identify seed
   *  pixel candidates. One global number for all detectors.
   *
   *  @param ClusterCut  This is the SNR threshold used to accept cluster
   *  candidates.
   *
   *  @param HistoInfoFileName This is the name of the XML file
   *  containing the histogram booking information.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelClusteringProcessor.h,v 1.12 2007-07-11 06:52:30 bulgheroni Exp $
   *
   */

  class EUTelClusteringProcessor :public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTelClusteringProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelClusteringProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelClusteringProcessor;
    }

    //! Default constructor 
    EUTelClusteringProcessor ();

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
    /*! It looks for clusters in the current event using the selected
     *  algorithm.
     *
     *  @see EUTelClusteringProcessor::fixedFrameClustering(LCEvent *)
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

    //! Reset the status map
    /*! This method is called at the beginning of the clustering
     *  procedure because it is possibly containing the position of
     *  the previous identified clusters. Hit pixels are identified by
     *  the value EUTELESCOPE::HITPIXEL; during the reset all of them
     *  are set to EUTELESCOPE::GOODPIXEL. This is not touching the
     *  bad pixels since them are marked with EUTELESCOPE::BADPIXEL.
     *
     *  @param status A pointer to the TrackerRawData with the status
     *  to be reset
     *
     *  @todo Consider the possibility to use instead of
     *  EUTELESCOPE::HITPIXEL, the clusterID to flag hit pixel. This
     *  is offering a very easy way to show on a 2D histograms where
     *  clusters have been found. It might be of any usefulness if we
     *  will try to write a piece of code to deconvolve merging
     *  clusters.
     */
    void resetStatus(IMPL::TrackerRawDataImpl * status);

#ifdef MARLIN_USE_AIDA
    //! Book histograms
    /*! This method is used to prepare the needed directory structure
     *  within the current ITree folder and books all required
     *  histograms. Histogram pointers are stored into
     *  EUTelHistogramMaker::_aidaHistoMap so that they can be
     *  recalled and filled from anywhere in the code.  
     */
    void bookHistos();

    //! Fill histograms
    /*! This method is called for each event and the cluster
     *  information are inserted into specific AIDA histograms.
     *
     *  @param evt The current event object
     */ 
    void fillHistos(LCEvent * evt);
#endif

  protected:
    
    //! Wrapper of the processEvent(LCEvent*) for fixed frame clustering
    /*! This method is called by the processEvent method in the case
     *  the user selected the EUTELESCOPE::FIXEDFRAME algorithm for
     *  clustering.
     *
     *  This algorithm is based on the reconstruction of clusters
     *  having a rectangular predefined maximum shape centered around
     *  the seed pixel. The seed pixel is defined as the one with the
     *  highest signal. Follows a brief description of the algorithm:
     *  
     *  \li The full data matrix is scanned searching for seed pixel
     *  candidates. A seed candidate is defined as a pixel with a
     *  signal to noise ratio in excess the
     *  EUTelClusteringProcessor::_seedPixelCut defined by the
     *  user. All candidates are added to a map
     *  (EUTelClusteringProcessor::_seedCandidateMap) where the first
     *  template element is the (float) pixel charge and the second is
     *  the pixel index. 
     *
     *  \li The use of a map has the advantage that the belonging
     *  elements are sorted according to the "key" value,
     *  i.e. according to the pixel signal. In this way, at the end of
     *  the matrix crossing, the last element of the map is the pixel
     *  seed candidate with the highest signal. A known limitation of
     *  this approach is the impossibility to store two pixels with
     *  exactly the same signal. Nevertheless, this approach is
     *  competitive compared to using the pixel index as map key. This
     *  second approach excludes the possibility to have two pixels
     *  with the same key, but requires a more complicated sorting
     *  algorithm. The seed candidate map has to be compulsory sorted,
     *  because the cluster building procedure has to start from a
     *  seed pixel.
     *
     *  \li Starting from the last entry of the seed candidate map
     *  (i.e. the pixel with the highest signal in the matrix), a
     *  candidate cluster is built around this seed. The clustering is
     *  done with two nested loops in way that the seed pixel is the
     *  center of the resulting cluster. Only pixels with a good
     *  status, effectively belonging to the matrix (1) and not yet
     *  belonging to the any other clusters are added to the current
     *  cluster candidate.
     *
     *  \li A cluster candidate is finally accepted as a good cluster
     *  if its SNR is passing the _clusterCur threshold. Each good
     *  cluster is added to the current event using a TrackerData
     *  class. The cellID encoding used is the
     *  EUTELESCOPE::CLUSTERDEFAULTENCODING where along with the
     *  detector number, also the cluster id, the seed pixel
     *  coordinates and the cluster sizes are stored so that the
     *  cluster can be reconstructed.
     *
     *  (1) The 2D coordinates of each pixel is determined using the
     *  pixel index information along with the size along X of the
     *  matrix. For this purpose, some utility methods have been
     *  defined; those methods don't perform any consistency check of
     *  the obtained results. This means that the use of such methods
     *  can result into pixel coordinates outside the actual valid
     *  range. For this reason, during the clustering, a check on the
     *  validity of the x, y pair is required.
     *
     *  @throw IncompatibleDataSetException in the case the two
     *  collections are found to be incompatible
     *
     *  @param evt The LCIO event has passed by processEvent(LCEvent*)
     */
    void fixedFrameClustering(LCEvent * evt);


    //! Input collection name.
    /*! The input collection is the calibrated data one coming from
     *  the EUTelCalibrateEventProcessor. It is, usually, called
     *  "data" and it is a collection of TrackerData
     */
    std::string _dataCollectionName;

    //! Noise collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _noiseCollectionName;

    //! Status collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _statusCollectionName;

    //! Pulse collection name.
    /*! This is the name used to store the output cluster
     *  collection. 
     */ 
    std::string _pulseCollectionName;

    
    //! Cluster collection name.
    /*! This is the name of the collection to store the in memory and
     *  on the disk the cluster. This will be a collection of
     *  TrackerData but it can be reimplemented depending on the
     *  clustering algorithm.
     */ 
    std::string _dummyCollectionName;


    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Clusterization algorithm
    /*! This string is used to select which clustering algorithm
     *  should be used. Follows a list of available algorithm:
     *
     *  \li <b>FixedFrame</b>: Selectable also using the
     *  EUTELESCOPE::FIXEDFRAME static constant, it allows to select
     *  clusters with a fixed size along x and y centered around the
     *  seed pixel. The seed pixel is identified as the one with the
     *  highest signal. The full matrix is scanned row wise both the
     *  index and the signal of pixel passing the seed pixel threshold
     *  are recorded into a map using the index as a key. The map is
     *  then sorted according to the signal and then, starting from
     *  the first entry in the map a new cluster is created around the
     *  seed pixel candidate. The resulting cluster is accepted if it
     *  passing also the cluster threshold. Once a pixel has been
     *  assigned to a cluster, to avoid double counting, it cannot be
     *  assigned to any other clusters and is removed from further
     *  operations.
     */
    std::string _clusteringAlgo;


    //! Cluster size along x in pixel
    /*! This parameter is used in the case the _clusteringAlgo is set
     *  to EUTELESCOPE::FIXEDFRAME and represents the maximum size a
     *  cluster can have along x. It has to be an odd number since the
     *  seed pixel has to lay in the cluster center.
     */
    int _xClusterSize;

    //! Cluster size along y in pixel
    /*! This parameter is used in the case the _clusteringAlgo is set
     *  to EUTELESCOPE::FIXEDFRAME and represents the maximum size a
     *  cluster can have along y. It has to be an odd number since the
     *  seed pixel has to lay in the cluster center.
     */
    int _yClusterSize;

    //! Threshold for seed pixel identification
    /*! This float number represents the threshold in SNR units for
     *  the identification of seed pixels. All pixels passing this
     *  threshold are considered seed pixel candidates and added to
     *  the candidate map. They are eventually removed from the map,
     *  if they are found to wing pixels surrounding a higher seed
     *  pixel.
     */
    float _seedPixelCut;

    //! Threshold for cluster identification
    /*! This float number represents the threshold in SNR units for
     *  the cluster identification. Once a cluster candidate is built
     *  centered around its seed, to be considered a real cluster the
     *  total SNR has to pass this cluster threshold.
     */
    float _clusterCut;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

    //! Fill histogram switch
    /*! This boolean is used to switch on and off the filling of
     *  histograms.
     */ 
    bool _fillHistos;

    //! The histogram information file
    /*! This string contain the name of the histogram information
     *  file. This is selected by the user in the steering file.
     * 
     *  @see eutelescope::EUTelHistogramManager
     *  @see eutelescope::EUTelHistogramInfo
     */ 
    std::string _histoInfoFileName;    


  private:

    //! The seed candidate pixel map. 
    /*! The map key is the pixel index in the matrix, while the float
     *  value the pixel signal
     */ 
    std::map<float, unsigned int> _seedCandidateMap;
    

    //! The current detector in the input collection
    /*! This int number represent the current detector being
     *  analyzed. It has been declared as a class data member because
     *  it is frequently used by several different methods.
     */
    int _iDetector;

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

    //! Total cluster found 
    /*! This array of int is used to store the total number of
     *  clusters found for each detector. The contents of this array
     *  is shown during the end()
     */ 
    IntVec _totCluster;

#ifdef MARLIN_USE_AIDA
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

  };

  //! A global instance of the processor
  EUTelClusteringProcessor gEUTelClusteringProcessor;      

}
#endif

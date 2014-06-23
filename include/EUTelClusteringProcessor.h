// Version: $Id$
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

// since v00-00-09 built only if GEAR is available
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelExceptions.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/EventModifier.h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// lcio includes <.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include <list>

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
   *  object but all of them should inherit from the
   *  EUTelVirtualCluster, that is the base cluster class in the
   *  EUTelescope framework. All these classes are working according
   *  to the Decorator Pattern, and they have at least a TrackerData
   *  pointer as data member.
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
   *  @since Since version v00-00-09, this processor requires GEAR to
   *  be initialized because the geometry information are no more
   *  taken from the input file Run Header but they are gathered from
   *  the GEAR description.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $$
   *
   */

  class EUTelClusteringProcessor :public marlin::Processor , public marlin::EventModifier {

  public:

    //a helper class for faster digital cluster finding on the sensors
    //(basically a two dimensional array)
    template <class T>
    class dim2array
    {
    public:
      dim2array() : array(NULL), size_x(0), size_y(0)
      {
        createarray();
      }
      dim2array(int x, int y) : array(NULL), size_x(x), size_y(y)
      {
        createarray();
      }
      dim2array(const unsigned int x, const unsigned int y, T value) : array(NULL), size_x(x), size_y(y)
      {
        createarray();
        for(unsigned int i =0; i < (size_x * size_y); ++i)
          array[i] = value;
      }
      dim2array(const dim2array &a) : array(NULL), size_x(a.size_x), size_y(a.size_y)
      {
        //size_x = a.size_x;
        //size_y = a.size_y;
        createarray();
        for(unsigned int i =0; i < size_x*size_y; ++i)
          array[i] = a.array[i];
      }
      dim2array& operator = (const dim2array &a)
      {
        size_x = a.size_x;
        size_y = a.size_y;
        delete [] array;
        createarray();
        for(unsigned int i =0; i < size_x*size_y; ++i)
          array[i] = a.array[i];
        return *this;  
      }
      T at(const int i, const int j) const
      {
        int index = j * size_x + i;
       //  if (index < 0 || index >= (size_x * size_y)) {
//           std::cout << "debug i x " << i << " " << j << std::endl;
//           abort();
//         }
        return array[index];
      }
      void pad(T v)
      {
        for(unsigned int i =0; i < size_x*size_y; ++i)
          array[i] = v;
      }
      unsigned int sizeX() const
      {
        return (size_x);
      }
      unsigned int sizeY() const
      {
        return (size_y);
      }
      void set(const unsigned int i, const unsigned int j, T value)
      {
        int index = j * size_x + i;
        array[index] = value;
      }
  
      ~dim2array()
      {
        delete [] array;
      }
  
    private:
      void createarray()
      {
        array = new T[size_x * size_y];
      }

      T *array;

      unsigned int size_x;
      unsigned int size_y;
    };
    //std::vector< dim2array<bool> > sensormatrix;


    class pixel
    {
    public:
      pixel() : x(0), y(0) {}
      pixel(unsigned int tmp_x, unsigned int tmp_y) : x(tmp_x), y(tmp_y)
      {
        x = tmp_x;
        y = tmp_y;
      }
      unsigned int x;
      unsigned int y;
    };
    
    //class of seed candidates
    class seed
    {
    public:
      seed(unsigned int tmp_x, unsigned int tmp_y, unsigned int tmp_nb, unsigned int cp) : x(tmp_x), y(tmp_y), neighbours(tmp_nb), p(cp)
      {
        x = tmp_x;
        y = tmp_y;
        neighbours = tmp_nb;
        p = cp;
      }
      //this operator is needed for the sort algorithm. the first
      //criteria is the number of neighbouring pixels and then the
      //second criteria is the number of fired pixel in a cluster
      //around the seed
      bool operator<(const seed& b) const 
      {
        //return (measuredZ < b.measuredZ);
        bool r = true;
        if(neighbours == b.neighbours)
          {   
            if(p < b.p)
              r = false;
          }
        else 
          if(neighbours < b.neighbours)
            r = false;
        return r;
      }
      unsigned int x; //x coordinate
      unsigned int y;//y coordinate
      unsigned int neighbours; //number of neighbours
      unsigned int p; //total number of fired pixel in the cluster formed by
      //this seed pixel candidate
    };



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

    virtual const std::string & name() const { return Processor::name() ; }

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

    //! Modify event method
    /*! Actually don't used
     *
     *  @param evt the current LCEvent event as passed by the ProcessMgr
     */
    virtual void modifyEvent( LCEvent * evt ) ;

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
     *  //todo Consider the possibility to use instead of
     *  EUTELESCOPE::HITPIXEL, the clusterID to flag hit pixel. This
     *  is offering a very easy way to show on a 2D histograms where
     *  clusters have been found. It might be of any usefulness if we
     *  will try to write a piece of code to deconvolve merging
     *  clusters.
     */
    void resetStatus(IMPL::TrackerRawDataImpl * status);

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! Book histograms
    /*! This method is used to prepare the needed directory structure
     *  within the current ITree folder and books all required
     *  histograms. Histogram pointers are stored into
     *  vectors for class-wide access
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


    //! Initialize the geometry information
    /*! This method is called to initialize the geometry information,
     *  namely the total number of sensors to be used in the cluster
     *  search and the boundaries for each sensors.
     *  In case of NZS data these information are contained also in
     *  the cell ID of the input collection, but in case of ZS only
     *  the sensor ID is available there.
     *
     *  @param evt The LCEvent to be used for geometry
     *  initialization.
     *
     *  @throw In case the @event does not contain all the needed
     *  information, a SkipEventException is thrown and the geometry
     *  will be initialize with the following event.
     */
    void initializeGeometry( LCEvent * evt ) throw ( marlin::SkipEventException );

    //! initialize statusCollection
    /*! values from hotpixel DB file are read into the statusCollection
     * hot pixels get status FiringPixel
     */
    void initializeStatusCollection();
   
    //! initialize HotPixelMapVec
    /*! values from hotpixel DB file are read into a vector of mapped
     * unique pixel index
     */
    void initializeHotPixelMapVec();
 
  protected:

    //! Method for fixed frame clustering
    /*! This method is called by the processEvent method in the case
     *  the user selected the EUTELESCOPE::FIXEDFRAME algorithm for
     *  clustering.
     *
     *  This algorithm is based on the reconstruction of clusters
     *  having a rectangular predefined maximum shape centered around
     *  the seed pixel. The seed pixel is defined as the one with the
     *  highest signal. Here follows a brief description of the algorithm:
     *
     *  \li The full data matrix is scanned searching for seed pixel
     *  candidates. A seed candidate is defined as a pixel with a
     *  signal to noise ratio in excess the
     *  EUTelClusteringProcessor::_seedPixelCut defined by the
     *  user. All candidates are added to a vector< pair< float, int > >
     *  (EUTelClusteringProcessor::_seedCandidateMap) where the first
     *  template element is the (float) pixel charge and the second is
     *  the pixel index.
     *
     *  \li This is the sorted using the std::algorithm library by the first element of the pair
     *  i.e. the pixel signal. In this way, at the end of
     *  the matrix crossing, the last element of the map is the pixel
     *  seed candidate with the highest signal. The seed candidate map has to be compulsory sorted,
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
     *  cluster candidate. //TODO: (Phillip Hamnett) CHECK THIS
     *
     *  \li A cluster candidate is finally accepted as a good cluster
     *  if its SNR is passing the _clusterCut threshold. Each good
     *  cluster is added to the current event using a TrackerData
     *  class. The cellID encoding used is the
     *  EUTELESCOPE::CLUSTERDEFAULTENCODING where along with the
     *  detector number, also the cluster id, the seed pixel
     *  coordinates and the cluster sizes are stored so that the
     *  cluster can be reconstructed.
     *
     *  (1) The 2D coordinates of each pixel are determined using the
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
     *  @param pulse The collection of pulses to append the found
     *  clusters.
     */
    void fixedFrameClustering(LCEvent * evt, LCCollectionVec * pulse);

    //! Method for zs Fixed Frame Clustering
    /*! This method is called by the processEvent method in the case
     */ 
    void zsFixedFrameClustering(LCEvent * evt, LCCollectionVec * pulse);

    //! Method for digital Fixed Frame Clustering
    /*! This method is called by the processEvent method in the case
     */ 
    void digitalFixedFrameClustering(LCEvent * evt, LCCollectionVec * pulse);


    //!HACK TAKI
    //! Methods for bricked pixel clustering
    //! zs Bricked Clustering
    /*! This method is called by the processEvent method in the case
     *  the user selected the EUTELESCOPE::BRICKEDCLUSTER algorithm for
     *  clustering.
     *
     *  This algorithm is based on the reconstruction of clusters
     *  surrounding the seed pixel in a bricked pixel structure.
     *  Similar to fixed frame clustering.
     *  At the moment(!) implemented for the nearest six neighbours
     *  only!
     *
     *  @throw IncompatibleDataSetException in the case the two
     *  collections are found to be incompatible
     *
     *  @param evt The LCIO event has passed by processEvent(LCEvent*)
     *  @param pulse The collection of pulses to append the found
     *  clusters.
     */
    void zsBrickedClustering(LCEvent * evt, LCCollectionVec * pulse);

    //!HACK TAKI
    //! Methods for bricked pixel clustering
    //! nzs Bricked Clustering
    void nzsBrickedClustering(LCEvent * evt, LCCollectionVec * pulse);

    //! TODO: Documentation
    void sparseClustering(LCEvent * evt, LCCollectionVec * pulse);


    //! Input collection name for NZS data
    /*! The input collection is the calibrated data one coming from
     *  the EUTelCalibrateEventProcessor. It is, usually, called
     *  "nzsdata" and it is a collection of TrackerData
     */
    std::string _nzsDataCollectionName;

    //! Input collection name for ZS data
    /*! The input collection is the calibrated data one coming from
     *  the EUTelCalibrateEventProcessor. It is, usually, called
     *  "zsdata" and it is a collection of TrackerData
     */
    std::string _zsDataCollectionName;

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

    //! Optional; HotPixel Collection Name.
    /*! This is the name used to store the output cluster
     *  collection.
     */
    std::string _hotPixelCollectionName;

 
    //! Pulse collection size
    size_t _initialPulseCollectionSize;

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

    //! Clusterization algorithm for full frames
    /*! This string is used to select which clustering algorithm
     *  should be used. Here follows a list of available algorithms:
     *
     *  \li <b>FixedFrame</b>: Selectable also using the
     *  EUTELESCOPE::FIXEDFRAME static constant, it allows to select
     *  clusters with a fixed size along x and y centered around the
     *  seed pixel. The seed pixel is identified as the one with the
     *  highest signal. The full matrix is scanned row wise both the
     *  index and the signal of pixel passing the seed pixel threshold
     *  are recorded into a vector< pair< float, int > >. This vector is
     *  then sorted according to the signal (float) and then, starting from
     *  the first entry in the map a new cluster is created around the
     *  seed pixel candidate. The resulting cluster is accepted if it
     *  is also passing the cluster threshold. Once a pixel has been
     *  assigned to a cluster, to avoid double counting, it cannot be
     *  assigned to any other clusters and is removed from further
     *  operations.
     *
     *  !HACK TAKI
     *  \li <b>Bricked</b>: Selectable also using the
     *  EUTELESCOPE::BRICKED static constant, it allows to
     *  select clusters within a bricked pixel structure.
     *  It will store the six surrounding pixels.
     */
    std::string _nzsClusteringAlgo;

    //! Clusterization algorithm for ZS frames
    /*! This string is used to select which clustering algorithm
     *  should be used. Follows a list of available algorithms:
     *
     *  \li <b>SparseCluster</b>: Selectable also using the
     *  EUTELESCOPE::SPARSECLUSTER static constant. Firstly the list
     *  of neighbor pixels is taken from the
     *  EUTelSparseDataImpl::findNeighborPixel, then they are added in
     *  a EUTelSparseClusterImpl cluster. The seed and cluster SNR
     *  cuts are then applied.
     *
     *  \li <b>SparseCluster2</b>:  Selectable also using the
     *  EUTELESCOPE::SPARSECLUSTER2 static constant. It works as
     *  SparseCluster but with improved performance.
     *
     *  !HACK TAKI
     *  \li <b>Bricked</b>: Selectable also using the
     *  EUTELESCOPE::BRICKED static constant, it allows to
     *  select clusters within a bricked pixel structure.
     *  For ZS data it will use as many pixels as it can
     *  find around the seed pixel (up to 6).
     */
    std::string _zsClusteringAlgo;

    //! Data format type options.
    /*! This string is used to optimise the clustering performance.
     *  The following types are known:
     *
     *  \li <b>Analog</b>: The data from the sensors is analog.
     *  Implies: any value in the range from Min to Max should be expected.
     *
     *  \li <b>Digital</b>: The data from the sensors is digital.
     *  Implies: the values in the range from Min to Max are expected to have
     *  some regular step.
     *
     *  \li <b>Binary</b>: The data from the sensors is binary.
     *  Implies: there are only two values of the signal possible - 0 and 1. 
     *  
     */
    std::string _dataFormatType;


    //! Cluster size along x in pixel
    /*! This parameter is used in the case the _clusteringAlgo is set
     *  to EUTELESCOPE::FIXEDFRAME and represents the maximum size a
     *  cluster can have along x. It has to be an odd number since the
     *  seed pixel has to lay in the cluster center.
     */
    int _ffXClusterSize;

    //! Cluster size along y in pixel
    /*! This parameter is used in the case the _clusteringAlgo is set
     *  to EUTELESCOPE::FIXEDFRAME and represents the maximum size a
     *  cluster can have along y. It has to be an odd number since the
     *  seed pixel has to lay in the cluster center.
     */
    int _ffYClusterSize;

    int _charge;

    //! Threshold for seed pixel identification
    /*! This float number represents the threshold in SNR units for
     *  the identification of seed pixels for Fixed Frame
     *  Algorithm. All pixels passing this threshold are considered
     *  seed pixel candidates and added to the candidate map. They are
     *  eventually removed from the map, if they are found to wing
     *  pixels surrounding a higher seed pixel.
     */
    float _ffSeedCut;

    //! Threshold for seed pixel in SparseCluster
    /*! The zero suppress reclustering algorithm may need a threshold
     *  for seed pixel SNR. To keep the framework as general as
     *  possible this value is different from the raw mode clustering
     *  threshold. Used only with EUTELESCOPE::SPARSECLUSTER and
     *  EUTELESCOPE::SPARSECLUSTER2.
     */
    float _sparseSeedCut;

    //! Threshold for cluster identification
    /*! This float number represents the threshold in SNR units for
     *  the cluster identification (Used only with Fixed Frame
     *  Algorithm). Once a cluster candidate is built centered around
     *  its seed, to be considered a real cluster the total SNR has to
     *  pass this cluster threshold.
     */
    float _ffClusterCut;

    //! Threshold for cluster SNR in SparseCluster
    /*! This value is used to accept clusters coming from ZS data. To
     *  keep the framework as general as possible this variable is
     *  kept different from the corresponding rawmode one. Used only
     *  with EUTELESCOPE::SPARSECLUSTER and
     *  EUTELESCOPE::SPARSECLUSTER2.
     */
    float _sparseClusterCut;

    int _sparseMinDistanceSquared;

    //! Minimum distance of neighbor pixels in SparseCluster
    /*! ZS pixel reclustering may need a minimum distance parameter to
     *  identify "close" pixels. Used only with
     *  EUTELESCOPE::SPARSECLUSTER and EUTELESCOPE::SPARSECLUSTER2.
     */
    float _sparseMinDistance;

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
  
  #ifndef DISALLOW_COPY_AND_ASSIGN
  //Following #define stops the accidental creation of a copy or assignment operator by causing a link error. Copy and Assignment operators not allowed because they are unnecessary and the cause of many bugs
  #define DISALLOW_COPY_AND_ASSIGN(EUTelClusteringProcessor) \
  EUTelClusteringProcessor(const EUTelClusteringProcessor&); \
  void operator=(const EUTelClusteringProcessor&);

  //Private Functions
  DISALLOW_COPY_AND_ASSIGN(EUTelClusteringProcessor)//See #define just above
  #endif

    //! read secondary collections
    /*!
     */
    void readCollections(LCEvent *evt);

    //! The seed candidate pixel map.
    /*! This is a vector which stores the seed index and the size of the signal. The signal is the floating point and the unsigned integer is the 
     */
    std::vector< std::pair<float,unsigned int> > _seedCandidateMap;    

    //! Total cluster found
    /*! This is a map correlating the sensorID number and the
     *  total number of clusters found on that sensor.
     *  The content of this map is show during end().
     */
    std::map< int, int > _totClusterMap;

    //! The number of detectors
    /*! The number of sensors in the telescope. This is retrieve from
     *  the run header
     */
    int _noOfDetector;
    
    //! List of excluded planes.
    /*! This vector contains a list of sensor ids for planes that have
     *   to be excluded from the clustering.
     */
    std::vector<int > _ExcludedPlanes;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
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

    //! Map for pointer to cluster signal histograms.
    std::map<int,AIDA::IBaseHistogram*> _clusterSignalHistos;

    //! Map for pointer to Cluster signal histogram (size along X).
    std::map<int,AIDA::IBaseHistogram*> _clusterSizeXHistos;

    //! Map for pointer to Cluster signal histogram (size along Y).
    std::map<int,AIDA::IBaseHistogram*> _clusterSizeYHistos;

     //! Map for pointer to Seed pixel signal histo 
    std::map<int,AIDA::IBaseHistogram*> _seedSignalHistos;

    //! Map for pointer to Hit map histogram 
     std::map<int,AIDA::IBaseHistogram*> _hitMapHistos;

    //! Map for pointer to Seed pixel SNR 
    std::map<int,AIDA::IBaseHistogram*> _seedSNRHistos;

    //! Map for pointer to Cluster noise histogram 
    std::map<int,AIDA::IBaseHistogram*> _clusterNoiseHistos;

    //! Map for pointer to Cluster SNR histogram 
    std::map<int,AIDA::IBaseHistogram*> _clusterSNRHistos;

    //! Map for pointer to Cluster vs Seed SNR histogram 
    std::map<int,AIDA::IBaseHistogram*> _cluster_vs_seedSNRHistos;

    //! Map for pointer to Event multiplicity histogram 
    std::map<int,AIDA::IBaseHistogram*> _eventMultiplicityHistos;

    //! Map (of maps) for pointers to histograms with cluster spectra with the X most significant pixels
    std::map<int, std::map<int,AIDA::IBaseHistogram*> > _clusterSignal_NHistos;

    //! Map (of maps) for pointers to histograms with cluster SRN spectra with the X most significant pixels
    std::map<int, std::map<int,AIDA::IBaseHistogram*> > _clusterSNR_NHistos;

    std::map<int, std::map<int,AIDA::IBaseHistogram*> > _clusterSignal_NxNHistos;
    std::map<int, std::map<int,AIDA::IBaseHistogram*> > _clusterSNR_NxNHistos;

#endif

    // gear stuff goes in here
    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    //! Geometry ready switch
    /*! This boolean reveals if the geometry has been properly
     *  initialized or not.
     */
    bool _isGeometryReady;

    //! Map relating telescope layer index and sensorID.
    /*! The first element is the sensor ID, while the second is the
     *  position of such a sensorID in the GEAR description.
     */
    std::map< int , int > _layerIndexMap;

    //! Map relating the DUT layer index and sensorID
    /*! The first element is the sensor ID, while the second is the
     *  position of such a sensorID in the DUT section of the GEAR
     *  description.
     *  So far, the DUT section had no more than one entry, but in the
     *  future this could be extended.
     *
     *  For the time being this map is created and properly filled but
     *  not yet used.
     */
    std::map< int, int > _dutLayerIndexMap;

    //! Map relating ancillary collection position and sensorID
    /*! The first element is the sensor ID, while the second is the
     *  position of such a sensorID in all the ancillary collections
     *  (noise, pedestal and status).
     */
    std::map< int, int > _ancillaryIndexMap;

    //! Inverse vector relation
    /*! This is the inverse relation with respect to the
     *  _ancillaryIndexMap. It contains the ordered list of sensorID
     */
    std::vector< int > _orderedSensorIDVec;

    //! SensorID vector
    /*! This is a vector of sensorID
     */
    std::vector< int > _sensorIDVec;

    //
    //! Zero Suppressed Data Collection
    LCCollectionVec *zsInputDataCollectionVec;
 
    //
    //! Non Zero Suppressed Data Collection
    LCCollectionVec *nzsInputDataCollectionVec;
    
    //
    //! pulse Collection 
    LCCollectionVec *pulseCollectionVec;
    
    //
    //! noise Collection 
    LCCollectionVec *noiseCollectionVec;
    
    //
    //! status Collection 
    LCCollectionVec *statusCollectionVec;



    //! Hot Pixel Collection 
    LCCollectionVec *hotPixelCollectionVec;

    //! keep what we have: Non ZS data (true/false)
    bool hasNZSData;
 
    //! keep what we have: ZS data (true/false)
    bool hasZSData;
  
    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  Key <int> is the sequential (counter) id of a hit,
     *  Value <int> - always status value = firing pixel
     */

    std::vector< std::map< int, int > > _hitIndexMapVec;          
    
  };

  //! A global instance of the processor
  EUTelClusteringProcessor gEUTelClusteringProcessor;

}
#endif // USE_GEAR
#endif

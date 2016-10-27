/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelProcessorGeometricClustering_H
#define EUTelProcessorGeometricClustering_H 1

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

namespace eutelescope {

  //! Geoemtric clustering processor for EUTelescope
  /*! This procssor used the Extended Geometry Framework (EGF) for a 
   *  correct spatial clustering. This means that via the EGF the
   *  geoemtrical positions as well as dimensions of the pixels are
   *  read in and used. 
   *
   *  The dimensions of the pixel are given by the surrounding rectangle.
   *  For simple rectengular pixels this is a correct description, 
   *  deviations of that pixel shape have to either adapt this processor
   *  or use it as an approximation.
   *
   *  Spatial proximity is simply defined as two rectangles touching.
   *  Since the EGF uses TGeo and floating point numbers are only stored
   *  in single precision, the code accounts for uncertainty by allowing
   *  a 1% deviation.
   *
   *  Given that the proximity is well defined, no additional arguments
   *  must be provided. If wanted, a time cut can be set. This will also
   *  require hits to be temporally in promximity. If not set not cut will
   *  be applied.
   *
   *  This clustering processor uses the @class EUTelGenericSparseClusterImpl 
   *  which derives from the new @class EUTelSimpleVirtualCluster base
   *  class.
   *
   *  The TrackerPulse cell id encoding is very similar to the
   *  EUTELESCOPE::CLUSTERDEFAULTENCODING but instead of having the
   *  quality, it has another field named ClusterType used to identify
   *  the class used to store the cluster information.
   *
   *  <h4>Input and Output collections</h4>
   *
   *  <b>Data Collection</b>: the input data TrackerData collection
   *  name. This collection is containing the zero-suppressed hits from
   *  previous analysis steps.
   *
   *  <b>Pulse Collection</b>: this is the TrackerPulse collection
   *  containing all clusters found in the event.
   *
   *  @param ZSDataCollectionName The name of the input data collection.
   *
   *  @param PulseCollectionName The name of the output TrackerPulse collection.
   *
   *  @param TCut This is the time cut value used to determine if hits are in
   *  temporal proximity. Values are in your detector specific time unit
   *
   *  @param HistoInfoFileName This is the name of the XML file
   *  containing the histogram booking information.
   *
   */

class EUTelProcessorGeometricClustering :public marlin::Processor , public marlin::EventModifier {

public:

    //! Returns a new instance of EUTelProcessorGeometricClustering
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelProcessorGeometricClustering.
     */
    virtual Processor* newProcessor() {
		return new EUTelProcessorGeometricClustering;
    }

    virtual const std::string & name() const { return Processor::name(); }

    //! Default constructor
    EUTelProcessorGeometricClustering ();

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
     *  @see EUTelProcessorGeometricClustering::fixedFrameClustering(LCEvent *)
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


	//TODO: Tobias
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
     *  search.
     *
     *  @param evt The LCEvent to be used for geometry
     *  initialization.
     *
     *  @throw In case the @event does not contain all the needed
     *  information, a SkipEventException is thrown and the geometry
     *  will be initialize with the following event.
     */
    void initializeGeometry( LCEvent * evt ) throw ( marlin::SkipEventException );

 
protected:
   
    //! Method for geometric clsutering
    /*! Algorithm which actually reads in teh collection of hit pixels
     *  and groups them together.
     *
     *  @param evt The LCIO event has passed by processEvent(LCEvent*)
     *  @param pulse The collection of pulses to append the found
     *  clusters.
     */
    void geometricClustering(LCEvent* evt, LCCollectionVec* pulse);

    //! Input collection name for ZS data
    /*! The input collection is the calibrated data one coming from
     *  the EUTelCalibrateEventProcessor. It is, usually, called
     *  "zsdata" and it is a collection of TrackerData
     */
    std::string _zsDataCollectionName;

    //! Pulse collection name.
    /*! This is the name used to store the output cluster
     *  collection.
     */
    std::string _pulseCollectionName;

    //! Pulse collection size
    size_t _initialPulseCollectionSize;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     *  events are counted from 0 and on a run base
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

	//! The time cut value as provided by the user.
	float _cutT;

private:
	DISALLOW_COPY_AND_ASSIGN(EUTelProcessorGeometricClustering)

    //! read secondary collections
    void readCollections(LCEvent *evt);

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

    //! Map for pointer to Hit map histogram 
     std::map<int,AIDA::IBaseHistogram*> _hitMapGeomHistos;

    //! Map for pointer to Cluster noise histogram 
    std::map<int,AIDA::IBaseHistogram*> _clusterNoiseHistos;

    //! Map for pointer to Event multiplicity histogram 
    std::map<int,AIDA::IBaseHistogram*> _eventMultiplicityHistos;

    //! Map for pointer to total cluster size histogram 
    std::map<int,AIDA::IBaseHistogram*>_clusterSizeTotalHistos;
#endif

    //! Geometry ready switch
    /*! This boolean reveals if the geometry has been properly
     *  initialized or not.
     */
    bool _isGeometryReady;

    //! SensorID vector
    /*! This is a vector of sensorID
     */
    std::vector< int > _sensorIDVec;

    //! Zero Suppressed Data Collection
    LCCollectionVec *_zsInputDataCollectionVec;
    
    //! pulse Collection 
    LCCollectionVec* _pulseCollectionVec;
  
};

//! A global instance of the processor
EUTelProcessorGeometricClustering gEUTelProcessorGeometricClustering;
}
#endif

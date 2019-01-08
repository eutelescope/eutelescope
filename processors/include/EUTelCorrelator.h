/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCORRELATOR_H
#define EUTELCORRELATOR_H

#if defined(USE_GEAR)
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRunHeader.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram2D.h>
#endif

// system includes <>
#include <map>
#include <string>
#include <vector>

namespace eutelescope {

  //! Hit and cluster correlator
  /*! This processor makes histograms that show the correlation
   *  between clusters of a detector and another.
   *  We study the correlation between centre of clusters' X and Y
   *  separately.
   *  At the end of our study we'll have (n*n - n) histograms for X
   *  and (n*n - n) for Y ( n is the number of our detectors ). We
   *  have this number of histograms because we don't create those
   *  histograms where there will be the correlation between a
   *  detector and himself.
   *  These histograms are put into different directories,
   *  one for X and one for Y, with cluster and hit correlation too.
   *  For each directory we have pairs of histograms that differ only
   *  in the order of DetectorID and these result with X and Y
   *  reversed, but otherwise they are the same.
   *  We do the same for hit correlation.
   *
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Cluster collection</b>: A collection with cluster
   *  produced by previous processors like EUTelClusteringProcessor or
   *  EUTelClusterFilter.
   */

  class EUTelCorrelator : public marlin::Processor {

  public:
    //! Returns a new instance of EUTelCorrelator
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelCorrelator.
     */
    virtual Processor *newProcessor() { return new EUTelCorrelator; }

    //! Default constructor
    EUTelCorrelator();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters
     */
    virtual void init();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader(LCRunHeader *run);

    //! Called every event
    /*! For each event, we loop over the input cluster collection and
     *  we fill in the correlation histograms using the x and y
     *  coordinate of the cluster center. Only clusters not belonging
     *  to the same sensors are used.
     *
     *  Histograms are booked in the first event.
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent(LCEvent *evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();

    //! Histogram booking
    /*! This method is used to book all the correlation histograms. It
     *  is called by processEvent when processing the first event.
     */
    void bookHistos();

    //! Internal function
    /*! Returns the ID of a plane selected as a reference
     *  plane for correlation plots
     */
    virtual int getFixedPlaneID() { return _fixedPlaneID; }

  protected:
  	//! Input collection name.
    /*! This is the name of the input hit collection.
     */
    std::string _inputHitCollectionName;

    //! Output collection name.
    /*! This is the name of the output hit collection.
     */
    std::string _outputHitCollectionName;
    
    //! set the plane to use as a reference/starting point for correlations
    int _fixedPlaneID;

    //! vector of correlation band cuts in X (upper limit)
    std::vector<float> _residualsXMax;
    //! vector of correlation band cuts in X (lower limit)
    std::vector<float> _residualsXMin;
    //! vector of correlation band cuts in Y (upper limit)
    std::vector<float> _residualsYMax;
    //! vector of correlation band cuts in Y (lower limit)
    std::vector<float> _residualsYMin;

	//! required number of correlated hits
    int _minNumberOfCorrelatedHits; 

    //! Cluster charge cut
    int _clusterChargeMin;

    //! How many events are needed to get reasonable correlation & offset values
    int _requiredEvents;

    //! Cluster collection list (EVENT::StringVec)
    EVENT::StringVec _clusterCollectionVec;

	//! Function for guessing the sensor offset
    std::vector<double> guessSensorOffset(int internalSensorID,
                                          int externalSensorID,
                                          std::vector<double> cluCenter);

  private:
    //! Initialization flag
    bool _isInitialize;

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! Correlation histogram matrix
    /*! This is used to store the pointers of each histogram
     */
    std::map<unsigned int, std::map<unsigned int, AIDA::IHistogram2D *>>
        _clusterXCorrelationMatrix;
    std::map<unsigned int, std::map<unsigned int, AIDA::IHistogram2D *>>
        _clusterYCorrelationMatrix;
    std::map<unsigned int, std::map<unsigned int, AIDA::IHistogram2D *>>
        _hitXCorrelationMatrix;
    std::map<unsigned int, std::map<unsigned int, AIDA::IHistogram2D *>>
        _hitYCorrelationMatrix;
    std::map<unsigned int, std::map<unsigned int, AIDA::IHistogram2D *>>
        _hitXCorrShiftMatrix;
    std::map<unsigned int, std::map<unsigned int, AIDA::IHistogram2D *>>
        _hitYCorrShiftMatrix;
#endif

    //! boolean to store if cluster/hit collection exists
    bool _hasClusterCollection;
    bool _hasHitCollection;

    //! vector of Sensor ID
    std::vector<int> _sensorIDVec;
    
    //! map of Sensor ID and z position
    std::map<int, int> _sensorIDtoZ;
  };
  
  //! A global instance of the processor
  EUTelCorrelator gEUTelCorrelator;
}

#endif // USE_GEAR
#endif

// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifdef EXPERIMENTAL
#ifndef EUTELCLUSTERFILTER_H
#define EUTELCLUSTERFILTER_H 1

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 

// system includes <>


namespace eutelescope {

  //! Cluster filter
  /*! This processor is used during the analysis chain to perform a
   *  further filtering on the available clusters. The main advantage
   *  of having this processor is that, especially when the user is
   *  getting used with a new detector, the cluster search
   *  (EUTelClusteringProcessor) should be executed with loose cuts,
   *  and then when the expected positions of peak signal is known,
   *  the previous results can be reprocessed very quickly just
   *  running this cluster filter.
   *
   *  Being born exactly during a beam test, this processor has to be
   *  considered always work in progress and it will take sometimes
   *  before it will reach its final and stable configuration.
   *
   *  The idea is to let the user select from as much as wide as
   *  possible list of cuts the selection criteria. For the time being
   *  it is particular focused on clusters (be careful still not hits,
   *  they are still in the detector frame of reference). 
   * 
   *  All cuts are declared on a detector basis. This means that if in
   *  the telescope there are 4 sensors and the user wants to cut on
   *  clusters having a total charge in excess than 100 ADC counts, so
   *  she/he has to provide the following in the steering file:
   *
   *  @code
   *  <parameter name="ClusterMinTotalCharge" type="IntVec"> 100 100 100 100 </parameter>
   *  @endcode
   *
   *  If instead she/he wants to cuts on the charge of a cluster but
   *  considering only the first N most significant pixels, the
   *  steering file should something like this:
   *   
   *  @code
   *  <parameter name="ClusterNMinCharge" type="IntVec"> 9 75 43 55 87 </parameter>
   *  @endcode
   *
   *  This last parameter will set a cut on the cluster charge
   *  (considering only 9 pixels) with a threshold of 75 ADC on the
   *  first sensor, 43 on the second and so on and so forth. Note that
   *  the first digit represents the number of significant pixels
   *
   *  The seed pixel can also be used as a cut, again in ADC value and
   *  one per sensor.
   *
   *  The user can select also upon the minimum and the maximum number
   *  of clusters per plane. So for example the following code will
   *  selects only events in each for each sensors there are between 1
   *  and 5 clusters.
   *  
   *  @code
   *  <parameter name="MinNoOfCluster" type="IntVec"> 1 1 1 1  </parameter>
   *  <parameter name="MaxNoOfCluster" type="IntVec"> 5 5 5 5  </parameter>
   *  @endcode
   *
   *  To switch off all those cuts the user can set them to zero or to
   *  a negative value. The system will then disable their
   *  functionality.
   *
   *  Another possible selection criterion is based on the cluster
   *  quality. In this case it is worth to remember that the cluster
   *  quality is defined via a specific enum ClusterQuality reported
   *  below:
   *
   *  @code
   *    enum ClusterQuality {
   *      kGoodCluster       = 0,
   *      kIncompleteCluster = 1L << 0,
   *      kBorderCluster     = 1L << 1,
   *      kMergedCluster     = 1L << 2
   *    };
   *  @endcode
   *
   *  This means that a good cluster corresponds to 0, an incomplete
   *  one to 1, etc. It is possible to sum up different cluster
   *  qualities, for example a merged cluster (4) that is also on the
   *  border (2) will have quality equal to 6.
   *  The idea behind is that the processor will select all clusters
   *  having the quality the users specify.

   *  Another selection can be done of the cluster position. For
   *  example the user may wants to keep only clusters having the
   *  central pixel within a certain range of pixels according to the
   *  following parametrization:
   *
   *  @code
   *  <parameter name="InsideRegion" type="IntVec"> detID xBotLeft  yBotLeft xTopRight yTopRight  </parameter> 
   *  <parameter name="OutsideRegion" type="IntVec"> detID xBotLeft  yBotLeft xTopRight yTopRight  </parameter> 
   *  @endcode
   *
   *  @param PulseCollectionName Name of the input cluster
   *  collection to be filtered
   *  
   *  @param FilteredPulseCollectionName Name of the filtered output
   *  collection. 
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelClusterFilter.h,v 1.1 2007-06-14 22:19:54 bulgheroni Exp $
   *
   *
   */

  class EUTelClusterFilter : public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTeleClusterSeparationProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTeleClusterSeparationProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelClusterFilter;
    }

    //! Default constructor 
    EUTelClusterFilter ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and reset all needed data
     *  members. 
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. 
     * 
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called every event
    /*  This is called for each event in the file. As a first thing,
     *  the system will check among all clusters found in the current
     *  event there are pairs of merging clusters on the same
     *  detector. If at least one pair of merging cluster is found,
     *  then the groupingMergingPairs(std::vector< pair<int, int> >,
     *  vector< set<int > > *) is called; otherwise it returns
     *  immediately.
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     *
     *  @throw UnknownDataTypeException if the cluster type stored in
     *  the TrackerPulse is unknown.
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

    //! This is the real method
    /*! This is the method where the real separation algorithm is
     *  written/called. 
     *
     *  @param setVector A STL vector of STL set containing the
     *  cluster of clusters to be separated.
     *
     *  @param collectionVec A pointer to the input cluster collection
     *
     *  @return true if the algorithm was successfully applied.
     */ 
    bool applySeparationAlgorithm(std::vector< std::set< int > > setVector, LCCollectionVec * collectionVec) const ;

    //! Groups together pairs of merging clusters.
    /*! The identification of merging clusters can very easily done on
     *  a pair basis. It means that after a double loop on clusters, a
     *  collection of pairs of merging clusters is available. To be
     *  more general, one should consider the possibility of having
     *  groups of merging clusters, a sort of cluster of
     *  clusters. This groupingMergingPairs method is exactly doing
     *  this, looking if there are pairs of merging clusters that can
     *  be grouped into a cluster of clusters.
     *  
     *  Of course this method is called if, and only if, at least a
     *  pair merging clusters have been found
     *
     *  @param pairVector A STL vector of STL pairs. For each pairs,
     *  the two integers represent the cluster indices within the
     *  clusterCollection
     *
     *  @param setVector A pointer to a STL vector of STL set. Each
     *  set has to be considered as the cluster of clusters defined
     *  above and the int value_type of the set represents the cluster
     *  index in the clusterCollection
     */ 
    void groupingMergingPairs(std::vector< std::pair<int , int> > pairVector, std::vector< std::set< int > > * setVector) const;

  protected:

    //! Input cluster collection name.
    /*! This is the name of the input cluster collection name
     */
    std::string _clusterCollectionName;

    //! Minimum distance
    /*! This is the minimum distance (in pixel units) allowed between
     *  two clusters. If their distance is below this number, then the
     *  two will be considered merging and the separation algorithm
     *  will be applied. If it is set to 0, only touching clusters
     *  will be considered merging.
     */
    float _minimumDistance;

    //! Separation algorithm name
    /*! This is the name of the algorithm used to divide merging
     *  clusters. @see EUTELESCOPE::FLAGONLY.
     */
    std::string _separationAlgo;

    //! Current run number.
    /*! This number is used to store the current run number
     */
    int _iRun;

    //! Current event number.
    /*! This number is used to store the current event number NOTE that
     * events are counted from 0 and on a run base
     */
    int _iEvt;

  };

  //! A global instance of the processor
  EUTelClusterFilter gEUTelClusterFilter;      

}
#endif
#endif

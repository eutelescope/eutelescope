// -*- C++ -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELCLUSTERSEPARATIONPROCESSOR_H
#define EUTELCLUSTERSEPARATIONPROCESSOR_H 1

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 

// system includes <>
#include <string>
#include <vector>
#include <set>

// forward declaration 


namespace eutelescope {

  //! Cluster separation processor.
  /*! This processor can be used to separate clusters touching or
   *  found so close that charge sharing cannot be avoided or
   *  excluded. Of course only clusters belonging to the same detector
   *  level, i.e. having the same detectorID are considered.
   *
   *  This is in principle not an easy task, because the way the two
   *  touching clusters are sharing the charge is not obvious and
   *  particular care should be used when dealing with this sort of
   *  problem.
   *
   *  First of all a criterium to define when two clusters are
   *  touching is needed. The safest way would be to scan each cluster
   *  and cross-check if a pixel is belonging to both clusters, or it
   *  is adjacent to a pixel belonging to another cluster. This is a
   *  very expensive procedure because it needs lots of
   *  operations. This processor is based on the following assumption:
   *  for each fixed frame cluster we can define an external circle
   *  having the following radius: <code> r = 0.5 * sqrt( pow(n,2) +
   *  pow(m.2) ) </code> where @c n and @c m are the cluster sizes in
   *  pixel units. We will consider as touching or overlapping, two
   *  clusters with a distance between the seed pixels lesser of equal
   *  to the sum of the two radii.
   *
   *  Here comes a list of all implemented
   *  separation algorithm:
   *
   *  \li <b>FlagOnly</b>: actually this is not a separation
   *  algorithm, because nothing will be done trying to separate the
   *  two components. The only effective result is all clusters having
   *  a distance below the minimum allowed one are flagged
   *  (eutelescope::ClusterQuality::kMergedCluster) so to exclude them
   *  from following analysis processors.
   *
   *  <h4>Input - Prerequisites</h4>
   *  <br><b>ClusterCollection</b>. This is the collection of cluster
   *  to be used. <br><b>MinimumDistance</b>. This floating point
   *  parameter is used to set which is the minimum distance between
   *  two separeted clusters. It is measured in pixel unit. If it set
   *  to 0, only touching clusters are considered
   *  merged. <br><b>SeparationAlgorithm</b>. This string is used to
   *  define different way to separate merging clusters.
   *  
   *  <h4>Output</h4>
   *  
   *  <br><b>ClusterCollection</b>. It outputs the same collection but
   *  now with separated pixels.
   *
   *  @param ClusterCollectionName Name of the input cluster
   *  collection 
   *  @param MinimumDistance Float number being the minimum
   *  distance allowed between clusters. 
   *  @param SeparationAlgorithm Name the algorithm to be used for
   *  cluster separation.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelClusterSeparationProcessor.h,v 1.1 2007-02-26 09:32:10 bulgheroni Exp $
   *
   *
   */

  class EUTelClusterSeparationProcessor : public marlin::Processor {

  public:

     
    //! Returns a new instance of EUTeleClusterSeparationProcessor
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTeleClusterSeparationProcessor.
     */
    virtual Processor * newProcessor() {
      return new EUTelClusterSeparationProcessor;
    }

    //! Default constructor 
    EUTelClusterSeparationProcessor ();

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
    /*  This is called for each event in the file. In the case this is
     *  the first event, then a cross-check is done on the
     *  compatibility of input raw data and pedestal collections. They
     *  have to contain the same number of detector planes and each of
     *  them should have exactly the same number of pixels. Of course
     *  this check is not exhaustive, but at least should avoid
     *  segmentation fault.
     *
     *  @throw IncompatibleDataSetException in the case the two
     *  collections are found to be incompatible
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

    //! Separation algorithm implementation
    /*! This function is called with the processEvent(LCEvent*) call
     *  back when two clusters have been found to be merging.
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
  EUTelClusterSeparationProcessor gEUTelClusterSeparationProcessor;      

}
#endif

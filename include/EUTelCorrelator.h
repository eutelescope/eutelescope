// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
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

// eutelescope includes ".h" 

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h> 
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>

// AIDA includes <.h>
#ifdef MARLIN_USE_AIDA
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram2D.h>
#endif


// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {

  //! Hit and cluster correlator
  /*! Please write here the description ... ;-)
   *
   *  
   *  <h4>Input collections</h4> 
   *  
   *  <b>Cluster collection</b>: A collection with cluster
   *  produced by previous processors like EUTelClusteringProcessor or
   *  EUTelClusterFilter. 
   *
   *
   *  @author Silvia Bonfanti, Uni. Insubria  <mailto:silviafisica@gmail.com>
   *  @author Loretta Negrini, Uni. Insubria  <mailto:loryneg@gmail.com>
   *  @version $Id: EUTelCorrelator.h,v 1.1 2008-07-11 15:35:30 bulgheroni Exp $
   *
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
    virtual Processor * newProcessor() {
      return new EUTelCorrelator;
    }

    //! Default constructor 
    EUTelCorrelator ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters
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
    virtual void processEvent (LCEvent * evt);


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


  protected:
    
    //! Input cluster collection name
    /*! This is the name of the collection containing the input clusters
     */
    std::string _inputCollectionName;
    
  private:
    
    //! Run number 
    int _iRun;
    
    //! Event number
    int _iEvt;

    //! Number of sensors
    int _noOfDetectors;

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


#ifdef MARLIN_USE_AIDA

    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names. 
     *  The histogram filling can proceed recalling an object through
     *  its name
     */ 
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    //! Correlation histogram matrix
    /*! This is used to store the pointers of each histogram
     */
    std::vector< std::vector< AIDA::IHistogram2D *  > > _clusterXCorrelationMatrix; 

    //! Base name of the correlation histogram
    static std::string _clusterXCorrelationHistoName;
  
#endif

 

  };

  //! A global instance of the processor
  EUTelCorrelator gEUTelCorrelator;      

}
#endif

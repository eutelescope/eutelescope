/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPROCESSORNOISYCLUSTERREMOVER_H
#define EUTELPROCESSORNOISYCLUSTERREMOVER_H

// eutelescope includes ".h"
#include "EUTelEventImpl.h"
#include "EUTelSimpleSparsePixel.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/LCCollectionVec.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes


namespace eutelescope {

  //! Processor to convert data to be compliant with EUTelGenericSparsePixel
  /*! EUTelescope stores data either in EUTelSimpleSparsePixel with 
   *  information (X,Y,signal) per hit or EUTelAPIXSparsePixel with
   *  (X,Y,signal, chip, time). 
   *  New processors should not have to deal with those different types,
   *  therefore there is the EUTelGenericSparsePixel which stores information
   *  (X,Y,signal,time).
   *  This processor takes zsdata in either of the above pixels and 
   *  stores them compliant to be interfaced with EUTelGenericSparsePixel
   *  in a TrackerDataImpl collection.
   *  If the original collection should be dumped, specify so in the 
   *  EUTelOutputProcessor ("Save")
   */

class EUTelProcessorNoisyClusterRemover : public marlin::Processor {

  public:

    //! Returns a new instance of EUTelProcessorNoisyClusterRemover
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelHotPixelKiller.
     */
    virtual Processor * newProcessor() {
      return new EUTelProcessorNoisyClusterRemover;
    }

    //! Default constructor
    EUTelProcessorNoisyClusterRemover();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and performs some asserts about
     *  the value of the provided parameters
     */
    virtual void init();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented.
     *
     *  @param run LCRunHeader of the this current run
     *
     *  @throw InvalidParameterException if a paramter is wrongly set
     */
    virtual void processRunHeader(LCRunHeader * run);

    //! Called every event
    /*! This is called for each event in the file. If the current @c
     *  evt is flagged to be used for update, then the selected
     *  algorithm wrapper is called
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     *
     *  @throw InvalidParameterException if information in the cellID
     *  are inconsistence
     */
    virtual void processEvent(LCEvent * evt);

    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished. Just printing a good bye message
     */
    virtual void end();


  protected:
	//! Input collection name for data	
	std::string _inputCollectionName;

	//! Output collection name for noise free data
	std::string _outputCollectionName; 
	
	//! Integer to track the collection size
	/*! Used to monitor if anything changed and the collection needs
	 *! to be updated 
         */
	unsigned int _initialOutputCollectionSize;

	//! Current run number.
	/*! This number is used to store the current run number */
	int _iRun;

	//! Current event number.
	/*! This number is used to store the current event number NOTE that
	 * events are counted from 0 and on a run base
         */
	int _iEvt;
};

//! A global instance of the processor
EUTelProcessorNoisyClusterRemover gEUTelProcessorNoisyClusterRemover;

}//namespace eutelescope

#endif //EUTelProcessorNoisyClusterRemover_H

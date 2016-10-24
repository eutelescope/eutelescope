/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPROCESSORNOISYPIXELREMOVER_H
#define EUTELPROCESSORNOISYPIXELREMOVER_H

// eutelescope includes ".h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <LCIOTypes.h>
#include <IMPL/LCCollectionVec.h>

// AIDA includes <.h>
//#include <AIDA/IBaseHistogram.h>

// system includes


namespace eutelescope {

  //! Processor to remove noisy clusters from a TrackerPulse collection
  /*! This processor clones the input TrackerPulse collection, but checks
   *  for any clusters if their quality flag: kNoisyCluster is set.
   *  If so, they are removed from the output collection. 
   */

class EUTelProcessorNoisyPixelRemover : public marlin::Processor {

  public:

    //! Returns a new instance of EUTelProcessorNoisyPixelRemover
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelProcessorNoisyPixelRemover.
     */
    virtual Processor* newProcessor() {
      return new EUTelProcessorNoisyPixelRemover;
    }

    //! Default constructor
    EUTelProcessorNoisyPixelRemover();

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
	//! Input collection name of the TrackerPulse collection	
	std::string _inputCollectionName;

	//! Output collection name for noise free pulses
	std::string _outputCollectionName; 
	
	//! Collection name for noisy pixel collection
	std::string _noisyPixelCollectionName; 
	
	std::map<int, std::vector<int>> _noisyPixelMap;
	bool _firstEvent = true;
};

//! A global instance of the processor
EUTelProcessorNoisyPixelRemover gEUTelProcessorNoisyPixelRemover;

}//namespace eutelescope

#endif //EUTelProcessorNoisyPixelRemover_H

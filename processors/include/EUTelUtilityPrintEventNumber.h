#ifndef EUTelUtilityPrintEventNumber_h
#define EUTelUtilityPrintEventNumber_h 1

// C++
#include <string>

// LCIO
#include "lcio.h"

// Marlin
#include "marlin/Processor.h"


namespace eutelescope {

  /**  Prints event numbers to the screen depending on the EveryNEvents variable.
   *
   *   @parameter EveryNEvents Print event number for every n-th event
   *
   *   @parameter printTimestamp Print the events timestamp as read from LCIO
   * 
   */
  class EUTelUtilityPrintEventNumber : public marlin::Processor {
	
  public:
	
    /* This method will be called by the marlin package
     * It returns a processor of the currend type
     */
    virtual Processor*  newProcessor() { 
      return new EUTelUtilityPrintEventNumber;
    }

    virtual const std::string & name() const { return Processor::name() ; }
	
    /* the default constructor
     * here the processor parameters are registered to the marlin package
     * other initialisation should be placed in the init method
     */
    EUTelUtilityPrintEventNumber() ;
	
    /* Called at the beginning of the job before anything is read.
     * Use to initialize the processor, e.g. book histograms
     */
    virtual void init() ;
	
    /* Called for every run.
     * in this processor it is not used
     */
    virtual void processRunHeader( lcio::LCRunHeader* run ) ;
	
    /* Called for every event - the working horse.
     * Here the real conversion take place
     */
    virtual void processEvent( lcio::LCEvent * evt ) ; 
	
    /* This method is only called if the check flag is set (default) in the main processor
     * (where the input files are given)
     * It should be used to create check plots which are not needed in a complete analysis
     * (where normaly the check flag is deactivated)
     */
    virtual void check( lcio::LCEvent * evt ) ; 
	
	
    /* Called after data processing for clean up.
     * e.g. saving control histogramms, delete created objects
     * (which are not stored in a collection/event), etc.
     */
    virtual void end() ;
	
	
  protected:
	
    /// variable how often to print the event number
    int _everyNEvents;

    // processor parameter to print the event timestamp as read from the LCIO file
    bool _printTimestamp;
    
    // count total events in case runs are concatenated:
    unsigned int totalevents;
    // count the total number of runs processed
    unsigned int totalruns;

  };
  
  //! A global instance of the processor
  EUTelUtilityPrintEventNumber gEUTelUtilityPrintEventNumber;

}

#endif
    
    


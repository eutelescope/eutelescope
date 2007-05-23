// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */


#ifndef EUTELOUTPUTPROCESSOR_H
#define EUTELOUTPUTPROCESSOR_H 1

// eutelescope includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/LCIOOutputProcessor.h"

// lcio includes <.h>
#include <lcio.h>
#include <IO/LCWriter.h>


namespace eutelescope {

  //! EUTelescope specific output processor
  /*! This processor is just a copy of the standard LCIOOutputProcess
   *  available with the standard version of Marlin, but it is used
   *  within the EUDET telescope to eventually add an End Of Run Event
   *  at the end.
   *
   *  Within the EUTelescope framework there are several processors
   *  actually needing to know which is the last event of a run,
   *  mainly because they need to rewind the data stream to the
   *  beginning. Since SIO is a serial I/O format, an easy workaround
   *  is to append at the end of every run an empty event tagged as
   *  the last one. 
   *
   *  The use of this processor is very helpful when for example an
   *  input LCIO file has been produced without the EORE and this
   *  should be appended in order to compatible with the full analysis
   *  chain. 
   *
   *  Another occasion in which the use of this processor is crucial
   *  is when the analysis of the input file is limited to a certain
   *  small range by the user with the MaxRecordNumber global
   *  parameter. Say for example that the input file have 1000 events
   *  (being the last one a EORE), if the analysis is limited to
   *  MaxRecordNumber = 100, the last processed event will be a DE and
   *  not the EORE. Consequently in the saved output file, the last
   *  event will not a EORE breaking the consistency of the analysis
   *  chain. The use of EUTelOutputProcessor instead of
   *  LCIOOutputProcessor allows to fix this problem because, if the
   *  last processed event was not a EORE, then a EORE is appended
   *  before closing the output file.
   *
   *  @see marlin::LCIOOutputProcessor
   *  @see eutelescope::EventType
   *  @see eutelescope::EUTelEventImpl
   *
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelOutputProcessor.h,v 1.1 2007-05-23 14:10:05 bulgheroni Exp $ 
   */

  class EUTelOutputProcessor : public marlin::LCIOOutputProcessor {
  
  public:  

    //*! Returns a new instance of EUTelOutputProcessor
    /*!  This method returns a new instance of this processor. It is
     *   called by Marlin execution framework and it should not be
     *   called/used by the final user.
     *
     *   @return a new EUTelOutputProcessor.
     */ 
    virtual Processor*  newProcessor() { 
      return new EUTelOutputProcessor ;
    }
    
    //! Default constructor
    EUTelOutputProcessor() ;

    //! Opening the LCIO output file 
    /*! This is executed only once in the whole execution. It opens
     *  the output file with the mode flag specified in the steering
     *  files and prints out the processor parameters.
     */
    virtual void init() ;

    //! Process the run header
    /*! This method processes each run header copying them in the output
     *  file.
     *
     *  @param run The LCRunHeader being processed.
     */
    virtual void processRunHeader( LCRunHeader* run) ;

    //! Process the event
    /*! This method processes the current event, removing the dropped
     *  collections and saving the other on disk. Before returning the
     *  current event type is stored in _eventType
     *
     *  @param evt The LCEvent being processed.
     */
    virtual void processEvent( LCEvent * evt ) ; 

    //! Close the output file
    /*! This is the method where the output file is closed. In the
     *  case the last processed event was 
     */
    virtual void end() ;


  protected:
    
    //! The current event type
    /*! This is actually the only reason for reimplementing the
     *  LCIOOutputProcessor. During the processEvent(LCEvent * evt) we
     *  check the current event type and store it in this
     *  variable. When the end() method is called, then accessing to
     *  _eventType we know if the previously processed event was an
     *  EORE or not. If not a new EORE is added.
     */ 
    EventType _eventType;

  } ;

  //! A global instance of EUTelOutputProcessor
  EUTelOutputProcessor  gEUTelOutputProcessor ;


} // end namespace eutelescope
#endif




// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELEVENTVIEWER_H
#define EUTELEVENTVIEWER_H 1

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <lcio.h>

// system includes <>
#include <vector>
#include <string>

namespace eutelescope {
  //! EUTelEventViewer
  /*! This processor is used to display the content of a EUTelescope
   *  file into a graphical window using the information stored in
   *  GEAR. This is still very experimental!!!
   *
   *  @param TrackerHitCollectionNameVec. This is a vector of string
   *  containing the names of all the TrackerHitCollection one would
   *  like to display
   *
   *  @param TrackCollectionNameVec. This is a vector of string
   *  containin the names of all the TrackCollection one would like to
   *  display.
   *
   *  @param LayerTrackerHit. The CED layer where hits are drawn. Set
   *  to -1 to disable it
   *
   *  @param LayerTrack. The CED layer where tracks are drawn.Set
   *  to -1 to disable it
   *  
   *  @param WaitForKeyboard. To stop the event loop at the end of
   *  each display waiting for any key from the keyboard to continue.
   *
   *  @param DetModel. This is the ID used to choose the proper
   *  geometry. For the telescope use 99999. The geometry information
   *  are then read from GEAR.
   *
   *  @Todo The universe described in the openGL window is well suited
   *  for describing a calorimeter that is much bigger than our small
   *  telescope. Drawing our telescope in the openGL window results in
   *  a so small telescope that can only hardly seen. The workaround
   *  currently used is based on a global parameter from the steering
   *  file. This is read by both the MarlinCED and by this
   *  processor. If the parameter is missing, the identity will be used.
   *
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id: EUTelEventViewer.h,v 1.1 2007-06-25 16:14:06 bulgheroni Exp $ 
   */
  class EUTelEventViewer : public marlin::Processor {
  
  public:
    
    //! Create a new processor instance
    /*! @return a new processor instance
     */ 
    virtual Processor*  newProcessor() { return new EUTelEventViewer ; }
    
    //! Default constructor;
    EUTelEventViewer() ;
    
    //! Init method;
    /*! It prints the parameters and initialize the MarlinCED
     *  interface.
     */
    virtual void init() ;
    
    //! Process the run header
    /*! Nothing to do.
     * 
     *  @param run The current run header
     */ 
    virtual void processRunHeader( LCRunHeader* run ) ;
    
    //! Process the current event
    /*! This is the real part of the class. The current event is taken
     *  and in case a EORE is found, the processor will immediately
     *  return. Otherwise, the current detector geometry is drawn in
     *  the graphical window. 
     *  After that each collections the user wants to display are
     *  drawn one by one on the corresponding layers.
     *
     *  In case the user wants to pause the processing, it is possible
     *  to use the wait for keyboard switch.
     *
     *  @param evt The current event
     */ 
    virtual void processEvent( LCEvent * evt ) ; 
    
    //! Check event
    /*! Nothing to do here
     *  
     *  @param evt The current event
     */ 
    virtual void check( LCEvent * evt ) ; 
  
    //! Finish processing
    /*! Nothing to do here
     */
    virtual void end() ;

    //! Find a color
    /*! This method can be used to get the proper color for the
     *  current collection.
     *
     *  @param counter Collection number
     *  @return A proper color value
     */ 
    int returnColor(int counter);
    
  protected:
    
    //! Tracker Hit collection names
    /*! This is the vector of tracker hit collection names. The user
     *  via the steering file can add as many collection she/he wish
     *  and they will be displayed using different colors. The reason
     *  for that is the possibility to compare different tracking /
     *  clustering algorithms.
     *
     */ 
    std::vector<std::string> _trackerHitCollectionNameVec;

    //! Track collection names
    /*! This is the vector of track collection names. The user
     *  via the steering file can add as many collection she/he wish
     *  and they will be displayed using different colors. The reason
     *  for that is the possibility to compare different tracking /
     *  clustering algorithms.
     *
     */ 
    std::vector<std::string> _trackCollectionNameVec;
    
    //! CED TrackerHit layer
    /*! This number corresponds to the CED layer where tracker hit
     *  will be drawn. Setting it to a negative number means switching
     *  it of.
     */ 
    int _layerTrackerHit;

    //! CED Track layer
    /*! This number corresponds to the CED layer where track will be
     *  drawn. Setting it to a negative number means switching it of.
     */ 
    int _layerTrack;

    //! Wait for keyboard switch
    /*! Set it to 1 to pause the event processing at the end of each
     *  event
     */
    int _waitForKeyboard; 

    //! The detector model
    int _detModel;
    

  } ;
  
}
#endif




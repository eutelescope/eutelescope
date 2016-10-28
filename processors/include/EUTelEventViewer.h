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
  //! Event viewer
  /*! This processor is used to display the content of a EUTelescope
   *  file into a graphical window using the information stored in
   *  GEAR.
   *
   *  The graphical environment is provided by CED (C Event Display)
   *  that is a server / client application. The server side that is
   *  implemented using openGL (glced) should be already running on
   *  the computer, before Marlin is executed.
   *
   *  This viewer is accessing MarlinCED, a singleton class playing
   *  the role of CED client, and forwarding there all the graphical
   *  objects. For every event, the graphical screen is reset and the
   *  geometry re-drawn along with the current hits and tracks.
   *
   *  The viewer can be used in a interactive way setting
   *  WaitForKeyboard to 1. In this way Marlin execution will be
   *  stopped as soon as the event is drawn waiting for the user to
   *  hit return.
   *
   *  Several objects can be displayed together, since for each
   *  collection type multiple selections are possible.
   *
   *  Starting from EUTelescope version 0.0.7 also tracks can be
   *  displayed as a spline connecting the fit hit. 
   *
   *  This processor is built only if GEAR and CED are available
   *  (-DUSE_GEAR and -DUSE_CED).
   *
   *  @image html CEDEvent.png "One event with a track candidate"
   *  @image html CEDTrack.png "One event with a reconstructed track"
   *
   *  <h4>Input collections</h4>
   *
   *  <b>Track collections</b>: this is a vector of all the track
   *  collections the user wants to display.
   *
   *  <b>Tracker hit collections</b>: this is a vector of all the
   *  tracker hit collections the user wants to display.
   *
   *  <b>Alignment Constants</b>: this is a LCGenericObject collection
   *  containing the results of the alignment procedure. This is just
   *  an optional collection and it is required only if the user
   *  wants to plot the real telescope geometry. See also the DetModel
   *  parameter. 
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
   *  @param AutoForwardDelay. In case you set the WaitForKeyBoard to
   *  false, then next event will be loaded automatically and the user
   *  can choose a delay in second in between.
   *
   *  @param DetModel. This is the ID used to choose the proper
   *  geometry. For the telescope use 99999. The geometry information
   *  are then read from GEAR. If the users wants to plot the planes
   *  taking into account the alignment corrections, than detector
   *  model has to be set to -1 and a suitable alignment collection
   *  has to be provided.
   *
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
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

    //! Alignment constant collection
    /*! This is the collection containing the alignment
     *  constants. Those numbers are applied to the planes so that
     *  they will show up in the event display exactly as they are in
     *  the reality. This is done only if _applyAlignmentToPlane is
     *  set to true.
     *
     *  When the planes are shifted, then the hits have to be
     *  corrected as well.
     */
    std::string _alignmentCollectionName;

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
    /*! Set it to true to pause the event processing at the end of each
     *  event
     */
    bool _waitForKeyboard;

    //! Auto forward delay
    /*! This is the time in second between two automatic event
     *  display. 
     *  This value is used only if _waitForKeyboard is false and has
     *  to be set to zero to go at full speed.
     */
    double _autoForwardDelay;
      
    //! The detector model
    int _detModel;

  } ;

  //! A global instance of the processor.
  EUTelEventViewer gEUTelEventViewer ;

}
#endif




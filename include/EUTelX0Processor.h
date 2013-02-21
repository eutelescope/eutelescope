// Version: $Id$
// Contact: Phillip Hamnett (phillip.hamnett@desy.de)
/*
 * This source code is part of the Eutelescope package of Marlin.
 * You are free to use this source files for your own development as
 * long as it stays in a public research context. You are not
 * allowed to use it for commercial purpose. You must put this
 * header with author names in all development based on this file.
 */

#ifndef EUTELX0PROCESSOR_H
#define EUTELX0PROCESSOR_H

//  Gear includes 
#ifdef USE_GEAR
#include <gear/SiPlanesLayerLayout.h>
#include <gear/SiPlanesParameters.h>
#endif
//  EUTelescope includes
#include "EUTelReferenceHit.h"
//  Marlin includes
#include <marlin/AIDAProcessor.h>
#include "marlin/Processor.h"
//  AIDA includes
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IFitResult.h>
#include <AIDA/IFitter.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IPlotter.h>
#include <AIDA/IPlotterFactory.h>
#include <AIDA/IPlotterRegion.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#include <RAIDA/IAnalysisFactoryROOT.h>
#include <RAIDA/IPlotterFactoryROOT.h>
#include <RAIDA/IPlotterROOT.h>
#endif
//  LCIO includes
#include <EVENT/LCEvent.h>
#include <EVENT/LCIO.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include "IMPL/TrackerHitImpl.h"
//  System includes
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
//  ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TVector3.h"

//  Forward Declarations
class AIDA::IAnalysisFactory;
class AIDA::IFitter;
class AIDA::IFitResult;
class AIDA::IPlotter;
class AIDA::IPlotterFactory;
class AIDA::IPlotterRegion;
class EVENT::TrackerHit;
class IMPL::LCCollectionVec;
class TMinuit;

namespace eutelescope {
//!EUTelX0Processor Class
/*!This class is used to work out the radiation length of a Device Under Test (DUT). What follows will be a list of each function and a brief description of its functionality, as well as each variable and its purpose:
*Public functions:
*-newProcessor() is necessary in all Marlin processors and is not used by the rest of the class.
*-EUTelX0Processor() is the default constructor, it does nothing and is not used, instead we use init() to initialise the class
*-init() is used to initialise the class, this fills in all the appropriate parameters from the Marlin steering file, and initialises all internal variables to 0 or a sensible default value
*-processRunHead(LCRunHeader *run) is called for every run, currently does nothing in this class
*-processEvent(LCEvent *evt) is called for every event in a run, and this currently deduces a track (if one exists) and records this tracks residual. This is used later to work out a kink angle which is necessary for the calculation of a radiation length.
*-end() is called after all the events are finished running, and this cleans up all the variables, before calculating the value of the material budget (X0) from all the events.
*-booksHistos() fills all the relevant histograms each event.
*-guessSensorID(const double *pos) is used to deduce which layer of the pixel telescope each hit was in, and it uses the GEAR file provided to do this (as well as the parameter refhit)
*Private functions (Should be accessed only from within this class):
*-basicFitter(LCCollection *alignedHitCollection) fits tracks based on relative X and Y position between two layers, if they are within a radius of each other determined by _cutValue1 then they are considered a track
*-calculateX0() is called in end() and is used to work out the value for the average radiation length in the material, this is the most important function in this processor
*-testtrack(LCCollection *trackCollection) is used if a fitting processor is used before teh X0 Processor and therefore no fitting is done in X0. This function just prints out all the values from the tracks as a test that they are working 
*-threePointResolution(LCCollection *alignedHitCollection) is a simple fitting algorithm that uses planes i and i+2 and takes the average x and y values to compare residuals against the actual values in plane i+1
*-createResiduals(LCCollection *trackCollection) this function creates the residuals from the basicFitter function
*/
class EUTelX0Processor : public marlin::Processor {

public:
  //! Returns a new instance of EUTelX0Processor
  /*! This method returns a new instance of this processor. It is
  *   called by Marlin execution framework and it shouldn't be
  *   called/used by the final user.
  *
  *   @return a new EUTelX0Processor.
  */
  virtual Processor * newProcessor() {return new EUTelX0Processor;}

  //! Default constructor sets all numeric values to -1, all pointers to NULL, all booleans to false, clears all STL containers and sets strings to ""
  EUTelX0Processor ();

  //! Called at the job beginning.
  /*! This is executed only once in the whole execution. It inputs
  *   the processor parameters and check that the GEAR
  *   environment is properly set up and accessible from Marlin.
  *   All the private members of the class are given their starting values here
  */
  virtual void init ();

  //! Called for every run.
  /*! It is called for every run, but does nothing in this code. 
  *   It must be left in the program due to the inheritence from Marlin
  *   @param run the LCRunHeader of the current run
  */
  virtual void processRunHeader (LCRunHeader * run);
  
  //! Called every event
  /*! This is called for each event in the file.
  *   This will, when completed, go through each event and 
  *   find the tracks (or accept the track from a track finding 
  *   processor). Then it will work out the scattering angle of
  *   of the track at the position of the DUT. Finally, when it
  *   has all the scattering angles it will use the values to
  *   produce an estimate of the radiation length of the material.
  *   Additionally, it may be possible to do a fine binning of how
  *   the radiation length changes with x and y position in the DUT.
  *   But this is for a later stage of the development.
  *   @param evt the current LCEvent event as passed by the
  */
  virtual void processEvent (LCEvent * evt);

  //! Called after data processing.
  /*! This method is called when the loop on events is
  *   finished. It cleans up code, and sets all values
  *   as per the constructor.
  */
  virtual void end();

  //! Histogram booking
  /*! Some control histograms are filled during this procedure in
  *   order to be able to perform easy check on the quality of the
  *   output hits and also to understand if the frame of reference
  *   conversion has been properly done. Of course this method is
  *   effectively doing something only in the case MARLIN_USE_AIDA.
  *   Currently this is not used, but that may change later.
  */
  void bookHistos();

  //!Guess Sensor ID
  /*!Works out which layer of the telescope the current element is in based on its position*/ 
  int guessSensorID(const double * hit );

private:
  #ifndef DISALLOW_COPY_AND_ASSIGN
  //Following #define stops the accidental creation of a copy or assignment operator by causing a link error. Copy and Assignment operators not allowed because they are unnecessary and the cause of many bugs
  #define DISALLOW_COPY_AND_ASSIGN(EUTelX0Processor) \
  EUTelX0Processor(const EUTelX0Processor&); \
  void operator=(const EUTelX0Processor&);
  
  //Private Functions
  DISALLOW_COPY_AND_ASSIGN(EUTelX0Processor)//See #define just above
  #endif
  
  //!Cut and Store Hits
  /*!This looks at the two layers and extrapolates a track line between them. It also allows for a cut to be made that stops tracks with angles greater than +/- anglecut.
  size1 and size2 refer to the specific sizes of the number of hits in the first and second layer.
  firstLayer and secondLayer inform the user in which two layers the extrapolation is being made. It is important that you get the order the right way around for the extrapolation. For instance, if going from layer 0 and 1 to layer 2 then the firstLayer = 0 and secondLayer = 1. But if going from layers 3 and 4 to 2 then the extrapolation is in the opposite direction and so firstLayer = 4 and secondLayer = 3.
  */
  void basicFitter(LCCollection *alignedHitCollection);
  //!Caculate X0
  /*!This function calculates the value of the material budget.
  */
  double calculateX0();
  void printTrackParameters(EVENT::Track* track);
  void printHitParameters(EVENT::TrackerHit* hit);
  void testtrack(LCCollection *trackCollection);
  void threePointResolution(EVENT::Track *track);
  void createResiduals(LCCollection *trackCollection);
  void kinkEstimate(EVENT::Track* track);

  //Private member values
  std::string _trackColName;
  double _cutValue1,_cutValue2;
  bool _debug;
  int _debugCount;  //TODO(Phillip Hamnett): Can this be deleted?
  int _eventNumber;
  std::map<std::string , std::vector< double > > _histoData;
  std::map<std::string , AIDA::IBaseHistogram * > _histoThing;  //This map contains all the histograms that are going to be plotted in this processor
  static std::string _histoResidualX;
  static std::string _histoResidualXPlane1;
  static std::string _histoResidualXPlane2;
  static std::string _histoResidualXPlane3;
  static std::string _histoResidualXPlane4;
  static std::string _histoResidualY;
  static std::string _histoResidualYPlane1;
  static std::string _histoResidualYPlane2;
  static std::string _histoResidualYPlane3;
  static std::string _histoResidualYPlane4;
  static std::string _histoResidualXY;
  static std::string _histoResidualXZ;
  static std::string _histoResidualYZ;
  std::map< int, std::vector< TVector3 > > _hitInfo;  //This stores the hit position in a TVector3. If there are multiple hits then they are all stored in the vector of TVector3's. The int key refers to the layer number
  IMPL::LCCollectionVec* _inputHitCollectionVec;  //Stores the information being brought in from the Marlin process which contains information about post-aligned hits
  IMPL::LCCollectionVec* _inputTrackCollectionVec;  //Stores the information being brought in from the Marlin process which contains information about post-aligned hits
  std::string _inputHitColName;
  std::string _inputHitCollectionName;  //Stores the name of the parameter to bring in from the Marlin steering file
  std::string _inputTrackColName;
  static const int _noLayers = 6;  //TODO(Phillip Hamnett): This is a hack just so that the computer knows how many layers there are, is this stored somewhere else or should it be made as a non-const variable in case at some future date more layers are added to the telescope?
  std::map< int, std::vector< TVector3 > > _projectedHits;  //int refers to 1 or 2, with 1 being the projection from the 01 plane and 2 being the projection from the 43 plane
  std::string _referenceHitCollectionName;  //Stores the name of the file to be brought in by the Marlin steering file that is used to determine layer number
  IMPL::LCCollectionVec* _referenceHitVec;   //Stores the information being brought in from the Marlin process containing information about the geometry of the setup for working out the layer number
  std::map<std::pair< Double_t, Double_t > , std::vector< TVector3 > > _residual;  //The pair of doubles refers to the x and y coordinates in layer 2. The vector of TVector3's refers to the positions of the projected hits
  std::map<std::pair< Double_t, Double_t > , std::vector< TVector3 > > _residualAngle;  //As above but the TVector3 doesn't contain xyz but instead theta, phi and alpha
  double _residualCut;  //This is the cut which determines the acceptable tracks
  std::map<std::pair< Double_t, Double_t > , std::vector< TVector3 > > _residualProfile; //TODO(Phillip Hamnett): Can this be joined with _residual? //Used as above but for created a profile histogram
  int _runNumber;
  std::string _trackCollectionName;
};
//! A global instance of the processor
EUTelX0Processor gEUTelX0Processor;
}
#endif
#endif

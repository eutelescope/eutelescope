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
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
//  ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TVector3.h"

//  Forward Declarations
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

private:
  #ifndef DISALLOW_COPY_AND_ASSIGN
  //Following #define stops the accidental creation of a copy or assignment operator by causing a link error. Copy and Assignment operators not allowed because they are unnecessary and the cause of many bugs
  #define DISALLOW_COPY_AND_ASSIGN(EUTelX0Processor) \
  EUTelX0Processor(const EUTelX0Processor&); \
  void operator=(const EUTelX0Processor&);
  DISALLOW_COPY_AND_ASSIGN(EUTelX0Processor)//See #define just above
  #endif
  
  //!Get the hits from the tracks
  /*!We need to get the hit positions from each track in order to form our own 'subtracks' between each plane and arm of the telescope. This will return a vector of TVector3's with each element of the vector being a hit*/
  std::vector< TVector3* > getHitsFromTrack(Track *track);

  //!Get Sigma
  /*!Works out the sigma value from a vector which contains scattering angles*/
  double getSigma(std::vector< double > angles);

  //!Print Track Parameters
  /*!This simply prints out all the LCIO information available about each track*/
  void printTrackParameters(EVENT::Track *track);

  //!Print Hit Parameters
  /*!This prints out the information about each hit from each track*/
  void printHitParameters(EVENT::TrackerHit *hit);

  //!Fills the Single Point Residual Plots
  /*!This fills plots with the residuals based on a straight line extrapolation between two planes, i.e. a straight line is drawn between planes 0 and 1, and then a histogram is filled with the difference between where that straight line hits plane 2 and where the actual hit is on plane 2*/
  void singlePointResolution(EVENT::Track *track);

  //!Fills the Triple Point Residual Plots
  /*!A straight line of best fit is drawn between the front three planes and the back three planes, and the difference between the track position at those points on the plane and the actual hits at those points is plotted as a residual*/
  void threePointResolution(EVENT::Track *track);

  //!Get Single Track Angles
  /*!This works out the angles of the tracks as they pass from plane to plane, in the pair of vectors which is returned from this you get the angles of the track in x and y respectively. The element of each vector is the plane that the angle passes between, e.g. element 0 contains the angle between plane 0 and 1*/
  std::pair< std::vector< double >, std::vector< double > > GetSingleTrackAngles(std::vector< TVector3* > hits);

  //!Get Triple Track Angles
  /*!This works out the angles of the tracks for the first three planes combined and the last three planes combined and stores them in a pair. With the first element of the pair being XZ angle and the second element being YZ angle. Each element of the vector is a different triplet, so each vector should just contain 2 elements, the first being the front three planes and the 2nd being the last three planes*/
  std::pair< std::vector< double >, std::vector< double > > GetTripleTrackAnglesStraightLines(std::vector< TVector3* > hits);
  std::pair< std::vector< double >, std::vector< double > > GetTripleTrackAnglesDoubleDafFitted(Track *frontthree, Track *backthree);
  void PlotTripleTrackAngleDoubleDafFitted(Track *track, bool front);

  //!Single Plane Track Scattering Angles
  /*!This works out the scattering angle between each plane, fills a histogram with this value and also stores data for use in a radiation length map later*/
  void SinglePlaneTrackScatteringAngles(std::vector< double > scatterx, std::vector< double > scattery, std::vector< TVector3* > hits);

  //!Triple Plane Track Scattering Angles
  /*!This works out the scattering angle between the front three planes and the back three planes and fills a histogram with the result. This also stores the data for use in a radiation length map later*/
  void TriplePlaneTrackScatteringAngles(std::vector< double > scatterx, std::vector< double > scattery, std::vector< TVector3 *> hits, bool doubledaf);

  //!Kink Estimate
  /*!This is the funciton which calls most of the other relevant functions in order to fill angle-related histograms*/
  void kinkEstimate(EVENT::Track *track);

  //!Conversion from X0 Map to Hitmap
  /*!This converts the integer values used for the binning of the radiation length map into global coordinates*/
  std::pair< double, double > ConversionX0mapToHitmap(int x, int y);

  //!Conversion from Hitmap to X0 Map
  /*!This converts the global coordinates into integer bins for use when filling the radiation length maps*/
  std::pair< int, int > ConversionHitmapToX0map(double x, double y);

  /***********************
  //Private member values*
  ***********************/

  //Beam energy
  double _beamEnergy;

  //Cut
  double _cut;

  //Double Daf Fitted
  bool _doubleDafFitted;

  //DUT position
  double _dutPosition;

  //Current event number
  int _eventNumber;

  //These two maps associate a name with a histogram for easy filling later
  std::map< std::string , TH1* > _histoThing;
  std::map< std::string , TH2* > _histoThing2D;
  
  //Stores the information being brought in from the Marlin process which contains information about post-aligned tracks
  IMPL::LCCollectionVec* _inputTrackCollectionVec;

  //The track collection name, as input from the steering file
  std::string _inputTrackColName;

  //The current run number
  int _runNumber;

  //The name of the track collection to use, the is input in the steering file
  std::string _trackCollectionName;
  std::string _trackCollectionName1;
  std::string _trackCollectionName2;

  //This is the number of bins to use in the residual plots
  int nobins;
  
  //This is the number of bins to use in the angle plots
  int nobinsangle;

  //The minimum bin in the residual plots
  double minbin;
  
  //The maximum bin in the residual plots
  double  maxbin;

  //The minimum bin in the angular plots
  double minbinangle;

  //The maximum bin in the angular plots
  double maxbinangle;

  //The number of bins in the radiation length map plots in the x direction
  int binsx;

  //The minimum on the spatial range of the radiation length map plots in x, this is measured in mm and anything below -11 is beyond the telescope sensor
  double minx;
  
  //The maximum on the spatial range of the radiation length map plots in x, this is measured in mm and anything above 11 is beyond the telescope sensor
  double maxx;

  //The number of bins in the radiation length map plots in the y direction
  int binsy;

  //The minimum on the spatial range of the radiation length map plots in y, this is measured in mm and anything below -6 is beyond the telescope sensor
  double miny;
  
  //The maximum on the spatial range of the radiation length map plots in y, this is measured in mm and anything above 6 is beyond the telescope sensor
  double maxy;

  //The granularity of the radiation length maps in x, measured in mm
  double binsizex;
  
  //The granularity of the radiation length maps in y, measured in mm
  double binsizey;

  //An attempt to make a folder for the histograms to live in
  TDirectory *X0ProcessorDirectory;

  TH1D *AngleXFrontThreePlanesDoubleDaf;
  TH1D *AngleYFrontThreePlanesDoubleDaf;
  TH1D *AngleXBackThreePlanesDoubleDaf;
  TH1D *AngleYBackThreePlanesDoubleDaf;
  TH2D *AngleXYFrontThreePlanesDoubleDaf;
  TH2D *AngleXYBackThreePlanesDoubleDaf;
  TH1D *ScatteringAngleXDoubleDaf;
  TH1D *ScatteringAngleYDoubleDaf;
  TH2D *ScatteringAngleXYDoubleDaf;
  TH2D *ScatteringAngleXTripleMapDoubleDaf;
  TH2D *ScatteringAngleYTripleMapDoubleDaf;
  TH2D *RadiationLengthMapDoubleDaf;
  TH1D *AngleXForwardTripleFirstThreePlanes;
  TH1D *AngleXForwardTripleLastThreePlanes;
  TH1D *AngleYForwardTripleFirstThreePlanes;
  TH1D *AngleYForwardTripleLastThreePlanes;
  TH2D *AngleXYForwardTripleFirstThreePlanes;
  TH2D *AngleXYForwardTripleLastThreePlanes;
  TH1D *ScatteringAngleXTriple;
  TH1D *ScatteringAngleYTriple;
  TH2D *ScatteringAngleXYTriple;
  TH2D *ScatteringAngleXTripleMap;
  TH2D *ScatteringAngleYTripleMap;
  TH2D *RadiationLengthTripleMap;
  TH1D *SinglePointResidualXPlane0;
  TH1D *SinglePointResidualXPlane1;
  TH1D *SinglePointResidualXPlane2;
  TH1D *SinglePointResidualXPlane3;
  TH1D *SinglePointResidualXPlane4;
  TH1D *SinglePointResidualXPlane5;
  TH1D *SinglePointResidualYPlane0;
  TH1D *SinglePointResidualYPlane1;
  TH1D *SinglePointResidualYPlane2;
  TH1D *SinglePointResidualYPlane3;
  TH1D *SinglePointResidualYPlane4;
  TH1D *SinglePointResidualYPlane5;
  TH1D *ThreePointResidualXPlane1;
  TH1D *ThreePointResidualXPlane2;
  TH1D *ThreePointResidualXPlane3;
  TH1D *ThreePointResidualXPlane4;
  TH1D *ThreePointResidualYPlane1;
  TH1D *ThreePointResidualYPlane2;
  TH1D *ThreePointResidualYPlane3;
  TH1D *ThreePointResidualYPlane4;
  TH1D *AngleXForwardPlane0;
  TH1D *AngleXForwardPlane1;
  TH1D *AngleXForwardPlane2;
  TH1D *AngleXForwardPlane3;
  TH1D *AngleXForwardPlane4;
  TH1D *AngleYForwardPlane0;
  TH1D *AngleYForwardPlane1;
  TH1D *AngleYForwardPlane2;
  TH1D *AngleYForwardPlane3;
  TH1D *AngleYForwardPlane4;
  TH2D *AngleXYForwardPlane0;
  TH2D *AngleXYForwardPlane1;
  TH2D *AngleXYForwardPlane2;
  TH2D *AngleXYForwardPlane3;
  TH2D *AngleXYForwardPlane4;
  TH1D *ScatteringAngleXPlane1;
  TH1D *ScatteringAngleXPlane2;
  TH1D *ScatteringAngleXPlane3;
  TH1D *ScatteringAngleXPlane4;
  TH1D *ScatteringAngleYPlane1;
  TH1D *ScatteringAngleYPlane2;
  TH1D *ScatteringAngleYPlane3;
  TH1D *ScatteringAngleYPlane4;
  TH2D *KinkAnglePlane1;
  TH2D *KinkAnglePlane2;
  TH2D *KinkAnglePlane3;
  TH2D *KinkAnglePlane4;
  std::map< int, TH2D* > ScatteringAngleXSingleMap;
  std::map< int, TH2D* > ScatteringAngleYSingleMap;
  std::map< int, TH2D* > RadiationLengthSingleMap;
//  TH2D *ScatteringAngleXPlane1Map;
//  TH2D *ScatteringAngleXPlane2Map;
//  TH2D *ScatteringAngleXPlane3Map;
//  TH2D *ScatteringAngleXPlane4Map;
//  TH2D *ScatteringAngleYPlane1Map;
//  TH2D *ScatteringAngleYPlane2Map;
//  TH2D *ScatteringAngleYPlane3Map;
//  TH2D *ScatteringAngleYPlane4Map;
//  TH2D *RadiationLengthPlane1Map;
//  TH2D *RadiationLengthPlane2Map;
//  TH2D *RadiationLengthPlane3Map;
//  TH2D *RadiationLengthPlane4Map;
  std::map< int, std::map< std::pair< int, int >, std::vector< double > > > ScatteringAngleXSingleMapData;
  std::map< int, std::map< std::pair< int, int >, std::vector< double > > > ScatteringAngleYSingleMapData;
  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleXTripleMapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleYTripleMapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleXTripleMapDataDaf; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleYTripleMapDataDaf; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleXPlane1MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleXPlane2MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleXPlane3MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleXPlane4MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleYPlane1MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleYPlane2MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleYPlane3MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
//  std::map< std::pair< int, int >, std::vector< double > > ScatteringAngleYPlane4MapData; //Pair gives the x and y bins of the track at the point of the DUT and the value of the double is the scattering angle
};
//! A global instance of the processor
EUTelX0Processor gEUTelX0Processor;
}
#endif
#endif

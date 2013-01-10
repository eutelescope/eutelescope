// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 * This source code is part of the Eutelescope package of Marlin.
 * You are free to use this source files for your own development as
 * long as it stays in a public research context. You are not
 * allowed to use it for commercial purpose. You must put this
 * header with author names in all development based on this file.
 * EUTelX0Processor written by Phillip Hamnett (phillip.hamnett@gmail.com)
 */

#ifndef EUTELX0PROCESSOR_H
#define EUTELX0PROCESSOR_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelReferenceHit.h"
#include "IMPL/TrackerHitImpl.h"
// marlin includes ".h"
#include "marlin/Processor.h"
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#include <AIDA/IPlotter.h>
#include <RAIDA/IAnalysisFactoryROOT.h>
#include <RAIDA/IPlotterFactoryROOT.h>
#include <RAIDA/IPlotterROOT.h>
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IPlotterFactory.h>
#include <AIDA/IPlotterRegion.h>
#include <AIDA/IFitter.h>
#include <AIDA/IFitResult.h>
// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCIO.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Track.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/IAxis.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <limits>
#include <sstream>
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMinuit.h"
#include "TVector3.h"
#include "TLine.h"
class TMinuit;
class EVENT::TrackerHit;
class IMPL::LCCollectionVec;
class AIDA::IAnalysisFactory;
class AIDA::IPlotter;
class AIDA::IPlotterFactory;
class AIDA::IPlotterRegion;
class AIDA::IFitter;
class AIDA::IFitResult;
#endif

namespace eutelescope {

class EUTelX0Processor : public marlin::Processor {

public:

 
	//! Returns a new instance of EUTelX0Processor
	/*! This method returns a new instance of this processor. It is
	* called by Marlin execution framework and it shouldn't be
	* called/used by the final user.
	*
	* @return a new EUTelX0Processor.
	*/
	virtual Processor * newProcessor() {return new EUTelX0Processor;}

	//! Default constructor
	EUTelX0Processor ();

	//! Called at the job beginning.
	/*! This is executed only once in the whole execution. It prints
	* out the processor parameters and check that the GEAR
	* environment is properly set up and accessible from Marlin.
	*/
	virtual void init ();

	//! Called for every run.
	/*! It is called for every run, and consequently the run counter
	* is incremented. The geometry ID of the file is compared with
	* the one provided by the GEAR geometry description. In case the
	* two are different, the user is asked to decide to quit or to
	* continue with a description that might be wrong.
	*
	* @param run the LCRunHeader of the this current run
	*/
	virtual void processRunHeader (LCRunHeader * run);

	//! Called every event
	/*! This is called for each event in the file. Each element of the
	* pulse collection is scanned and the center of the cluster is
	* translated into the external frame of reference thanks to the
	* GEAR geometry description.
	*
	* @throw UnknownDataTypeException if the cluster type is unknown
	*
	* @param evt the current LCEvent event as passed by the
	* ProcessMgr
	*/
	virtual void processEvent (LCEvent * evt);

	//! Called after data processing.
	/*! This method is called when the loop on events is
	* finished.
	*/
	virtual void end();

	//! Histogram booking
	/*! Some control histograms are filled during this procedure in
	* order to be able to perform easy check on the quality of the
	* output hits and also to understand if the frame of reference
	* conversion has been properly done. Of course this method is
	* effectively doing something only in the case MARLIN_USE_AIDA.
	*/
	void bookHistos();
 

	//!Guess Sensor ID
	/*!Works out which layer of the telescope the current element is in based on its position
	*/ 
	int guessSensorID(const double * hit );

private:
//Following #define stops the accidental creation of a copy or assignment operator by causing a link error. Copy and Assignment operators not allowed because they are unnecessary and the cause of many bugs
#define DISALLOW_COPY_AND_ASSIGN(EUTelX0Processor) \
	EUTelX0Processor(const EUTelX0Processor&); \
	void operator=(const EUTelX0Processor&);
	
	//Private Functions
	DISALLOW_COPY_AND_ASSIGN(EUTelX0Processor)//See #define just above

	//!Cut and Store Hits
	/*!This looks at the two layers and extrapolates a track line between them. It also allows for a cut to be made that stops tracks with angles greater than +/- anglecut.
	   size1 and size2 refer to the specific sizes of the number of hits in the first and second layer.
	   firstLayer and secondLayer inform the user in which two layers the extrapolation is being made. It is important that you get the order the right way around for the extrapolation. For instance, if going from layer 0 and 1 to layer 2 then the firstLayer = 0 and secondLayer = 1. But if going from layers 3 and 4 to 2 then the extrapolation is in the opposite direction and so firstLayer = 4 and secondLayer = 3.
	*/
	void basicFitter(LCCollection *alignedHitCollection);
	void cutAndStoreHits(size_t size1, size_t size2, int firstLayer, int secondLayer, double anglecut);
	//!Caculate X0
	/*!This function calculates the value of the material budget.
	*/
	double calculateX0();
        void testtrack(LCCollection *trackCollection);
	void threePointResolution(LCCollection *alignedHitCollection);
	void createResiduals(LCCollection *trackCollection);

        std::string _correctedHitColName;
        double _cutValue1,_cutValue2;
	bool _debug;
	int _debugCount;  //TODO(Phillip Hamnett): Can this be deleted?
	int _noTracks;
	int _noEvents;
	std::string _fileName;
	std::map<std::string , AIDA::IBaseHistogram * > _histoThing;  //This map contains all the histograms that are going to be plotted in this processor
	std::map<std::string , std::vector< double > > _histoData;
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
        static std::string _histoResidualXZ;
        static std::string _histoResidualYZ;
        static std::string _histoResidualXY;
	std::map< int, std::vector< TVector3 > > _hitInfo;  //This stores the hit position in a TVector3. If there are multiple hits then they are all stored in the vector of TVector3's. The int key refers to the layer number
        ofstream _hitsfiles;
	IMPL::LCCollectionVec* _inputCollectionVec;  //Stores the information being brought in from the Marlin process which contains information about post-aligned hits
        std::string _inputHitColName;
	std::string _inputHitCollectionName;  //Stores the name of the parameter to bring in from the Marlin steering file
        std::string _inputTrackColName;
	static const int _noLayers = 5;  //TODO(Phillip Hamnett): This is a hack just so that the computer knows how many layers there are, is this stored somewhere else or should it be made as a non-const variable in case at some future date more layers are added to the telescope?
	//! Event number
	int _iEvt;
	//! Run number
	int _iRun;
	std::map< int, std::vector< TVector3 > > _projectedHits;  //int refers to 1 or 2, with 1 being the projection from the 01 plane and 2 being the projection from the 43 plane
	std::string _referenceHitCollectionName;  //Stores the name of the file to be brought in by the Marlin steering file that is used to determine layer number
        std::map<std::pair< Double_t, Double_t > , std::vector< TVector3 > > _residual;  //The pair of doubles refers to the x and y coordinates in layer 2. The vector of TVector3's refers to the positions of the projected hits
	std::map<std::pair< Double_t, Double_t > , std::vector< TVector3 > > _residualAngle;  //As above but the TVector3 doesn't contain xyz but instead theta, phi and alpha
	std::map<std::pair <Double_t, Double_t > , std::vector< TVector3 > > _residualProfile; //TODO(Phillip Hamnett): Can this be joined with _residual? //Used as above but for created a profile histogram
	IMPL::LCCollectionVec* _referenceHitVec;   //Stores the information being brought in from the Marlin process containing information about the geometry of the setup for working out the layer number
	
};
	//! A global instance of the processor
	EUTelX0Processor gEUTelX0Processor;
}
#endif
#endif

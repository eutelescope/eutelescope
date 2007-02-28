// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelCalculateEtaProcessor.cc,v 1.2 2007-02-28 08:15:59 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelCalculateEtaProcessor.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelFFClusterImpl.h"

// marlin includes ".h"
#include "marlin/Exceptions.h"
#include "marlin/Processor.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"
#include "PseudoHistogram.h"

#ifdef MARLIN_USE_AIDA
// aida includes <.h>
#include "marlin/AIDAProcessor.h"
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>
#include <AIDA/IAxis.h>
#endif

#ifdef MARLIN_USE_ROOT
#include "ROOTProcessor.h"
#include <TH1.h>
#endif 

// lcio includes <.h> 
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>


// system includes <>
#include <iostream>
#include <string>
#include <sstream>


using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;

#ifdef MARLIN_USE_HISTOGRAM
string EUTelCalculateEtaProcessor::_cogHistogramXName = "CoG-X";
string EUTelCalculateEtaProcessor::_cogHistogramYName = "CoG-Y";
string EUTelCalculateEtaProcessor::_cogIntegralXName  = "Integral-CoG-X";
string EUTelCalculateEtaProcessor::_cogIntegralYName  = "Integral-CoG-Y";
string EUTelCalculateEtaProcessor::_etaHistoXName     = "EtaProfile-X";
string EUTelCalculateEtaProcessor::_etaHistoYName     = "EtaProfile-Y";
string EUTelCalculateEtaProcessor::_cogHisto2DName    = "CoG-Histo2D";
#endif 

EUTelCalculateEtaProcessor::EUTelCalculateEtaProcessor () : Processor("EUTelCalculateEtaProcessor") {

  // modify processor description
  _description =
    "EUTelCalculateEtaProcessor calculates the eta function for a given set of clusters";

  _isEtaCalculationFinished = false;

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERDATA, "ClusterCollectionName",
			   "Input cluster collection",
			   _clusterCollectionName, string ("cluster"));

  registerProcessorParameter("EventNumber",
			     "Write here how many events you want to use for eta calculation (-1 for all)",
			     _nEvent, static_cast<int> ( -1 ));

  IntVec noOfBinExample;
  noOfBinExample.push_back(1000);
  noOfBinExample.push_back(1000);
  registerProcessorParameter("NumberOfBins",
			     "Write here in how many bins the seed pixel should be divided (x and y)",
			     _noOfBin, noOfBinExample, noOfBinExample.size());

  registerProcessorParameter("ClusterQualitySelection",
			     "To use only kGoodQuality write 0 here",
			     _clusterQuality, static_cast<int> ( 0 ));
  registerProcessorParameter("ClusterTypeSelection",
			     "Write FULL: full cluster, NxMPixel: for a NxM sub-cluster, NPixel: to use only N pixel",
			     _clusterTypeSelection, string("FULL"));

  IntVec xyCluSizeExample;
  xyCluSizeExample.push_back(3);
  xyCluSizeExample.push_back(3);
  registerProcessorParameter("NxMPixelClusterSize",
			     "The size along x and y of the subcluster (only for NxMPixel)",
			     _xyCluSize, xyCluSizeExample, xyCluSizeExample.size());
 
  registerProcessorParameter("NPixelSize",
			     "The number of pixel with the highest signal (only for NPixel)",
			     _nPixel, static_cast<int> ( 5 ));
 
}


void EUTelCalculateEtaProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // reset stl vectors
  _cogHistogramX.clear();
  _cogHistogramY.clear();
  _integralHistoX.clear();
  _integralHistoY.clear();
  
}

void EUTelCalculateEtaProcessor::processRunHeader (LCRunHeader * rdr) {

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  EUTelRunHeaderImpl *  runHeader = static_cast<EUTelRunHeaderImpl*>(rdr);

  _noOfDetector = runHeader->getNoOfDetector();

  int tempEvent = min( runHeader->getNoOfEvent(),
		       Global::parameters->getIntVal("MaxRecordNumber") ) - 1;
  
  if ( ( _nEvent == -1 ) || ( _nEvent >= tempEvent ) ) {
    _nEvent = tempEvent;
  }
  
  
  if ( !_isEtaCalculationFinished ) {

    const double min = -0.5;
    const double max = +0.5;
    int xNoOfBin = _noOfBin[0];
    int yNoOfBin = _noOfBin[1];
    
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      PseudoHistogram * histoX = new PseudoHistogram(xNoOfBin, min, max);
      _cogHistogramX.push_back(histoX);
      
      PseudoHistogram * histoY = new PseudoHistogram(yNoOfBin, min, max);
      _cogHistogramY.push_back(histoY);
      
      PseudoHistogram * integralX = new PseudoHistogram(xNoOfBin, min, max);
      _integralHistoX.push_back(integralX);
      
      PseudoHistogram * integralY = new PseudoHistogram(yNoOfBin, min, max);
      _integralHistoY.push_back(integralY);
      
      
#ifdef MARLIN_USE_ROOT
      {
	string name, title, path;
	{
	  stringstream sname;
	  sname  << _cogHistogramXName << "-" << iDetector ;
	  name = sname.str();
	  
	  stringstream stitle;
	  stitle << "CoG shift along X on detector " << iDetector;
	  title = stitle.str();
	  
	  stringstream spath;
	  spath << "detector-" << iDetector << "/" ;
	  path = spath.str();
	}
	TH1D * cogHistoX = new TH1D(name.c_str(), title.c_str(), xNoOfBin, min, max);
	ROOTProcessor::addTObject(this, cogHistoX);
	
	// ** TO BE COMPLETED ** //
	//
	// As soon as the ROOTProcess will be in an advanced
	// developement stage, all histogramming functions will be
	// made both for AIDA and also for ROOT

      }
#endif 

#ifdef MARLIN_USE_AIDA
      {
	string name, title, path;
	{
	  stringstream sname;
	  sname  << _cogHistogramXName << "-" << iDetector ;
	  name = sname.str();
	  
	  stringstream stitle;
	  stitle << "CoG shift along X on detector " << iDetector;
	  title = stitle.str();
	  
	  stringstream spath;
	  spath << "detector-" << iDetector << "/" ;
	  path = spath.str();
	}
	
	AIDAProcessor::tree(this)->mkdir(path.c_str());
	
	AIDA::IHistogram1D * cogHistoX = AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), xNoOfBin, 
												   min, max);
	cogHistoX->setTitle(title.c_str());
	_aidaHistoMap.insert( make_pair(name, cogHistoX) );
	
	{
	  stringstream sname;
	  sname  << _cogHistogramYName << "-" << iDetector ;
	  name = sname.str();
	  
	  stringstream stitle;
	  stitle << "CoG shift along Y on detector " << iDetector;
	  title = stitle.str();
	}
	
	AIDA::IHistogram1D * cogHistoY = AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), yNoOfBin, 
												   min, max);
	cogHistoY->setTitle(title.c_str());
	_aidaHistoMap.insert( make_pair(name, cogHistoY) );
	
	{
	  stringstream sname;
	  sname  << _cogIntegralXName << "-" << iDetector ;
	  name = sname.str();
	  
	  stringstream stitle;
	  stitle << "Integral CoG (x) shift histogram on " << iDetector;
	  title = stitle.str();
	}
	AIDA::IHistogram1D * cogIntegralHistoX = 
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), xNoOfBin, min, max);
	cogIntegralHistoX->setTitle(title.c_str());
	_aidaHistoMap.insert( make_pair(name, cogIntegralHistoX) );
	
	{
	  stringstream sname;
	  sname  << _cogIntegralYName << "-" << iDetector ;
	  name = sname.str();
	  
	  stringstream stitle;
	  stitle << "Integral CoG (y) shift histogram on " << iDetector;
	  title = stitle.str();
	}
	AIDA::IHistogram1D * cogIntegralHistoY = 
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (path + name).c_str(), yNoOfBin, min, max);
	cogIntegralHistoY->setTitle(title.c_str());
	_aidaHistoMap.insert( make_pair(name, cogIntegralHistoY) );
	
	{
	  stringstream sname;
	  sname << _etaHistoXName << "-" << iDetector ;
	  name = sname.str();
	  
	  stringstream stitle;
	  stitle << "Eta profile x for detector " << iDetector;
	  title = stitle.str();
	}
	AIDA::IProfile1D * etaHistoX = AIDAProcessor::histogramFactory(this)
	  ->createProfile1D( (path + name).c_str(), xNoOfBin, min, max, min, max);
	etaHistoX->setTitle(title.c_str());
	_aidaHistoMap.insert( make_pair(name, etaHistoX) );
	
	{
	  stringstream sname;
	  sname << _etaHistoYName << "-" << iDetector ;
	  name = sname.str();
	  
	  stringstream stitle;
	  stitle << "Eta profile y for detector " << iDetector;
	  title = stitle.str();
	}
	AIDA::IProfile1D * etaHistoY = AIDAProcessor::histogramFactory(this)
	  ->createProfile1D( (path + name).c_str(), yNoOfBin, min, max, min, max);
	etaHistoY->setTitle(title.c_str());
	_aidaHistoMap.insert( make_pair(name, etaHistoY) );

	{
	  stringstream sname;
	  sname << _cogHisto2DName << "-" << iDetector ;
	  name = sname.str();

	  stringstream stitle;
	  stitle << "2D Histo with the CoG within the seed pixel for detector " << iDetector;
	  title = stitle.str();
	}
	AIDA::IHistogram2D * cogHisto2D = AIDAProcessor::histogramFactory(this)
	  ->createHistogram2D( (path + name).c_str(), xNoOfBin, min, max, yNoOfBin, min, max);
	cogHisto2D->setTitle(title.c_str());
	_aidaHistoMap.insert( make_pair(name, cogHisto2D) );
      }
#endif
    }
  }
  
  // increment the run counter
  ++_iRun;
  
}


void EUTelCalculateEtaProcessor::processEvent (LCEvent * evt) {


  if ( !_isEtaCalculationFinished ) {

    if (_iEvt % 10 == 0) 
      cout << "[" << name() << "] Filling CoG histograms event " << _iEvt << endl;
    
    LCCollectionVec * clusterCollectionVec    = dynamic_cast < LCCollectionVec * > (evt->getCollection(_clusterCollectionName));
    CellIDDecoder<TrackerDataImpl> cellDecoder(clusterCollectionVec);
    
    for (int iCluster = 0; iCluster < clusterCollectionVec->getNumberOfElements() ; iCluster++) {
      
      EUTelFFClusterImpl * cluster = static_cast< EUTelFFClusterImpl * > ( clusterCollectionVec->getElementAt(iCluster) );
      int detectorID = cluster->getDetectorID();
      float xShift, yShift;
      
      if ( cluster->getClusterQuality() == static_cast<ClusterQuality> (_clusterQuality) ) {
	
	if ( _clusterTypeSelection == "FULL" ) {
	  cluster->getCenterOfGravityShift(xShift, yShift);
	} else if ( _clusterTypeSelection == "NxMPixel" ) {
	  cluster->getCenterOfGravityShift(xShift, yShift, _xyCluSize[0], _xyCluSize[1]);
	} else if ( _clusterTypeSelection == "NPixel" ) {
	  cluster->getCenterOfGravityShift(xShift, yShift, _nPixel);
	}
	
	_cogHistogramX[detectorID]->fill(static_cast<double>(xShift), 1.0);
	_cogHistogramY[detectorID]->fill(static_cast<double>(yShift), 1.0);
	
#ifdef MARLIN_USE_AIDA
	{
	  string name;
	  {
	    stringstream ss;
	    ss << _cogHistogramXName << "-" << detectorID;
	    name = ss.str();
	  }
	  (dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]))->fill(xShift);
	  
	  {
	    stringstream ss;
	    ss << _cogHistogramYName << "-" << detectorID;
	    name = ss.str();
	  }
	  (dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]))->fill(yShift);

	  {
	    stringstream ss;
	    ss << _cogHisto2DName << "-" << detectorID ;
	    name = ss.str();
	  }
	  (dynamic_cast< AIDA::IHistogram2D* > (_aidaHistoMap[name]))->fill(xShift, yShift);
	}
#endif
      
#ifdef MARLIN_USE_ROOT
	{
	  string name;
	  {
	    stringstream ss;
	    ss << _cogHistogramXName << "-" << detectorID;
	    name = ss.str();
	  }
	  (dynamic_cast<TH1D*> (ROOTProcessor::getTObject(this, name.c_str())))->Fill(xShift);
	}
#endif
      }
    }
    
    ++_iEvt;
    
    if ( isLastEvent() )  finishCalculation();
  }

  setReturnValue( "isEtaCalculationFinished" , _isEtaCalculationFinished);
  
}



void EUTelCalculateEtaProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelCalculateEtaProcessor::finishCalculation() {
  
  double integral;
  vector<vector<double > > xEtaBinCenterVec;
  vector<vector<double > > xEtaValueVec;

  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    
    for (int iBin = 1; iBin < _cogHistogramX[iDetector]->getNumberOfBins(); iBin++ ) {
      double x = _cogHistogramX[iDetector]->getBinCenter(iBin);
      integral = _cogHistogramX[iDetector]->integral(1, iBin);
      _integralHistoX[iDetector]->fill(x, integral);
      
    }
    
    vector<double > xEtaBinCenter;
    vector<double > xEtaBinValue;
    
    for (int iBin = 1; iBin < _integralHistoX[iDetector]->getNumberOfBins(); iBin++) {
      xEtaBinCenter.push_back( _integralHistoX[iDetector]->getBinCenter(iBin) );
      xEtaBinValue.push_back(  _integralHistoX[iDetector]->getBinContent(iBin) / integral );
    }
    xEtaBinCenterVec.push_back(xEtaBinCenter);
    xEtaValueVec.push_back(xEtaBinValue);
 
    for (int iBin = 1; iBin < _cogHistogramY[iDetector]->getNumberOfBins(); iBin++ ) {
      double y = _cogHistogramY[iDetector]->getBinCenter(iBin);
      integral = _cogHistogramY[iDetector]->integral(1, iBin);
      _integralHistoY[iDetector]->fill(y, integral);
    }




#ifdef MARLIN_USE_AIDA
    integral = 0;
    string name;
    {
      stringstream ss;
      ss << _cogIntegralXName << "-" << iDetector;
      name = ss.str();
    }
    AIDA::IHistogram1D * integralHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    {
      stringstream ss;
      ss << _cogHistogramXName << "-" << iDetector;
      name = ss.str();
    }
    AIDA::IHistogram1D * cogHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    integral = 0;
    for (int iBin = 0; iBin < cogHisto->axis().bins(); iBin++) {
      double x = cogHisto->axis().binLowerEdge(iBin) + 0.5 * cogHisto->axis().binWidth(iBin);
      integral += cogHisto->binHeight(iBin);
      integralHisto->fill(x, integral);
    }

    {
      stringstream ss;
      ss << _etaHistoXName << "-" << iDetector ;
      name = ss.str();
    }
    AIDA::IProfile1D * etaXHisto = dynamic_cast< AIDA::IProfile1D* > (_aidaHistoMap[name]);

    for (int iBin = 0; iBin < integralHisto->axis().bins(); iBin++) {
      double x = integralHisto->axis().binLowerEdge(iBin) + 0.5 * integralHisto->axis().binWidth(iBin);
      double v = integralHisto->binHeight(iBin);
      etaXHisto->fill(x, v/integral - 0.5);
    }

    {
      stringstream ss;
      ss << _cogIntegralYName << "-" << iDetector;
      name = ss.str();
    }
    integralHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    {
      stringstream ss;
      ss << _cogHistogramYName << "-" << iDetector;
      name = ss.str();
    }
    cogHisto = dynamic_cast< AIDA::IHistogram1D* > (_aidaHistoMap[name]);

    integral = 0;
    for (int iBin = 0; iBin < cogHisto->axis().bins(); iBin++) {
      double y = cogHisto->axis().binLowerEdge(iBin) + 0.5 * cogHisto->axis().binWidth(iBin);
      integral += cogHisto->binHeight(iBin);
      integralHisto->fill(y, integral);
    }

    {
      stringstream ss;
      ss << _etaHistoYName << "-" << iDetector ;
      name = ss.str();
    }
    AIDA::IProfile1D * etaYHisto = dynamic_cast< AIDA::IProfile1D* > (_aidaHistoMap[name]);

    for (int iBin = 0; iBin < integralHisto->axis().bins(); iBin++) {
      double y = integralHisto->axis().binLowerEdge(iBin) + 0.5 * integralHisto->axis().binWidth(iBin);
      double v = integralHisto->binHeight(iBin);
      etaYHisto->fill(y, v/integral - 0.5);
    }

#endif

  } 

  _isEtaCalculationFinished = true;
  throw RewindDataFilesException(this);

}


void EUTelCalculateEtaProcessor::end() {
  cout << "[" << name() <<"] Successfully finished" << endl;



}


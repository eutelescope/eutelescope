/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELPREALIGNER_H
#define EUTELPREALIGNER_H

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/TrackerHitImpl.h>

// gear includes <.h>
#include <gear/SiPlanesLayerLayout.h>
#include <gear/SiPlanesParameters.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#endif

// system includes <>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace eutelescope {

  //class implementation of PreAligner
  class PreAligner {
  
  private:
    float pitchX, pitchY;
    std::vector<int> histoX, histoY;
    float minX, maxX;
    float range;
    float zPos;
    int iden;
    
    //! function to get maximum bin of given histogram
    float getMaxBin(std::vector<int> &histo) {
      int maxBin(0), maxVal(0);
      //loop over bins
      for(size_t ibin = 0; ibin < histo.size(); ibin++) {
	if(histo.at(ibin) > maxVal) {
	  maxBin = ibin;
	  maxVal = histo.at(ibin);
	}
      }
      
      //check if maxBin is not at the edges
      if(maxBin == 0 || maxBin == static_cast<int>(histo.size())) {
        streamlog_out(WARNING3)
	  << "At least one sensor frame might be empty or heavily "
	  "misaligned. Please check the GEAR file!"
	  << " MaxBin: " << maxBin << " histo.size(): " << histo.size()
	  << std::endl;
        return static_cast<float>(maxBin);
      }
      	
      //get weighted position from three neighboring bins:
      float weight(0.0);
      double pos1(0.0);
      double pos2(0.0);
      double pos3(0.0);
      try {
	//use logarithms to be safe even with large number of bin entries
	pos1 = log(maxBin - 1) + log(histo.at(maxBin - 1));
	pos2 = log(maxBin) + log(histo.at(maxBin));
	pos3 = log(maxBin + 1) + log(histo.at(maxBin + 1));
	weight = log((histo.at(maxBin - 1)) + (histo.at(maxBin)) +
                     (histo.at(maxBin + 1)));
      } catch(...) {
	streamlog_out(ERROR) 
	  << "Could not execute prealignment bin content retrieval. The "
	  "sensor frame might be empty or heavily misaligned. Please "
	  "check the GEAR file!"
	  << std::endl;
      }
      return (exp(pos1-weight) + exp(pos2-weight) + exp(pos3-weight));
    }
    
  public:
  PreAligner(float pitchX, float pitchY, float zPos, int iden)
    : pitchX(pitchX), pitchY(pitchY), minX(-40.0), maxX(40),
      range(maxX - minX), zPos(zPos), iden(iden) {
      
      histoX.assign(int(range / pitchX), 0);
      histoY.assign(int(range / pitchY), 0);
    }
    
    void *current() { return this; }
    
    float getZPos() const { 
      return (zPos); 
    }
    
    int getIden() const { 
      return (iden); 
    }
    
    //add point if within bounds, throw away data that is out of bounds
    void addPoint(float x, float y) {
      try {
	histoX.at(static_cast<int>((x - minX) / pitchX)) += 1;
      } catch(std::out_of_range &e) {;}
      try {
	histoY.at(static_cast<int>((y - minX) / pitchY)) += 1;
      } catch(std::out_of_range &e) {;}
    }	
    
    float getPeakX() { 
      return ((getMaxBin(histoX) * pitchX) + minX); 
    }
    
    float getPeakY() { 
      return ((getMaxBin(histoY) * pitchY) + minX); 
    }
  };
  
  
  class EUTelPreAligner : public marlin::Processor {
    
  public:
    //! Returns a new instance of EUTelPreAligner
    virtual Processor *newProcessor() { 
      return new EUTelPreAligner; 
    }
    
    //! default constructor
    EUTelPreAligner();
    
    virtual void init();
    
    virtual void processRunHeader(LCRunHeader *run);
    
    virtual void processEvent(LCEvent *evt);
    
    virtual void end();
    
    //! Histogram booking
    void bookHistos();
    
  private:
    //! How many events are needed to get reasonable plots (correlation & offset)
    int _requiredEvents;
    
    //! Sensor ID vector
    std::vector<int> _sensorIDVec;
    
    //! map for sensor ID to position along Z id
    std::map<int, int> _sensorIDtoZOrderMap;
    
    //! fixed plane for prealignment
    int _fixedID;

    //! Residual cuts [relative to the first upstream plane!]
    //! vector of correlation band cuts in X (upper limit)
    std::vector<float> _residualsXMax;
    //! vector of correlation band cuts in X (lower limit)
    std::vector<float> _residualsXMin;
    //! vector of correlation band cuts in Y (upper limit)
    std::vector<float> _residualsYMax;
    //! vector of correlation band cuts in Y (lower limit)
    std::vector<float> _residualsYMin;

    //! Minimal number of correlated hits
    int _minNumberOfCorrelatedHits;

    //! Boolean for turning histogram creation on and off
    bool _histogramSwitch;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    std::map<unsigned int, AIDA::IBaseHistogram *> _hitXCorr;
    std::map<unsigned int, AIDA::IBaseHistogram *> _hitYCorr;
#endif

  protected:
    int _iRun;
    int _iEvt;
    std::string _inputHitCollectionName;
    std::string _alignmentConstantLCIOFile;
    std::string _GEARFileSuffix;
    bool _dumpGEAR;
    std::vector<PreAligner> _preAligners;
    std::vector<int> _ExcludedPlanesXCoord;
    std::vector<int> _ExcludedPlanesYCoord;
    std::vector<int> _ExcludedPlanes;
  };
  
  //! A global instance of the processor
  EUTelPreAligner gEUTelPreAligner;
}
#endif

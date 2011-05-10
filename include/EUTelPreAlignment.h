// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// built only if USE_GEAR
#if defined(USE_GEAR)
#ifndef EUTELPREALIGNMENT_H
#define EUTELPREALIGNMENT_H

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/TrackerHitImpl.h>

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif


// system includes <>
#include <iostream>
#include <string>
#include <map>
#include <cstdio>

namespace eutelescope {
class PreAligner{
private:
  float pitchX, pitchY;
  std::vector<int> histoX, histoY;
  float minX, maxX;
  float range;
  float zPos;
  int iden;
  float getMaxBin(std::vector<int>& histo){
    int maxBin(0), maxVal(0);
    int sum(0);
    for(size_t ii = 0; ii < histo.size(); ii++){
      sum += histo.at(ii);
      if(histo.at(ii) > maxVal){ 
	maxBin = ii; 
	maxVal = histo.at(ii);
      }
    }
    //Get weighted position from 3 neighboring bins
    float weight(0.0);
    float pos(0.0);
    for(int ii = -1; ii <2; ii++){
      int tmpBin = maxBin + ii;
      if(maxBin < 0 or maxBin >= (int)histo.size()){ continue;}
try
{
      pos += tmpBin * histo.at(tmpBin);
      weight += histo.at(tmpBin);
}
catch(...)
{
      std::cout << "Could not execute bin content retrieval. The sensor frame empty or heavily misalgined" << std::endl; 
}
    }
    return(pos/weight);
  }
public:
  PreAligner(float pitchX, float pitchY, float zPos, int iden): 
    pitchX(pitchX), pitchY(pitchY), 
    minX(-20.0), maxX(20), range(maxX - minX),
    zPos(zPos), iden(iden){
    histoX.assign( int( range / pitchX ), 0);
    histoY.assign( int( range / pitchY ), 0);
  }
  void* current(){return this; } 
  float getZPos() const { return(zPos); }
  int getIden() const { return(iden); }
  void addPoint(float x, float y){
    //Add to histo if within bounds, throw away data that is out of bounds
    try{
      histoX.at( static_cast<int> ( (x - minX)/pitchX) ) += 1; 
    } catch (std::out_of_range& e) {;}
    try{
      histoY.at( static_cast<int> ( (y - minX)/pitchY) ) += 1; 
    } catch (std::out_of_range& e) {;}
  }
  float getPeakX(){
    return( (getMaxBin(histoX) * pitchX) + minX) ;
  }
  float getPeakY(){
    return( (getMaxBin(histoY) * pitchY) + minX) ;
  }


};
  



  class EUTelPreAlign:public marlin::Processor {

  public:
    virtual Processor * newProcessor() {
      return new EUTelPreAlign;
    }
    //! Default constructor
    EUTelPreAlign ();
    virtual void init ();
    virtual void processRunHeader (LCRunHeader * run);
    virtual void processEvent (LCEvent * evt);
    virtual void end();
    virtual bool hitContainsHotPixels( TrackerHitImpl   * hit) ;

    //! Called for first event per run
    /*! Reads hotpixel information from hotPixelCollection into hotPixelMap
     * to be used in the sensor exclusion area logic 
     */
    virtual void  FillHotPixelMap(LCEvent *event);

  private:
    //! Hot pixel collection name.
    /*! 
     * this collection is saved in a db file to be used at the clustering level
     */
    std::string _hotPixelCollectionName;
     
    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  first level key   sensor unique 
     *              value sensor map
     *  sensor map key    unique row number
     *             value  vector of column numbers.
     */
    
    std::map<std::string, bool> _hotPixelMap;
    
    //! bool tag if PreAlign should run anyway or not;
    /*! default 0
     */   
    bool _UsefullHotPixelCollectionFound;

  protected:
    std::string _inputHitCollectionName;
    std::string _alignmentConstantLCIOFile;
    int _iRun;
    int _iEvt;
    int _fixedID;
    float _fixedZ;
    gear::SiPlanesParameters * _siPlanesParameters;
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;
    std::vector<PreAligner> _preAligners;
  };
  //! A global instance of the processor
  EUTelPreAlign gEUTelPreAlign;

}
#endif
#endif // GEAR

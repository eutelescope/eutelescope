#ifndef DeadColumnFinder_h
#define DeadColumnFinder_h 1

#ifdef USE_GEAR

#include "marlin/Processor.h"
#include "lcio.h"
#include <algorithm>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
//#include "IMPL/TrackerHitImpl.h"
//#include "EUTelSparseDataImpl.h"
//#include <EVENT/LCRunHeader.h>
//#include <EVENT/LCEvent.h> 
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"    
#include <UTIL/LCTime.h>
#include <UTIL/CellIDEncoder.h>

#include "TH1.h"
#include "TH2.h"

class DeadColumnFinder : public marlin::Processor { 
public: 
  virtual Processor* newProcessor() {return new DeadColumnFinder;} 
  DeadColumnFinder();
  virtual void init() ;
  virtual void processEvent( LCEvent * evt ) ;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  void bookHistos();
#endif

  virtual void end();

protected:
  std::string _zsDataCollectionName;
  LCCollectionVec *zsInputDataCollectionVec;
  std::string _deadColumnFile;
  std::string _deadColumnCollectionName;
  bool _fillHistos;
private:
  bool _isFirstEvent;
  int _nLayer;
  std::vector<int> _xPixel;
  std::vector<int> _yPixel;
  gear::SiPlanesParameters * _siPlanesParameters;
  gear::SiPlanesLayerLayout * _siPlanesLayerLayout;
  int _nEvent;
  std::vector<std::vector<bool> > isDead;
  std::map<int,TH2I*> hitMap;
};
#endif
#endif

#ifndef EUTelProcessorDeadColumnFinder_h
#define EUTelProcessorDeadColumnFinder_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <algorithm>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include <UTIL/LCTime.h>
#include <UTIL/CellIDEncoder.h>

#include "TH1.h"
#include "TH2.h"

class EUTelProcessorDeadColumnFinder : public marlin::Processor {
public:
  virtual Processor* newProcessor() {return new EUTelProcessorDeadColumnFinder;}
  EUTelProcessorDeadColumnFinder();
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
  int _nEvent;
  std::vector<std::vector<bool> > isDead;
  std::map<int,TH2I*> hitMap;
};
#endif

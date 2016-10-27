#include "EUTelProcessorDeadColumnFinder.h"
#include "EUTELESCOPE.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "marlin/Global.h"

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace eutelescope;
using namespace gear;

EUTelProcessorDeadColumnFinder aEUTelProcessorDeadColumnFinder;

EUTelProcessorDeadColumnFinder::EUTelProcessorDeadColumnFinder()
: Processor("EUTelProcessorDeadColumnFinder"),
  _zsDataCollectionName(""),
  _fillHistos(false),
  _nLayer(0),
  _xPixel(),
  _yPixel(),
  _nEvent(0),
  isDead(0)
  {
    _description="Search of dead columns in the chip";
    registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                           "Input of Zero Suppressed data",
                           _zsDataCollectionName, string ("zsdata") );
    registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling",                     _fillHistos, static_cast< bool > ( true ) );
    registerOptionalParameter("DeadColumnFileName","This is the name of the LCIO file containing the pixels belonging to a dead column",
                           _deadColumnFile, static_cast< string > ( "dead.slcio" ) );
    registerOptionalParameter("DeadColumnCollectionName", "This is the name of the dead column collection",
                           _deadColumnCollectionName, static_cast< string > ( "deadColumn") );
    _isFirstEvent = true;
  }

void EUTelProcessorDeadColumnFinder::init() {
  _nLayer = geo::gGeometry().nPlanes();
  vector<int> tmp(4,0);
  for (int iLayer=0; iLayer<_nLayer; iLayer++ )
  {
    _xPixel.push_back(geo::gGeometry().siPlaneXNpixels(iLayer));
    _yPixel.push_back(geo::gGeometry().siPlaneYNpixels(iLayer));
    vector<bool> isDeadTmp(_xPixel[iLayer],false);
    isDead.push_back(isDeadTmp);
//    cerr << iLayer << "\t" << nFiredPixel[0][iLayer] << endl;
  }
}

void EUTelProcessorDeadColumnFinder::processEvent(LCEvent *evt)
{
  _nEvent++;
//  cerr << evt->getEventNumber() << endl;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  if ( _fillHistos && _isFirstEvent )
  {
    bookHistos();
    _isFirstEvent = false;
  }
#endif
  try
  {
    zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _zsDataCollectionName ) ) ;
  } catch ( lcio::DataNotAvailableException )
  {
//    cerr << "_zsDataCollectionName " << _zsDataCollectionName.c_str() << " not found " << endl;
    return;
  }
  for ( size_t iDetector = 0 ; iDetector < zsInputDataCollectionVec->size(); iDetector++ )
  {
    TrackerDataImpl* zsData = dynamic_cast<TrackerDataImpl*>(zsInputDataCollectionVec->getElementAt(iDetector));
    auto sparseData = EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(zsData);
    for( size_t iPixel = 0; iPixel < sparseData.size(); iPixel++ )
    {
      auto& sparsePixel = sparseData.at(iPixel);
      hitMap[iDetector]->Fill(sparsePixel.getXCoord(), sparsePixel.getYCoord());
      if (iPixel != sparseData.size()-1)
      {
        auto& sparsePixel2 = sparseData.at( iPixel+1 );
        if (sparsePixel.getXCoord() == sparsePixel2.getXCoord() && sparsePixel.getYCoord() == sparsePixel2.getYCoord())
        {
          isDead[iDetector][sparsePixel.getXCoord()] = true;
          if (sparsePixel.getXCoord()%2 == 0) isDead[iDetector][sparsePixel.getXCoord()+1] = true;
          else isDead[iDetector][sparsePixel.getXCoord()-1] = true;
        }
//          cerr << "Same pixel (" << sparsePixel->getXCoord() << ", " << sparsePixel->getYCoord() << ") appearing twice in event " << evt->getEventNumber() << endl;
      }
    }
  }
}

void EUTelProcessorDeadColumnFinder::bookHistos()
{
  for (int iLayer=0; iLayer<_nLayer; iLayer++ )
  {
    hitMap[iLayer] = new TH2I(Form("hitMap_%d",iLayer),Form("Hit map of layer %d",iLayer),_xPixel[iLayer],0,_xPixel[iLayer],_yPixel[iLayer]-1,0,_yPixel[iLayer]-1);
  }
}

void EUTelProcessorDeadColumnFinder::end()
{
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try
  {
    lcWriter->open( _deadColumnFile, LCIO::WRITE_NEW );
  }
  catch ( IOException& e )
  {
    streamlog_out ( ERROR4 ) << e.what() << endl << "Sorry for quitting. " << endl;
    exit(-1);
  }
  LCRunHeaderImpl * lcHeader  = new LCRunHeaderImpl;
  lcHeader->setRunNumber( 0 );

  lcWriter->writeRunHeader(lcHeader);

  delete lcHeader;

  LCEventImpl * event = new LCEventImpl;
  event->setRunNumber( 0 );
  event->setEventNumber( 0 );

  LCTime * now = new LCTime;
  event->setTimeStamp( now->timeStamp() );
  delete now;

  LCCollectionVec * deadColumnCollection = new LCCollectionVec( lcio::LCIO::TRACKERDATA );

  streamlog_out ( MESSAGE5 ) << "Average number of hits per event:" << endl;
  for (int iLayer=0; iLayer<_nLayer; iLayer++ )
    streamlog_out ( MESSAGE5 ) << "Layer " << iLayer << "\t" << (double)hitMap[iLayer]->GetEntries()/_nEvent << endl;
  for (int iLayer=0; iLayer<_nLayer; iLayer++ )
  {
    CellIDEncoder< TrackerDataImpl > deadColumnEncoder  ( eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, deadColumnCollection);
    deadColumnEncoder["sensorID"] = iLayer;
    deadColumnEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;
    auto currentFrame = std::make_unique<lcio::TrackerDataImpl>();
    deadColumnEncoder.setCellID( currentFrame.get() );
    auto sparseFrame = std::make_unique<eutelescope::EUTelTrackerDataInterfacerImpl<eutelescope::EUTelGenericSparsePixel>>(currentFrame.get());
    vector<int> hitPixels(hitMap[iLayer]->GetNbinsX(),0) ;
//    for (int x=0; x<hitMap[iLayer]->GetNbinsX(); x++)
//      hitPixels[x] = 0;
    for (int x=1; x<hitMap[iLayer]->GetNbinsX()+1; x++)
    {
//      hitPixels[x-1] = 0;
      for (int y=1; y<hitMap[iLayer]->GetNbinsY()+1; y++)
      {
        if (hitMap[iLayer]->GetBinContent(x,y) != 0)
        {
//          if ( x == 1) cerr << iLayer << "\t" << y << "\t" << hitMap[iLayer]->GetBinContent(x,y) << endl;
          hitPixels[x-1]++;
//          if (hitPixels > 1) deadColumn = false;
//          cerr << iLayer << "\t" << x << "\t" << y << endl;
//          return;
        }
      }
    }
    for (int x=0; x<hitMap[iLayer]->GetNbinsX(); x++)
    {
//      bool deadColumn = false;
      if (hitPixels[x] <=1 && hitPixels[x+1] <=1)
      {
        isDead[iLayer][x] = true;
        isDead[iLayer][x+1] = true;
//      deadColumn = true;
      }
      else if (x > 0 && x < hitMap[iLayer]->GetNbinsX()-2 && hitPixels[x-1] > 100 && hitPixels[x+2] > 100 && (double)hitPixels[x]/hitPixels[x-1] < 0.7 && (double)hitPixels[x]/hitPixels[x+2] < 0.7 && (double)hitPixels[x+1]/hitPixels[x-1] < 0.7 && (double)hitPixels[x+1]/hitPixels[x+2] < 0.7)
      {
        isDead[iLayer][x] = true;
        isDead[iLayer][x+1] = true;
//      deadColumn = true;
      }
/*      if (deadColumn)
      {
        cerr << "Dead double column found in layer " << iLayer << " at X=" << x << " and " << x+1 << endl;
//        cerr << x-1 << "\t" << x << "\t" <<  x+1 << "\t" << x+2 << endl;
//        cerr << hitPixels[x-1] << "\t" <<  hitPixels[x] << "\t" << hitPixels[x+1] << "\t" <<  hitPixels[x+2] << endl;
      }
*/    }
    deadColumnCollection->push_back( currentFrame.release() );
    for (int x=0; x<_xPixel[iLayer];x++)
    {
      if (isDead[iLayer][x])
      {
        for (int y=0; y<_yPixel[iLayer]; y++)
        {
          EUTelGenericSparsePixel sparsePixel;
          sparsePixel.setXCoord(x);
          sparsePixel.setYCoord(y);
          sparsePixel.setSignal(1);
          sparseFrame->push_back(sparsePixel);
         }
         streamlog_out ( MESSAGE5 ) << "Dead double column found in layer " << iLayer << " at X=" << x << endl;
       }
    }
  }
  event->addCollection( deadColumnCollection, _deadColumnCollectionName);
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
}

#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelPreAlignment.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelVirtualCluster.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelSparseCluster2Impl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes
#include <AIDA/IHistogram3D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <marlin/AIDAProcessor.h>
#endif

// gear includes <.h>
#include "marlin/Global.h"
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// system includes <>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <memory>
#include <cstdio>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gear;

EUTelPreAlign::EUTelPreAlign () :Processor("EUTelPreAlign") {
  _description = "Apply alignment constants to hit collection";

  registerInputCollection (LCIO::TRACKERHIT, "InputHitCollectionName",
                           "The name of the input hit collection",
                           _inputHitCollectionName, string ("hit"));
  registerOptionalParameter ("FixedPlane", "SensorID of fixed plane", _fixedID, 0);
  registerOptionalParameter("AlignmentConstantLCIOFile","Name of LCIO db file where alignment constantds will be stored", 
			    _alignmentConstantLCIOFile, std::string( "alignment.slcio" ) );

  registerOptionalParameter("HotPixelCollectionName", "This is the name of the hot pixel collection to be saved into the output slcio file",
                             _hotPixelCollectionName, static_cast< string > ( "hotpixel_apix" ));

}


void EUTelPreAlign::init () {
  // this method is called only once even when the rewind is active
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;  _iEvt = 0;

  _UsefullHotPixelCollectionFound = 0; 

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));
  for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) {
    if(_siPlanesLayerLayout->getID(iPlane) == _fixedID){ 
      //Get Zpos of ref plane
      _fixedZ = _siPlanesLayerLayout->getSensitivePositionZ(iPlane); 
    } else {
      //Get 
      _preAligners.push_back( PreAligner(  _siPlanesLayerLayout->getSensitivePitchX(iPlane),
					   _siPlanesLayerLayout->getSensitivePitchY(iPlane),
					   _siPlanesLayerLayout->getSensitivePositionZ(iPlane),
					   _siPlanesLayerLayout->getID(iPlane)) );
      printf("iPlane:%5d  %8.3f  %8.3f  %8.3f  %5d  \n",iPlane,
                                           _siPlanesLayerLayout->getSensitivePitchX(iPlane),
					   _siPlanesLayerLayout->getSensitivePitchY(iPlane),
					   _siPlanesLayerLayout->getSensitivePositionZ(iPlane),
					   _siPlanesLayerLayout->getID(iPlane) 
);
    }
  }
}

void EUTelPreAlign::processRunHeader (LCRunHeader * rdr) {
  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl( rdr ) ) ;
  runHeader->addProcessor( type() );
  ++_iRun;
}


void  EUTelPreAlign::FillHotPixelMap(LCEvent *event)
{
    LCCollectionVec *hotPixelCollectionVec = 0;
    try 
    {
      hotPixelCollectionVec = static_cast< LCCollectionVec* > ( event->getCollection( _hotPixelCollectionName  ) );
      streamlog_out ( MESSAGE ) << "_hotPixelCollectionName " << _hotPixelCollectionName.c_str() << " found" << endl; 
    }
    catch (...)
    {
      streamlog_out ( MESSAGE ) << "_hotPixelCollectionName " << _hotPixelCollectionName.c_str() << " not found" << endl; 
      return;
    }

        CellIDDecoder<TrackerDataImpl> cellDecoder( hotPixelCollectionVec );
// 	CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputCollectionVec );
	
        for(int i=0; i<  hotPixelCollectionVec->getNumberOfElements(); i++)
        {
           TrackerDataImpl* hotPixelData = dynamic_cast< TrackerDataImpl *> ( hotPixelCollectionVec->getElementAt( i ) );
	   SparsePixelType  type         = static_cast<SparsePixelType> (static_cast<int> (cellDecoder( hotPixelData )["sparsePixelType"]));
//		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputCollectionVec->getElementAt( i ) );
//		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );


	   int sensorID              = static_cast<int > ( cellDecoder( hotPixelData )["sensorID"] );
           streamlog_out ( MESSAGE ) << "sensorID: " << sensorID << " type " << kEUTelAPIXSparsePixel << " ?= " << type << endl; 

           if( type  == kEUTelAPIXSparsePixel)
           {  
           auto_ptr<EUTelSparseDataImpl<EUTelAPIXSparsePixel > > apixData(new EUTelSparseDataImpl<EUTelAPIXSparsePixel> ( hotPixelData ));
           //auto_ptr<EUTelAPIXSparsePixel> apixPixel( new EUTelAPIXSparsePixel );
	   EUTelAPIXSparsePixel apixPixel;
           //Push all single Pixels of one plane in the apixPixelVec
           for ( unsigned int iPixel = 0; iPixel < apixData->size(); iPixel++ ) 
           {
              std::vector<int> apixColVec();
              apixData->getSparsePixelAt( iPixel, &apixPixel);
//	       apixPixelVec.push_back(new EUTelAPIXSparsePixel(apixPixel));
              streamlog_out ( MESSAGE ) << iPixel << " of " << apixData->size() << " HotPixelInfo:  " << apixPixel.getXCoord() << " " << apixPixel.getYCoord() << " " << apixPixel.getSignal() << " " << apixPixel.getChip() << " " << apixPixel.getTime()<< endl;
              try
              {
                 char ix[100];
                 sprintf(ix, "%d,%d,%d", sensorID, apixPixel.getXCoord(), apixPixel.getYCoord() ); 
                 _hotPixelMap[ix] = true;             
              }
              catch(...)
              {
                 std::cout << "can not add pixel " << std::endl;
//               std::cout << sensorID << " " << apixPixel.getXCoord() << " " << apixPixel.getYCoord() << " " << std::endl;   
              }
           }

           
           }
            else
           {
              _UsefullHotPixelCollectionFound = 0;
           }  	
       }
}

void EUTelPreAlign::processEvent (LCEvent * event) {

  if ( isFirstEvent() )
  {
    FillHotPixelMap(event);
  }

  ++_iEvt;
  if ( _iEvt % 10000 == 0 )
    streamlog_out ( MESSAGE4 ) << "Processing event "
                               << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                               << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber()
                               << setfill(' ')
                               << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);
  
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } else if ( evt->getEventType() == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
                               << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  try {
    LCCollectionVec * inputCollectionVec         = dynamic_cast < LCCollectionVec * > (evt->getCollection(_inputHitCollectionName));
    //Loop over hits in fixed plane
    for (size_t ref = 0; ref < inputCollectionVec->size(); ref++) {
      TrackerHitImpl   * refHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( ref ) ) ;
      const double * refPos = refHit->getPosition();
      if( std::fabs(refPos[2] - _fixedZ) > 2.5) { continue; }
      
      for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) {
	TrackerHitImpl   * hit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
	if( hitContainsHotPixels(hit) ) continue;

        const double * pos = hit->getPosition();
	if( std::fabs(pos[2] - _fixedZ) < 2.5) { continue; }
	bool gotIt(false);
	for(size_t ii = 0; ii < _preAligners.size(); ii++){
	  PreAligner& pa = _preAligners.at(ii);
	  if( std::fabs(pa.getZPos() - pos[2]) > 2.5  ) { continue; }
	  gotIt = true;
	  pa.addPoint(refPos[0] - pos[0], refPos[1] - pos[1]);
//printf("%5d %5d [%5d] %p %5d \n", iHit, ii, _preAligners.size(), pa.current(), pa.getIden()); //, pa.getPeakX(), pa.getPeakY() );
	  break;
	}
	if(not gotIt) {
	  streamlog_out ( ERROR ) << "Mismatched hit at " << pos[2] << endl;
	}
      }
    }
  }
  catch (DataNotAvailableException& e) { 
    streamlog_out  ( WARNING2 ) <<  "No input collection " << _inputHitCollectionName << " found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

  if ( isFirstEvent() ) _isFirstEvent = false;
}

bool EUTelPreAlign::hitContainsHotPixels( TrackerHitImpl   * hit) 
{
        bool skipHit = 0;
//printf("EUTelPreAlign::hitContainsHotPixels \n");
 
        EUTelVirtualCluster * cluster = 0;

        try
        {
            LCObjectVec clusterVector = hit->getRawHits();


            if ( hit->getType() == kEUTelBrickedClusterImpl ) {

               // fixed cluster implementation. Remember it
               //  can come from
               //  both RAW and ZS data
   
                cluster = new EUTelBrickedClusterImpl(static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
                
            } else if ( hit->getType() == kEUTelDFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelDFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } else if ( hit->getType() == kEUTelFFClusterImpl ) {
              
              // fixed cluster implementation. Remember it can come from
              // both RAW and ZS data
              cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl *> ( clusterVector[0] ) );
            } 
            else if ( hit->getType() == kEUTelAPIXClusterImpl ) 
            {
              
//              cluster = new EUTelSparseClusterImpl< EUTelAPIXSparsePixel >
//                 ( static_cast<TrackerDataImpl *> ( clusterVector[ 0 ]  ) );

                // streamlog_out(MESSAGE4) << "Type is kEUTelAPIXClusterImpl" << endl;
                TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );
                // streamlog_out(MESSAGE4) << "Vector size is: " << clusterVector.size() << endl;

                cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);
	      
	        // CellIDDecoder<TrackerDataImpl> cellDecoder(clusterFrame);
                eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel > *apixCluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelAPIXSparsePixel >(clusterFrame);
                
                int sensorID = apixCluster->getDetectorID();//static_cast<int> ( cellDecoder(clusterFrame)["sensorID"] );
//                cout << "Pixels at sensor " << sensorID << ": ";

                for (int iPixel = 0; iPixel < apixCluster->size(); ++iPixel) 
                {
                    int pixelX, pixelY;
                    EUTelAPIXSparsePixel apixPixel;
                    apixCluster->getSparsePixelAt(iPixel, &apixPixel);
                    pixelX = apixPixel.getXCoord();
                    pixelY = apixPixel.getYCoord();
//                    cout << "(" << pixelX << "|" << pixelY << ") ";
//                    cout << endl;

                    try
                    {                       
//                       printf("pixel %3d %3d was found in the _hotPixelMap = %1d (0/1) \n", pixelX, pixelY, _hotPixelMap[sensorID][pixelX][pixelY]  );                       
                       char ix[100];
                       sprintf(ix, "%d,%d,%d", sensorID, apixPixel.getXCoord(), apixPixel.getYCoord() ); 
                       if( _hotPixelMap[ix]  )
                       { 
                          skipHit = true; 	      
//                          streamlog_out(MESSAGE4) << "Skipping hit due to hot pixel content." << endl;
//                          printf("pixel %3d %3d was found in the _hotPixelMap \n", pixelX, pixelY  );
                          delete cluster;                        
                          return true; // if TRUE  this hit will be skipped
                       }
                       else
                       { 
                          skipHit = false; 	      
                       } 
                    } 
                    catch (...)
                    {
//                       printf("pixel %3d %3d was NOT found in the _hotPixelMap \n", pixelX, pixelY  );
                    }
                    if(skipHit ) 
                    {
//                       printf("pixel %3d %3d was found in the _hotPixelMap \n", pixelX, pixelY  );
                    }
                    else
                    { 
//                       printf("pixel %3d %3d was NOT found in the _hotPixelMap \n", pixelX, pixelY  );
                    }

                }

                if (skipHit) 
                {
//                 streamlog_out(MESSAGE4) << "Skipping hit due to hot pixel content." << endl;
//                    continue;
                }
                else
                {
//                  streamlog_out(MESSAGE4) << "Cluster/hit is fine for preAlignment!" << endl;
                }
                
                delete cluster; 
                return skipHit; // if TRUE  this hit will be skipped
            } 
            
       }
       catch(...)
       { 
          // if anything went wrong in the above return FALSE, meaning do not skip this hit
          printf("something went wrong in EUTelPreAlign::hitContainsHotPixels \n"); 
          return 0;
       }

//       if(cluster != 0) delete cluster; 
 
       // if none of the above worked return FALSE, meaning do not skip this hit
       return 0;

}

void EUTelPreAlign::end() {
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
    exit(-1);
  }

  std::cout << "Writing to " << _alignmentConstantLCIOFile << std::endl;

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

  LCCollectionVec * constantsCollection = new LCCollectionVec( LCIO::LCGENERICOBJECT );
  for(size_t ii = 0 ; ii < _preAligners.size(); ii++){
  printf("ii %5d %5d \n", ii, _preAligners.size()  );
  printf("id %5d              \n", _preAligners.at(ii).getIden() );
  printf("id %5d %8.3f        \n", _preAligners.at(ii).getIden(), _preAligners.at(ii).getPeakX()  );
  printf("id %5d %8.3f %8.3f  \n", _preAligners.at(ii).getIden(), _preAligners.at(ii).getPeakX(), _preAligners.at(ii).getPeakY() );



  EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
    constant->setXOffset( -1.0 * _preAligners.at(ii).getPeakX());
    constant->setYOffset( -1.0 * _preAligners.at(ii).getPeakY());
    constant->setSensorID( _preAligners.at(ii).getIden() );
    constantsCollection->push_back( constant );
    streamlog_out ( MESSAGE ) << (*constant) << endl;
    // std::cout << "Iden: " << _preAligners.at(ii).getIden()
    // 	      << " Aligned x: "
    // 	      << _preAligners.at(ii).getPeakX()
    // 	      << " Aligned y:" 
    // 	      << _preAligners.at(ii).getPeakY() << std::endl;
  }
  event->addCollection( constantsCollection, "alignment" );
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
}
#endif

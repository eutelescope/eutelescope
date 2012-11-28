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
#include "marlin/AIDAProcessor.h"

// lcio includes <.h>
#include <IO/LCWriter.h>
#include <UTIL/LCTime.h>
#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitImpl.h>
//#include <TrackerHitImpl2.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// aida includes <.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
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

  registerProcessorParameter ("Events",
                              "How many events are needed to get reasonable approximation to the X,Y shift (pre-alignment)? (default=1000)",
                              _events, static_cast <int> (1000) );

  registerOptionalParameter("ReferenceCollection","reference hit collection name ", _referenceHitCollectionName, static_cast <string> ("reference") );
 
  registerOptionalParameter("ApplyToReferenceCollection","Do you want the reference hit collection to be corrected by the shifts and tilts from the alignment collection? (default - false )",  _applyToReferenceHitCollection, static_cast< bool   > ( false ));
 
  registerOptionalParameter("ResidualsXMin","Minimal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMin, std::vector<float > (static_cast<float > (6), -10.) );

  registerOptionalParameter("ResidualsYMin","Minimal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMin, std::vector<float > (static_cast<float > (6), -10.) );

  registerOptionalParameter("ResidualsXMax","Maximal values of the hit residuals in the X direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsXMax, std::vector<float > (static_cast<float > (6),  10.) );

  registerOptionalParameter("ResidualsYMax","Maximal values of the hit residuals in the Y direction for a correlation band. Note: these numbers are ordered according to the z position of the sensors and NOT according to the sensor id.",_residualsYMax, std::vector<float > (static_cast<float > (6),  10.) );

  registerOptionalParameter ("MinNumberOfCorrelatedHits",
                              "If there are more then this number of correlated hits (planes->track candidate) (default=5)",
                              _minNumberOfCorrelatedHits, static_cast <int> (5) );

}


void EUTelPreAlign::init () {
  // this method is called only once even when the rewind is active
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;  _iEvt = 0;

  _UsefullHotPixelCollectionFound = 0; 

  _referenceHitVec = 0;

  // clear the sensor ID vector
  _sensorIDVec.clear();

  // clear the sensor ID map
  _sensorIDVecMap.clear();
  _sensorIDtoZOrderMap.clear();


  // clear the sensor ID vector (z-axis order)
  _sensorIDVecZOrder.clear();

//
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
//       printf("iPlane:%5d  %8.3f  %8.3f  %8.3f  %5d  \n",iPlane,
//                                            _siPlanesLayerLayout->getSensitivePitchX(iPlane),
// 					   _siPlanesLayerLayout->getSensitivePitchY(iPlane),
// 					   _siPlanesLayerLayout->getSensitivePositionZ(iPlane),
// 					   _siPlanesLayerLayout->getID(iPlane) 
// );
    }
  }


   _siPlaneZPosition = new double[ _siPlanesLayerLayout->getNLayers() ];
   for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); iPlane++ ) 
   {
    _siPlaneZPosition[ iPlane ] = _siPlanesLayerLayout->getLayerPositionZ(iPlane);
    int sensorID = _siPlanesLayerLayout->getID( iPlane );

    _sensorIDVec.push_back( sensorID );
    _sensorIDVecMap.insert( make_pair( sensorID, iPlane ) );

    // count number of the sensors to the left of the current one:
    int _sensors_to_the_left = 0;
    for ( int jPlane = 0 ; jPlane < _siPlanesLayerLayout->getNLayers(); jPlane++ ) 
    {
        if( _siPlanesLayerLayout->getLayerPositionZ(jPlane) + 1e-06 < _siPlaneZPosition[ iPlane ] )
        {
            _sensors_to_the_left++;
        }
    }
    streamlog_out (DEBUG) << " iPlane " << iPlane << " sensor_#_along_Z_axis " << _sensors_to_the_left << "[z= " << setprecision (3) << _siPlaneZPosition[iPlane] << " ] [sensorID " << sensorID << " ]  " << endl;

    _sensorIDVecZOrder.push_back( _sensors_to_the_left );
    _sensorIDtoZOrderMap.insert(make_pair( sensorID, _sensors_to_the_left));
   }

   string tempHistoName = "";
   for(int i=0;i<_sensorIDVecZOrder.size(); i++)
   { 
       tempHistoName =  "PreAlign/hitXCorr_0_to_" + to_string( i ) ;
       AIDA::IHistogram1D* histo1Da  = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(), 100 , -10., 10.) ;
       _hitXCorr.insert( make_pair( i, histo1Da) );
       //printf("XCorr  %5d %p \n", i, histo1Da );

       tempHistoName =  "PreAlign/hitYCorr_0_to_" + to_string( i ) ;
       AIDA::IHistogram1D* histo1Db = AIDAProcessor::histogramFactory(this)->createHistogram1D( tempHistoName.c_str(), 100 , -10., 10.) ;
       _hitYCorr.insert( make_pair( i, histo1Db) );
       //printf("YCorr  %5d %p \n", i, histo1Db );

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
//	        apixPixelVec.push_back(new EUTelAPIXSparsePixel(apixPixel));
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
//                std::cout << sensorID << " " << apixPixel.getXCoord() << " " << apixPixel.getYCoord() << " " << std::endl;   
               }
             }
           }
           else if( type  ==  kEUTelSimpleSparsePixel )
           {  
              auto_ptr<EUTelSparseClusterImpl< EUTelSimpleSparsePixel > > m26Data( new EUTelSparseClusterImpl< EUTelSimpleSparsePixel >   ( hotPixelData ) );

              std::vector<EUTelSimpleSparsePixel*> m26PixelVec;
	      EUTelSimpleSparsePixel m26Pixel;
  	      //Push all single Pixels of one plane in the m26PixelVec

             for ( unsigned int iPixel = 0; iPixel < m26Data->size(); iPixel++ ) 
             {
              std::vector<int> m26ColVec();
              m26Data->getSparsePixelAt( iPixel, &m26Pixel);
//	       apixPixelVec.push_back(new EUTelAPIXSparsePixel(apixPixel));
              streamlog_out ( MESSAGE ) << iPixel << " of " << m26Data->size() << " HotPixelInfo:  " << m26Pixel.getXCoord() << " " << m26Pixel.getYCoord() << " " << m26Pixel.getSignal() << endl;
              try
              {
                 char ix[100];
                 sprintf(ix, "%d,%d,%d", sensorID, m26Pixel.getXCoord(), m26Pixel.getYCoord() ); 
                 _hotPixelMap[ix] = true;             
              }
              catch(...)
              {
                 std::cout << "can not add pixel " << std::endl;
                 std::cout << sensorID << " " << m26Pixel.getXCoord() << " " << m26Pixel.getYCoord() << " " << std::endl;   
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
    if ( _applyToReferenceHitCollection ) 
    {
       try{
       _referenceHitVec = dynamic_cast < LCCollectionVec * > (event->getCollection( _referenceHitCollectionName));
       }
       catch(...)
       {
         _referenceHitVec = 0;
         _applyToReferenceHitCollection = false;
       }
    }
  }

  ++_iEvt;

  if(_iEvt > _events) return;

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

    std::vector<float>   hitX(0.0f);
    std::vector<float>   hitY(0.0f);
    std::vector<PreAligner*>   prealign;


    //Loop over hits in fixed plane
    for (size_t ref = 0; ref < inputCollectionVec->size(); ref++) {
      TrackerHitImpl   * refHit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( ref ) ) ;
      const double * refPos = refHit->getPosition();
      // identify fixed plane
      if( guessSensorID(refPos) != _fixedID ) continue;
      hitX.clear();
      hitY.clear();
      prealign.clear();

      for (size_t iHit = 0; iHit < inputCollectionVec->size(); iHit++) 
      {
	TrackerHitImpl   * hit   = dynamic_cast< TrackerHitImpl * >  ( inputCollectionVec->getElementAt( iHit ) ) ;
	if( hitContainsHotPixels(hit) ) continue;

        const double * pos = hit->getPosition();
        int iHitID = guessSensorID(pos);  
        if( iHitID == _fixedID ) continue;
        bool gotIt(false);
	for(size_t ii = 0; ii < _preAligners.size(); ii++)
        {
	  PreAligner& pa = _preAligners.at(ii);
          if( pa.getIden() != iHitID  ) { continue; }
	  gotIt = true;

          double correlationX =  refPos[0] - pos[0] ;
          double correlationY =  refPos[1] - pos[1] ;

          int idZ = _sensorIDtoZOrderMap[ iHitID ];
          if( 
              (_residualsXMin[idZ] < correlationX ) && ( correlationX < _residualsXMax[idZ]) &&
              (_residualsYMin[idZ] < correlationY ) && ( correlationY < _residualsYMax[idZ]) 
            )
          {
            hitX.push_back( correlationX );
            hitY.push_back( correlationY );
            prealign.push_back(&pa);
          }
	  break;
	}
	if(not gotIt) 
        {
	  streamlog_out ( ERROR ) << "Mismatched hit at " << pos[2] << endl;
	}
      }

      if( prealign.size() > _minNumberOfCorrelatedHits && hitX.size() == hitY.size() )
      {
         for( int ii=0;ii<prealign.size();ii++)
         {  
            prealign[ii]->addPoint( hitX[ii], hitY[ii] );
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
        bool skipHit = false;
 
        EUTelVirtualCluster * cluster = 0;

        try
        {
            LCObjectVec clusterVector = hit->getRawHits();


            if ( hit->getType() == kEUTelSparseClusterImpl ) 
            {
      
              TrackerDataImpl * clusterFrame = static_cast<TrackerDataImpl*> ( clusterVector[0] );
 //              cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel >(clusterFrame);
              eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel > *cluster = new eutelescope::EUTelSparseClusterImpl< eutelescope::EUTelSimpleSparsePixel >(clusterFrame);
              int sensorID = cluster->getDetectorID();
 
              for ( unsigned int iPixel = 0; iPixel < cluster->size(); iPixel++ ) 
              {
                EUTelSimpleSparsePixel m26Pixel;
                cluster->getSparsePixelAt( iPixel, &m26Pixel);
                int pixelX, pixelY;
                pixelX = m26Pixel.getXCoord();
                pixelY = m26Pixel.getYCoord();
//                streamlog_out ( MESSAGE ) << iPixel << " of " << cluster->size() << "real Pixel:  " << pixelX << " " << pixelY << " " <<  endl;
 
                try
                {                       
                       char ix[100];
                       sprintf(ix, "%d,%d,%d", sensorID, pixelX, pixelY ); 
//                       printf("real pixel at %s \n", ix);
//
//  std::map< std::string, bool >::iterator iter = _hotPixelMap.begin();
//  while ( iter != _hotPixelMap.end() ) {
//    std::cout << "Found hot pixel " << iter->second << " at " << iter->first << std::endl;
//    ++iter;
//  }
   
//                       std::cout << _hotPixelMap.first() << " " << _hotPixelMap.second() << std::endl; 
                       std::map<std::string, bool >::const_iterator z = _hotPixelMap.find(ix);
                       if(z!=_hotPixelMap.end() && _hotPixelMap[ix] == true  )
                       { 
                          skipHit = true; 	      
//                          streamlog_out(MESSAGE4) << "Skipping hit due to hot pixel content." << endl;
//                          printf("pixel %3d %3d was found in the _hotPixelMap \n", pixelX, pixelY  );
                          delete cluster;                        			  
                          return true; // if TRUE  this hit will be skipped
                       }
                    } 
                    catch (...)
                    {
//                       printf("pixel %3d %3d was NOT found in the _hotPixelMap \n", pixelX, pixelY  );
                    }
              }
            } else if ( hit->getType() == kEUTelBrickedClusterImpl ) {

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
	 streamlog_out(ERROR4) << "something went wrong in EUTelPreAlign::hitContainsHotPixels " << endl;
          return 0;
       }

       if(cluster) delete cluster; 
 
       // if none of the above worked return FALSE, meaning do not skip this hit
       return 0;

}


int EUTelPreAlign::guessSensorID(const double * hit ) 
{

  int sensorID = -1;
  double minDistance =  numeric_limits< double >::max() ;

  if( _referenceHitVec == 0 || _applyToReferenceHitCollection == false )
  {
    // use z information of planes instead of reference vector
    for ( int iPlane = 0 ; iPlane < _siPlanesLayerLayout->getNLayers(); ++iPlane ) {
      double distance = std::abs( hit[2] - _siPlaneZPosition[ iPlane ] );
      if ( distance < minDistance ) {
	minDistance = distance;
	sensorID = _siPlanesLayerLayout->getID( iPlane );
      }
    }
    if ( minDistance > 30  ) {
      // advice the user that the guessing wasn't successful 
      streamlog_out( WARNING3 ) << "A hit was found " << minDistance << " mm far from the nearest plane\n"
	"Please check the consistency of the data with the GEAR file: hitPosition[2]=" << hit[2] <<       endl;
    }
    
    return sensorID;
  }

      for(size_t ii = 0 ; ii <  _referenceHitVec->getNumberOfElements(); ii++)
      {
        EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
        
        TVector3 hit3d( hit[0], hit[1], hit[2] );
        TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
        TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
 
        double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
//        printf("iPlane %5d   hitPos:  [%8.3f;%8.3f%8.3f]  distance: %8.3f \n", refhit->getSensorID(), hitPosition[0],hitPosition[1],hitPosition[2], distance  );
        if ( distance < minDistance ) 
        {
           minDistance = distance;
           sensorID = refhit->getSensorID();
//           printf("sensorID: %5d \n", sensorID );
        }    

      }

  return sensorID;
}



      
void EUTelPreAlign::end() {
  LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
  try {
    lcWriter->open( _alignmentConstantLCIOFile, LCIO::WRITE_NEW    );
  } catch ( IOException& e ) {
    streamlog_out ( ERROR4 ) << e.what() << endl;
    exit(-1);
  }

  streamlog_out ( MESSAGE ) << "Writing to " << _alignmentConstantLCIOFile << endl;

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


  for(size_t ii=0; ii<_sensorIDVec.size(); ii++)
  {
    streamlog_out(DEBUG) << " sensor " << _sensorIDVec[ii] << endl;
    bool ifound = false;
    for(size_t jj=0; jj< _preAligners.size(); jj++)
    {
      int sensorID = _preAligners.at(jj).getIden();
      if( _sensorIDVec[ii] == sensorID ) { ifound = true; break; }
    }
    if( ifound == false)
    {
      streamlog_out(DEBUG) << " not found ... so, adding " << endl;
      EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
      constant->setXOffset( 0.0 );
      constant->setYOffset( 0.0 );
      constant->setSensorID( _sensorIDVec[ii] );
      constantsCollection->push_back( constant );
      streamlog_out ( MESSAGE ) << (*constant) << endl;
      continue; 
    }
    streamlog_out(DEBUG) << " found in  usedID " << endl;
  }


  for(size_t ii = 0 ; ii < _preAligners.size(); ii++){
   EUTelAlignmentConstant* constant = new EUTelAlignmentConstant();
   if( abs( _preAligners.at(ii).getPeakX() ) <1000. )
      constant->setXOffset( -1.0 * _preAligners.at(ii).getPeakX());
   else
      constant->setXOffset( -1.0 * 0.0                           );
 
   if( abs( _preAligners.at(ii).getPeakY() ) <1000. )
      constant->setYOffset( -1.0 * _preAligners.at(ii).getPeakY());
   else
      constant->setYOffset( -1.0 * 0.0                           );
 
    int sensorID = _preAligners.at(ii).getIden();
    constant->setSensorID( sensorID );
    constantsCollection->push_back( constant );
    streamlog_out ( MESSAGE ) << (*constant) << endl;
    // std::cout << "Iden: " << _preAligners.at(ii).getIden()
    // 	      << " Aligned x: "
    // 	      << _preAligners.at(ii).getPeakX()
    // 	      << " Aligned y:" 
    // 	      << _preAligners.at(ii).getPeakY() << std::endl;
  }


  streamlog_out(DEBUG) << " adding Collection " << "alignment " << endl;
 
  event->addCollection( constantsCollection, "alignment" );
  lcWriter->writeEvent( event );
  delete event;
  lcWriter->close();
}
#endif

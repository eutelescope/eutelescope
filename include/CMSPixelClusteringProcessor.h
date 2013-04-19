// Version: $Id$
/*========================================================================*/
/*          CMSPixel ClusteringProcessor (clustering of zs data)          */
/*          Author: Simon Spannagel                                       */
/*                (simon.spannagel@student.kit.edu or s.spannagel@cern.ch)*/
/*          Created       14 mar 2012                                     */
/*          Last modified 24 apr 2012                                     */
/*========================================================================*/

// based on EUTelAPIXClusteringProcessor:
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $ $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef CMSPIXELCLUSTERINGPROCESSOR_H
#define CMSPIXELCLUSTERINGPROCESSOR_H 1

#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelExceptions.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/EventModifier.h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// lcio includes <.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include <list>

namespace eutelescope {


  class CMSPixelClusteringProcessor : public marlin::Processor , public marlin::EventModifier {

  public:
    
    virtual Processor * newProcessor() {
      return new CMSPixelClusteringProcessor;
    }
    
    CMSPixelClusteringProcessor();
    ~CMSPixelClusteringProcessor() { };
    CMSPixelClusteringProcessor(const CMSPixelClusteringProcessor&);

    void operator=(const CMSPixelClusteringProcessor&);

    virtual const std::string & name() const { return Processor::name() ; }
    virtual void init ();
    virtual void processRunHeader (LCRunHeader * run);
    virtual void processEvent (LCEvent * evt);
    virtual void modifyEvent( LCEvent * evt ) ;
    virtual void check (LCEvent * evt);
    virtual void end();
    
    //! initialize HotPixelMapVec
    /*! values from hotpixel DB file are read into a vector of mapped
     * unique pixel index
     */
    void initializeHotPixelMapVec();

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	void bookHistos();
	void fillHistos(LCEvent * evt);
	static std::string _clusterSignalHistoName;
	static std::string _hitMapHistoName;
	static std::string _chargeMapHistoName;
	static std::string _pixelPerEventHistoName;
	static std::string _clusterPerEventHistoName;
	static std::string _pixelSignalHistoName;
	static std::string _clusterSizeName;
	static std::string _clusterSizeVsChargeName;
	static std::string _clusterSizeXName;
	static std::string _clusterSizeYName;
	static std::string _cluster1pxHistoName;
	static std::string _cluster2pxHistoName;
	static std::string _cluster3pxHistoName;
	static std::string _cluster4pxHistoName;	
	static std::string _clusterMorepxHistoName;
	static std::string _cluster1_2pxHistoName;			
#endif
    void initializeGeometry( LCEvent * event ) throw ( marlin::SkipEventException );
    
    protected:
	void Clustering(LCEvent * evt, LCCollectionVec * pulse);
	std::string _zsDataCollectionName;
	std::string _clusterCollectionName;
    
    int _iRun;
	int _iEvt; 
	bool _isFirstEvent;
	int _iClusters;
	std::vector<int> _iPlaneClusters;
	unsigned int _initialClusterCollectionSize;
	int _minNPixels;
	int _minXDistance;
	int _minYDistance;
	int _minDiagDistance;
	int _minCharge;
	bool _fillHistos;
	
    //! Hot Pixel Collection 
    LCCollectionVec *hotPixelCollectionVec;
    
    //! Vector of map arrays, keeps record of hit pixels 
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created. 
     *  Key <int> is the sequential (counter) id of a hit,
     *  Value <int> - always status value = firing pixel
     */

    std::vector< std::map< int, int > > _hitIndexMapVec;          
    
	unsigned int _noOfDetector;
	bool _isGeometryReady;
	std::vector<int> _sensorIDVec;
	gear::SiPlanesParameters* _siPlanesParameters;
	gear::SiPlanesLayerLayout*_siPlanesLayerLayout;
	std::map< int , int > _layerIndexMap;
	std::vector< int > _orderedSensorIDVec;
	   
	
	std::string _histoInfoFileName;
	std::string _hotPixelCollectionName;
	std::vector<int > _clusterSpectraNVector;
	std::vector<int > _clusterSpectraNxNVector;
	std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;


  };

  //! A global instance of the processor
  CMSPixelClusteringProcessor gCMSPixelClusteringProcessor;

}
#endif // USE_GEAR
#endif

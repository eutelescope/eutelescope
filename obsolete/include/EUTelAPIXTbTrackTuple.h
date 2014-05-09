// Version: $Id$
#ifndef EUTelFitTuple_h
#define EUTelFitTuple_h 1

#include "marlin/Processor.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelReferenceHit.h"


// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCFlagImpl.h>


// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/ITuple.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

#include <gsl/gsl_matrix_double.h>
#include <TFile.h>
#include <TTree.h>
namespace eutelescope {
  class EUTelAPIXTbTrackTuple : public marlin::Processor {
    
  public:
    virtual Processor*  newProcessor() { return new EUTelAPIXTbTrackTuple; }
    
    EUTelAPIXTbTrackTuple() ;
    virtual void init() ;
    virtual void processRunHeader( LCRunHeader* run ) ;
    virtual void processEvent( LCEvent * evt ) ;
    virtual void check( LCEvent * /*evt*/ ){;} ;
    virtual void end() ;
    
  protected:
    //TbTrack additions
    void prepareTree();
    void readAlignment(LCEvent * event);
    void invertAlignment(EUTelAlignmentConstant * alignment);
    void invertGear();
    void reverseAlign(double& x, double& y, double& z, int iden, double nomZPos);
    int  guessSensorID( TrackerHit *hit);
    void clear();
    gsl_matrix* invertLU(int dim, gsl_matrix* m);
    int readZsHits(std::string colName, LCEvent* event);
    int readZsHitsFromClusters(std::string colName, LCEvent* event);
    int readTracks(LCEvent* event);
    int readClusters( std::string colName, LCEvent* event);
    int readHits( std::string hitColName, LCEvent* event );
    void setClusterIdInHits();
    void getDUTRot(EUTelAlignmentConstant * alignment);

    bool _foundAllign;
    bool _doScales;
    std::vector<std::string> _alignColNames;
    std::map<int, bool> _rotationstored;
    int _countrotstored;

    std::map<int, std::vector< std::vector<double> > > _alignShift;
    std::map<int, std::vector< std::vector<double> > > _alignRotations;
    //std::map< int, std::vector< std::vector<double> > >
    std::map<int, std::vector<gsl_matrix*> > _alignRot;
    std::map<int, std::vector<double> > _gearShift;
    std::map<int, std::vector<double> > _gearOffset;
    std::map<int, std::vector<double> > _gearSize;
    std::map<int, gsl_matrix* > _gearRot;
    std::map<int, std::pair<double, double> > _gearPitch;
    std::map<int, std::vector<double> > _gearEulerRot;

    gear::SiPlanesParameters * _siPlanesParameters;
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    std::string _inputTrackColName ;
    std::string _inputTrackerHitColName ;
    std::string _inputTelPulseCollectionName;
    std::string _inputDutPulseCollectionName;
    std::string _telZsColName;
    std::string _dutZsColName;
    bool _clusterBased;

    std::string _path2file;

    //! reference HitCollection name 
    /*!
     */
    std::string _referenceHitCollectionName;
    LCCollectionVec* _referenceHitVec;    
 
    // Internal processor variables
    // ----------------------------

    int _nRun ;
    int _nEvt ;
    int _runNr;
    int _evtNr;

    bool _isFirstEvent;
    
    TFile* _file;

    TTree* _eutracks;
    int _nTrackParams;
    std::vector<double> *_xPos;
    std::vector<double> *_yPos;
    std::vector<double> *_dxdz;
    std::vector<double> *_dydz;
    std::vector<int>    *_trackIden;
    std::vector<int>    *_trackNum;
    std::vector<double> *_chi2;
    std::vector<double> *_ndof;    

    TTree* _zstree;
    int _nPixHits;
    std::vector<int> *p_col;
    std::vector<int> *p_row;
    std::vector<int> *p_tot;
    std::vector<int> *p_iden;
    std::vector<int> *p_lv1;
    std::vector<int> *p_chip;
    std::vector<int> *p_clusterId;

    TTree* _clutree;
    std::vector<int> *_clusize;
    std::vector<int> *_clusizeX;
    std::vector<int> *_clusizeY;
    std::vector<int> *_cluposX;
    std::vector<int> *_cluposY;
    std::vector<int> *_clucharge;
    std::vector<int> *_cluSensorId; 
    std::vector<int> *_cluClusterId;    
    std::map<IMPL::TrackerDataImpl*, int> *_cluPointer;
    //std::vector< >   *_cluPointerToPixHits;
    
    TTree* _euhits;
    int _nHits;
    std::vector<double> *_hitXPos;
    std::vector<double> *_hitYPos;
    std::vector<double> *_hitZPos;
    std::vector<int>    *_hitClusterId;
    std::vector<int>    *_hitSensorId;
    std::vector<IMPL::TrackerDataImpl*>    *_hitPointerToCluster;

    TTree* _rottree;
    std::vector<int> *_rotDUTId;  
    std::vector<double> *_alpha;
    std::vector<double> *_beta; 
    std::vector<double> *_gamma;  
    std::vector<double> *_rotZY;
    std::vector<double> *_rotZX; 
    std::vector<double> *_rotXY;
    std::vector<double> *_rotZYerr;
    std::vector<double> *_rotZXerr; 
    std::vector<double> *_rotXYerr;  
  };


  //! A global instance of the processor.
  EUTelAPIXTbTrackTuple aEUTelAPIXTbTrackTuple;
}

#endif




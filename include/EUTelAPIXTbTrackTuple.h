#ifndef EUTelFitTuple_h
#define EUTelFitTuple_h 1

#include "marlin/Processor.h"
#include "EUTelAlignmentConstant.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"

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
    virtual void check( LCEvent * evt ){;} ;
    virtual void end() ;
    
  protected:
    //TbTrack additions
    void prepareTree();
    void readAlignment(LCEvent * event);
    void invertAlignment(EUTelAlignmentConstant * alignment);
    void invertGear();
    void reverseAlign(double& x, double& y, double& z, int iden);
    void clear();
    gsl_matrix* invertLU(int dim, gsl_matrix* m);
    int readZsHits(std::string colName, LCEvent* event);
    int readTracks(LCEvent* event);
    int readClusters( std::string colName, LCEvent* event);

    bool _foundAllign;
    std::string _dutAlignmentCollectionName;
    std::string _telAlignmentCollectionName;

    std::map<int, std::vector<double> > _alignShift;
    std::map<int, gsl_matrix* > _alignRot;
    std::map<int, std::vector<double> > _gearShift;
    std::map<int, gsl_matrix* > _gearRot;
    std::map<int, std::pair<double, double> > _gearPitch;

    gear::SiPlanesParameters * _siPlanesParameters;
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    std::string _inputTrackColName ;
    std::string _inputTelPulseCollectionName;
    std::string _inputDutPulseCollectionName;
    std::string _telZsColName;
    std::string _dutZsColName;

    std::string _path2file;

    // Internal processor variables
    // ----------------------------

    int _nRun ;
    int _nEvt ;
    int _runNr;
    int _evtNr;
    
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

    TTree* _clutree;
    std::vector<int> *_clusize;
    std::vector<int> *_clusizeX;
    std::vector<int> *_clusizeY;
    std::vector<int> *_cluposX;
    std::vector<int> *_cluposY;
    std::vector<int> *_clucharge;
    std::vector<int> *_cluid;        
  };


  //! A global instance of the processor.
  EUTelAPIXTbTrackTuple aEUTelAPIXTbTrackTuple;
}

#endif




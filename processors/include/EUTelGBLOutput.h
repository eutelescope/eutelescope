/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELGBLOUTPUT_H
#define EUTELGBLOUTPUT_H

// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// system includes <>
#include <string>
#include <vector>
#include <limits>

// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TVectorT.h>

namespace eutelescope {

  class EUTelGBLOutput : public marlin::Processor {

  public:
    virtual Processor *newProcessor() { return new EUTelGBLOutput; }

    EUTelGBLOutput();
    virtual void init();
    virtual void processRunHeader(LCRunHeader *run);
    virtual void processEvent(LCEvent *evt);
    virtual void check(LCEvent * /*evt*/) { ; };
    virtual void end();

  protected:
    void clear();

    std::vector<std::string> _inputHitCollections;
    std::vector<std::string> _inputZsCollections;
    std::string _path2file;

    bool _onlyWithTracks;
    bool _tracksLocalSystem;
    bool _dumpHeader;
    bool _applyCenterShift;
    std::vector<int> _selectedPlanes;
    std::map<int, float> _xShift;
    std::map<int, float> _yShift;

    int _nRun;
    int _nEvt;
    int _runNr;
    int _evtNr;

    TFile *_file;

    TTree *_eutracks;
    int _nTrackParams;
    std::vector<int> *_planeID;
    std::vector<int> *_trackID;
    int _triggerID;
    int _timestamp;
    std::vector<double> *_xPos;
    std::vector<double> *_yPos;
    std::vector<double> *_omega;
    std::vector<double> *_phi;
    std::vector<double> *_kinkx;
    std::vector<double> *_kinky;
    std::vector<double> *_chi2;
    std::vector<int> *_ndof;
    TTree *_euhits;
    int _nHits;
    std::vector<int> *_hitSensorID;
    std::vector<double> *_hitXPos;
    std::vector<double> *_hitYPos;
    std::vector<double> *_hitZPos;
    
    TTree *_zstree;
    int _nPixHits;
    std::vector<int> *_zsID;
    std::vector<int> *_zsX;
    std::vector<int> *_zsY;
    std::vector<double> *_zsSignal;
    std::vector<int> *_zsTime;

    TTree *_versionTree;
    std::vector<double> *_versionNo;
  };

  //! A global instance of the processor.
  EUTelGBLOutput gEUTelGBLOutput;
}
#endif

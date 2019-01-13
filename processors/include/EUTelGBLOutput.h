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
    std::vector<int> _SelectedPlanes;

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
    std::vector<int> *_hitSensorId;
    std::vector<double> *_hitXPos;
    std::vector<double> *_hitYPos;
    std::vector<double> *_hitZPos;
    
    TTree *_zstree;
    std::vector<int> *zs_id;
    std::vector<int> *zs_x;
    std::vector<int> *zs_y;
    std::vector<double> *zs_signal;
    std::vector<int> *zs_time;
  };

  //! A global instance of the processor.
  EUTelGBLOutput gEUTelGBLOutput;
}
#endif

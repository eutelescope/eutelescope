/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
 
#ifndef EUTELAPIXTBTRACKTUPLE_H
#define EUTELAPIXTBTRACKTUPLE_H

// marlin includes ".h"
#include "marlin/Processor.h"

// system includes <>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TVectorT.h>

namespace eutelescope {

  class EUTelAPIXTbTrackTuple : public marlin::Processor {

  public:
    virtual Processor *newProcessor() { return new EUTelAPIXTbTrackTuple; }

    EUTelAPIXTbTrackTuple();
    virtual void init();
    virtual void processRunHeader(LCRunHeader *run);
    virtual void processEvent(LCEvent *evt);
    virtual void check(LCEvent * /*evt*/) { ; };
    virtual void end();

  protected:
    // TbTrack additions
    void prepareTree();
    void clear();

    bool readZsHits(std::string colName, LCEvent *event);
    bool readTracks(LCEvent *event);
    bool readHits(std::string hitColName, LCEvent *event);

    std::string _inputTrackColName;
    std::string _inputTrackerHitColName;
    std::string _inputDutPulseCollectionName;
    std::string _dutZsColName;

    std::string _path2file;

    std::vector<int> _DUTIDs;
    std::map<int, float> _xSensSize;
    std::map<int, float> _ySensSize;

    // Internal processor variables
    // ----------------------------
    int _nRun;
    int _nEvt;
    int _runNr;
    int _evtNr;

    bool _isFirstEvent;

    TFile *_file;

    TTree *_eutracks;
    int _nTrackParams;
    std::vector<double> *_xPos;
    std::vector<double> *_yPos;
    std::vector<double> *_dxdz;
    std::vector<double> *_dydz;
    std::vector<int> *_trackIden;
    std::vector<int> *_trackNum;
    std::vector<double> *_chi2;
    std::vector<double> *_ndof;

    TTree *_zstree;
    int _nPixHits;
    std::vector<int> *p_col;
    std::vector<int> *p_row;
    std::vector<int> *p_tot;
    std::vector<int> *p_iden;
    std::vector<int> *p_lv1;
    std::vector<int> *p_chip;
    std::vector<int> *p_hitTime;
    std::vector<double> *p_frameTime;

    TTree *_euhits;
    int _nHits;
    std::vector<double> *_hitXPos;
    std::vector<double> *_hitYPos;
    std::vector<double> *_hitZPos;
    std::vector<int> *_hitSensorId;

    TTree *_versionTree;
    std::vector<double> *_versionNo;
  };

  //! A global instance of the processor.
  EUTelAPIXTbTrackTuple aEUTelAPIXTbTrackTuple;
}
#endif

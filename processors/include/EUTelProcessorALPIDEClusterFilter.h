/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTelProcessorALPIDEClusterFilter_H
#define EUTelProcessorALPIDEClusterFilter_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

namespace eutelescope {

  class EUTelProcessorALPIDEClusterFilter : public marlin::Processor {

  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelProcessorALPIDEClusterFilter)
    bool _clusterAvailable;
    LCCollectionVec *zsInputDataCollectionVec;
    std::string _zsDataCollectionName;
    std::vector<std::vector<std::vector<std::vector<int>>>>PixelsOfEvents;
    std::vector<std::vector<std::vector<std::vector<int>>>>PixelsOfEventsAfterShift;
    unsigned int _nDeep;
    float _Range;


  public:
    virtual Processor * newProcessor() {
      return new EUTelProcessorALPIDEClusterFilter;
    }

    EUTelProcessorALPIDEClusterFilter ();

  
    virtual void init ();
    virtual bool SameCluster(int iEvent, int iCluster, int jEvent,int jCluster);
    virtual void AddCluster(int iEvent, int iCluster, int jEvent,int jCluster);
    virtual void DeleteCluster(int jEvents, int jCluster);
    virtual void processEvent (LCEvent * evt);
    virtual void end();
    virtual void readCollections(LCCollectionVec * zsInputDataCollectionVec);
    virtual void writeCollection(LCCollectionVec * sparseClusterCollectionVec, LCCollectionVec * pulseCollection);
    virtual void filter();
    
  protected:
    //! Pulse collection size
    size_t _initialPulseCollectionSize;

    //! Pulse collection name.
    /*! This is the name used to store the output cluster
     *  collection.
     */
    std::string _pulseCollectionName;

    std::string _sparseClusterCollectionName;

    int ID;
    bool _shiftedRun;
    int _nLayers;
    std::map< int, int > _totClusterMap;
  };

  EUTelProcessorALPIDEClusterFilter gEUTelProcessorALPIDEClusterFilter;

}
#endif
#endif

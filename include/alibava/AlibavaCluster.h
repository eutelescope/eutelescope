/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVACLUSTER_H_
#define ALIBAVACLUSTER_H_ 1

// alibava includes ".h"
#include "ALIBAVA.h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <vector>

namespace alibava
{
   const float _unrealisticSignal = -1000000;
    class AlibavaCluster {
        std::vector<int> _channums;
        std::vector<float> _signals;
        int _chipNum;
        int _seedChanNum;
        int _clusterID;
        bool _isSensitiveAxisX;
        double _signalPolarity;
    public:
        AlibavaCluster();
        AlibavaCluster(lcio::TrackerDataImpl* trkdata);
        ~AlibavaCluster();
        
        int getChanNum(int imember);
        float getSignal(int imember);
        float getTotalSignal();
        float getTotalSNR(EVENT::FloatVec noiseVec);
        std::vector<float> getSNRs(EVENT::FloatVec noiseVec);
        
        void add(int achannum, float asignal);
        int getClusterSize();
        bool has_seed();
        void print();
       
	float getSignalOnChannel(int channelnum); 
        float getCenterOfGravity();
        
        void createTrackerData(lcio::TrackerDataImpl * alibavaCluster);
        
        std::vector<int> getChanNums();
        std::vector<float> getSignals();
        
        
        float getEta();

        ///////////////////////
        // Setters - Getters //
        ///////////////////////
        
        
        // setter / getter for _chipNum
        int getChipNum();
        void setChipNum(int chipnum);
        
        // setter / getter for _seedChanNum
        int getSeedChanNum();
        void setSeedChanNum(int seedChanNum);
        
        // setter / getter for _clusterID
        int getClusterID();
        void setClusterID(int clusterID);
        
        // setter / getter for _isSensitiveAxisX
        void setIsSensitiveAxisX(bool isSensitiveAxisX);
        bool getIsSensitiveAxisX();
        
        // setter / getter for _signalPolarity
        double getSignalPolarity();
        void setSignalPolarity(double signalPolarity);
        
    protected:
        
    };
    
} // end of alibava namespace

#endif /* ALIBAVACLUSTER_H_ */

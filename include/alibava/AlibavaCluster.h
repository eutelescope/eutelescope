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

using namespace std;
using namespace lcio;

namespace alibava
{
	
	class AlibavaCluster {
		vector<int> _channums;
		vector<float> _signals;
		float _eta;
		int _chipNum;
		int _seedChanNum;
		int _clusterID;
		bool _isSensitiveAxisX;
		int _signalPolarity;
		
	public:
		AlibavaCluster();
		AlibavaCluster(TrackerDataImpl* trkdata);
		~AlibavaCluster();
		
		int getChanNum(int imember);
		float getSignal(int imember);
		float getTotalSignal();
		float getTotalSNR(EVENT::FloatVec noiseVec);
		
		
		void add(int achannum, float asignal);
		int getClusterSize();
		bool has_seed();
		void print();
		
		float getCenterOfGravity();
		
		void createTrackerData(TrackerDataImpl * alibavaCluster);
		
		///////////////////////
		// Setters - Getters //
		///////////////////////

		// setter / getter for _eta
		float getEta();
		void setEta(float eta);
		
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
		int getSignalPolarity();
		void setSignalPolarity(int signalPolarity);

/*
		// setter / getter for _sensorIDOffset
		int getSensorIDOffset();
		void setSensorIDOffset(int sensorIDOffset);
		
		// setter / getter for _missingCorrdinateValue
		int getMissingCorrdinateValue();
		void setMissingCorrdinateValue(int missingCorrdinateValue);
*/
	protected:
		
	};
	
} // end of alibava namespace

#endif /* ALIBAVACLUSTER_H_ */

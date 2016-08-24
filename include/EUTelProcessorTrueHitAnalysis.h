#ifndef EUtelProcessorTrueHitAnalysis_H
#define EUTelProcessorTrueHitAnalysis_H 1

#include "EUTELESCOPE.h"

#include "marlin/Processor.h"

#include <AIDA/IBaseHistogram.h>

#include <IMPL/LCCollectionVec.h>

#include <cmath>
//#include <utility>
#include <string>
#include <vector>
#include <array>
#include <map>

namespace eutelescope {

class EUTelProcessorTrueHitAnalysis : public marlin::Processor {

public:

	virtual Processor* newProcessor() { return new EUTelProcessorTrueHitAnalysis;}

	virtual const std::string & name() const { return Processor::name();}

	EUTelProcessorTrueHitAnalysis();

	virtual void init();

	virtual void processRunHeader(LCRunHeader* rhdr);

	virtual void processEvent(LCEvent* evt);

	virtual void check(LCEvent* evt);

	virtual void end();

	//std::vector<std::pair<double const*, double const*>> pairHits(std::vector<double const*> in_a, std::vector<double const*> in_b);

	int findPairIndex(double const* a, std::vector<double const*> vec);

	void bookHistos();
	void fillHistos(LCEvent* event);

protected:

	std::string _trueHitCollectionName;
	std::string _reconstructedHitCollectionName;
	std::string _rawDataCollectionName;

	int _iRun;
	int _iEvt;
	bool _isFirstEvent;

	std::string _histoInfoFileName;

private:

	int _noOfDetector;
	bool isGeometryReady;
	std::vector<int> _sensorIDVec;

	LCCollectionVec* _trueHitCollectionVec;
	LCCollectionVec* _reconstructedHitCollectionVec;
	LCCollectionVec* _rawDataCollectionVec;

	std::array<std::map<int, AIDA::IBaseHistogram*>, 8> _1DHistos;
	std::array<std::map<int, AIDA::IBaseHistogram*>, 4> _2DHistos;
	/*std::map<int, AIDA::IBaseHistogram*> _xDiffClusterSize1Histos;
	std::map<int, AIDA::IBaseHistogram*> _yDiffClusterSize1Histos;
	std::map<int, AIDA::IBaseHistogram*> _xDiffClusterSize2Histos;
	std::map<int, AIDA::IBaseHistogram*> _yDiffClusterSize2Histos;
	std::map<int, AIDA::IBaseHistogram*> _xDiffClusterSize3Histos;
	std::map<int, AIDA::IBaseHistogram*> _yDiffClusterSize3Histos;
	std::map<int, AIDA::IBaseHistogram*> _xDiffClusterSize4upHistos;
	std::map<int, AIDA::IBaseHistogram*> _yDiffClusterSize4upHistos;*/

	void readCollections(LCEvent* event);
};

EUTelProcessorTrueHitAnalysis gEUTelProcessorTrueHitAnalysis;
}

#endif

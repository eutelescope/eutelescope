#ifndef EUtelProcessorSpuriousClusterFinder_H
#define EUTelProcessorSpuriousClusterFinder_H 1

#include "EUTELESCOPE.h"

#include "marlin/Processor.h"

#include <AIDA/IBaseHistogram.h>

#include <IMPL/LCCollectionVec.h>

#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <map>

namespace eutelescope {

class EUTelProcessorSpuriousClusterFinder : public marlin::Processor {

public:

	virtual Processor* newProcessor() { return new EUTelProcessorSpuriousClusterFinder;}

	virtual const std::string & name() const { return Processor::name();}

	EUTelProcessorSpuriousClusterFinder();

	virtual void init();

	virtual void processRunHeader(LCRunHeader* rhdr);

	virtual void processEvent(LCEvent* evt);

	virtual void check(LCEvent* evt);

	virtual void end();

	int findPairIndex(double const* a, std::vector<double const*> vec);

	void bookHistos();
	void fillHistos(LCEvent* event);

protected:

	std::string _trueHitCollectionName;
	std::string _reconstructedHitCollectionName;

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

	std::array<std::map<int, AIDA::IBaseHistogram*>, 13> _1DHistos;
	std::array<std::map<int, AIDA::IBaseHistogram*>, 9> _2DHistos;

	void readCollections(LCEvent* event);
};

EUTelProcessorSpuriousClusterFinder gEUTelProcessorSpuriousClusterFinder;
}

#endif

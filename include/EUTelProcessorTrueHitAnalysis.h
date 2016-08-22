#ifndef EUtelProcessorTrueHitAnalysis_H
#define EUTelProcessorTrueHitAnalysis_H 1

#include "EUTELESCOPE.h"

#include "marlin/Processor.h"

#include <AIDA/IBaseHistogram.h>

#include <IMPL/LCCollectionVec.h>

#include <cmath>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {

class EUTelProcessorTrueHitAnalysis : public marlin::Processor {

public:

	virtual Processor* newProcessor() { return new EUTelProcessorTrueHitAnalysis;}

	virtual const std::string & name() const { return Processor::name();}

	EUTelProcessorTrueHitAnalysis ();

	virtual void init();

	virtual void processRunHeader(LCRunHeader* rhdr);

	virtual void initializeGeometry(LCEvent* evt);

	virtual void processEvent(LCEvent* evt);

	virtual void check(LCEvent* evt);

	virtual void end();

	std::vector<std::pair<double*, double*>> pairHits(std::vector<double*>& in_a, std::vector<double*>& in_b);

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

	std::map<int, AIDA::IBaseHistogram*> _xDiffHistos;
	std::map<int, AIDA::IBaseHistogram*> _yDiffHistos;

	void readCollections(LCEvent* event);
};

EUTelProcessorTrueHitAnalysis gEUTelProcessorTrueHitAnalysis;
}

#endif

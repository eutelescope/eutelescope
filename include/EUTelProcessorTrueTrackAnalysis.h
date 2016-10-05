#ifndef EUtelProcessorTrueTrackAnalysis_H
#define EUTelProcessorTrueTrackAnalysis_H 1

#include "EUTELESCOPE.h"

#include "marlin/Processor.h"

#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>

#include <IMPL/LCCollectionVec.h>

#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <map>

namespace eutelescope {

class EUTelProcessorTrueTrackAnalysis : public marlin::Processor {

public:

	virtual Processor* newProcessor() { return new EUTelProcessorTrueTrackAnalysis;}

	virtual const std::string & name() const { return Processor::name();}

	EUTelProcessorTrueTrackAnalysis();

	virtual void init();

	virtual void processRunHeader(LCRunHeader* rhdr);

	virtual void processEvent(LCEvent* evt);

	virtual void check(LCEvent* evt);

	virtual void end();

	int findPairIndex(const double* a, std::vector<const double*> vec);

	void bookHistos();
	void fillHistos(LCEvent* event);

protected:

	std::string _trueHitCollectionName;
	std::string _trueFitpointCollectionName;
	std::string _recoFitpointCollectionName;

	int _iRun;
	int _iEvt;

private:

	std::vector<int> _sensorIDVec;

	bool _collectionsAreInitiated;

	std::map<int, std::vector<const double*>> _trueHitMap;
	LCCollectionVec* _trueFitpointCollectionVec;
	LCCollectionVec* _recoFitpointCollectionVec;

	std::array<std::map<int, AIDA::IHistogram1D*>, 4> _1DHistos;
	std::array<std::map<int, AIDA::IHistogram2D*>, 2> _2DHistos;
	std::array<std::map<int, AIDA::IProfile1D*>, 8> _1DProfileHistos;

	void readCollections(LCEvent* event);
};

EUTelProcessorTrueTrackAnalysis gEUTelProcessorTrueTrackAnalysis;
}

#endif

#include "EUTelProcessorTrueHitAnalysis.h"

#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelHistogramManager.h"
#include "EUTelUtility.h"

#include "marlin/Processor.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"

#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>

#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>

#include <cmath>
#include <utility>
#include <string>
#include <vector>
#include <array>
#include <memory>
#include <map>

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelProcessorTrueHitAnalysis::EUTelProcessorTrueHitAnalysis() :
	Processor("EUTelProcessorTrueHitAnalysis"),
	_trueHitCollectionName(""),
	_reconstructedHitCollectionName(""),
	_iRun(0),
	_iEvt(0),
	_isFirstEvent(true),
	_histoInfoFileName(""),
	_noOfDetector(0),
	_sensorIDVec(),
	_trueHitCollectionVec(nullptr),
	_reconstructedHitCollectionVec(nullptr),
	_xDiffHistos(),
	_yDiffHistos()
{

	_description = "EUTelProcessorTrueHitAnalysis compares the true simulated hits with the reconstructed hits.";

	registerInputCollection(LCIO::TRACKERHIT, "TrueHitCollectionName", "Input of True Hit data", _trueHitCollectionName, std::string("true_hits"));

	registerInputCollection(LCIO::TRACKERHIT, "ReconstructedHitCollectionName", "Input of Reconstructed Hit data", _reconstructedHitCollectionName, std::string("hits"));
}

void EUTelProcessorTrueHitAnalysis::init() {

	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);
	_sensorIDVec = geo::gGeometry().sensorIDsVec();

	printParameters();

	_iRun = 0;
	_iEvt = 0;
}

void EUTelProcessorTrueHitAnalysis::processRunHeader(LCRunHeader* rhdr) {

	std::unique_ptr<EUTelRunHeaderImpl> runHeader(new EUTelRunHeaderImpl(rhdr));
	runHeader->addProcessor(type());

	++_iRun;
}

void EUTelProcessorTrueHitAnalysis::readCollections(LCEvent* event) {

	try {

		_trueHitCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_trueHitCollectionName));
		streamlog_out(DEBUG4) << "_trueHitCollectionVec: " << _trueHitCollectionName.c_str() << " found\n";
	}
	catch (lcio::DataNotAvailableException& e) {

		streamlog_out(DEBUG4) << "_trueHitCollectionVec: " << _trueHitCollectionName.c_str() << " not found in event " << event->getEventNumber() << std::endl;
	}

	try {

		_reconstructedHitCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_reconstructedHitCollectionName));
		streamlog_out(DEBUG4) << "_reconstructedHitCollectionVec: " << _reconstructedHitCollectionName.c_str() << "found\n";
	}
	catch (lcio::DataNotAvailableException& e) {

		streamlog_out(DEBUG4) << "_reconstructedHitCollectionVec: " << _reconstructedHitCollectionName.c_str() << " not found in event " << event->getEventNumber() << std::endl;
	}
}

	void EUTelProcessorTrueHitAnalysis::processEvent(LCEvent* event) {

		++_iEvt;

		readCollections(event);

		if (isFirstEvent()) bookHistos();

		EUTelEventImpl* eutelEvent = static_cast<EUTelEventImpl*>(event);
		if (eutelEvent->getEventType() == kEORE) {

			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		}
		else if (eutelEvent->getEventType() == kUNKNOWN) {

			streamlog_out(WARNING2) << "Event number " << eutelEvent->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}

		fillHistos(event);

		_isFirstEvent = false;
	}

	void EUTelProcessorTrueHitAnalysis::check (LCEvent* event) {
	}

	void EUTelProcessorTrueHitAnalysis::end() {

		streamlog_out(MESSAGE4) << "Successfully finished\n";
	}

	std::vector<std::pair<double*, double*>> EUTelProcessorTrueHitAnalysis::pairHits(std::vector<double*>&  in_a, std::vector<double*>&  in_b) {

		std::vector<std::pair<double*, double*>> out;

		for (size_t i = 0; i < in_a.size(); i++) {

			int min_dist_index = 0;
			double min_distance = sqrt(pow(in_a[i][0] - in_b[0][0], 2) + pow(in_a[i][1] - in_b[0][1], 2));

				for (size_t j = 1; j < in_b.size(); j++) {

					if (sqrt(pow(in_a[i][0] - in_b[j][0], 2) + pow(in_a[i][1] - in_b[j][1], 2)) < min_distance) min_dist_index = j;
				}

				out.push_back(std::pair<double*, double*>(in_a[i], in_b[min_dist_index]));

				in_b.erase(in_b.begin()+min_dist_index);
		}

		return std::move(out);
	}

	void EUTelProcessorTrueHitAnalysis::fillHistos(LCEvent* event) {

		EUTelEventImpl* eutelEvent = static_cast<EUTelEventImpl*>(event);
		EventType type = eutelEvent->getEventType();

		if (type == kEORE) {

			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		}
		else if (type == kUNKNOWN) {
			// if it is unknown we had already issued a warning to the user at
			// the beginning of the processEvent. If we get to here, it means
			// that the assumption that the event was a data event was
			// correct, so no harm to continue...
		}

		try {

			LCCollectionVec* _trueHitCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_trueHitCollectionName));
			CellIDDecoder<TrackerHitImpl> trueHitDecoder(_trueHitCollectionVec);

			LCCollectionVec* _reconstructedHitCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_reconstructedHitCollectionName));
			CellIDDecoder<TrackerHitImpl> reconstructedHitDecoder(_reconstructedHitCollectionVec);

			std::map<int, std::vector<double const*>> trueHitMap;
			std::map<int, std::vector<double const*>> reconstructedHitMap;

		for (int i = 0; i < _trueHitCollectionVec->getNumberOfElements(); i++) {

			TrackerHitImpl* trueHit = dynamic_cast<TrackerHitImpl*>(_trueHitCollectionVec->getElementAt(i));
			int detectorID = static_cast<int>(trueHitDecoder(trueHit)["sensorID"]);

			trueHitMap[detectorID].push_back(trueHit->getPosition());
		}

		for (int i = 0; i < _reconstructedHitCollectionVec->getNumberOfElements(); i++) {

			TrackerHitImpl* reconstructedHit = dynamic_cast<TrackerHitImpl*>(_reconstructedHitCollectionVec->getElementAt(i));
			int detectorID = static_cast<int>(reconstructedHitDecoder(reconstructedHit)["sensorID"]);

			reconstructedHitMap.at(detectorID).push_back(reconstructedHit->getPosition());
		}

		for (auto& detectorID: trueHitMap) {

			std::vector<std::pair<double*, double*>> hitPairs = pairHits(trueHitMap.at(detectorID.first), reconstructedHitMap.at(detectorID.first));

			for (size_t i = 0; i < hitPairs.size(); i++) {

				double diff_x = hitPairs[i].first[0] - hitPairs[i].second[0];
				double diff_y = hitPairs[i].first[1] - hitPairs[i].second[1];

				(dynamic_cast<AIDA::IHistogram1D*>(_xDiffHistos.at(detectorID.first)))->fill(diff_x);
				(dynamic_cast<AIDA::IHistogram1D*>(_yDiffHistos.at(detectorID.first)))->fill(diff_y);
			}
		}
	}
	catch (lcio::DataNotAvailableException& e) {

		return;
	}

	return;
}

void EUTelProcessorTrueHitAnalysis::bookHistos() {

	streamlog_out(DEBUG5) << "Booking histograms\n";

	std::unique_ptr<EUTelHistogramManager> histoMng = std>>make_unique<EUTelHistogramManager>(_histoInfoFileName);
	EUTelHistogramInfo* histoInfo;
	bool isHistoManagerAvailable;

	try {
		isHistoManagerAvailable = histoMng->init();
	}
	catch ( std::ios::failure& e) {
		streamlog_out ( WARNING2 ) << "I/O problem with " << _histoInfoFileName << "\n"
                               << "Continuing without histogram manager"  << std::endl;
		isHistoManagerAvailable = false;
	}
	catch ( ParseException& e ) {
		streamlog_out ( WARNING2 ) << e.what() << "\n"
                               << "Continuing without histogram manager" << std::endl;
		isHistoManagerAvailable = false;
	}

	std::string _xDiffHistoName = "xPositionDifference";
	std::string _yDiffHistoName = "yPositionDifference";

	std::string tempHistoName;
	std::string basePath;

	for (size_t i = 0; i < _sensorIDVec.size(); i++) {

		int sensorID = _sensorIDVec[i];
		basePath = "detector_" + to_string(sensorID);

		int xDiffNBin = 501;
		double xDiffMin = -0.002;
		double xDiffMax = 0.1;
		std::string xDiffTitle = "difference in x position between true simulated hits and reconstructed hits";
		if (isHistoManagerAvailable) {

			histoInfo = histoMng->getHistogramInfo(_xDiffHistoName);

			if (histoInfo) {

				streamlog_out(DEBUG2) << (* histoInfo) << std::endl;
				xDiffNBin = histoInfo->_xBin;
				xDiffMin = histoInfo->_xMin;
				xDiffMax = histoInfo->_xMax;

				if (histoInfo->title != "") xDiffTitle = histoInfo->_title;
			}
		}

		int yDiffNBin = 501;
		double yDiffMin = -0.002;
		double yDiffMax = 0.1;
		std::string yDiffTitle = "difference in y position between true simulated hits and reconstructed hits";
		if (isHistoManagerAvailable) {

			histoInfo = histoMng->getHistogramInfo(_yDiffHistoName);

			if (histoInfo) {

				streamlog_out(DEBUG2) << (* histoInfo) << std::endl;
				yDiffNBin = histoInfo->_xBin;
				yDiffMin = histoInfo->_xMin;
				yDiffMax = histoInfo->_xMax;

				if (histoInfo->title != "") yDiffTitle = histoInfo->_title;
			}
		}

		tempHitsoName = _xDiffHistoName + "_d" + to_string(sensorID);
		_xDiffHistos.insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+tempHistoName).c_str(), xDiffNBin, xDiffMin, xDiffMax)));
		_xDiffHistos[sensorID]->setTitle(xDiffTitle.c_str());

		tempHitsoName = _yDiffHistoName + "_d" + to_string(sensorID);
		_yDiffHistos.insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+tempHistoName).c_str(), yDiffNBin, yDiffMin, yDiffMax)));
		_yDiffHistos[sensorID]->setTitle(yDiffTitle.c_str());
	}

	streamlog_out(DEBUG5) << "end of Booking histograms\n";
}

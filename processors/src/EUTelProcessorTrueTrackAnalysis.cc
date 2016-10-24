#include "EUTelProcessorTrueTrackAnalysis.h"

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
#include <AIDA/IProfile1D.h>
#include <AIDA/ITree.h>

#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/LCCollectionVec.h>

#include <cmath>
#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <array>
#include <memory>
#include <map>

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelProcessorTrueTrackAnalysis::EUTelProcessorTrueTrackAnalysis() :
	Processor("EUTelProcessorTrueTrackAnalysis"),
	_trueHitCollectionName(""),
	_trueFitpointCollectionName(""),
	_recoFitpointCollectionName(""),
	_iRun(0),
	_iEvt(0),
	_sensorIDVec(),
	_collectionsAreInitiated(true),
	_trueHitMap(),
	_trueFitpointCollectionVec(nullptr),
	_recoFitpointCollectionVec(nullptr),
	_1DHistos(),
	_2DHistos(),
	_1DProfileHistos()
{

	_description = "EUTelProcessorTrueTrackAnalysis compares the reconstructed track fitpoints with the true hits and true track fitpoints";

	registerInputCollection(LCIO::TRACKERHIT, "TrueHitCollectionName", "Input of true hit data", _trueHitCollectionName, std::string("true_hits"));

	registerInputCollection(LCIO::TRACKERHIT, "TrueTrackFitpointCollectionName", "Input of true track fitpoints", _trueFitpointCollectionName, std::string("true_track_fitpoints"));

	registerInputCollection(LCIO::TRACKERHIT, "ReconstructedTrackFitpointCollectionName", "Input of reconstructed track fitpoints", _recoFitpointCollectionName, std::string("fitpoints"));
}

void EUTelProcessorTrueTrackAnalysis::init() {

	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);
	_sensorIDVec = geo::gGeometry().sensorIDsVec();

	for (size_t i = 0; i < _sensorIDVec.size(); i++) {

		_trueHitMap.insert(std::make_pair(_sensorIDVec[i], std::vector<const double*>()));
	}

	printParameters();

	bookHistos();

	_iRun = 0;
	_iEvt = 0;
}

void EUTelProcessorTrueTrackAnalysis::processRunHeader(LCRunHeader* rhdr) {

	std::unique_ptr<EUTelRunHeaderImpl> runHeader(new EUTelRunHeaderImpl(rhdr));
	runHeader->addProcessor(type());

	++_iRun;
}

void EUTelProcessorTrueTrackAnalysis::readCollections(LCEvent* event) {

	//monitor whether all input collections are found
	_collectionsAreInitiated = true;

	//read in the true hits
	try {

		LCCollectionVec* trueHitCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_trueHitCollectionName));
		CellIDDecoder<TrackerHitImpl> trueHitDecoder(trueHitCollectionVec);

		//fill the map with the true hits
		for (int i = 0; i < trueHitCollectionVec->getNumberOfElements(); i++) {

			TrackerHitImpl* trueHit = dynamic_cast<TrackerHitImpl*>(trueHitCollectionVec->getElementAt(i));

			if (trueHit->getQuality() == 1) {

				int sensorID = static_cast<int>(trueHitDecoder(trueHit)["sensorID"]);
				_trueHitMap.at(sensorID).push_back(trueHit->getPosition());
			}
		}

		streamlog_out(DEBUG4) << "trueHitCollectionVec: " << _trueHitCollectionName.c_str() << " found\n";
	}
	catch (lcio::DataNotAvailableException& e) {

		streamlog_out(WARNING4) << "trueHitCollectionVec: " << _trueHitCollectionName.c_str() << " not found in event " << event->getEventNumber() << std::endl;

		_collectionsAreInitiated = false;
	}

	//read in true track fitpoints
	try {

		_trueFitpointCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_trueFitpointCollectionName));
		streamlog_out(DEBUG4) << "_trueFitpointCollectionVec: " << _trueFitpointCollectionName.c_str() << "found\n";
	}
	catch (lcio::DataNotAvailableException& e) {

		streamlog_out(WARNING4) << "_trueFitpointCollectionVec: " << _trueFitpointCollectionName.c_str() << " not found in event " << event->getEventNumber() << std::endl;

		_collectionsAreInitiated = false;
	}

	//read in reconstructed track fitpoints
	try {

		_recoFitpointCollectionVec = dynamic_cast<LCCollectionVec*>(event->getCollection(_recoFitpointCollectionName));
		streamlog_out(DEBUG4) << "_recoFitpointCollectionVec: " << _recoFitpointCollectionName.c_str() << "found\n";
	}
	catch (lcio::DataNotAvailableException& e) {

		streamlog_out(WARNING4) << "_recoFitpointCollectionVec: " << _recoFitpointCollectionName.c_str() << " not found in event " << event->getEventNumber() << std::endl;

		_collectionsAreInitiated = false;
	}
}

void EUTelProcessorTrueTrackAnalysis::processEvent(LCEvent* event) {

	++_iEvt;

	readCollections(event);

	EUTelEventImpl* eutelEvent = static_cast<EUTelEventImpl*>(event);
	if (eutelEvent->getEventType() == kEORE) {

		streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
		return;
	}
	else if (eutelEvent->getEventType() == kUNKNOWN) {

		streamlog_out(WARNING2) << "Event number " << eutelEvent->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	if (_collectionsAreInitiated) fillHistos(event);

	for (size_t i = 0; i < _sensorIDVec.size(); i++) {

		_trueHitMap.at(_sensorIDVec[i]).clear();
	}
}

void EUTelProcessorTrueTrackAnalysis::check (LCEvent*) {
}

void EUTelProcessorTrueTrackAnalysis::end() {

	streamlog_out(MESSAGE4) << "Successfully finished\n";
}

int EUTelProcessorTrueTrackAnalysis::findPairIndex(const double* a, std::vector<const double*> vec) {

	if (vec.size() == 1) return 0;
	else {

		int min_dist_index = 0;
		double min_distance = sqrt(pow(a[0] - vec[0][0], 2) + pow(a[1] - vec[0][1], 2));

		for (size_t i = 1; i < vec.size(); i++) {

			double current_distance = sqrt(pow(a[0] - vec[i][0], 2) + pow(a[1] - vec[i][1], 2));

			if (current_distance < min_distance) {

				min_dist_index = i;
				min_distance = current_distance;
			}
		}

		return min_dist_index;
	}
}

void EUTelProcessorTrueTrackAnalysis::fillHistos(LCEvent* event) {

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

	CellIDDecoder<TrackerHitImpl> fitpointDecoder(EUTELESCOPE::HITENCODING);
	//CellIDDecoder<TrackerHitImpl> recoFitpointDecoder(_recoFitpointCollectionVec);
	//CellIDDecoder<TrackerDataImpl> rawDataDecoder("sensorID:7,sparsePixelType:5");

	//load reconstructed fitpoints and find the closest true hit from the _trueHitMap
	for (int i = 0; i < _recoFitpointCollectionVec->getNumberOfElements(); i++) {

		TrackerHitImpl* recoFitpoint = dynamic_cast<TrackerHitImpl*>(_recoFitpointCollectionVec->getElementAt(i));
		int sensorID = static_cast<int>(fitpointDecoder(recoFitpoint)["sensorID"]);

		if (_trueHitMap.at(sensorID).size() == 0) {

			streamlog_out(WARNING2) << "found an unpaired hit at event " << event->getEventNumber() << " and at plane " << sensorID << std::endl;

			continue;
		}

		/*TrackerDataImpl* zsData = dynamic_cast<TrackerDataImpl*>(recoFitpoint->getRawHits().front());

		int pixelType = static_cast<int>(rawDataDecoder(zsData)["sparsePixelType"]);
		int pixelEntries;
		if (pixelType == 1) pixelEntries = 3;//EUTelSimpleSparsePixel
		else if (pixelType == 2) pixelEntries = 4;//EUTelGenericSparsePixel
		else if (pixelType == 3) pixelEntries = 8;//EUTelGeometricPixel
		else if (pixelType == 4) pixelEntries = 7;//EUTelMuPixel
		else {

			streamlog_out(ERROR3) << "Unrecognized pixel type in reconstructed hit raw data at event " << event->getEventNumber() << std::endl;

			return;
		}

		std::vector<float> chargeValues = zsData->getChargeValues();
		std::vector<float> XHitPixels, YHitPixels;
		XHitPixels.push_back(chargeValues[0]);
		YHitPixels.push_back(chargeValues[1]);
		for (size_t j = pixelEntries; j < chargeValues.size(); j += pixelEntries) {

			//add the hit pixel index to the appropriate vector if it is not yet stored there
			if (std::find(XHitPixels.begin(), XHitPixels.end(), chargeValues[j]) == XHitPixels.end()) XHitPixels.push_back(chargeValues[j]);
			if (std::find(YHitPixels.begin(), YHitPixels.end(), chargeValues[j+1]) == YHitPixels.end()) YHitPixels.push_back(chargeValues[j+1]);
		}
		int clusterSizeX = static_cast<int>(XHitPixels.size());
		int clusterSizeY = static_cast<int>(YHitPixels.size());
		if (clusterSizeX > 4) clusterSizeX = 4;
		if (clusterSizeY > 4) clusterSizeY = 4;*/

		const double* fitpointPos = recoFitpoint->getPosition();
		int pairIndex = findPairIndex(fitpointPos, static_cast<std::vector<const double*>>(_trueHitMap.at(sensorID)));
		const double* pair = (static_cast<std::vector<const double*>>(_trueHitMap.at(sensorID)))[pairIndex];

		double diff_x = (fitpointPos[0] - pair[0])*1000;
		double diff_y = (fitpointPos[1] - pair[1])*1000;

		//fill histograms for all sensors
		_1DHistos[0].at(sensorID)->fill(diff_x);
		_1DHistos[1].at(sensorID)->fill(diff_y);
		_1DProfileHistos[0].at(sensorID)->fill(fitpointPos[0], diff_x);
		_1DProfileHistos[1].at(sensorID)->fill(fitpointPos[0], diff_y);
		_1DProfileHistos[2].at(sensorID)->fill(fitpointPos[1], diff_x);
		_1DProfileHistos[3].at(sensorID)->fill(fitpointPos[1], diff_y);
		_2DHistos[0].at(sensorID)->fill(diff_x, diff_y, 1.0);

		_trueHitMap.at(sensorID).erase(_trueHitMap.at(sensorID).begin()+pairIndex);
	}

	//load true track fitpoints into map
	std::map<int, std::vector<const double*>> trueFitpointMap;
	for (size_t i = 0; i < _sensorIDVec.size(); i++) {

		trueFitpointMap.insert(std::make_pair(_sensorIDVec[i], std::vector<const double*>()));
	}

	for (int i = 0; i < _trueFitpointCollectionVec->getNumberOfElements(); i++) {

		TrackerHitImpl* trueFitpoint = dynamic_cast<TrackerHitImpl*>(_trueFitpointCollectionVec->getElementAt(i));
		int sensorID = static_cast<int>(fitpointDecoder(trueFitpoint)["sensorID"]);

		trueFitpointMap.at(sensorID).push_back(trueFitpoint->getPosition());
	}

	//load reconstructed fitpoints and compare them to the closest true fitpoint
	for (int i = 0; i < _recoFitpointCollectionVec->getNumberOfElements(); i++) {

		TrackerHitImpl* recoFitpoint = dynamic_cast<TrackerHitImpl*>(_recoFitpointCollectionVec->getElementAt(i));
		int sensorID = static_cast<int>(fitpointDecoder(recoFitpoint)["sensorID"]);

		if (trueFitpointMap.at(sensorID).size() == 0) {

			streamlog_out(WARNING2) << "found an unpaired fitpoint at event " << event->getEventNumber() << " and at plane " << sensorID << std::endl;

			continue;
		}

		const double* fitpointPos = recoFitpoint->getPosition();
		int pairIndex = findPairIndex(fitpointPos, static_cast<std::vector<const double*>>(trueFitpointMap.at(sensorID)));
		const double* pair = (static_cast<std::vector<const double*>>(trueFitpointMap.at(sensorID)))[pairIndex];

		double diff_x = (fitpointPos[0] - pair[0])*1000;
		double diff_y = (fitpointPos[1] - pair[1])*1000;

		//fill histograms for all sensors
		_1DHistos[2].at(sensorID)->fill(diff_x);
		_1DHistos[3].at(sensorID)->fill(diff_y);
		_1DProfileHistos[4].at(sensorID)->fill(fitpointPos[0], diff_x);
		_1DProfileHistos[5].at(sensorID)->fill(fitpointPos[0], diff_y);
		_1DProfileHistos[6].at(sensorID)->fill(fitpointPos[1], diff_x);
		_1DProfileHistos[7].at(sensorID)->fill(fitpointPos[1], diff_y);
		_2DHistos[1].at(sensorID)->fill(diff_x, diff_y, 1.0);

		trueFitpointMap.at(sensorID).erase(trueFitpointMap.at(sensorID).begin()+pairIndex);
	}

	return;
}

void EUTelProcessorTrueTrackAnalysis::bookHistos() {

	streamlog_out(DEBUG5) << "Booking histograms\n";

	for (size_t i = 0; i < _sensorIDVec.size(); i++) {

		int sensorID = _sensorIDVec[i];
		std::string basePath = "detector_" + to_string(sensorID);
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		//set up true hit residual histograms 
		//set up 1D histograms
		for (size_t j = 0; j < _1DHistos.size()/2; j++) {

			std::string coordinate = (j%2 == 0)?"x":"y";

			std::string histoName = coordinate + "TrueHitResidual_d" + to_string(sensorID);

			int histoNBin = 201;
			double histoMin = -50;
			double histoMax = 50;
			std::string histoTitle = "difference in " + coordinate + " position between true simulated hits and reconstructed fitpoints from the reconstructed track;#Delta" + coordinate + " /#mum;Entries";

			_1DHistos[j].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histoNBin, histoMin, histoMax)));
			_1DHistos[j].at(sensorID)->setTitle(histoTitle.c_str());
		}
		
		//set up 2D histograms
		for (size_t j = 0; j < _2DHistos.size()/2; j++) {

			std::string histoName = "2DTrueHitResidual_d" + to_string(sensorID);

			int histoXNBin = 201;
			double histoXMin = -50;
			double histoXMax = 50;

			int histoYNBin = 201;
			double histoYMin = -50;
			double histoYMax = 50;
			std::string histoTitle = "difference in position between true simulated hits and reconstructed fitpoints from reconstructed tracks;#Deltax /#mum;#Deltay /#mum;Entries";

			_2DHistos[j].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histoXNBin, histoXMin, histoXMax, histoYNBin, histoYMin, histoYMax)));
			_2DHistos[j].at(sensorID)->setTitle(histoTitle.c_str());
		}

		//set up 1D profile histograms
		for (size_t j = 0; j < _1DProfileHistos.size()/2; j++) {

			std::string coordinate_1 = (j%2 == 0)?"x":"y";
			std::string coordinate_2 = (j<_1DProfileHistos.size()/4)?"x":"y";

			std::string histoName = "TrueHitProfile_" + coordinate_1 + "Resid_vs_" + coordinate_2 + "Pos_d" + to_string(sensorID);

			int profileXNBin = 201;
			double profileXMin = -10;
			double profileXMax = 10;

			double profileYMin = -100;
			double profileYMax = 100;
			std::string histoTitle = "difference in " + coordinate_1 + " position between true simulated hits and reconstructed fitpoints from reconstructed tracks against " + coordinate_2 + " position of reconstructed fitpoints;" + coordinate_2 + " position /#mum;#Delta" + coordinate_1 + " /#mum";

			_1DProfileHistos[j].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createProfile1D((basePath+histoName).c_str(), profileXNBin, profileXMin, profileXMax, profileYMin, profileYMax)));
			_1DProfileHistos[j].at(sensorID)->setTitle(histoTitle.c_str());
		}

		//set up histograms fitpoints residual histograms
		//set up 1D histograms
		for (size_t j = _1DHistos.size()/2; j < _1DHistos.size(); j++) {

			std::string coordinate = (j%2 == 0)?"x":"y";

			std::string histoName = coordinate + "FitpointResidual_d" + to_string(sensorID);

			int histoNBin = 201;
			double histoMin = -16;
			double histoMax = 16;
			std::string histoTitle = "difference in " + coordinate + " position between true fitpoints from the true track and reconstructed fitpoints from the reconstructed track;#Delta" + coordinate + " /#mum;Entries";

			_1DHistos[j].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histoNBin, histoMin, histoMax)));
			_1DHistos[j].at(sensorID)->setTitle(histoTitle.c_str());
		}

		//set up 2D histograms
		for (size_t j = _2DHistos.size()/2; j < _2DHistos.size(); j++) {

			std::string histoName = "2DFitpointResidual_d" + to_string(sensorID);

			int histoXNBin = 201;
			double histoXMin = -16;
			double histoXMax = 16;

			int histoYNBin = 201;
			double histoYMin = -16;
			double histoYMax = 16;
			std::string histoTitle = "difference in position between true fitpoints from true tracks and reconstructed fitpoints from reconstructed tracks;#Deltax /#mum;#Deltay /#mum;Entries";

			_2DHistos[j].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histoXNBin, histoXMin, histoXMax, histoYNBin, histoYMin, histoYMax)));
			_2DHistos[j].at(sensorID)->setTitle(histoTitle.c_str());
		}

		//set up 1D profile histograms
		for (size_t j = _1DProfileHistos.size()/2; j < _1DProfileHistos.size(); j++) {

			std::string coordinate_1 = (j%2 == 0)?"x":"y";
			std::string coordinate_2 = (j<3*_1DProfileHistos.size()/4)?"x":"y";

			std::string histoName = "FitpointProfile_" + coordinate_1 + "Resid_vs_" + coordinate_2 + "Pos_d" + to_string(sensorID);

			int profileXNBin = 201;
			double profileXMin = -10;
			double profileXMax = 10;

			double profileYMin = -100;
			double profileYMax = 100;
			std::string histoTitle = "difference in " + coordinate_1 + " position between true fitpoints from true tracks and reconstructed fitpoints from reconstructed tracks against " + coordinate_2 + " position of reconstructed fitpoints;" + coordinate_2 + " position /#mum;#Delta" + coordinate_1 + " /#mum";

			_1DProfileHistos[j].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createProfile1D((basePath+histoName).c_str(), profileXNBin, profileXMin, profileXMax, profileYMin, profileYMax)));
			_1DProfileHistos[j].at(sensorID)->setTitle(histoTitle.c_str());
		}
	}

	streamlog_out(DEBUG5) << "end of Booking histograms\n";
}

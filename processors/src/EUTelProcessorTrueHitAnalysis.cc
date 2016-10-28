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
#include <algorithm>
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
	_sensorIDVec(),
	_trueHitCollectionVec(nullptr),
	_reconstructedHitCollectionVec(nullptr),
	_1DHistos(),
	_2DHistos(),
	_xClustSize2Histos(),
	_yClustSize3Histos()
{

	_description = "EUTelProcessorTrueHitAnalysis compares the true simulated hits with the reconstructed hits";

	registerInputCollection(LCIO::TRACKERHIT, "TrueHitCollectionName", "Input of True Hit data", _trueHitCollectionName, std::string("true_hits"));

	registerInputCollection(LCIO::TRACKERHIT, "ReconstructedHitCollectionName", "Input of Reconstructed Hit data", _reconstructedHitCollectionName, std::string("hit"));
}

void EUTelProcessorTrueHitAnalysis::init() {

	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);
	_sensorIDVec = geo::gGeometry().sensorIDsVec();

	printParameters();

	_iRun = 0;
	_iEvt = 0;

	bookHistos();
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

	EUTelEventImpl* eutelEvent = static_cast<EUTelEventImpl*>(event);
	if (eutelEvent->getEventType() == kEORE) {

		streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
		return;
	}
	else if (eutelEvent->getEventType() == kUNKNOWN) {

		streamlog_out(WARNING2) << "Event number " << eutelEvent->getEventNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
	}

	fillHistos(event);
}

void EUTelProcessorTrueHitAnalysis::check (LCEvent*) {
}

void EUTelProcessorTrueHitAnalysis::end() {

	streamlog_out(MESSAGE4) << "Successfully finished\n";
}

int EUTelProcessorTrueHitAnalysis::findPairIndex(double const* a, std::vector<double const*> vec) {

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

void EUTelProcessorTrueHitAnalysis::fillHistos(LCEvent* event) {

	EUTelEventImpl* eutelEvent = static_cast<EUTelEventImpl*>(event);
	EventType type = eutelEvent->getEventType();

	if (type == kEORE) {

		streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
		return;
	}
	else if (type == kUNKNOWN) {
		//comment from EUTelProcessorGeometricClustering.cc:
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
		CellIDDecoder<TrackerDataImpl> rawDataDecoder("sensorID:7,sparsePixelType:5");

		std::map<int, std::vector<double const*>> trueHitMap;

		for (size_t i = 0; i < _sensorIDVec.size(); i++) {

			trueHitMap.insert(std::make_pair(_sensorIDVec[i], std::vector<const double*>()));
		}

		//fill the map with the true hits
		for (int i = 0; i < _trueHitCollectionVec->getNumberOfElements(); i++) {

			TrackerHitImpl* trueHit = dynamic_cast<TrackerHitImpl*>(_trueHitCollectionVec->getElementAt(i));
			int detectorID = static_cast<int>(trueHitDecoder(trueHit)["sensorID"]);

			trueHitMap.at(detectorID).push_back(trueHit->getPosition());
		}

		//load reconstructed hits and find the closest true hit from the trueHitMap
		for (int i = 0; i < _reconstructedHitCollectionVec->getNumberOfElements(); i++) {

			TrackerHitImpl* reconstructedHit = dynamic_cast<TrackerHitImpl*>(_reconstructedHitCollectionVec->getElementAt(i));
			int detectorID = static_cast<int>(reconstructedHitDecoder(reconstructedHit)["sensorID"]);
			if (trueHitMap.at(detectorID).size() == 0) {

				streamlog_out(WARNING2) << "found an unpaired hit at event " << event->getEventNumber() << std::endl;
				continue;
			}

			TrackerDataImpl* zsData = dynamic_cast<TrackerDataImpl*>(reconstructedHit->getRawHits().front());

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
				if (std::find(XHitPixels.begin(), XHitPixels.end(), chargeValues[j]) == XHitPixels.end()) {

					XHitPixels.push_back(chargeValues[j]);
				}
				if (std::find(YHitPixels.begin(), YHitPixels.end(), chargeValues[j+1]) == YHitPixels.end()) {

					YHitPixels.push_back(chargeValues[j+1]);
				}
			}
			int clusterSizeX = static_cast<int>(XHitPixels.size());
			int clusterSizeY = static_cast<int>(YHitPixels.size());
			if (clusterSizeX > 4) clusterSizeX = 4;
			if (clusterSizeY > 4) clusterSizeY = 4;

			double const* hitPos = reconstructedHit->getPosition();
			int pairIndex = findPairIndex(hitPos, static_cast<std::vector<double const*>>(trueHitMap.at(detectorID)));
			double const* pair = (static_cast<std::vector<double const*>>(trueHitMap.at(detectorID)))[pairIndex];

			double diff_x = (hitPos[0] - pair[0])*1000;
			double diff_y = (hitPos[1] - pair[1])*1000;

			//fill histograms for all sensors
			_1DHistos[2*(clusterSizeX-1)].at(detectorID)->fill(diff_x);
			_1DHistos[2*clusterSizeY-1].at(detectorID)->fill(diff_y);
			_2DHistos[2*(clusterSizeX-1)].at(detectorID)->fill(diff_x, diff_y, 1.0);
			_2DHistos[2*clusterSizeY-1].at(detectorID)->fill(diff_x, diff_y, 1.0);

			//fill histograms for x cluster size 2 for DUT
			if ((detectorID >= 10) && (clusterSizeX == 2)) {

				//fill seperate histograms for the case where there are only two hit pixels
				if (clusterSizeY == 1) {

					int TOTDiff = abs(chargeValues[2] - chargeValues[2+pixelEntries]);
					if (TOTDiff > 4) TOTDiff = 4;
					_xClustSize2Histos[2*TOTDiff]->fill(diff_x);
				}

				//in for general y cluster sizes, take the TOT difference as the
				//difference between the total TOT value of all pixels with the same x-index
				int leftTOT = chargeValues[2];
				int rightTOT = 0;
				for (size_t j = pixelEntries; j < chargeValues.size(); j += pixelEntries) {

					if (chargeValues[j] == chargeValues[0]) leftTOT += chargeValues[j+2];
					else rightTOT += chargeValues[j+2];
				}

				int TOTDiff = abs(leftTOT - rightTOT);
				if (TOTDiff > 4) TOTDiff = 4;
				_xClustSize2Histos[2*TOTDiff+1]->fill(diff_x);
			}

			if ((detectorID == 0) && (clusterSizeY == 3)) {

				_yClustSize3Histos[clusterSizeX-1]->fill(diff_y);
			}

			trueHitMap.at(detectorID).erase(trueHitMap.at(detectorID).begin()+pairIndex);
		}
	}
	catch (lcio::DataNotAvailableException& e) {

		return;
	}

	return;
}

void EUTelProcessorTrueHitAnalysis::bookHistos() {

	streamlog_out(DEBUG5) << "Booking histograms\n";

	for (size_t i = 0; i < _sensorIDVec.size(); i++) {

		int sensorID = _sensorIDVec[i];
		std::string basePath = "detector_" + to_string(sensorID);
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		//set up 1D histograms
		for (size_t i = 0; i < _1DHistos.size(); i++) {

			std::string coordinate = (i%2 == 0)?"x":"y";
			std::string clusterSize = to_string(i/2+1) + ((i+3>_1DHistos.size())?"+":"");

			std::string histoName = coordinate + "PositionDifference_" + coordinate + "ClusterSize" + clusterSize + "_d" + to_string(sensorID);

			int histoNBin = 161;
			double histoMin = -10;
			double histoMax = 10;
			std::string histoTitle = "difference in " + coordinate + " position between true simulated hits and reconstructed hits for " + coordinate + " cluster size " + clusterSize + ";#Delta" + coordinate + " /#mum;Entries";

			if (sensorID >= 10) {//sensor is a DUT, with much large pixel pitch
				if (i%2 == 0) {//adjust histogram for x-coordinate

					histoMin *= 20;
					histoMax *= 20;
				}
				else {//adjust histogram for y-coordinate

					histoMin *= 6;
					histoMax *= 6;
				}
			}

			_1DHistos[i].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histoNBin, histoMin, histoMax)));
			_1DHistos[i].at(sensorID)->setTitle(histoTitle.c_str());
		}

		//set up 2D histograms
		for (size_t i = 0; i < _2DHistos.size(); i++) {

			std::string coordinate = (i%2 == 0)?"x":"y";
			std::string clusterSize = to_string(i/2+1) + ((i+3>_2DHistos.size())?"+":"");

			std::string histoName = "2DPositionDifference_" + coordinate + "ClusterSize" + clusterSize + "_d" + to_string(sensorID);

			int histoXNBin = 161;
			double histoXMin = -10;
			double histoXMax = 10;

			int histoYNBin = 161;
			double histoYMin = -10;
			double histoYMax = 10;

			std::string histoTitle = "difference in position between true simulated hits and reconstructed hits for " + coordinate + " cluster size " + clusterSize + ";#Deltax /#mum;#Deltay /#mum;Entries";

			if (sensorID >= 10) {

				histoXMin *= 20;
				histoXMax *= 20;
				histoYMin *= 20;
				histoYMax *= 20;
			}

			_2DHistos[i].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histoXNBin, histoXMin, histoXMax, histoYNBin, histoYMin, histoYMax)));
			_2DHistos[i].at(sensorID)->setTitle(histoTitle.c_str());
		}

		//set up histograms for DUT with x cluster size 2
		if (sensorID >= 10) {

			std::string tempPath = basePath + "x_ClusterSize_2";
			AIDAProcessor::tree(this)->mkdir(tempPath.c_str());
			tempPath.append("/");

			for (size_t i = 0; i < 10; i++) {

				std::string modifier = (i%2==0)?"_yClustSize1":"";
				std::string TOTDiff = to_string(i/2) + ((i+3>10)?"+":"");

				std::string histoName = "xPositionDifference" + modifier + "_TOTDiff" + TOTDiff + "_d" + to_string(sensorID);

				int histoNBin = 161;
				double histoMin = -140;
				double histoMax = 140;
				std::string histoTitle = "difference in x position between true simulated hits and reconstructed hits for x cluster size 2 with TOT difference " + TOTDiff + ";#Deltax /#mum;Entries";

				_xClustSize2Histos.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D((tempPath+histoName).c_str(), histoNBin, histoMin, histoMax));
				_xClustSize2Histos[i]->setTitle(histoTitle.c_str());
			}
		}

		//set up histograms for Mimosa26 with y cluster size 3
		if (sensorID == 0) {

			std::string tempPath = basePath + "y_ClusterSize_3";
			AIDAProcessor::tree(this)->mkdir(tempPath.c_str());
			tempPath.append("/");

			for (size_t i = 0; i < 4; i++) {

				std::string xClustSize = to_string(i+1) + ((i+2>4)?"+":"");

				std::string histoName = "yPositionDifference_xClusterSize" + xClustSize + "_d" + to_string(sensorID);

				int histoNBin = 161;
				double histoMin = -20;
				double histoMax = 20;
				std::string histoTitle = "difference in y position between true simulated hits and reconstructed hits for y cluster size 3 with x cluster size " + xClustSize + ";#Deltay /#mum;Entries";

				_yClustSize3Histos.push_back(AIDAProcessor::histogramFactory(this)->createHistogram1D((tempPath+histoName).c_str(), histoNBin, histoMin, histoMax));
				_yClustSize3Histos[i]->setTitle(histoTitle.c_str());
			}
		}
	}

	streamlog_out(DEBUG5) << "end of Booking histograms\n";
}

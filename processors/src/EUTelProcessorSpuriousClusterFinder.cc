#include "EUTelProcessorSpuriousClusterFinder.h"

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

EUTelProcessorSpuriousClusterFinder::EUTelProcessorSpuriousClusterFinder() :
	Processor("EUTelProcessorSpuriousClusterFinder"),
	_trueHitCollectionName(""),
	_reconstructedHitCollectionName(""),
	_iRun(0),
	_iEvt(0),
	_sensorIDVec(),
	_trueHitCollectionVec(nullptr),
	_reconstructedHitCollectionVec(nullptr),
	_1DHistos(),
	_2DHistos()
{

	_description = "EUTelProcessorSpuriousClusterFinder finds and anlalyses spurious clusters (x cluster size 2, y cluster size 3, with a particular cluster shape)";

	registerInputCollection(LCIO::TRACKERHIT, "TrueHitCollectionName", "Input of True Hit data", _trueHitCollectionName, std::string("true_hits"));

	registerInputCollection(LCIO::TRACKERHIT, "ReconstructedHitCollectionName", "Input of Reconstructed Hit data", _reconstructedHitCollectionName, std::string("hit"));
}

void EUTelProcessorSpuriousClusterFinder::init() {

	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);
	_sensorIDVec = geo::gGeometry().sensorIDsVec();

	printParameters();

	_iRun = 0;
	_iEvt = 0;

	bookHistos();
}

void EUTelProcessorSpuriousClusterFinder::processRunHeader(LCRunHeader* rhdr) {

	std::unique_ptr<EUTelRunHeaderImpl> runHeader(new EUTelRunHeaderImpl(rhdr));
	runHeader->addProcessor(type());

	++_iRun;
}

void EUTelProcessorSpuriousClusterFinder::readCollections(LCEvent* event) {

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

void EUTelProcessorSpuriousClusterFinder::processEvent(LCEvent* event) {

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

void EUTelProcessorSpuriousClusterFinder::check (LCEvent*) {
}

void EUTelProcessorSpuriousClusterFinder::end() {

	streamlog_out(MESSAGE4) << "Successfully finished\n";
}

int EUTelProcessorSpuriousClusterFinder::findPairIndex(double const* a, std::vector<double const*> vec) {

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

void EUTelProcessorSpuriousClusterFinder::fillHistos(LCEvent* event) {

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
		std::map<int, std::vector<double>> trueHitEDepMap;

		for (size_t i = 0; i < _sensorIDVec.size(); i++) {

			trueHitMap.insert(std::make_pair(_sensorIDVec[i], std::vector<double const*>()));
			trueHitEDepMap.insert(std::make_pair(_sensorIDVec[i], std::vector<double>()));
		}

		//fill the map with the true hits
		for (int i = 0; i < _trueHitCollectionVec->getNumberOfElements(); i++) {

			TrackerHitImpl* trueHit = dynamic_cast<TrackerHitImpl*>(_trueHitCollectionVec->getElementAt(i));
			int detectorID = static_cast<int>(trueHitDecoder(trueHit)["sensorID"]);

			trueHitMap.at(detectorID).push_back(trueHit->getPosition());
			trueHitEDepMap.at(detectorID).push_back(trueHit->getEDepError());
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
				if (std::find(XHitPixels.begin(), XHitPixels.end(), chargeValues[j]) == XHitPixels.end()) XHitPixels.push_back(chargeValues[j]);
				if (std::find(YHitPixels.begin(), YHitPixels.end(), chargeValues[j+1]) == YHitPixels.end()) YHitPixels.push_back(chargeValues[j+1]);
			}
			int clusterSizeX = static_cast<int>(XHitPixels.size());
			int clusterSizeY = static_cast<int>(YHitPixels.size());
			if (clusterSizeX > 4) clusterSizeX = 4;
			if (clusterSizeY > 4) clusterSizeY = 4;

			double const* hitPos = reconstructedHit->getPosition();
			int pairIndex = findPairIndex(hitPos, static_cast<std::vector<double const*>>(trueHitMap.at(detectorID)));
			double const* pair = (static_cast<std::vector<double const*>>(trueHitMap.at(detectorID)))[pairIndex];
			double edepTotal = (static_cast<std::vector<double>>(trueHitEDepMap.at(detectorID)))[pairIndex];

			double diff_x = (hitPos[0] - pair[0])*1000;
			double diff_y = (hitPos[1] - pair[1])*1000;

			if ((clusterSizeX == 2) && (clusterSizeY == 3)) {//potentially a spurious cluster

				std::vector<int> shape(2, 0);
				for (size_t j = pixelEntries; j < chargeValues.size(); j+= pixelEntries) {

					shape.push_back(static_cast<int>(chargeValues[j]-chargeValues[0]));
					shape.push_back(static_cast<int>(chargeValues[j+1]-chargeValues[1]));
				}

				bool isSpuriousCluster = false;
				if ((static_cast<int>(shape.size()) == 8) && (shape[2] == -1) && (shape[3] == 1) && (shape[4] == 0) && (shape[5] == 1) && (shape[6] == -1) && (shape[7] && 2)) isSpuriousCluster = true;

				if (isSpuriousCluster) {

					//fill histograms for spurious clusters
					_1DHistos[0].at(detectorID)->fill(diff_x);
					_1DHistos[1].at(detectorID)->fill(diff_y);
					_1DHistos[2].at(detectorID)->fill(hitPos[0]);
					_1DHistos[3].at(detectorID)->fill(hitPos[1]);
					_1DHistos[4].at(detectorID)->fill(pair[0]);
					_1DHistos[5].at(detectorID)->fill(pair[1]);
					_1DHistos[6].at(detectorID)->fill(edepTotal);
					_1DHistos[13].at(detectorID)->fill(diff_x);
					_1DHistos[14].at(detectorID)->fill(diff_y+9.2);
					_1DHistos[17].at(detectorID)->fill(diff_x);
					_1DHistos[18].at(detectorID)->fill(diff_y+9.2);
					_2DHistos[0].at(detectorID)->fill(diff_x, diff_y, 1.0);
					_2DHistos[1].at(detectorID)->fill(hitPos[0], hitPos[1], 1.0);
					_2DHistos[2].at(detectorID)->fill(pair[0], pair[1], 1.0);
					_2DHistos[9].at(detectorID)->fill(diff_x, diff_y+9.2, 1.0);
					_2DHistos[12].at(detectorID)->fill(diff_x, diff_y+9.2, 1.0);
					for (size_t j = 0; j < shape.size(); j+=2) {

						_2DHistos[3].at(detectorID)->fill(shape[j], shape[j+1], 1.0);
					}
				}

				if (isSpuriousCluster == false) {

					//fill histograms for non spurious clusters
					_1DHistos[7].at(detectorID)->fill(diff_x);
					_1DHistos[8].at(detectorID)->fill(diff_y);
					_1DHistos[9].at(detectorID)->fill(hitPos[0]);
					_1DHistos[10].at(detectorID)->fill(hitPos[1]);
					_1DHistos[11].at(detectorID)->fill(pair[0]);
					_1DHistos[12].at(detectorID)->fill(pair[1]);
					_2DHistos[5].at(detectorID)->fill(diff_x, diff_y, 1.0);
					_2DHistos[6].at(detectorID)->fill(hitPos[0], hitPos[1], 1.0);
					_2DHistos[7].at(detectorID)->fill(pair[0], pair[1], 1.0);
					for (size_t j = 0; j < shape.size(); j+=2) {

						_2DHistos[8].at(detectorID)->fill(shape[j], shape[j+1], 1.0);
					}
				}

				for (size_t j = 0; j < shape.size(); j+=2) {

					_2DHistos[4].at(detectorID)->fill(shape[j], shape[j+1], 1.0);
				}
			}

			if ((clusterSizeX == 2) && (clusterSizeY == 2)) {

				std::vector<int> shape(2, 0);
				for (size_t j = pixelEntries; j < chargeValues.size(); j+= pixelEntries) {

					shape.push_back(static_cast<int>(chargeValues[j]-chargeValues[0]));
					shape.push_back(static_cast<int>(chargeValues[j+1]-chargeValues[1]));
				}

				bool is2x2Cluster = false;
				if (static_cast<int>(shape.size()) == 8) is2x2Cluster = true;

				if (is2x2Cluster) {

					_1DHistos[15].at(detectorID)->fill(diff_x);
					_1DHistos[16].at(detectorID)->fill(diff_y);
					_1DHistos[17].at(detectorID)->fill(diff_x);
					_1DHistos[18].at(detectorID)->fill(diff_y);
					_2DHistos[10].at(detectorID)->fill(diff_x, diff_y, 1.0);
					_2DHistos[12].at(detectorID)->fill(diff_x, diff_y, 1.0);
					for (size_t j = 0; j < shape.size(); j+=2) {

						_2DHistos[11].at(detectorID)->fill(shape[j], shape[j+1], 1.0);
					}
				}
				else {

					_1DHistos[19].at(detectorID)->fill(diff_x);
					_1DHistos[20].at(detectorID)->fill(diff_y);
					_2DHistos[13].at(detectorID)->fill(diff_x, diff_y, 1.0);
					for (size_t j = 0; j < shape.size(); j+=2) {

						_2DHistos[14].at(detectorID)->fill(shape[j], shape[j+1], 1.0);
					}
				}
			}

			trueHitMap.at(detectorID).erase(trueHitMap.at(detectorID).begin()+pairIndex);
		}
	}
	catch (lcio::DataNotAvailableException& e) {

		return;
	}

	return;
}

void EUTelProcessorSpuriousClusterFinder::bookHistos() {

	streamlog_out(DEBUG5) << "Booking histograms\n";

	for (size_t i = 0; i < _sensorIDVec.size(); i++) {

		int sensorID = _sensorIDVec[i];
		std::string basePath = "detector_" + to_string(sensorID);
		AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		//sets the number of bins for all 1D histos
		int histo1DNBin = 161;

		//sets the histogram parameters for all 1D x position difference histos
		double xPosDiffMin = -12;
		double xPosDiffMax = 12;

		//sets the histogram parameters for all 1D y position difference histos
		double yPosDiffMin = -12;
		double yPosDiffMax = 12;

		//sets the histogram parameters for all 1D hit x position histograms
		double xHitPosMin = -11;
		double xHitPosMax = 11;

		//sets the histogram parameters for all 1D hit y position histograms
		double yHitPosMin = -6;
		double yHitPosMax = 6;

		//sets the number of bins for all 2D non-shape histograms
		double histo2DXNBin = 161;
		double histo2DYNBin = 161;

		//sets the histogram parameters for all 2D position difference histograms
		double xPosDiff2DMin = -12;
		double xPosDiff2DMax = 12;
		double yPosDiff2DMin = -12;
		double yPosDiff2DMax = 12;

		//sets the histogram parameters for all 2D hit position histograms
		double xHitPos2DMin = -11;
		double xHitPos2DMax = 11;
		double yHitPos2DMin = -11;
		double yHitPos2DMax = 11;

		//sets the histogram parameters for all 2D shape histograms
		double shapeXNBin = 6;
		double shapeXMin = -2;
		double shapeXMax = 4;

		double shapeYNBin = 6;
		double shapeYMin = -2;
		double shapeYMax = 4;

		//DUTs have larger pixel pitch, so histogram limits must be adjusted
		//to account for larger spread in position difference
		if (sensorID >= 10) {

			xPosDiffMin = -200;
			xPosDiffMax = 200;

			yPosDiffMin = -80;
			yPosDiffMax = 80;

			xPosDiff2DMin = -200;
			xPosDiff2DMax = 200;
			yPosDiff2DMin = -200;
			yPosDiff2DMax = 200;
		}

		std::string histoName, histoTitle;

		//set up 1D histograms for spurious clusters
		//set up x position difference histogram
		histoName = "xPositionDifference_d" + to_string(sensorID);
		histoTitle = "difference in x position between true simulated hits and reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltax /#mum;Entries";
		_1DHistos[0].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xPosDiffMin, xPosDiffMax)));
		_1DHistos[0].at(sensorID)->setTitle(histoTitle.c_str());

		//set up y position difference histogram
		histoName = "yPositionDifference_d" + to_string(sensorID);
		histoTitle = "difference in y position between true simulated hits and reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltay /#mum;Entries";
		_1DHistos[1].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yPosDiffMin, yPosDiffMax)));
		_1DHistos[1].at(sensorID)->setTitle(histoTitle.c_str());

		//set up x position difference histogram
		histoName = "xPositionDifference_shifted_d" + to_string(sensorID);
		histoTitle = "difference in x position between true simulated hits and reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3) with y difference offset by 9.2 #mum;#Deltax /#mum;Entries";
		_1DHistos[13].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xPosDiffMin, xPosDiffMax)));
		_1DHistos[13].at(sensorID)->setTitle(histoTitle.c_str());

		//set up y position difference histogram
		histoName = "yPositionDifference_shifted_d" + to_string(sensorID);
		histoTitle = "difference in y position between true simulated hits and reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3) with y difference offset by 9.2 #mum;#Deltay /#mum;Entries";
		_1DHistos[14].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yPosDiffMin, yPosDiffMax)));
		_1DHistos[14].at(sensorID)->setTitle(histoTitle.c_str());

		//set up reconstructed hit x position histogram
		histoName = "xReconstructedHitPosition_d" + to_string(sensorID);
		histoTitle = "x position of reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltay /#mum;Entries";
		_1DHistos[2].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xHitPosMin, xHitPosMax)));
		_1DHistos[2].at(sensorID)->setTitle(histoTitle.c_str());

		//set up reconstructed hit y position histogram
		histoName = "yReconstrutedHitPosition_d" + to_string(sensorID);
		histoTitle = "y position of reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltay /#mum;Entries";
		_1DHistos[3].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yHitPosMin, yHitPosMax)));
		_1DHistos[3].at(sensorID)->setTitle(histoTitle.c_str());

		//set up true hit x position histogram
		histoName = "xTrueHitPosition_d" + to_string(sensorID);
		histoTitle = "x position of true simulated hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltay /#mum;Entries";
		_1DHistos[4].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xHitPosMin, xHitPosMax)));
		_1DHistos[4].at(sensorID)->setTitle(histoTitle.c_str());

		//set up true hit y position histogram
		histoName = "yTrueHitPosition_d" + to_string(sensorID);
		histoTitle = "y position for true simulated hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltay /#mum;Entries";
		_1DHistos[5].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yHitPosMin, yHitPosMax)));
		_1DHistos[5].at(sensorID)->setTitle(histoTitle.c_str());

		//set up true hit total energy deposition histo
		histoName = "TrueHitEDepTotal_d" + to_string(sensorID);
		histoTitle = "total deposition energy for true simulated hits for spurious clusters (x cluster size 2, y cluster size 3);Total deposition energy;Entries";
		_1DHistos[6].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), 160, 0, 80)));
		_1DHistos[6].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 1D histograms for non-spurious clusters
		//set up x position difference histogram
		histoName = "xPositionDifference_NonSpurious_d" + to_string(sensorID);
		histoTitle = "difference in x position between true simulated hits and reconstructed hits for non spurious clusters;#Deltax /#mum;Entries";
		_1DHistos[7].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xPosDiffMin, xPosDiffMax)));
		_1DHistos[7].at(sensorID)->setTitle(histoTitle.c_str());

		//set up y position difference histogram
		histoName = "yPositionDifference_NonSpurious_d" + to_string(sensorID);
		histoTitle = "difference in y position between true simulated hits and reconstructed hits for non spurious clusters;#Deltay /#mum;Entries";
		_1DHistos[8].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yPosDiffMin, yPosDiffMax)));
		_1DHistos[8].at(sensorID)->setTitle(histoTitle.c_str());

		//set up reconstructed hit x position histogram
		histoName = "xReconstructedHitPosition_NonSpurious_d" + to_string(sensorID);
		histoTitle = "x position of reconstructed hits for non spurious clusters;#Deltay /#mum;Entries";
		_1DHistos[9].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xHitPosMin, xHitPosMax)));
		_1DHistos[9].at(sensorID)->setTitle(histoTitle.c_str());

		//set up reconstructed hit y position histogram
		histoName = "yReconstrutedHitPosition_NonSpurious_d" + to_string(sensorID);
		histoTitle = "y position of reconstructed hits for non spurious clusters;#Deltay /#mum;Entries";
		_1DHistos[10].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yHitPosMin, yHitPosMax)));
		_1DHistos[10].at(sensorID)->setTitle(histoTitle.c_str());

		//set up true hit x position histogram
		histoName = "xTrueHitPosition_NonSpurious_d" + to_string(sensorID);
		histoTitle = "x position of true simulated hits for non spurious clusters;#Deltay /#mum;Entries";
		_1DHistos[11].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xHitPosMin, xHitPosMax)));
		_1DHistos[11].at(sensorID)->setTitle(histoTitle.c_str());

		//set up true hit y position histogram
		histoName = "yTrueHitPosition_NonSpurious_d" + to_string(sensorID);
		histoTitle = "y position for true simulated hits for non spurious clusters;#Deltay /#mum;Entries";
		_1DHistos[12].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yHitPosMin, yHitPosMax)));
		_1DHistos[12].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 1D histograms for x cluster size 2 y cluster size 2 non 2x2 clusters
		//set up x position difference histogram
		histoName = "xPositionDifference_Non2x2_d" + to_string(sensorID);
		histoTitle = "difference in x position between true simulated hits and reconstructed hits for x cluster size 2 y cluster size 2 non 2x2;#Deltax /#mum;Entries";
		_1DHistos[19].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xPosDiffMin, xPosDiffMax)));
		_1DHistos[19].at(sensorID)->setTitle(histoTitle.c_str());

		//set up y position difference histogram
		histoName = "yPositionDifference_Non2x2_d" + to_string(sensorID);
		histoTitle = "difference in y position between true simulated hits and reconstructed hits for x cluster size 2 y cluster size 2 non 2x2;#Deltay /#mum;Entries";
		_1DHistos[20].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yPosDiffMin, yPosDiffMax)));
		_1DHistos[20].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 1D histograms for 2x2 clusters
		//set up x position difference histogram
		histoName = "xPositionDifference_2x2Cluster_d" + to_string(sensorID);
		histoTitle = "difference in x position between true simulated hits and reconstructed hits for 2x2 clusters;#Deltax /#mum;Entries";
		_1DHistos[15].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xPosDiffMin, xPosDiffMax)));
		_1DHistos[15].at(sensorID)->setTitle(histoTitle.c_str());

		//set up y position difference histogram
		histoName = "yPositionDifference_2x2Cluster_d" + to_string(sensorID);
		histoTitle = "difference in y position between true simulated hits and reconstructed hits for 2x2 clusters;#Deltay /#mum;Entries";
		_1DHistos[16].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yPosDiffMin, yPosDiffMax)));
		_1DHistos[16].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 1D histograms for 2x2 and spurious clusters
		//set up x position difference histogram
		histoName = "xPositionDifference_2x2+spurious_d" + to_string(sensorID);
		histoTitle = "difference in x position between true simulated hits and reconstructed hits for 2x2 and spurious clusters;#Deltax /#mum;Entries";
		_1DHistos[17].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, xPosDiffMin, xPosDiffMax)));
		_1DHistos[17].at(sensorID)->setTitle(histoTitle.c_str());

		//set up y position difference histogram
		histoName = "yPositionDifference_2x2+spurious_d" + to_string(sensorID);
		histoTitle = "difference in y position between true simulated hits and reconstructed hits for 2x2 and spurious clusters;#Deltay /#mum;Entries";
		_1DHistos[18].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram1D((basePath+histoName).c_str(), histo1DNBin, yPosDiffMin, yPosDiffMax)));
		_1DHistos[18].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 2D histograms for spurious clusters
		//set up position difference histogram
		histoName = "2DPositionDifference_d" + to_string(sensorID);
		histoTitle = "difference in position between true simulated hits and reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[0].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xPosDiff2DMin, xPosDiff2DMax, histo2DYNBin, yPosDiff2DMin, yPosDiff2DMax)));
		_2DHistos[0].at(sensorID)->setTitle(histoTitle.c_str());

		//set up position difference shifted histogram
		histoName = "2DPositionDifference_shifted_d" + to_string(sensorID);
		histoTitle = "difference in position between true simulated hits and reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3) with y difference offset by 9.2 #mum;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[9].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xPosDiff2DMin, xPosDiff2DMax, histo2DYNBin, yPosDiff2DMin, yPosDiff2DMax)));
		_2DHistos[9].at(sensorID)->setTitle(histoTitle.c_str());

		//set up reconstructed hit position histogram
		histoName = "2DReconstructedHitPosition_d" + to_string(sensorID);
		histoTitle = "position of reconstructed hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[1].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xHitPos2DMin, xHitPos2DMax, histo2DYNBin, yHitPos2DMin, yHitPos2DMax)));
		_2DHistos[1].at(sensorID)->setTitle(histoTitle.c_str());

		//set up true hit position histogram
		histoName = "2DTrueHitPosition_d" + to_string(sensorID);
		histoTitle = "position of true simulated hits for spurious clusters (x cluster size 2, y cluster size 3);#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[2].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xHitPos2DMin, xHitPos2DMax, histo2DYNBin, yHitPos2DMin, yHitPos2DMax)));
		_2DHistos[2].at(sensorID)->setTitle(histoTitle.c_str());

		//set up shape histogram for spurious clusters
		histoName = "2DClusterShape_SpuriousClusters_d" + to_string(sensorID);
		histoTitle = "shape of clusters for spurious clusters (x cluster size 2, y cluster size 3);#Deltax /#mum;#Deltay /#mum;Entries";

		_2DHistos[3].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), shapeXNBin, shapeXMin, shapeXMax, shapeYNBin, shapeYMin, shapeYMax)));
		_2DHistos[3].at(sensorID)->setTitle(histoTitle.c_str());

		//set up shape histogram for all clusters with x cluster size 2 and y cluster size 3
		histoName = "2DClusterShape_d" + to_string(sensorID);
		histoTitle = "shape of all clusters of x cluster size 2 and y cluster size 3;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[4].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), shapeXNBin, shapeXMin, shapeXMax, shapeYNBin, shapeYMin, shapeYMax)));
		_2DHistos[4].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 2D histograms for non spurious clusters
		//set up position difference histogram
		histoName = "2DPositionDifference_NonSpurious_d" + to_string(sensorID);
		histoTitle = "difference in position between true simulated hits and reconstructed hits for non spurious clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[5].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xPosDiff2DMin, xPosDiff2DMax, histo2DYNBin, yPosDiff2DMin, yPosDiff2DMax)));
		_2DHistos[5].at(sensorID)->setTitle(histoTitle.c_str());

		//set up reconstructed hit position histogram
		histoName = "2DReconstructedHitPosition_NonSpurious_d" + to_string(sensorID);
		histoTitle = "position of reconstructed hits for non spurious clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[6].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xHitPos2DMin, xHitPos2DMax, histo2DYNBin, yHitPos2DMin, yHitPos2DMax)));
		_2DHistos[6].at(sensorID)->setTitle(histoTitle.c_str());

		//set up true hit position histogram
		histoName = "2DTrueHitPosition_NonSpurious_d" + to_string(sensorID);
		histoTitle = "position of true simulated hits for non spurious clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[7].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xHitPos2DMin, xHitPos2DMax, histo2DYNBin, yHitPos2DMin, yHitPos2DMax)));
		_2DHistos[7].at(sensorID)->setTitle(histoTitle.c_str());

		//set up shape histogram for non spurious clusters
		histoName = "2DClusterShape_NonSpurious_d" + to_string(sensorID);
		histoTitle = "shape of clusters for non spurious clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[8].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), shapeXNBin, shapeXMin, shapeXMax, shapeYNBin, shapeYMin, shapeYMax)));
		_2DHistos[8].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 2D histograms for x cluster size 2 y cluster size 2 non 2x2 clusters
		//set up position difference histogram
		histoName = "2DPositionDifference_Non2x2_d" + to_string(sensorID);
		histoTitle = "difference in position between true simulated hits and reconstructed hits for x cluster size 2 y cluster size 2 non 2x2;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[13].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xPosDiff2DMin, xPosDiff2DMax, histo2DYNBin, yPosDiff2DMin, yPosDiff2DMax)));
		_2DHistos[13].at(sensorID)->setTitle(histoTitle.c_str());

		//set up shape histogram for non 2x2 clusters
		histoName = "2DClusterShape_Non2x2_d" + to_string(sensorID);
		histoTitle = "shape of clusters for x cluster size 2 y cluster size 2 non 2x2 clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[14].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), shapeXNBin, shapeXMin, shapeXMax, shapeYNBin, shapeYMin, shapeYMax)));
		_2DHistos[14].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 2D histograms for 2x2 clusters
		//set up position difference histogram
		histoName = "2DPositionDifference_2x2Clusters_d" + to_string(sensorID);
		histoTitle = "difference in position between true simulated hits and reconstructed hits for 2x2 clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[10].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xPosDiff2DMin, xPosDiff2DMax, histo2DYNBin, yPosDiff2DMin, yPosDiff2DMax)));
		_2DHistos[10].at(sensorID)->setTitle(histoTitle.c_str());

		//set up shape histogram for 2x2 clusters
		histoName = "2DClusterShape_2x2Clusters_d" + to_string(sensorID);
		histoTitle = "shape of clusters for 2x2 clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[11].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), shapeXNBin, shapeXMin, shapeXMax, shapeYNBin, shapeYMin, shapeYMax)));
		_2DHistos[11].at(sensorID)->setTitle(histoTitle.c_str());

		//set up 2D histograms for 2x2 and spurious clusters
		//set up position difference histogram
		histoName = "2DPositionDifference_2x2+spurious_d" + to_string(sensorID);
		histoTitle = "difference in position between true simulated hits and reconstructed hits for 2x2 and spurious clusters;#Deltax /#mum;#Deltay /#mum;Entries";
		_2DHistos[12].insert(std::make_pair(sensorID, AIDAProcessor::histogramFactory(this)->createHistogram2D((basePath+histoName).c_str(), histo2DXNBin, xPosDiff2DMin, xPosDiff2DMax, histo2DYNBin, yPosDiff2DMin, yPosDiff2DMax)));
		_2DHistos[12].at(sensorID)->setTitle(histoTitle.c_str());
	}

	streamlog_out(DEBUG5) << "end of Booking histograms\n";
}

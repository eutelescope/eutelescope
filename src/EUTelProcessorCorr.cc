/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelProcessorCorr.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "CellIDReencoder.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/AIDAProcessor.h"

//Standard C++ libraries 
#include <vector>

// lcio includes <.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>

using namespace eutelescope;

struct hPnt {
	hPnt(double x, double y): x(x), y(y){};
	double x, y;
};

EUTelProcessorCorr::EUTelProcessorCorr():
Processor("EUTelProcessorCorr"),
_hitCollectionNameInput(),
_sensorIDVec() 
{
		_description ="EUTelProcessorCorr fills correlation histograms. It is possible to fill correlations for all different permutations of axis flips and changes.";
		registerInputCollection(LCIO::TRACKERHIT, "hitCollectionNameInput", "Input hit collection name", _hitCollectionNameInput, std::string("inputhit"));
		registerOptionalParameter("distCut", "Set distance cut to discard correlations further away than this", _distCut,  static_cast<float>(0)); //Not used right now
}

void EUTelProcessorCorr::init() {
	geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME, EUTELESCOPE::DUMPGEOROOT);
	_sensorIDVec = geo::gGeometry().sensorIDsVec();
	bookHistos();	
}

void EUTelProcessorCorr::processEvent(LCEvent* event) {

		EUTelEventImpl* evt	= static_cast<EUTelEventImpl*>(event);				
		if( evt->getEventType() == kEORE ) {
				streamlog_out( MESSAGE5 ) << "EORE found: nothing else to do." << std::endl;
				return;
		} else if( evt->getEventType() == kUNKNOWN ) {
				streamlog_out( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
						<< " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}

		//Opens collection for input.
		LCCollection* inputCollection = nullptr;
		try {
				inputCollection = evt->getCollection(_hitCollectionNameInput);
		} catch (DataNotAvailableException& e) {
				streamlog_out( WARNING2 ) << _hitCollectionNameInput << " collection not available" << std::endl;
				return;
		}

		std::string encoding = inputCollection->getParameters().getStringVal( LCIO::CellIDEncoding );
		if(encoding.empty()) {
			encoding = EUTELESCOPE::HITENCODING;
		}
		lcio::CellIDDecoder<TrackerHitImpl> hitDecoder ( encoding );

		std::vector<hPnt> firstPlanePoints;
		std::map<int, std::array<std::vector<hPnt>, perm::last>> otherPoints;
		
		//We loop over all hits in the event and store them via the hPnt helper point struct in containers,
		//only once this is done we process them
		for(int iHit = 0; iHit < inputCollection->getNumberOfElements(); ++iHit) {  
			TrackerHitImpl*	inputHit = static_cast<TrackerHitImpl*>(inputCollection->getElementAt(iHit));

			int sensorID = hitDecoder(inputHit)["sensorID"];
			const double* inputPos = inputHit->getPosition();
	
			if(sensorID == _sensorIDVec.front()) {
				firstPlanePoints.emplace_back(inputPos[0], inputPos[1]);
			} else {
				otherPoints[sensorID].at(perm::x_y).emplace_back(inputPos[0], inputPos[1]);
				otherPoints[sensorID].at(perm::mx_y).emplace_back(-inputPos[0], inputPos[1]);
				otherPoints[sensorID].at(perm::x_my).emplace_back(inputPos[0], -inputPos[1]);
				otherPoints[sensorID].at(perm::mx_my).emplace_back(-inputPos[0], -inputPos[1]);
				otherPoints[sensorID].at(perm::y_x).emplace_back(inputPos[1], inputPos[0]);
				otherPoints[sensorID].at(perm::my_x).emplace_back(-inputPos[1], inputPos[0]);
				otherPoints[sensorID].at(perm::y_mx).emplace_back(inputPos[1], -inputPos[0]);
				otherPoints[sensorID].at(perm::my_mx).emplace_back(-inputPos[1], -inputPos[0]);
			}
		}

		auto calcDistAndFillHisto = [&](perm p){
			for (auto& val: otherPoints) {
				auto sensorID = val.first;
				auto perms = val.second;
				
				double dist = 9999999999;

				hPnt* p1 = nullptr;
				hPnt* p2 = nullptr;

				for(auto& xPlanePnt: perms.at(p)){
					for(auto& fstPlanePnt: firstPlanePoints) {
						double tDist = pow(xPlanePnt.x-fstPlanePnt.x, 2)+pow(xPlanePnt.y-fstPlanePnt.y, 2);
						if( tDist < dist) {
							p1 = &xPlanePnt;
							p2 = &fstPlanePnt;
							dist = tDist;
						}
					}
					//here min pnt
					if(p1 && p2){
						auto p1d = *p1;
						auto p2d = *p2;
						_histoMap[sensorID].at(p).xx->fill(p1d.x, p2d.x, 1);
						_histoMap[sensorID].at(p).xy->fill(p1d.x, p2d.y, 1);
						_histoMap[sensorID].at(p).yx->fill(p1d.y, p2d.x, 1);
						_histoMap[sensorID].at(p).yy->fill(p1d.y, p2d.y, 1);
					}
				}
			}	
		};
		
		calcDistAndFillHisto(perm::x_y);
		calcDistAndFillHisto(perm::mx_y);
		calcDistAndFillHisto(perm::x_my);
		calcDistAndFillHisto(perm::mx_my);
		calcDistAndFillHisto(perm::y_x);
		calcDistAndFillHisto(perm::my_x);
		calcDistAndFillHisto(perm::y_mx);
		calcDistAndFillHisto(perm::my_mx);
}

void EUTelProcessorCorr::end() {
	streamlog_out(MESSAGE4) << "Successfully finished" << std::endl;
}

void EUTelProcessorCorr::bookHistos() {
	// histograms are grouped in loops and detectors
	streamlog_out ( DEBUG5 )  << "Booking histograms " << std::endl;

	int firstSensorID = _sensorIDVec.front();
	float firstSensorSizeX, firstSensorSizeY, firstSensorSizeZ;
	geo::EUTelGenericPixGeoDescr* firstSensorGeoDescr = geo::gGeometry().getPixGeoDescr( firstSensorID );
	firstSensorGeoDescr->getSensitiveSize(firstSensorSizeX, firstSensorSizeY, firstSensorSizeZ);

	int firstBinX = 2*ceil(firstSensorSizeX);
	int firstBinY = 2*ceil(firstSensorSizeY);

	float firstBoundX = firstSensorSizeX/2+2;
	float firstBoundY = firstSensorSizeY/2+2;

	std::map<perm, std::string> permTitleMap;

	permTitleMap[perm::x_y] = "plusXplusY";
	permTitleMap[perm::mx_y] = "minusXplusY";
	permTitleMap[perm::x_my] = "plusXminusY";
	permTitleMap[perm::mx_my] = "minusXminusY";
	permTitleMap[perm::y_x] = "plusYplusX";
	permTitleMap[perm::my_x] = "minusYplusX";
	permTitleMap[perm::y_mx] = "plusYminusX";
	permTitleMap[perm::my_mx] = "minusYminusX";

	for (auto mIt = _sensorIDVec.begin()+1; mIt != _sensorIDVec.end(); ++mIt) {
		int sensorID = *mIt;
		geo::EUTelGenericPixGeoDescr* thisSensorGeoDescr = geo::gGeometry().getPixGeoDescr( sensorID );

		//Sizes and positions are in mm!
		float sizeX, sizeY, sizeZ;
		thisSensorGeoDescr->getSensitiveSize(sizeX, sizeY, sizeZ);

		//Trying to figure out a smart way of setting the histo bins
		int thisBinX = 10*ceil(sizeX);
		int thisBinY = 10*ceil(sizeY);

		float thisBoundX = sizeX/2+2;
		float thisBoundY = sizeY/2+2;

		std::string basePath = "detector_" + to_string( sensorID );
		marlin::AIDAProcessor::tree(this)->mkdir(basePath.c_str());
		basePath.append("/");

		//In the 'normal' case the XY axis are not flipped, this is relevant for the histo boundaries
		auto subCreateHistoNorm = [&](perm p){
			std::string newBasePath = basePath;
			newBasePath.append(permTitleMap.at(p));
			marlin::AIDAProcessor::tree(this)->mkdir(newBasePath.c_str());	
			newBasePath.append("/");
	
			_histoMap[sensorID].at(p).xx = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_XXcorr").c_str(), thisBinX, -thisBoundX, thisBoundX, firstBinX, -firstBoundX, firstBoundX 
			);
			 _histoMap[sensorID].at(p).xy = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_XYcorr").c_str(), thisBinX, -thisBoundX, thisBoundX, firstBinY, -firstBoundY, firstBoundY 
			);
			 _histoMap[sensorID].at(p).yx = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_YXcorr").c_str(), thisBinY, -thisBoundY, thisBoundY, firstBinX, -firstBoundX, firstBoundX 
			);
			 _histoMap[sensorID].at(p).yy = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_YYcorr").c_str(), thisBinY, -thisBoundY, thisBoundY, firstBinY, -firstBoundY, firstBoundY 
			);
		};

		//The permutations with flipped axis need to set the boundaries differently
		auto subCreateHistoInv = [&](perm p){
			std::string newBasePath = basePath;
			newBasePath.append(permTitleMap.at(p));
			marlin::AIDAProcessor::tree(this)->mkdir(newBasePath.c_str());	
			newBasePath.append("/");
	
			_histoMap[sensorID].at(p).xx = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_XXcorr").c_str(), thisBinY, -thisBoundY, thisBoundY, firstBinX, -firstBoundX, firstBoundX 
			);
			 _histoMap[sensorID].at(p).xy = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_XYcorr").c_str(), thisBinY, -thisBoundY, thisBoundY, firstBinY, -firstBoundY, firstBoundY 
			);
			 _histoMap[sensorID].at(p).yx = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_YXcorr").c_str(), thisBinX, -thisBoundX, thisBoundX, firstBinX, -firstBoundX, firstBoundX 
			);
			 _histoMap[sensorID].at(p).yy = marlin::AIDAProcessor::histogramFactory(this)->createHistogram2D( 
				(newBasePath + "_YYcorr").c_str(), thisBinX, -thisBoundX, thisBoundX, firstBinY, -firstBoundY, firstBoundY 
			);
		};

		//Call the lambda corresponding to if the axis have been flipped or not
		subCreateHistoNorm(perm::x_y);
		subCreateHistoNorm(perm::mx_y);
		subCreateHistoNorm(perm::x_my);
		subCreateHistoNorm(perm::mx_my);

		subCreateHistoInv(perm::y_x);
		subCreateHistoInv(perm::my_x);
		subCreateHistoInv(perm::y_mx);
		subCreateHistoInv(perm::my_mx);
	}
}

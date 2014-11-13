#include "EUTelProcessorTrackAnalysis.h"

using namespace eutelescope;

EUTelProcessorTrackAnalysis::EUTelProcessorTrackAnalysis() :
Processor("EUTelProcessorTrackAnalysis"){
	
	registerInputCollection(LCIO::TRACK, "TrackInputCollectionName", "Input track collection name",_trackInputCollectionName,std::string("TrackCandidatesCollection"));
	registerOptionalParameter("SensorIDs", "A vector of the sensor IDs to plot",_sensorIDs,IntVec());

	
	
}


void EUTelProcessorTrackAnalysis::init(){
	initialiseResidualVsPositionHistograms();
	EUTelTrackAnalysis*	analysis = new EUTelTrackAnalysis(_mapFromSensorIDToHistogramX,_mapFromSensorIDToHistogramY,_mapFromSensorIDToKinkXZ,_mapFromSensorIDToKinkYZ) ;
	_analysis = analysis;
}

void EUTelProcessorTrackAnalysis::processRunHeader(LCRunHeader * run) {}

void EUTelProcessorTrackAnalysis::check(LCEvent * evt){}

void EUTelProcessorTrackAnalysis::processEvent(LCEvent * evt){

	EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

	if (event->getEventType() == kEORE) {
		streamlog_out(DEBUG4) << "EORE found: nothing else to do." << endl;
		return;
	}else if (event->getEventType() == kUNKNOWN) {
		streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << endl;
	}
	LCCollection* eventCollection = NULL;
	try {
		eventCollection = evt->getCollection(_trackInputCollectionName);
		streamlog_out(DEBUG1) << "collection : " << _trackInputCollectionName << " retrieved" << std::endl;
	}catch (DataNotAvailableException e) {
		streamlog_out(MESSAGE0) << _trackInputCollectionName << " collection not available" << std::endl;
		throw marlin::SkipEventException(this);
	}
	if (eventCollection != NULL) {
		streamlog_out(DEBUG2) << "Collection contains data! Continue!" << endl;
		for (int iTrack = 0; iTrack < eventCollection->getNumberOfElements(); ++iTrack){
			EUTelTrack track = *(static_cast<EUTelTrack*> (eventCollection->getElementAt(iTrack)));
			_analysis->plotResidualVsPosition(track);	
			_analysis->plotIncidenceAngles(track);
		}

	}	
	
	
	
	
}

void EUTelProcessorTrackAnalysis::end(){}

void	EUTelProcessorTrackAnalysis::initialiseResidualVsPositionHistograms(){
	int NBinX;
	double MinX;
	double MaxX;
	int NBinY;
	double MinY;
	double MaxY;
	double MinZ;
	double MaxZ;

	std::string _histoInfoFileName="histoInfo.xml";
	auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
	bool                    isHistoManagerAvailable;

	try {
			isHistoManagerAvailable = histoMgr->init( );
	} catch ( ios::failure& e ) {
			streamlog_out( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
							<< "Continuing without histogram manager using default settings"    << endl;
			isHistoManagerAvailable = false;
	} catch ( ParseException& e ) {
			streamlog_out( ERROR5 ) << e.what( ) << "\n"
							<< "Continuing without histogram manager using default settings" << endl;
			isHistoManagerAvailable = false;
	}

	std::stringstream sstm;
	std::string residGblFitHistName;
	std::string histTitle;
	for (int i = 0; i < _sensorIDs.size() ; ++i){
		/////////////////////////////////////////////////////////////////////////////XY residual plots with position
		sstm << "ResidualX" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "ResidualsX. Plane " <<  _sensorIDs.at(i) << ";X direction; Y direction";
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 40;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 20;
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -5;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 5;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : -20;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 20;
		AIDA::IProfile2D *  residGblFitX =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (residGblFitX) {
				residGblFitX->setTitle(histTitle);
				_mapFromSensorIDToHistogramX.insert(std::make_pair(_sensorIDs.at(i), residGblFitX));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	for (int i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "ResidualY" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "ResidualsY. Plane " <<  _sensorIDs.at(i) << ";X direction; Y direction";
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 40;//every 500 micron there is a bin
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 20;//every 500 micron there is a bin
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -5;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 5;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : -20;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 20;
		AIDA::IProfile2D *  residGblFitY =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (residGblFitY) {
				residGblFitY->setTitle(histTitle);
				_mapFromSensorIDToHistogramY.insert(std::make_pair(_sensorIDs.at(i), residGblFitY));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////The incidence angles for each plane
	for (int i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "IncidenceXZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Incidence local Tx (XZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 40000;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.005;
		AIDA::IHistogram1D * incidenceGblFitXZ = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX); 

		if (incidenceGblFitXZ){
				incidenceGblFitXZ->setTitle(histTitle);
				_mapFromSensorIDToKinkXZ.insert(std::make_pair(_sensorIDs.at(i), incidenceGblFitXZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	for (int i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "IncidenceYZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Incidence local Ty (YZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 40;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.001 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.001;
		AIDA::IHistogram1D * incidenceGblFitYZ = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D(residGblFitHistName, NBinX, MinX, MaxX); 

		if (incidenceGblFitYZ) {
				incidenceGblFitYZ->setTitle(histTitle);
				_mapFromSensorIDToKinkYZ.insert(std::make_pair(_sensorIDs.at(i), incidenceGblFitYZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	/////////////////////////////////////////////////////////////////////////////////////// 
}

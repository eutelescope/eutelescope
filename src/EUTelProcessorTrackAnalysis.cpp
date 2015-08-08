/* To create a new analysis follow the series of steps below:
 * 1) Create a new set of histograms in this processor and link the histogram's pointers to the state's locations; i.e the plane number.
 * 		Histograms are made using the input to the processor and then matched using the input from the states location.
 * 2) Create a new input member function to set this in the class (EUTelTrackAnalysis).
 * 3) Now pass the track from this processor to EUTelTrackAnalysis via a function as shown in processEvent below.
 * 4)You now have the trackand histogram. Do the analysis and output to that histogram or anyone oyu want.   */
#include "EUTelProcessorTrackAnalysis.h"

using namespace eutelescope;

EUTelProcessorTrackAnalysis::EUTelProcessorTrackAnalysis() :
Processor("EUTelProcessorTrackAnalysis"){
	
	registerInputCollection(LCIO::TRACK, "TrackInputCollectionName", "Input track collection name",_trackInputCollectionName,std::string("TrackCandidatesCollection"));
	registerOptionalParameter("SensorIDs", "A vector of the sensor IDs to plot",_sensorIDs,IntVec());
    registerOptionalParameter("HistogramInfoFilename", "Name of histogram info xml file", _histoInfoFileName, std::string("histoinfo.xml"));

	
	
}


void EUTelProcessorTrackAnalysis::init(){
	try{
		initialiseResidualVsPositionHistograms();
		//Some initialised in the constructor in part 2.
		EUTelTrackAnalysis*	analysis = new EUTelTrackAnalysis(_mapFromSensorIDToHistogramX,_mapFromSensorIDToHistogramY,_mapFromSensorIDToKinkXZ,_mapFromSensorIDToKinkYZ,_beamEnergy); 

		//Others here.
		analysis->setSensorIDTo2DPValuesWithPosition(_mapFromSensorIDToPValueHisto);
		analysis->setSensorIDToPValuesVsIncidenceAngleYZ(_mapFromSensorIDToPValuesVsIncidenceYZ);
		analysis->setSensorIDToPValuesVsIncidenceAngleXZ(_mapFromSensorIDToPValuesVsIncidenceXZ);
		analysis->setPValueBeamEnergy(_pValueVsBeamEnergy);
		_analysis = analysis;
	}catch(...){	
		streamlog_out(MESSAGE9)<<"There is an unknown error in EUTelProcessorTrackAnalysis-init()" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}

}

void EUTelProcessorTrackAnalysis::processEvent(LCEvent * evt){
	try{
		EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		}else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}
        streamlog_out(DEBUG2) << "Collection contains data! Continue!" << std::endl;
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        std::vector<EUTelTrack> tracks = reader.getTracks(evt, _trackInputCollectionName);
        for (size_t iTrack = 0; iTrack < tracks.size(); ++iTrack){
            EUTelTrack track = tracks.at(iTrack); 
            _analysis->plotResidualVsPosition(track);	
            _analysis->plotIncidenceAngles(track);
            if(track.getChi2()/track.getNdf() < 5.0){
                _analysis->plotBeamEnergy(track);
                _analysis->plotPValueVsBeamEnergy(track);
            }

            _analysis->plotPValueWithPosition(track);
            _analysis->plotPValueWithIncidenceAngles(track);
        }

	}catch(...){	
		streamlog_out(MESSAGE9)<<"There is an unknown error in EUTelProcessorTrackAnalysis-processEvent" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
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

	std::auto_ptr<EUTelHistogramManager> histoMgr( new EUTelHistogramManager( _histoInfoFileName ));
	EUTelHistogramInfo    * histoInfo;
    bool isHistoManagerAvailable;
	try {
			isHistoManagerAvailable = histoMgr->init( );
	} catch ( std::ios::failure& e ) {
			streamlog_out( ERROR5 ) << "I/O problem with " << _histoInfoFileName << "\n"
							<< "Continuing without histogram manager using default settings"    << std::endl;
			isHistoManagerAvailable = false;
	} catch ( marlin::ParseException& e ) {
			streamlog_out( ERROR5 ) << e.what( ) << "\n"
							<< "Continuing without histogram manager using default settings" << std::endl;
			isHistoManagerAvailable = false;
	}
	isHistoManagerAvailable = false;


	std::stringstream sstm;
	std::string residGblFitHistName;
	std::string histTitle;
	for (size_t i = 0; i < _sensorIDs.size() ; ++i){
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
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
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
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "IncidenceXZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Incidence local Tx (XZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 60;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
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
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "IncidenceYZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Incidence local Ty (YZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 60;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
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
	/////////////////////////////////////////////////////////////////////////////////////The profile of the p-values with incidence.
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << " Profile of p-values Vs IncidenceXZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Profile of p-values Vs Incidence Angle, local Tx (XZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 60;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
		AIDA::IProfile1D *pValueVsIncidenceXZ = marlin::AIDAProcessor::histogramFactory(this)->createProfile1D(residGblFitHistName, NBinX, MinX, MaxX, 0, 1); 

		if (pValueVsIncidenceXZ){
				pValueVsIncidenceXZ->setTitle(histTitle);
				_mapFromSensorIDToPValuesVsIncidenceXZ.insert(std::make_pair(_sensorIDs.at(i), pValueVsIncidenceXZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	for(size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << " Profile of p-values Vs IncidenceYZ" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "Profile of p-values Vs Incidence Angle, local Ty (YZ plane). Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 60;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-0.05 ;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 0.01;
		AIDA::IProfile1D * pValueVsIncidenceYZ = marlin::AIDAProcessor::histogramFactory(this)->createProfile1D(residGblFitHistName, NBinX, MinX, MaxX, 0,1); 

		if (pValueVsIncidenceYZ) {
				pValueVsIncidenceYZ->setTitle(histTitle);
			_mapFromSensorIDToPValuesVsIncidenceYZ.insert(std::make_pair(_sensorIDs.at(i), pValueVsIncidenceYZ));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
	/////////////////////////////////////////////////////////////////////////////////////// 

	/////////////////////////////////////////////////////////////////////////////////////////p-value with position
	for (size_t i = 0; i < _sensorIDs.size() ; ++i){
		sstm << "P-value vs Position" << _sensorIDs.at(i);
		residGblFitHistName = sstm.str();
		sstm.str(std::string());
		sstm << "P-value. Plane " <<  _sensorIDs.at(i);
		histTitle = sstm.str();
		sstm.str(std::string(""));
		histoInfo = histoMgr->getHistogramInfo(residGblFitHistName);
		NBinX = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xBin : 20;
		MinX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMin :-10.6;
		MaxX =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_xMax : 10.6;
		NBinY = ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yBin : 10;
		MinY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMin : -5.3;
		MaxY =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_yMax : 5.3;
		MinZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMin : 0;
		MaxZ =  ( isHistoManagerAvailable && histoInfo ) ? histoInfo->_zMax : 1;
		AIDA::IProfile2D *  pValueHisto =	marlin::AIDAProcessor::histogramFactory(this)->createProfile2D(residGblFitHistName,  NBinX, MinX, MaxX, NBinY, MinY, MaxY, MinZ,MaxZ);
		if (pValueHisto) {
				pValueHisto->setTitle(histTitle);
				_mapFromSensorIDToPValueHisto.insert(std::make_pair(_sensorIDs.at(i), pValueHisto));
		} else {
				streamlog_out(ERROR2) << "Problem booking the " << (residGblFitHistName) << std::endl;
				streamlog_out(ERROR2) << "Very likely a problem with path name. Switching off histogramming and continue w/o" << std::endl;
		}
		sstm.str(std::string(""));
	}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////Beam Energy
	_beamEnergy  = marlin::AIDAProcessor::histogramFactory(this)->createHistogram1D("BeamEnergy", 1000, 0, 6); 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////P-value with energy
	_pValueVsBeamEnergy = marlin::AIDAProcessor::histogramFactory(this)->createProfile1D("p-Value Vs Beam Energy", 50, 0, 6, 0,1); 

}

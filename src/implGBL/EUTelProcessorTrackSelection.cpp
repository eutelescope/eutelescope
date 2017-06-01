#include "EUTelProcessorTrackSelection.h"

using namespace eutelescope;
/// Combined input can create a particular type of track, with hits on certain planes, and a particular chi2.
/**
 * \param [in] TrackInputCollectionName Name of collection containing GBL tracks.
 * \param [in] TrackOutputCollectionName Name of collection containing GBL tracks to output.
 * \param [in] sensorsIDsMustHaveHit If plane specified then this plane must have a hit.
 * \param [in] sensorsIDsMustNotHaveHit If plane specified then this plane must not have a hit.
 * \param [in] chi2NormCut Chi2 cut to apply to all tracks 
 */
EUTelProcessorTrackSelection::EUTelProcessorTrackSelection() :
Processor("EUTelProcessorTrackSelection"){
	
	registerInputCollection(LCIO::TRACK, "TrackInputCollectionName", "Input track collection name",_trackInputCollectionName,std::string("TrackCandidatesCollection"));
    registerOutputCollection(LCIO::TRACK,"TracksOutputCollectionName","Output tracks collection name",_tracksOutputCollectionName,std::string("TrackCollection"));
	registerOptionalParameter("sensorsIDsMustHaveHit", "These sensors must have a hit to pass selection",_mustHave,std::vector<int> () );
	registerOptionalParameter("sensorsIDsMustNotHaveHit", "The sensors must have a missing hit on this plane to pass.",_mustNotHave,std::vector<int> () );
	registerOptionalParameter("chi2NormCut", "The tracks must have chi2/ndf lower than this value ",_chi2NormCut,double(5) );


}


void EUTelProcessorTrackSelection::init(){
	try{
        _selector = new EUTelTrackSelection();
		std::string name("test.root");
		geo::gGeometry().initializeTGeoDescription(name,false);

	}catch(...){	
        streamlog_out(MESSAGE9)<<"There is an unknown error in EUTelProcessorTrackSelection-init" <<std::endl;
		throw marlin::StopProcessingException( this ) ;

	}

}

void EUTelProcessorTrackSelection::processEvent(LCEvent * evt){
	try{
		EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		}else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}
        std::vector<EUTelTrack> tracksOut;
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        std::vector<EUTelTrack> tracks = reader.getTracks(evt, _trackInputCollectionName);
        for (size_t i = 0; i < tracks.size(); ++i){
            EUTelTrack track = tracks.at(i); 
            streamlog_out(DEBUG5)<<"Chi2/ndf " << track.getChi2()/track.getNdf() <<" Cut: " << _chi2NormCut <<std::endl;
            if(track.getChi2()/track.getNdf() < _chi2NormCut){
//            bool pass = _selector->removeTracksWithHitsOnPlane(track,_mustNotHave);
                bool pass = true;
                if(pass){
                    tracksOut.push_back(track);
                }
            }
        }
        outputLCIO(evt,tracksOut);
    }	
    catch (DataNotAvailableException e) {
//		streamlog_out(MESSAGE9) << "Data not avaliable skip event. " << std::endl;
		throw marlin::SkipEventException(this);
	}
	catch(std::string &e){
//		streamlog_out(MESSAGE9) << e << std::endl;
		throw marlin::SkipEventException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(MESSAGE9) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
	catch(...){
		streamlog_out(MESSAGE9)<<"Unknown exception in processEvent function of EUTelProcessorGBLTrackSelection" <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}

	
}

void EUTelProcessorTrackSelection::end(){}

void EUTelProcessorTrackSelection::outputLCIO(LCEvent* evt, std::vector<EUTelTrack>  tracks){
    if(!tracks.empty()){
        for(unsigned int i=0 ; i< tracks.size(); i++){
            streamlog_out(DEBUG1)<<"Found "<<tracks.size()<<" track for event " << evt->getEventNumber() <<".  Track number  " << i <<std::endl;
            tracks.at(i).print();
        }
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        reader.getColVec(tracks, evt, _tracksOutputCollectionName);
    }

}


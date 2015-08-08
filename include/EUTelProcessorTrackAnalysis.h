// LCIO
#include <EVENT/LCCollection.h>
#include "lcio.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

//EUTelescope
#include "EUTelUtility.h"
#include "EUTelEventImpl.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelTrackAnalysis.h"
// AIDA
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IProfile2D.h>
#endif // MARLIN_USE_AIDA
#include "EUTelReaderGenericLCIO.h"

namespace eutelescope {

	class  EUTelProcessorTrackAnalysis : public marlin::Processor {

  	private:
  	DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackAnalysis)
        
    public:

    virtual Processor* newProcessor() {
    	return new  EUTelProcessorTrackAnalysis;
    }

    EUTelProcessorTrackAnalysis();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every event - the working horse.*/
    virtual void processEvent(LCEvent * evt);

   	/** Called after data processing for clean up. **/
		virtual void end();
		void initialiseResidualVsPositionHistograms();
		
		std::string _trackInputCollectionName;
		EUTelTrackAnalysis* _analysis;
		IntVec _sensorIDs;
		std::map< int,  AIDA::IProfile2D* > _mapFromSensorIDToHistogramX;
		std::map< int,  AIDA::IProfile2D* > _mapFromSensorIDToHistogramY;
		std::map< int,  AIDA::IProfile2D* > _mapFromSensorIDToPValueHisto;
		std::map< int,   AIDA::IHistogram1D *> _mapFromSensorIDToKinkXZ;
		std::map< int,  AIDA::IHistogram1D * > _mapFromSensorIDToKinkYZ;
		std::map< int,  AIDA::IProfile1D* > _mapFromSensorIDToPValuesVsIncidenceXZ;
		std::map< int,  AIDA::IProfile1D* > _mapFromSensorIDToPValuesVsIncidenceYZ;
		AIDA::IHistogram1D * _beamEnergy;
		AIDA::IProfile1D *_pValueVsBeamEnergy;
        std::string _histoInfoFileName;
	};

    EUTelProcessorTrackAnalysis gEUTelProcessorTrackAnalysis;

}

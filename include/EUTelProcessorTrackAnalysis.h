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



using namespace lcio;
using namespace marlin;
using namespace std;


namespace eutelescope {

	class  EUTelProcessorTrackAnalysis : public Processor {

  	private:
  	DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackAnalysis)      // prevent users from making (default) copies of processors
        
    public:

    virtual Processor* newProcessor() {
    	return new  EUTelProcessorTrackAnalysis;
    }

    EUTelProcessorTrackAnalysis();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every run.*/
    virtual void processRunHeader(LCRunHeader* run);

    /** Called for every event - the working horse.*/
    virtual void processEvent(LCEvent * evt);

    virtual void check(LCEvent * evt);

   	/** Called after data processing for clean up. **/
		virtual void end();
		void initialiseResidualVsPositionHistograms();
		
		std::string _trackInputCollectionName;
		EUTelTrackAnalysis* _analysis;
		IntVec _sensorIDs;
		map< int,  AIDA::IProfile2D* > _mapFromSensorIDToHistogramX;
		map< int,  AIDA::IProfile2D* > _mapFromSensorIDToHistogramY;
		map< int,   AIDA::IHistogram1D *> _mapFromSensorIDToKinkXZ;
		map< int,  AIDA::IHistogram1D * > _mapFromSensorIDToKinkYZ;

	};

    EUTelProcessorTrackAnalysis gEUTelProcessorTrackAnalysis;

}

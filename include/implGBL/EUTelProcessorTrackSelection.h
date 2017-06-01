/*
 * EUTelProcessorTrackSelection.h 
 * 
 * Created on: April 19th 2015 
 *     author:Alexander Morton 
 * 
 *          This processor will take in GBLTrack objects and remove tracks not meeting all cut requirements.  
 *           
 *          
 *          
 *          
 * 
 * 
 *
 */




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
#include "EUTelTrackSelection.h"
#include "EUTelGBLFitter.h"
#include "EUTelReaderGenericLCIO.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

namespace eutelescope {

	class  EUTelProcessorTrackSelection : public marlin::Processor {

  	private:
  	DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackSelection)
        
    public:

    virtual Processor* newProcessor() {
    	return new  EUTelProcessorTrackSelection;
    }

    EUTelProcessorTrackSelection();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every event - the working horse.*/
    virtual void processEvent(LCEvent * evt);

    /** Called after data processing for clean up. **/
		virtual void end();

    void outputLCIO(LCEvent* evt, std::vector<EUTelTrack>  track);

    EUTelTrackSelection* _selector;
    std::string _trackInputCollectionName;
    std::string _tracksOutputCollectionName;
    double _chi2NormCut;
    std::vector<int> _mustHave;
    std::vector<int> _mustNotHave;

	};

    EUTelProcessorTrackSelection gEUTelProcessorTrackSelection;

}

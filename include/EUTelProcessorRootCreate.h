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
// AIDA
#include <marlin/AIDAProcessor.h>
#include "EUTelReaderGenericLCIO.h"
#include "EUTelGeometryTelescopeGeoDescription.h"


namespace eutelescope {

	class  EUTelProcessorRootCreate : public marlin::Processor {

  	private:
  	DISALLOW_COPY_AND_ASSIGN(EUTelProcessorRootCreate)
        
    public:

    virtual Processor* newProcessor() {
    	return new  EUTelProcessorRootCreate;
    }

    EUTelProcessorRootCreate();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every event - the working horse.*/
    virtual void processEvent(LCEvent * evt);

   	/** Called after data processing for clean up. **/
		virtual void end();
		std::string _trackInputCollectionName;
        std::string _histogramName;

	};

    EUTelProcessorRootCreate gEUTelProcessorRootCreate;

}

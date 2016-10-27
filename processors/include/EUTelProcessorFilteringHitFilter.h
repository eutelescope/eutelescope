// Version: $Id$
#ifndef EUTelProcessorFilteringHitFilter_h
#define EUTelProcessorFilteringHitFilter_h 1

// LCIO
#include "lcio.h"

#include "marlin/Processor.h"

#include "IMPL/TrackerHitImpl.h"
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>

// C++
#include <string>


namespace eutelescope {

  class EUTelProcessorFilteringHitFilter : public marlin::Processor {
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorFilteringHitFilter;
        }

        EUTelProcessorFilteringHitFilter();

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init();

        /** Called for every run.
         */
        virtual void processRunHeader(lcio::LCRunHeader* run);

        /** Called for every event - the working horse.
         */
        virtual void processEvent(lcio::LCEvent * evt);

        virtual void check(lcio::LCEvent * evt);

        /** Called after data processing for clean up.
         */
        virtual void end();


    protected:

        //! Input TrackerHit collection name
        std::string _hitInputCollectionName;
        
        //! Output TrackerHit collection name
        std::string _hitOutputCollectionName;

        //! Input HotPixel collection name
        std::string _hotpixelCollectionName;

        IntVec _wantPlaneID;
        
        int _nProcessedRuns;
        int _nProcessedEvents;

        // treat hits with hotpixels
        std::map<std::string, bool > _hotPixelMap;
 
    };

    //! A global instance of the processor
    EUTelProcessorFilteringHitFilter gEUTelProcessorFilteringHitFilter;
    
} // eutelescope

#endif // EUTelProcessorFilteringHitFilter_h

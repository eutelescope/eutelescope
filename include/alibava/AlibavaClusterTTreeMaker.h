/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVACLUSTERTTREEMAKER_H
#define ALIBAVACLUSTERTTREEMAKER_H 1

// alibava includes ".h"
#include "AlibavaBaseHistogramMaker.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"
#include "TTree.h"

// system includes <>
#include <string>
#include <list>
#include <vector>

namespace alibava {
    
    class AlibavaClusterTTreeMaker :  public alibava::AlibavaBaseProcessor  {
        
    public:
        
        
        //! Returns a new instance of AlibavaClusterTTreeMaker
        /*! This method returns an new instance of the this processor.  It
         *  is called by Marlin execution framework and it shouldn't be
         *  called/used by the final user.
         *
         *  @return a new AlibavaClusterTTreeMaker.
         */
        virtual Processor * newProcessor () {
            return new AlibavaClusterTTreeMaker;
        }
        
        //! Default constructor
        AlibavaClusterTTreeMaker ();
        
        //! Called at the job beginning.
        /*! This is executed only once in the whole execution. It prints
         *  out the processor parameters and reset all needed data
         *  members. In principle this could also be the right place to
         *  book histograms, but since those are also grouped according to
         *  the detector numbers we need to have this parameter available.
         */
        virtual void init ();
        
        //! Called for every run.
        /*! At the beginning of every run the run header is read and
         *  processed by this method. As a first thing, the input run
         *  header is dynamically re-casted to a EUTelRunHeaderImpl and
         *  then important things like the number of detectors and the
         *  pixel detector boundaries are dumped from the file. After that
         *  the EUTelPedestalNoiseProcess::bookHistos() is called.
         *
         *  @param run the LCRunHeader of the this current run
         */
        virtual void processRunHeader (LCRunHeader * run);
        
        //! Called every event
        /*! Since the behavior of the PedestalNoise processor is different
         *  if this is the first or one of the following loop, this method
         *  is just calling
         *  AlibavaClusterTTreeMaker::firstLoop(LCEvent*) or
         *  AlibavaClusterTTreeMaker::otherLoop(LCEvent*)
         *
         *  @param evt the current LCEvent event as passed by the
         *  ProcessMgr
         */
        virtual void processEvent (LCEvent * evt);
        
        
        //! Check event method
        /*! This method is called by the Marlin execution framework as
         *  soon as the processEvent is over. It can be used to fill check
         *  plots. For the time being there is nothing to check and do in
         *  this slot.
         *
         *  @param evt The LCEvent event as passed by the ProcessMgr
         *
         */
        virtual void check (LCEvent * evt);
        
        
        //! Book histograms
        /*! This method is used to prepare the needed directory structure
         *  within the current ITree folder and books all required
         *  histograms. Histogram pointers are stored into
         *  AlibavaBaseProcessor::_rootObjectMap so that they can be
         *  recalled and filled from anywhere in the code.  Apart from the
         *  histograms listed in AlibavaClusterTTreeMaker::fillHistos()
         *  there is also a common mode histo described here below:
         *
         *  @see AlibavaClusterTTreeMaker::fillHistos() for the todos
         */
        virtual void bookHistos();
        
        //! Fill histograms
        /*! This method is used to fill in histograms for each event.
         *
         */
        virtual void fillHistos(TrackerDataImpl * trkdata);
        
        
        //! Called after data processing.
        /*! This method is called when the loop on events is finished. It
         *  is checking whether the calculation is properly finished or
         *  not.
         *  A very common error occurs when the file finished without a
         *  EORE or when the MaxRecordNumber was set to low to loop over
         *  all the needed record. To check this is very easy because we
         *  just have to crosscheck if _iLoop is equal to noOfCMIterations.
         */
        virtual void end();
        
        
    protected:
        
        
        ///////////////
        // Tree Name //
        ///////////////
        
        // Name of the tree
        std::string _treeName;
        
        // TTree that will be created
        TTree *_tree;
        
        // Tree branches
        int _runnumber;
        int _eventnumber;
        int _chipNum;
        int _clusterID;
        int _clusterSize;
        float _totalSignal;
        float _totalSNR;
        int _seedChanNum;
        std::vector<int> *_channums;
        std::vector<double> *_signals;
        std::vector<double> *_snrs;
        bool _isSensitiveAxisX;
        float _signalPolarity;
        
    };
    
    //! A global instance of the processor
    AlibavaClusterTTreeMaker gAlibavaClusterTTreeMaker;
    
}

#endif

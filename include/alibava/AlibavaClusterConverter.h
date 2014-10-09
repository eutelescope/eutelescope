/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVACLUSTERCONVERTER_H
#define ALIBAVACLUSTERCONVERTER_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>


namespace alibava {
	
	//! Example Alibava processor for Marlin.

	
	class AlibavaClusterConverter:public alibava::AlibavaBaseProcessor   {
		
	public:
		
		
		//! Returns a new instance of AlibavaClusterConverter
		/*! This method returns an new instance of the this processor.  It
		 *  is called by Marlin execution framework and it shouldn't be
		 *  called/used by the final user.
		 *
		 *  @return a new AlibavaClusterConverter.
		 */
		virtual Processor * newProcessor () {
			return new AlibavaClusterConverter;
		}
		
		//! Default constructor
		AlibavaClusterConverter ();
		
		//! Called at the job beginning.
		virtual void init ();
		
		//! Called for every run.
		virtual void processRunHeader (LCRunHeader * run);
		
		//! Called every event
		virtual void processEvent (LCEvent * evt);
				
		//! Check event method
		virtual void check (LCEvent * evt);
		
		//! Book histograms
		/*! This method is used to book histograms
		 */
		void bookHistos();
		
		//! Fill histograms
		/*! This method is used to fill in histograms for each channel. 
		 *
		 */
		void fillHistos(TrackerDataImpl * trkdata);
		
		//! Called after data processing.
		/*! This method is called when the loop on events is finished. It
		 *  is checking whether the calculation is properly finished or
		 *  not.
		 */
		virtual void end();
		
		// cluster collection names for EUTel
		// The collection name of cluster pulse
		string _pulseCollectionName;
		// The collection name of sparse cluster
		string _sparseCollectionName;
		
		
		// SensorID
		int _sensorIDStartsFrom;
		
		// missing coordinate value
		// The value that should be stored in missing coordinate.
		// This number has to be integer since it will be used as channel number of the missing coordinate
		int _missingCorrdinateValue;
		
		
	protected:
			
	
		
	};
	
	//! A global instance of the processor
	AlibavaClusterConverter gAlibavaClusterConverter;
	
}

#endif

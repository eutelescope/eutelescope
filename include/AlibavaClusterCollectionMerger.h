/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


#ifndef ALIBAVACLUSTERCOLLECTIONMERGER_H
#define ALIBAVACLUSTERCOLLECTIONMERGER_H 1

// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaRunHeaderImpl.h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>
#include <vector>

namespace alibava {
	
	
	class AlibavaClusterCollectionMerger : public marlin::DataSourceProcessor    {
		
	public:
		
		//! Default constructor
		AlibavaClusterCollectionMerger ();
		
		//! New processor
		/*! Return a new instance of a AlibavaClusterCollectionMerger. It is
		 *  called by the Marlin execution framework and shouldn't be used
		 *  by the final user.
		 */
		virtual AlibavaClusterCollectionMerger * newProcessor ();
		
		//! Creates events from the eudaq software
		virtual void readDataSource (int numEvents);
		
		//! Init method
		/*! It is called at the beginning of the cycle and it prints out
		 *  the parameters.
		 */
		virtual void init ();
		
		//! End method
		/*! It prints out a good bye message
		 */
		virtual void end ();
		
		
	protected:
		
		///////////////
		// Telescope //
		///////////////

		//! This is the input file name that telescope cluster collections stored
		std::string _telescopeFileName;
		
		//! Name of the cluster pulse collection of telescope data
		std::string _telescopePulseCollectionName;
		
		//! Name of the sparse cluster collection of telescope data
		std::string _telescopeSparseCollectionName;

		
		/////////////
		// Alibava //
		/////////////
		
		//! This is the input file name that alibava cluster collections stored
		std::string _alibavaFileName;
		
		//! Name of the cluster pulse collection of alibava data
		std::string _alibavaPulseCollectionName;
		
		//! Name of the sparse cluster collection of alibava data
		std::string _alibavaSparseCollectionName;

		
		////////////
		// Output //
		////////////
		
		//! Name of the merged/output cluster pulse collection
		std::string _outputPulseCollectionName;
		
		//! Name of the merged/output sparse cluster collection.
		/*! This parameter is hard coded in other EUTelescope
		 *  processors. It has to be "original_zsdata"
		 */
		std::string _outputSparseCollectionName;
		
		//! Event ID difference
		int _eventIDDiff;
	
	private:
		
		void copyClustersInCollection(LCCollectionVec * outputPulseColVec, LCCollectionVec * outputSparseColVec, LCCollectionVec * inputPulseColVec, LCCollectionVec * inputSparseColVec );

	};
	
	//! A global instance of the processor
	AlibavaClusterCollectionMerger gAlibavaClusterCollectionMerger;
	
} // end of alibava namespace
#endif


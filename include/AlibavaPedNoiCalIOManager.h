/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVAPEDNOICALIOMANAGER_H
#define ALIBAVAPEDNOICALIOMANAGER_H 1
// alibava includes ".h"
#include "ALIBAVA.h"
#include "AlibavaRunHeaderImpl.h"

// lcio includes <.h>
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/LCCollectionVec.h>

// system includes <>
#include <string>

using namespace std;
using namespace lcio;

namespace alibava {
	
	class AlibavaPedNoiCalIOManager{
		
	public:
		AlibavaPedNoiCalIOManager();
		~AlibavaPedNoiCalIOManager();
		
		void createFile(string filename, IMPL::LCRunHeaderImpl* runHeader);
		
		void addToFile(string filename, string collectionName, int chipnum, EVENT::FloatVec datavec);
		
		EVENT::FloatVec getPedNoiCalForChip(string filename, string collectionName, unsigned int chipnum);
		
	private:
		// gets data vector from an event
		EVENT::FloatVec getDataFromEventForChip(LCEvent* evt, string collectionName, unsigned int chipnum);
		
		// returns true if the collection exists in the event
		bool doesCollectionExist(LCEvent* evt, string collectionName);
		
		// returns the element number of the TrackerData with the chosen chipnum in collection col.
		// returns -1 if it doesn't exist.
		int getElementNumberOfChip(LCCollectionVec* col, int chipnum);
		
		// checks if the file exists
		bool doesFileExist(string astring);
		
		// to create a file with an empty runheader
		void createFile(string filename);
		
		// returns a copy of runHeader of existing AlibavaPedNoiCalFile
		LCRunHeaderImpl* getRunHeader(string filename);
		
		// returns a copy of lcevent of existing AlibavaPedNoiCalFile
		LCEventImpl* getEvent(string filename);
		
		
		
	};
}


#endif

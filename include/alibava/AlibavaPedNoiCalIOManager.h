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

namespace alibava {
	
	class AlibavaPedNoiCalIOManager{
		
	public:
		AlibavaPedNoiCalIOManager();
		~AlibavaPedNoiCalIOManager();
		
		void createFile(std::string filename, lcio::LCRunHeaderImpl* runHeader);
		
		void addToFile(std::string filename, std::string collectionName, int chipnum, lcio::FloatVec datavec);
		
		lcio::FloatVec getPedNoiCalForChip(std::string filename, std::string collectionName, unsigned int chipnum);
		
	private:
		// gets data vector from an event
		lcio::FloatVec getDataFromEventForChip(lcio::LCEvent* evt, std::string collectionName, unsigned int chipnum);
		
		// returns true if the collection exists in the event
		bool doesCollectionExist(lcio::LCEvent* evt, std::string collectionName);
		
		// returns the element number of the TrackerData with the chosen chipnum in collection col.
		// returns -1 if it doesn't exist.
		int getElementNumberOfChip(lcio::LCCollectionVec* col, int chipnum);
		
		// checks if the file exists
		bool doesFileExist(std::string astring);
		
		// to create a file with an empty runheader
		void createFile(std::string filename);
		
		// returns a copy of runHeader of existing AlibavaPedNoiCalFile
		lcio::LCRunHeaderImpl* getRunHeader(std::string filename);
		
		// returns a copy of lcevent of existing AlibavaPedNoiCalFile
		lcio::LCEventImpl* getEvent(std::string filename);
		
		
		
	};
}


#endif

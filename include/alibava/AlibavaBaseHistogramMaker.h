/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */

#ifndef ALIBAVABASEHISTOGRAMMAKER_H
#define ALIBAVABASEHISTOGRAMMAKER_H 1

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
#include <vector>

namespace alibava {
	
	// Helper classes
	class ReplaceXMLValues {
	public:
		ReplaceXMLValues();
		std::string _histoName; // name of the histogram this will apply
		std::string _variableName; // Hasto be one of these: xBin, xMin, xMax, yBin, yMin, yMax
		float _newValue; // for xBin and yBin will changed to interger!
	};

	class InsertToXMLTitle {
	public:
		InsertToXMLTitle();
		std::string _histoName; // name of the histogram this will apply
		std::string _titleName; // Hasto be one of these: title, labelX, labelY
		std::string _whichSide; // insert to left or right
		std::string _stringToAdd; // what to insert
	};
	
	class AlibavaBaseHistogramMaker : public alibava::AlibavaBaseProcessor   {
		
	public:
		
		
		//! Default constructor
		AlibavaBaseHistogramMaker (std::string processorName);
		
		// Inherited from marlin::Processor
		/* init()
		 * processRunHeader(LCRunHeader * run)
		 * processEvent (LCEvent * evt)
		 * end()
		 * bookHistos()
		 * fillHistos(TrackerData* )
		 */
		
		////////////////////
		// Histo XML File //
		////////////////////
		// Path of the XML file containin the list of histograms
		std::string _histoXMLFileName;
		// Top tag in the XML file
		std::string _topTag;
		// The tag (inside the top tag) that should be processed by this histogram maker
		std::string _tagToProcess;
		
		// Reads and creates the histograms defined in _histoXMLFileName
		void processHistoXMLFile();

		// checks if the histogram list in XML file is complete
		bool checkListOfHistosCreatedByXMLFile();

		////////////////////////
		// Pedestal and Noise //
		////////////////////////
		// Chooses if pedestal and noise should be plotted
		bool _plotPedestalAndNoise;
		// function to book and plot pedestal and noise histograms
		void plotPedestalAndNoise();
		
		
		///////////
		// Event //
		///////////
		// Stores the event numbers that user wants to plot
		IntVec _plotEvents;
		// Plot random events
		float _plotXPercentOfEvents;
		
		//! The function that returns name of the event histogram for each chip
		std::string getEventHistoName(int eventnum, int ichip);
		
		// Books histogram for the event
		virtual void bookEventHisto(int eventnum) = 0;

		// plots histogram for the event
		virtual void fillEventHisto(int eventnum, TrackerDataImpl * trkdata)=0;
		
		////////////
		// Others //
		////////////
		virtual void fillListOfHistos()=0;
		
		// if you want to change any variable defined in HistoXMLFile use this function
		virtual void createRulesToChangeXMLValues() = 0;
		/////////////////////////////
		// To change XML variables //
		/////////////////////////////
		
		// changeXMLVariable method is used to replace xBin, xMin, xMax, yBin, yMin, yMax
		// so variableName has to be one of these: xBin, xMin, xMax, yBin, yMin, yMax (case sensitive)
		// for xBin and yBin newValue will changed to integer
		void changeXMLVariable(std::string histoName, std::string variableName, float newValue);
		
		// addToXMLTitle method is used to change title, labelX and labelY
		// so titleName has to be one of these: title, labelX and labelY
		// std::stringToAdd is the std::string that will added to the defined title
		// whichSide defines which side of the title you want to add the "stringToAdd",
		//          if you set it "left":  title = std::stringToAdd + title
		//          if you set it "right": title = title + std::stringToAdd
		//          if you set it "all":   title = std::stringToAdd
		
		void addToXMLTitle(std::string histoName, std::string titleName, std::string whichSide, std::string stringToAdd);
		
		// adds the histo name to the histogram check list
		void addToHistoCheckList(std::string histoName);

		// adds the histo name per chip to the histogram check list
		void addToHistoCheckList_PerChip(std::string histoName);
		
	protected:

		////////////////////////
		// Pedestal and Noise //
		////////////////////////
		
		//! Name of the Pedestal histogram
		std::string _pedestalHistoName;
		//! Name of the Noise histogram
		std::string _noiseHistoName;
		//! Name of the directory to store Pedestal and Noise histograms
		std::string _pedestalAndNoiseHistoDir;
		//! The function that returns name of the pedestal histogram for each chip
		std::string getPedestalHistoName(int ichip);
		//! The function that returns name of the noise histogram for each chip
		std::string getNoiseHistoName(int ichip);

		///////////
		// Event //
		///////////

		//! Name of the Event histogram
		std::string _eventHistoName;
		//! Name of the directory to store Event histograms
		std::string _eventHistoDir;
		//! Name of the masked Events histogram
		std::string _maskedEventsHistoName;
		// checks if eventnum is in _plotEvents vector
		bool isEventToBePlotted(int eventnum);

		/////////////////
		// Calibration //
		/////////////////
		//! Name of the calibration charge histogram
		std::string _calChargeHistoName;
		//! Name of the delay histogram
		std::string _delayHistoName;
		
		////////////
		// Others //
		////////////
		// Histogram lists that has to be defined in _histoXMLFileName
		std::vector<std::string> _histoList;

		
		//! If one wants to multiply signal by some number
		float _multiplySignalby;
		
		// total number of event
		int _totalNumberOfEvents;
		
		// just adding "_chip" and chipnumber to the end of std::string
		std::string getHistoNameForChip(std::string histoName, int ichip);

		// To change XML variables
		std::vector<InsertToXMLTitle> _insertToXMLTitle;
		std::vector<ReplaceXMLValues> _replaceXMLValues;

		void updateXMLVariables(std::string histoName, int* xBin, float* xMin, float* xMax, int* yBin, float* yMin, float* yMax, std::string* title, std::string* labelX, std::string* labelY);
	};
		
}

#endif

/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


#ifndef ALIBAVABASEPROCESSOR_H
#define ALIBAVABASEPROCESSOR_H 1

// alibava includes ".h"
#include "ALIBAVA.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <EVENT/LCParameters.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>

using namespace std;

namespace alibava {
	
	//! Pedestal and noise  processor for Marlin.
	
	
	class AlibavaBaseProcessor:public marlin::Processor   {
		
	public:
		
		//! Default constructor
		AlibavaBaseProcessor (std::string processorName);
	
		
		//! Root Object Map
		/*! The object might be needed in different
		 *  places, while it is usually a good reason to keep the booking
		 *  procedure in one place only. To recall an object pointer
		 *  from a different place in the code they are stored within this
		 *  object map. The second object of the pair has to a
		 *  TObjecy since this is the base class of all
		 *  different kind of root objects. It has also
		 *  to be a pointer since, this is a pure virtual class and we
		 *  want to use <code>dynamic_cast</code> to convert them back to
		 *  their original cast.
		 */
		std::map<std::string , TObject * > _rootObjectMap;
		
		// checks if the root object exists in _rootObjectMap
		bool doesRootObjectExists(std::string aHistoName);
		
		
		///////////////////////////
		// Input/Output Collection
		///////////////////////////

		// getter and setter for _inputCollectionName
		void setInputCollectionName(std::string inputCollectionName);
	   std::string getInputCollectionName();

		// getter and setter for _outputCollectionName
		void setOutputCollectionName(std::string outputCollectionName);
	   std::string getOutputCollectionName();

		///////////////////////////
		// Pedestal & Noise
		///////////////////////////
		
		// Sets pedestal and noise values
		// Note that chip selection, pedestalfile name and pedestal and noise collection names have to be set before calling this method.
		void setPedestals();
		
		// Checks the size and existance of the pedestal and noise values of each selected chips
		void checkPedestals();
		
		///////////////////////////
		// Pedestal
		///////////////////////////

		// getter and setter for _pedestalCollectionName
		void setPedestalCollectionName(std::string pedestalCollectionName);
		std::string getPedestalCollectionName();

		// to access the pedestal values of a chip
		EVENT::FloatVec getPedestalOfChip(int chipnum);
		
		// to access the pedestal value of a channel
		float getPedestalAtChannel(int chipnum, int channum);

		///////////////////////////
		// Noise
		///////////////////////////
		// getter and setter for _noiseCollectionName
		void setNoiseCollectionName(std::string noiseCollectionName);
		std::string getNoiseCollectionName();
				
		// to access the noise values of a chip
		EVENT::FloatVec getNoiseOfChip(int chipnum);
		
		// to access the noise value of a channel
		float getNoiseAtChannel(int chipnum, int channum);

		///////////////////////////
		// Charge Calibration
		///////////////////////////
	
		// Sets charge calibration values
		// Note that chip selection, calibrationFile name and charge calibration collection names have to be set before calling this method.
		void setCalibration();
		
		// Checks the size and existance of the charge calibration values of each selected chips
		void checkCalibration();
		
		// getter and setter for _chargeCalCollectionName
		void setChargeCalCollectionName(std::string chargeCalCollectionName);
		std::string getChargeCalCollectionName();

		// to access the chargeCal values of a chip
		EVENT::FloatVec getChargeCalOfChip(int chipnum);
		
		// to access the chargeCal value of a channel
		float getChargeCalAtChannel(int chipnum, int channum);

		
		///////////////////////////
		// Others
		///////////////////////////
		// getter and setter for _nChips
		void setChipSelection(EVENT::IntVec chipselection);
		EVENT::IntVec getChipSelection();

		// get number of selected chips
		unsigned int getNumberOfChips();
		
		// to access the mask value of a channel		
		bool isMasked(int chipnum, int ichan);
		
		//! Applies _channelsToBeUsed parameter
		/*! Make sure you set _channelsToBeUsed parameter
		 *  and _nChips before using this function
		 */
		void setChannelsToBeUsed();
		
		// decodes channel masking string
		void decodeMaskingString(string istring, int *onchip, int *fromchannel, int *tochannel );
		// checks if the decoded channel masking makes sense
		bool isMaskingValid(int onchip, int fromchannel, int tochannel );
		// prints channel masking
		void printChannelMasking();

		// returns chip number of the TrackerDataImpl using ALIBAVA::ALIBAVADATA_ENCODE and ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM
		int getChipNum(TrackerDataImpl * trkdata);
		
		//returns the min chip number
		unsigned int getMinChipNumber();
		//returns the max chip number
		unsigned int getMaxChipNumber();
		
		// returns true if the chip is in the list of selected chip
		bool isChipValid(int ichip);
		// returns true if the channel number is valid
		bool isChannelValid(int ichan);
		
		// !!!
		// All the parameters that will be registered should be a public member!!!
		// !!!
		
		
		//! Input collection name.
		/*! A string containing input collection name
		 */
		std::string _inputCollectionName;

		//! Output collection name
		/*! A string containing output collection name
		 */
		std::string _outputCollectionName;
		
		//! Pedestal file path.
		/*! Path of the alibava pedestal file 
		 */
		std::string _pedestalFile;
		
		//! Charge calibration file path.
		/*! Path of the alibava charge calibration file
		 */
		std::string _calibrationFile;
		
		//! Pedestal collection name.
		/*! We have three different output collections that are going to be
		 *  saved in the file within the first event possibly so that the
		 *  next processor will have access to the pedestal, noise and
		 *  status immediately.  We have a collection for the pedestal, one
		 *  for the noise and another for the status.  The first two are
		 *  self explaining, while few more words are needed for the status
		 *  one.  This one is used, for example, to set a bad pixel mask or
		 *  to flag pixel as already belonging to an already reconstructed
		 *  cluster.
		 */
		std::string _pedestalCollectionName;
		
		//! Noise collection name.
		/*! See _pedestalCollectionName for the detailed description
		 */
		std::string _noiseCollectionName;

		//! Charge calibration collection name.
		/*! See _pedestalCollectionName for the detailed description
		 */
		std::string _chargeCalCollectionName;

		//! StringVec parameter for the processor to set the channels to be used
		/*! The format of _channelsToBeUsed string parameter shoul be like $ChipNumber:StartChannel-EndChannel$
		 *   ex. $0:5-20$ $0:30-100$ $1:50-70$ means from chip 0 channels between 5-20 and 30-100, from chip 1 channels between 50-70
		 *    will be used (all numbers included).
		 *   the rest will be masked and not used
		 *   Note that the numbers should be in ascending order
		 *   and there should be no space between two $ character
		 */

		EVENT::StringVec _channelsToBeUsed;
	
		////////////////////////////
		// Skipping Masked Events //
		////////////////////////////
		
		// choose if you want to skip masked event
		bool _skipMaskedEvents;
		
		// number of skipped Events
		int _numberOfSkippedEvents;

		// returns if pedestal is set and valid
		bool isPedestalValid();

		// returns if noise is set and valid
		bool isNoiseValid();
		
	protected:
		//! Selected chips
		/*! This is the selected chip numbers in the input collection
		 *
		 */
		EVENT::IntVec _chipSelection;
		
		//! Array to store mask information for all channels
		/*!
		 *  Channel ichan will not be used if _isMasked[ichan]=true
		 */
		bool _isMasked[ALIBAVA::NOOFCHIPS][ALIBAVA::NOOFCHANNELS];
		
		void setAllMasksTo(bool abool);
		
		
		
		// a map to store pedestal values for chips
		std::map<int , EVENT::FloatVec > _pedestalMap;
		
		// a map to store noise values for chips
		std::map<int , EVENT::FloatVec > _noiseMap;
		
		// a map to store charge calibration values for chips
		std::map<int , EVENT::FloatVec > _chargeCalMap;

		
		bool _isPedestalValid;
		bool _isNoiseValid;
		
		bool _isCalibrationValid;

	};
	
	//! A global instance of the processor
	//AlibavaBaseProcessor gAlibavaBaseProcessor;
	
	

	
}

#endif


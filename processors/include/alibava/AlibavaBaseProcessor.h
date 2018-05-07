/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
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

namespace alibava
{
    // Pedestal and noise  processor for Marlin.
    class AlibavaBaseProcessor : public marlin::Processor
    {

	public:

	    // Default constructor
	    AlibavaBaseProcessor ( std::string processorName );

	    std::map < std::string , TObject * > _rootObjectMap;

	    // checks if the root object exists in _rootObjectMap
	    bool doesRootObjectExists ( std::string aHistoName );

	    // Input/Output Collection
	    // getter and setter for _inputCollectionName
	    void setInputCollectionName ( std::string inputCollectionName );
	    std::string getInputCollectionName ( );

	    // getter and setter for _outputCollectionName
	    void setOutputCollectionName ( std::string outputCollectionName );
	    std::string getOutputCollectionName ( );


	    // Pedestal & Noise
	    // Sets pedestal and noise values
	    // Note that chip selection, pedestalfile name and pedestal and noise collection names have to be set before calling this method.
	    void setPedestals ( );

	    // Checks the size and existance of the pedestal and noise values of each selected chips
	    void checkPedestals ( );

	    // Pedestal
	    // getter and setter for _pedestalCollectionName
	    void setPedestalCollectionName ( std::string pedestalCollectionName );
	    std::string getPedestalCollectionName ( );

	    // to access the pedestal values of a chip
	    EVENT::FloatVec getPedestalOfChip ( int chipnum );

	    // to access the pedestal value of a channel
	    float getPedestalAtChannel ( int chipnum, int channum );

	    // Noise
	    // getter and setter for _noiseCollectionName
	    void setNoiseCollectionName ( std::string noiseCollectionName );
	    std::string getNoiseCollectionName ( );

	    // to access the noise values of a chip
	    EVENT::FloatVec getNoiseOfChip ( int chipnum );

	    // to access the noise value of a channel
	    float getNoiseAtChannel ( int chipnum, int channum );

	    // Charge Calibration
	    // Sets charge calibration values
	    // Note that chip selection, calibrationFile name and charge calibration collection names have to be set before calling this method.
	    void setCalibration ( );

	    // Checks the size and existance of the charge calibration values of each selected chips
	    void checkCalibration ( );

	    // getter and setter for _chargeCalCollectionName
	    void setChargeCalCollectionName ( std::string chargeCalCollectionName );
	    std::string getChargeCalCollectionName ( );

	    // to access the chargeCal values of a chip
	    EVENT::FloatVec getChargeCalOfChip ( int chipnum );

	    // to access the chargeCal value of a channel
	    float getChargeCalAtChannel ( int chipnum, int channum );

	    // Others
	    // getter and setter for _nChips
	    void setChipSelection ( EVENT::IntVec chipselection );
	    EVENT::IntVec getChipSelection ( );

	    // get number of selected chips
	    unsigned int getNumberOfChips ( );

	    // to access the mask value of a channel
	    bool isMasked ( int chipnum, int ichan );

	    //! Applies _channelsToBeUsed parameter
	    void setChannelsToBeUsed ( );

	    // decodes channel masking std::string
	    void decodeMaskingString ( std::string istring, int *onchip, int *fromchannel, int *tochannel );

	    // checks if the decoded channel masking makes sense
	    bool isMaskingValid ( int onchip, int fromchannel, int tochannel );

	    // prints channel masking
	    void printChannelMasking ( );

	    // returns chip number of the TrackerDataImpl using ALIBAVA::ALIBAVADATA_ENCODE and ALIBAVA::ALIBAVADATA_ENCODE_CHIPNUM
	    int getChipNum ( TrackerDataImpl * trkdata );

	    // returns the min chip number
	    unsigned int getMinChipNumber ( );

	    //returns the max chip number
	    unsigned int getMaxChipNumber ( );

	    // returns true if the chip is in the list of selected chip
	    bool isChipValid ( int ichip );

	    // returns true if the channel number is valid
	    bool isChannelValid ( int ichan );

	    //! Input collection name.
	    std::string _inputCollectionName;

	    // Output collection name
	    std::string _outputCollectionName;

	    // Pedestal file path.
	    std::string _pedestalFile;

	    // Charge calibration file path.
	    std::string _calibrationFile;

	    // Pedestal collection name.
	    std::string _pedestalCollectionName;

	    // Noise collection name.
	    std::string _noiseCollectionName;

	    // Charge calibration collection name.
	    std::string _chargeCalCollectionName;

	    // StringVec parameter for the processor to set the channels to be used
	    /*! The format of _channelsToBeUsed std::string parameter shoul be like $ChipNumber:StartChannel-EndChannel$
	     *   ex. $0:5-20$ $0:30-100$ $1:50-70$ means from chip 0 channels between 5-20 and 30-100, from chip 1 channels between 50-70
	     *    will be used (all numbers included).
	     *   the rest will be masked and not used
	     *   Note that the numbers should be in ascending order
	     *   and there should be no space between two $ character
	     */
	    EVENT::StringVec _channelsToBeUsed;

	    // Skipping Masked Events
	    bool _skipMaskedEvents;

	    // number of skipped Events
	    int _numberOfSkippedEvents;

	    // returns if pedestal is set and valid
	    bool isPedestalValid ( );

	    // returns if noise is set and valid
	    bool isNoiseValid ( );

	protected:

	    // Selected chips
	    EVENT::IntVec _chipSelection;

	    // Array to store mask information for all channels
	    bool _isMasked[ALIBAVA::NOOFCHIPS][ALIBAVA::NOOFCHANNELS];

	    void setAllMasksTo ( bool abool );

	    // a map to store pedestal values for chips
	    std::map < int , EVENT::FloatVec > _pedestalMap;

	    // a map to store noise values for chips
	    std::map < int , EVENT::FloatVec > _noiseMap;

	    // a map to store charge calibration values for chips
	    std::map < int , EVENT::FloatVec > _chargeCalMap;

	    bool _isPedestalValid;
	    bool _isNoiseValid;

	    bool _isCalibrationValid;

    };

    //! A global instance of the processor
    //AlibavaBaseProcessor gAlibavaBaseProcessor;
}

#endif

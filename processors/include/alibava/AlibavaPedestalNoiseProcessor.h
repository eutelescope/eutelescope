/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
 */

#ifndef ALIBAVAPEDESTALNOISEPROCESSOR_H
#define ALIBAVAPEDESTALNOISEPROCESSOR_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerRawDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>

namespace alibava
{
    //! Pedestal and noise  processor for Marlin.
    class AlibavaPedestalNoiseProcessor : public alibava::AlibavaBaseProcessor
    {
	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaPedestalNoiseProcessor;
	    }

	    AlibavaPedestalNoiseProcessor ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( TrackerDataImpl * trkdata );

	    virtual void end ( );

	protected:

	    // Name of the Pedestal histogram 
	    std::string _pedestalHistoName;
	    // Name of the Noise histogram 
	    std::string _noiseHistoName;

	    // Name of the Temperature histogram
	    std::string _temperatureHistoName;

	    // The name of the histogram used to calculate pedestal and noise
	    std::string _chanDataHistoName;

	    //! The name of the fits used to calculate pedestal and noise
	    std::string _chanDataFitName;

	    //! The function that returns name of the histogram for each channel
	    std::string getChanDataHistoName ( unsigned int ichip, unsigned int ichan );

	    //! The function that returns name of the fit for each channel
	    std::string getChanDataFitName ( unsigned int ichip, unsigned int ichan );

	    //! The function that returns name of the pedestal histogram for each chip
	    std::string getPedestalHistoName ( unsigned int ichip );

	    //! The function that returns name of the noise histogram for each chip
	    std::string getNoiseHistoName ( unsigned int ichip );

	    //! Calculates and saves pedestal and noise values
	    void calculatePedestalNoise ( );

    };

    //! A global instance of the processor
    AlibavaPedestalNoiseProcessor gAlibavaPedestalNoiseProcessor;
}

#endif

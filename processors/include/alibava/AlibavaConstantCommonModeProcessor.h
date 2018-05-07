/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 *
 *  modified by: Eda Yildirim eda.yildirim@cern.ch
 */

#ifndef ALIBAVACONSTANTCOMMONMODEPROCESSOR_H
#define ALIBAVACONSTANTCOMMONMODEPROCESSOR_H 1

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

namespace alibava
{
    //! Common mode processor for Marlin.

    class AlibavaConstantCommonModeProcessor : public alibava::AlibavaBaseProcessor
    {
	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaConstantCommonModeProcessor;
	    }

	    AlibavaConstantCommonModeProcessor ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( TrackerDataImpl * trkdata, int event );

	    virtual void end ( );

	    std::string _commonmodeCollectionName;

	    std::string _commonmodeerrorCollectionName;

	    int _Niteration;

	    float _NoiseDeviation;

	    std::string _commonmodeMethod;

	    void setCommonModeCollectionName ( std::string CommonModeCollectionName );
	    std::string getCommonModeCollectionName ( );

	    void setCommonModeErrorCollectionName ( std::string CommonModeErrorCollectionName );
	    std::string getCommonModeErrorCollectionName ( );

	    void setCommonModeVec ( EVENT::FloatVec common );
	    EVENT::FloatVec getCommonModeVec ( );

	    void setCommonModeErrorVec ( EVENT::FloatVec commonerror );
	    EVENT::FloatVec getCommonModeErrorVec ( );

	protected:

	    std::string _commonmodeHistoName;

	    std::string _commonmodeerrorHistoName;

	    void calculateConstantCommonMode ( TrackerDataImpl * trkdata );

	    std::string getCommonCorrectionName ( );

	    EVENT::FloatVec _commonmode;

	    EVENT::FloatVec _commonmodeerror;

    };

    AlibavaConstantCommonModeProcessor gAlibavaConstantCommonModeProcessor;
}

#endif

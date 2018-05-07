/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef ALIBAVACLUSTERING_H
#define ALIBAVACLUSTERING_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"
#include "AlibavaEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>

// system includes <>
#include <string>
#include <list>

using namespace std;

    //! This is used for the landau gaus function and fits
    //Details see: http://root.cern.ch/root/html/tutorials/fit/langaus.C.html

    Double_t langaufun ( Double_t *x, Double_t *par );
    TF1 *langaufit ( TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF );
    Int_t langaupro ( Double_t *params, Double_t &maxx, Double_t &FWHM );

namespace alibava
{

    //! Example Alibava processor for Marlin.
    class AlibavaClustering:public alibava::AlibavaBaseProcessor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaClustering;
	    }

	    AlibavaClustering ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( );

	    void fillClusterHisto ( int clusize );

	    void fillEtaHisto ( float etaratio );

	    void fillEtaHisto2 ( float etaratio );

	    void fillEtaHisto2TDC ( float etaratio, float tdc );

	    void fillEtaHistoPos ( float etaratio, int ichan );

	    void fillChargeDistHisto ( float a, float b, float c, float d, float e, float f, float g );

	    void fillSNRHisto ( float signal );

	    void fillSignalHisto ( float signal );

	    void fillHitmapHisto ( int ichan, int negclustersize, int posclustersize );

	    void fillSeedHisto ( int ichan );

	    void fillSeedChargeHisto ( float signal );

	    void fillCogHisto ( float cog );

	    void findSeedClusters  (TrackerDataImpl * trkdata, LCCollectionVec * clusterCollection, LCCollectionVec * sparseClusterCollectionVec, AlibavaEventImpl * alibavaEvent );

	    virtual void end ( );

	    void setClusterCollectionName ( std::string clusterCollectionName );

	    std::string getClusterCollectionName ( );
	    std::string _clusterCollectionName;		

	    float _seedcut;

	    float _clustercut;

	    float _clustercharge[5];

	    int _clustercount;

	    std::string _sparseclusterCollectionName;

	    int _polarity;

	    std::string _nonsensitiveaxis;

	    void dolandaugausfit ( string tempHistoName );

	    int _clusterminsize;
	    int _clustermaxsize;

	    bool _usefir;
	    bool _writecoefficients;
	    bool _writezero;
	    bool _readcoefficients;

	    double _readcoefficient1;
	    double _readcoefficient2;
	    double _initcoefficient1;
	    double _initcoefficient2;

	    LCCollectionVec * filteredcollectionVec;

	    string _filteredCollectionName;

	    string _filterFileName;

	    void fillclusterspereventhisto ( int a, int b );

	protected:

	    IMPL::LCRunHeaderImpl* _runHeader;

    };

    //! A global instance of the processor
    AlibavaClustering gAlibavaClustering;

}

#endif

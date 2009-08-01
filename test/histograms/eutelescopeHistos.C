// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#include "eutelescopeHistos.h"

using namespace std;


// implementations

void showPedeNoisePlot( const char * filename, const char *  detector  ) {

  // check if this file is already open
  TFile * inputFile = closeAndReopenFile( filename );

  // typical names for the processor folder
  vector< string > folderNames;
  folderNames.push_back( pedeProcessorFolderName.Data() );
  folderNames.push_back( "Pedestal" );

  // DUT related folder names
  string dutFolderName = pedeProcessorFolderName.Data() + toString( "_dut" );
  folderNames.push_back( dutFolderName );
  dutFolderName = pedeProcessorFolderName.Data() + toString( "_dep" );
  folderNames.push_back( dutFolderName );
  dutFolderName = pedeProcessorFolderName.Data() + toString( "_taki" );
  folderNames.push_back( dutFolderName );
  dutFolderName = pedeProcessorFolderName.Data() + toString( "_tel" );
  folderNames.push_back( dutFolderName );

  // this function will create the following canvases:
  //
  // -> 1 general canvas with the pede map and pede histo (3 detectors per canvas)
  // -> 1 general canvas with the noise map and noise histo for each detectors (3 detectors per canvas)
  // -> 1 general canvas with the status map for each destector. (6 detectors per canvas )

  // PEDE MAP and DIST
  UInt_t nDetPerCanvas =  3;
  TString pedeCanvasBaseName    = "PedeCanvas";

  // close all canvases with these names
  closeCanvases( pedeCanvasBaseName );

  // look into the input file for a folder named
  TDirectoryFile * pedeProcessorFolder = checkFolder( folderNames, inputFile );
  if ( pedeProcessorFolder == 0x0 ) {
    return;
  }

  Int_t loop = 1;
  vector< string > loopFolderNames;
  loopFolderNames.push_back( "loop-" + toString( loop ) );
  loopFolderNames.push_back( "loop_" + toString( loop ) ) ;


  TDirectoryFile * loopFolder = NULL;
  for ( size_t iName = 0; iName < loopFolderNames.size() ; ++iName ) {

    loopFolder = (TDirectoryFile *) pedeProcessorFolder->Get( loopFolderNames.at( iName ).c_str() );

    if ( loopFolder != NULL ) break;
  }

  if ( loopFolder == NULL ) {
    cerr << "None of the loop folder possibilities:" << endl;
    for ( size_t iName = 0; iName < loopFolderNames.size() ; ++iName ) {
      cerr << "\t" << loopFolderNames.at( iName ) << endl;
    }
    cerr << "was found in " << filename << endl;
    return;
  }

  // guess the number of sensors from the number of subfolder in the loopfolder
  UInt_t nDetector = loopFolder->GetListOfKeys()->GetSize();
  vector< string > detectorFolderNames;
  vector< int >    sensorIDs;
  string separator = "-_";

  for ( size_t i = 0 ; i < nDetector; ++i ) {
    string name( loopFolder->GetListOfKeys()->At( i )->GetName() );
    detectorFolderNames.push_back( name );
    sensorIDs.push_back( atoi(name.substr( name.find_last_of( separator ) + 1,  name.length() ).c_str()) );
  }
  UInt_t nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  vector<TCanvas * > canvasVec;
  vector<TPad * >    padVec;

  Double_t titleHeight = 0.10;
  Int_t canvasWidth  = 800;
  Int_t canvasHeight = 800;

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string pedeCanvasName  = string(pedeCanvasBaseName.Data()) + "_" + toString(iCanvas);
    string pedeCanvasTitle = string(runName) + " - Pedestal histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * pedeCanvas = new TCanvas( pedeCanvasName.c_str(), pedeCanvasTitle.c_str(), canvasWidth, canvasHeight);
    pedeCanvas->Range(0,0,1,1);
    pedeCanvas->SetBorderSize(0);
    pedeCanvas->SetFrameFillColor(0);
    pedeCanvas->SetBorderMode(0);
    canvasVec.push_back( pedeCanvas );

    // title pad
    TPad * titlePad = new TPad("pedeTitle","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( pedeCanvasTitle.c_str() );
    title->Draw();
    pedeCanvas->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() <  2 * nDetector ) {
        padVec.push_back( smallPad );
      }
    }
  }

  gStyle->SetTitleFontSize( 0.05 );
  gStyle->SetTitleFillColor( kCyan - 9 );
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.35);
  gStyle->SetPalette( 1, 0 );
  gStyle->SetOptStat("emr");

  // now plot the histograms in the right pad, again starting from the pedestal
  UInt_t iPad = 0;
  string newTitle;

  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) loopFolder->Get( detectorFolderNames.at(iDetector).c_str() );
    vector< string > mapNames;
    mapNames.push_back( "PedeMap-d" + toString( sensorIDs.at( iDetector ) ) + "-l" + toString( loop ) );
    mapNames.push_back( "PedeMap_d" + toString( sensorIDs.at( iDetector ) ) + "_l" + toString( loop ) );

    TH2D * map = NULL;
    for ( size_t iName = 0 ; iName < mapNames.size(); ++iName ) {
      map = (TH2D* ) detectorFolder->Get( mapNames.at( iName ).c_str() );
      if ( map != NULL ) break;
    }
    if ( map == NULL ) {
      cerr << "None of the map histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < mapNames.size() ; ++iName ) {
        cerr << "\t" << mapNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string( map->GetTitle() ) + " - Detector " + toString( sensorIDs.at( iDetector ) );
    map->SetTitle( newTitle.c_str() );
    map->SetXTitle("x [pixel]");
    map->SetYTitle("y [pixel]");
    map->SetStats( false );
    padVec[iPad++]->cd();
    map->SetContour(99);
    map->Draw("colz");


    vector< string > histoNames;
    histoNames.push_back( "PedeDist-d" + toString( sensorIDs.at( iDetector ) ) + "-l" + toString( loop ) );
    histoNames.push_back( "PedeDist_d" + toString( sensorIDs.at( iDetector ) ) + "_l" + toString( loop ) );

    TH1D * histo = NULL;
    for ( size_t iName = 0 ; iName < histoNames.size(); ++iName ) {
      histo  = (TH1D*) detectorFolder->Get( histoNames.at( iName ).c_str() );
      if ( histo != NULL ) break;
    }
    if ( histo == NULL ) {
      cerr << "None of the histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < histoNames.size() ; ++iName ) {
        cerr << "\t" << histoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string( histo->GetTitle() ) + " - Detector " + toString( sensorIDs.at( iDetector ) );
    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle("Pedestal [ADC]");
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad++]->cd();
    histo->Draw();
  }

  // now let's move to the noise histos and map. Again nDetPerCanvas = 3;
  nDetPerCanvas = 3;
  nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  // reset the padVec
  padVec.clear();

  TString noiseCanvasBaseName   = "NoiseCanvas";
  closeCanvases( noiseCanvasBaseName );

  // prepare the new canvases
  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string noiseCanvasName  = string(noiseCanvasBaseName.Data()) + "_" + toString(iCanvas);
    string noiseCanvasTitle = string(runName) + " - Noise histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * noiseCanvas = new TCanvas( noiseCanvasName.c_str(), noiseCanvasTitle.c_str(), canvasWidth, canvasHeight);
    noiseCanvas->Range(0,0,1,1);
    noiseCanvas->SetBorderSize(0);
    noiseCanvas->SetFrameFillColor(0);
    noiseCanvas->SetBorderMode(0);
    canvasVec.push_back( noiseCanvas );

    // title pad
    TPad * titlePad = new TPad("noiseTitle","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( noiseCanvasTitle.c_str() );
    title->Draw();
    noiseCanvas->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() <  2 * nDetector ) {
        padVec.push_back( smallPad );
      }
    }
  }

  iPad = 0;
  // now the same as above for the noise
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) loopFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > mapNames;
    mapNames.push_back(  "NoiseMap-d" + toString( sensorIDs.at( iDetector ) ) + "-l" + toString( loop ) );
    mapNames.push_back(  "NoiseMap_d" + toString( sensorIDs.at( iDetector ) ) + "_l" + toString( loop ) );
    TH2D *  map = NULL;
    for ( size_t iName = 0 ; iName < mapNames.size(); ++iName ) {
      map = (TH2D* ) detectorFolder->Get( mapNames.at( iName) .c_str() );
      if ( map != NULL ) break;
    }
    if ( map == NULL ) {
      cerr << "None of the map histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < mapNames.size() ; ++iName ) {
        cerr << "\t" << mapNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }


    newTitle = string( map->GetTitle() ) + " - Detector " + toString( sensorIDs.at( iDetector ) );
    map->SetTitle( newTitle.c_str() );
    map->SetXTitle("x [pixel]");
    map->SetYTitle("y [pixel]");
    map->SetStats( false );
    padVec[iPad++]->cd();
    map->SetContour(99);
    map->Draw("colz");

    vector< string > histoNames;
    histoNames.push_back( "NoiseDist-d" + toString( sensorIDs.at( iDetector ) ) + "-l" + toString( loop ) );
    histoNames.push_back( "NoiseDist_d" + toString( sensorIDs.at( iDetector ) ) + "_l" + toString( loop ) );

    TH1D * histo = NULL;
    for ( size_t iName = 0 ; iName < histoNames.size(); ++iName ) {
      histo  = (TH1D*) detectorFolder->Get( histoNames.at( iName ).c_str() );
      if ( histo != NULL ) break;
    }
    if ( histo == NULL ) {
      cerr << "None of the histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < histoNames.size() ; ++iName ) {
        cerr << "\t" << histoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string( histo->GetTitle() ) + " - Detector " + toString( sensorIDs.at( iDetector ) );
    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle("Noise [ADC]");
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad++]->cd();
    histo->Draw();
  }


  // now the status canvas. This has nDetPerCanvas = 6;
  nDetPerCanvas = 6;
  nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }
  padVec.clear();
  TString statusCanvasBaseName  = "StatusCanvas";
  closeCanvases( statusCanvasBaseName );

  // prepare all the needed canvases for the status
  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {
    string statusCanvasName  = string( statusCanvasBaseName.Data() ) + "_" + toString( iCanvas ) ;
    string statusCanvasTitle = string( runName ) + " - Status map " + toString (iCanvas + 1) + " / " + toString( nCanvas );

    TCanvas * statusCanvas  = new TCanvas( statusCanvasName.c_str(), statusCanvasTitle.c_str(), canvasWidth, canvasHeight);
    statusCanvas->Range(0,0,1,1);
    statusCanvas->SetBorderSize(0);
    statusCanvas->SetFrameFillColor(0);
    statusCanvas->SetBorderMode(0);
    canvasVec.push_back( statusCanvas );

    // title pad
    TPad * titlePad = new TPad("pedeTitle","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( statusCanvasTitle.c_str() );
    title->Draw();
    statusCanvas->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() <  nDetector ) {
        padVec.push_back( smallPad );
      }
    }
  }

  iPad = 0;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) loopFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > mapNames;
    mapNames.push_back( "StatusMap-d" + toString( sensorIDs.at( iDetector ) ) + "-l" + toString( loop ) );
    mapNames.push_back( "StatusMap_d" + toString( sensorIDs.at( iDetector ) ) + "_l" + toString( loop ) );

    TH2D * map = NULL;
    for ( size_t iName = 0 ; iName < mapNames.size(); ++iName ) {
      map = (TH2D* ) detectorFolder->Get( mapNames.at( iName) .c_str() );
      if ( map != NULL ) break;
    }
    if ( map == NULL ) {
      cerr << "None of the map histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < mapNames.size() ; ++iName ) {
        cerr << "\t" << mapNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string( map->GetTitle() ) + " - Detector " + toString( sensorIDs.at( iDetector ) );
    map->SetTitle( newTitle.c_str() );
    map->SetXTitle("x [pixel]");
    map->SetYTitle("y [pixel]");
    map->SetStats( false );
    padVec[iPad++]->cd();
    map->SetContour(99);
    map->Draw("colz");

  }

  string detectorString( detector );
  if ( detectorString == "mimotel" ) {

    // a few global constants
    // const UInt_t kXPixel = 256;
    const UInt_t kYPixel = 256;

    const UInt_t kNChan  = 4;
    UInt_t xLimit[ kNChan + 1  ] = { 0, 64, 128 , 192, 256 };
    Double_t channel[ kNChan ] ;

    Double_t * meanNoise = new Double_t[ kNChan * nDetector ];
    Double_t * rmsNoise = new Double_t[ kNChan * nDetector ];
    const UInt_t kNColor = 10;

    int kColor[kNColor] = {
      kRed,
      kBlue,
      kBlack,
      kGreen,
      kTeal,
      kOrange,
      kPink
    };

    // in case the detector is a mimotel then we can make also the
    // channel by channel analysis.
    nDetPerCanvas =  3;
    TString canvasBaseName = "AdvancedNoiseCanvas";

    // close all canvases with these names;
    closeCanvases( canvasBaseName );

    nCanvas   = nDetector / nDetPerCanvas;
    if ( nDetector % nDetPerCanvas != 0 ) {
      ++nCanvas;
    }

    padVec.clear();

    for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

      string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
      string canvasTitle = string(runName) + " - Advanced noise histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

      TCanvas * canvas = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
      canvas->Range(0,0,1,1);
      canvas->SetBorderSize(0);
      canvas->SetFrameFillColor(0);
      canvas->SetBorderMode(0);
      canvasVec.push_back( canvas );

      // title pad
      TPad * titlePad = new TPad("pedeTitle","title",0, 1 - titleHeight,1,1);
      titlePad->Draw();
      titlePad->SetBorderMode(0);
      titlePad->SetBorderSize(0);
      titlePad->SetFrameFillColor(0);
      titlePad->cd();
      TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
      title->SetBorderSize(1);
      title->SetLabel( canvasTitle.c_str() );
      title->Draw();
      canvas->cd();

      // big pad for the rest
      TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
      bigPad->Draw();
      bigPad->cd();
      bigPad->SetBorderMode(0);
      bigPad->SetBorderSize(0);
      bigPad->SetFrameFillColor(0);

      // divide the bigPad in 2 x 3 TPad and add them to the subpad list
      Int_t nX = 4, nY = 3;
      bigPad->Divide(nX, nY);

      for ( Int_t i = 0; i < nX * nY; i++ ) {
        TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
        smallPad->SetBorderMode(0);
        smallPad->SetBorderSize(0);
        smallPad->SetFrameFillColor(0);
        if ( padVec.size() <  4 * nDetector ) {
          padVec.push_back( smallPad );
        }
      }
    }

    iPad = 0;
    for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

      TDirectoryFile * detectorFolder = (TDirectoryFile*) loopFolder->Get( detectorFolderNames.at(iDetector).c_str() );

      vector< string > mapNames;
      mapNames.push_back(  "NoiseMap-d" + toString( sensorIDs.at( iDetector ) ) + "-l" + toString( loop ) );
      mapNames.push_back(  "NoiseMap_d" + toString( sensorIDs.at( iDetector ) ) + "_l" + toString( loop ) );
      TH2D *  map = NULL;
      for ( size_t iName = 0 ; iName < mapNames.size(); ++iName ) {
        map = (TH2D* ) detectorFolder->Get( mapNames.at( iName) .c_str() );
        if ( map != NULL ) break;
      }
      if ( map == NULL ) {
        cerr << "None of the map histo name possibilities:" << endl;
        for ( size_t iName = 0; iName < mapNames.size() ; ++iName ) {
          cerr << "\t" << mapNames.at( iName ) << endl;
        }
        cerr << "was found in " << filename << endl;
        return;
      }

      for ( UInt_t iChan = 0 ; iChan < kNChan ; ++iChan ) {
        if ( iDetector == 0 ) {
          channel[iChan] = iChan - 0.5;
        }

        string tempName  = "Noise_Dist_d" + toString( sensorIDs.at( iDetector ) ) + "_ch" + toString( iChan ) ;
        string tempTitle = "Noise Det. " + toString( sensorIDs.at( iDetector ) ) + " - Ch. " + toString( iChan ) ;

        TH1D * noiseDistCh = new TH1D( tempName.c_str(), tempTitle.c_str(), 50, 0., 10. );
        noiseDistCh->SetXTitle( "Noise [ADC}");
        noiseDistCh->SetFillColor( kColor[iDetector % kNColor]  );

        // let's start looping on pixels now
        for ( size_t yPixel = 1 ; yPixel <= kYPixel ; ++yPixel ) {
          for ( size_t xPixel = xLimit[ iChan ] + 1; xPixel <= xLimit[ iChan +1 ] ; ++xPixel ) {
            double noise = map->GetBinContent( xPixel , yPixel );
            noiseDistCh->Fill( noise );
          }
        }
        padVec[iPad]->cd();
        noiseDistCh->Draw();
        padVec[iPad]->Update();
        TPaveText * padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
        padTitle->SetY1NDC(0.90);
        padTitle->SetX2NDC(0.80);
        padVec[iPad]->Modified( true );
        meanNoise[kNChan * iDetector + iChan] = noiseDistCh->GetMean() ;
        rmsNoise[kNChan * iDetector + iChan]  = noiseDistCh->GetRMS() ;

        padVec[iPad]->Update();
        TPaveStats * st = (TPaveStats*) noiseDistCh->GetListOfFunctions()->FindObject("stats");
        st->SetX1NDC( 0.55 );
        st->SetX2NDC( 0.98 );
        st->SetY1NDC( 0.60 );
        st->SetY2NDC( 0.85 );
        padVec[iPad]->Modified();
        ++iPad;
      }
    }


    // create one more canvas for the comparison among different
    // channels on different sensors
    padVec.clear();

    string canvasName  = "NoiseComparison";
    string canvasTitle = "Noise comparison";
    Int_t  canvasWidth1  = 1000;
    Int_t  canvasHeight1 =  400;


    TCanvas * canvas = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth1, canvasHeight1);
    canvas->Range(0,0,1,1);
    canvas->SetBorderSize(0);
    canvas->SetFrameFillColor(0);
    canvas->SetBorderMode(0);
    canvasVec.push_back( canvas );

    // title pad
    TPad * titlePad = new TPad("pedeTitle","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    canvas->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    padVec.push_back( bigPad );

    iPad = 0;

    TMultiGraph * mg     = new TMultiGraph();
    TLegend     * legend = new TLegend(0.6358, 0.6207, 0.7956, 0.9705);
    mg->SetTitle("Channel by channel noise comparison");
    Double_t * x = new Double_t[ kNChan ];
    Double_t * errX = new Double_t[ kNChan ];
    for ( UInt_t iDetector = 0; iDetector < nDetector; ++iDetector ) {
      if ( iDetector == 0 ) {
        for ( UInt_t iChan = 0; iChan < kNChan ; ++iChan ) {
          x[iChan] =  +0.5 + iChan  ;
          errX[iChan] = 0.5;
        }
      }

      TGraphErrors * gr = new TGraphErrors( kNChan, x, &meanNoise[ iDetector * kNChan ], NULL, &rmsNoise[ iDetector * kNChan ] );
      string grName   = "Detector_d" + toString( sensorIDs.at( iDetector ) );
      string grTitle  = "Detector " + toString( sensorIDs.at( iDetector ) );
      gr->GetHistogram()->SetBit( TH1::kNoTitle, true ) ;
      gr->SetName( grName.c_str() );
      gr->SetTitle( grTitle.c_str() );
      gr->SetMarkerColor( kColor[ iDetector % kNColor ] );
      gr->SetMarkerStyle( 20 + ( (iDetector + 1 ) % kNColor ) );
      gr->GetXaxis()->SetNdivisions( 0 );
      gr->GetXaxis()->SetRangeUser(-0.5, kNChan + 0.5 );
      gr->GetYaxis()->SetTitle("Noise [ADC]");
      legend->AddEntry( gr, gr->GetTitle(), "P" );
      if ( iDetector == 0 )  mg->Add( gr, "AP1" );
      else mg->Add( gr, "P1");
      //gr->Draw("ALP");

    }
    padVec[iPad]->cd();
    mg->SetBit( TH1::kNoTitle, true ) ;
    mg->Draw(  );
    legend->Draw();
    Double_t step = 1.0 / kNChan;
    Double_t xPos[kNChan] = { 0.17, 0.39, 0.63, 0.83 } ;
    TPaveText * padTitle = new TPaveText(0.073, 0.8539, 0.599, 0.984,"NDC");
    padTitle->AddText("Noise comparison channel per channel");
    padTitle->Draw();
    for ( UInt_t iChan = 0 ; iChan < kNChan ; ++iChan ) {
      //
      if ( iChan != 0 ) {
        TLine * line = new TLine ;
        TLine * line2 = line->DrawLineNDC(iChan * step ,0.10,iChan * step,0.90);
        line2->SetLineStyle( 2 );
      }

      TPaveLabel * label = new TPaveLabel;
      label->SetX1NDC( xPos[iChan] - 0.03 );
      label->SetX2NDC( xPos[iChan] + 0.03 );
      label->SetY1NDC( 0.01 );
      label->SetY2NDC( 0.07);
      label->SetOption("");
      string l = "Ch. " + toString( iChan );
      label->SetLabel( l.c_str() );
      label->SetTextSize( 1.0) ;
      label->Draw();
    }

  }


  string path( prepareOutputFolder( "PedeNoise" ));
  for ( UInt_t iCanvas = 0 ; iCanvas < canvasVec.size(); ++iCanvas ) {
    string figName = path + canvasVec[iCanvas]->GetName() + pictureOutputFormat.Data();
    canvasVec[iCanvas]->SaveAs( figName.c_str() );
  }

}

void showCorrelationPlot( const char * filename ) {

  // check if this file is already open
  TFile * inputFile = closeAndReopenFile( filename );

  if ( inputFile == 0x0 ) {
    cerr << "Problems opening file " << filename << endl;
    return;
  }

  // typical names for the processor folder
  vector< string > folderNames;
  folderNames.push_back( correlatorFolderName.Data() );
  folderNames.push_back( "Correlation" ) ;
  folderNames.push_back( "CorrelationAfterFilter" ) ;
  folderNames.push_back( "CorrelationBeforeFilter" ) ;

  string dutFolderName = correlatorFolderName.Data() + toString( "_dut" );
  folderNames.push_back( dutFolderName );
  dutFolderName = correlatorFolderName.Data() + toString( "_dep" );
  folderNames.push_back( dutFolderName );
  dutFolderName = correlatorFolderName.Data() + toString( "_taki" );
  folderNames.push_back( dutFolderName );

  // this function will create the following canvases:
  //
  // --> ClusterXCanvas
  // --> ClusterYCanvas
  // --> HiXCanvas
  // --> HitYCanvas

  UInt_t nDetPerCanvas = 3;
  UInt_t nDetector = 0;
  UInt_t nCanvas = 0;
  UInt_t iPad;
  vector<TCanvas * > canvasVec;
  vector<TPad *    > padVec;

  // look into the input file for a folder named correlator
  TDirectoryFile * correlatorFolder = checkFolder( folderNames, inputFile ) ;
  if ( correlatorFolder == 0x0 ) {
    return;
  }

  // this folder can contain up to 4 subfolder, one for x cluster, one
  // for y cluster, one for x hit and last one for y hit.
  TDirectoryFile * clusterXFolder = (TDirectoryFile*) correlatorFolder->Get("ClusterX");
  if ( clusterXFolder != 0x0 ) {

    TString canvasBaseName    = "ClusterXCanvas";
    closeCanvases( canvasBaseName );

    // ok there is a cluster x folder, let's do it!
    padVec.clear();

    UInt_t nHisto = clusterXFolder->GetListOfKeys()->GetSize();
    set< int > sensorIDSet;
    string separator = "-_";

    for ( size_t i = 0 ; i < nHisto; ++i ) {
      string name  = clusterXFolder->GetListOfKeys()->At( i )->GetName();
      size_t last  = name.find_last_not_of( separator );
      size_t first = name.find_last_not_of( separator, last );
      sensorIDSet.insert( atoi( name.substr( first, last ).c_str() ) );
    }

    nDetector = sensorIDSet.size();
    nCanvas = nDetector / nDetPerCanvas;
    if ( nDetector % nDetPerCanvas != 0 ) {
      ++nCanvas;
    }

    Double_t titleHeight = 0.10;
    Int_t canvasWidth  = 800;
    Int_t canvasHeight = 800;
    for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

      string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
      string canvasTitle = string(runName) + " - Cluster correlation (x) " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

      TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
      c->Range(0,0,1,1);
      c->SetBorderSize(0);
      c->SetFrameFillColor(0);
      c->SetBorderMode(0);
      canvasVec.push_back( c );

      // title pad
      TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
      titlePad->Draw();
      titlePad->SetBorderMode(0);
      titlePad->SetBorderSize(0);
      titlePad->SetFrameFillColor(0);
      titlePad->cd();
      TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
      title->SetBorderSize(1);
      title->SetLabel( canvasTitle.c_str() );
      title->Draw();
      c->cd();

      // big pad for the rest
      TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
      bigPad->Draw();
      bigPad->cd();
      bigPad->SetBorderMode(0);
      bigPad->SetBorderSize(0);
      bigPad->SetFrameFillColor(0);

      // divide the bigPad in 2 x 3 TPad and add them to the subpad list
      Int_t nX = 2, nY = 3;
      bigPad->Divide(nX, nY);


      for ( Int_t i = 0; i < nX * nY; i++ ) {
        TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
        smallPad->SetBorderMode(0);
        smallPad->SetBorderSize(0);
        smallPad->SetFrameFillColor(0);
        if ( padVec.size() < 2 * (nDetector - 1) ) {
          padVec.push_back( smallPad);
        }
      }
    }

    gStyle->SetPalette(1,0);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(11);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.35);
    gStyle->SetTitleFontSize( 0.05 );
    gStyle->SetTitleFillColor( kCyan - 9 );

    iPad = 0;
    string newTitle;
    set< int >::iterator iter = sensorIDSet.begin();
    for ( size_t i = 0 ; i < sensorIDSet.size() - 1 ; ++i, ++iter ) {
      int iDetector = *iter;

      string histoName = "ClusterXCorrelationHisto_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TH2D * histo = (TH2D*) clusterXFolder->Get(histoName.c_str());

      newTitle = "Cluster X: " + toString( iDetector ) + " vs " + toString( iDetector + 1 );
      histo->SetTitle( newTitle.c_str() );

      newTitle = "Detector " + toString( iDetector ) + " [pixel] ";
      histo->SetXTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      histo->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      histo->Draw("colz");
      padVec[iPad]->Update();
      TPaveText * padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      padVec[iPad]->Modified( true );

      ++iPad;

      // find xMin, xMax, yMin, yMax
      Double_t xMin = histo->GetXaxis()->GetXmin();
      Double_t xMax = histo->GetXaxis()->GetXmax();
      Double_t yMin = histo->GetYaxis()->GetXmin();
      Double_t yMax = histo->GetYaxis()->GetXmax();

      // now make a slice of the histo in a 1 / 3 of the height
      Int_t yBin = histo->GetYaxis()->FindBin( ( yMax - yMin ) / 3. ) ;
      cout << "yBin =" << yBin  << endl;
      string profName = string(histo->GetName()) + "_p1";
      TH1D * profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      Int_t xBin = profile->GetMaximumBin();

      Double_t xA = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yA = (yMax - yMin) / 3;

      // now at 2 / 3
      yBin = histo->GetYaxis()->FindBin( 2 * ( yMax - yMin ) / 3. ) ;
      profName = string(histo->GetName()) + "_p2";
      profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      xBin = profile->GetMaximumBin();

      Double_t xB = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yB = 2 * ( yMax - yMin ) / 3.;

      // cout << "A = (" << xA << ", " << yA << ") and B = (" << xB << ", " << yB << ")" << endl;

      // this is the candidate for the correlation line
      Double_t p1 = (yA - yB) / (xA - xB);
      Double_t p0 = yA - p1 * xA;

      // shift some pixels above
      Double_t shiftUp = 15;
      Double_t p1_up = p1;
      Double_t p0_up = p0 + shiftUp;

      // Point C:
      Double_t xC, yC = 0;
      xC = ( yMin - p0_up) / p1_up;
      if ( xC < xMin ) {
        xC = xMin;
        yC = p1_up * xMin + p0_up;
      } else if ( xC > xMin && xC <= xMax ) {
        yC = yMin;
      } else if ( xC > xMax ) {
        xC = xMax;
        yC = p1_up * xMax + p0_up;
      }

      // Point D
      Double_t xD, yD = 0;
      xD = (yMax - p0_up) / p1_up;
      if ( xD > xMax ) {
        xD = xMax;
        yD = p1_up * xMax + p0_up;
      } else if ( xD > xMin && xD <= xMax ) {
        yD = yMax;
      } else if ( xD < xMin ) {
        xD = xMin;
        yD = p1_up * xMin + p0_up;
      }

      // cout << "C = (" << xC << ", " << yC << ") and D = (" << xD << ", " << yD << ")" << endl;

      // shift some pixels above
      Double_t shiftDown = 15;
      Double_t p1_down = p1;
      Double_t p0_down = p0 - shiftDown;

      // Point E
      Double_t xE, yE = 0;
      xE = ( yMax - p0_down ) / p1_down ;
      if ( xE > xMax ) {
        xE = xMax;
        yE = p1_down * xMax + p0_down;
      } else if ( xE > xMin && xE < xMax ) {
        yE = yMax;
      } else if ( xE < xMin ) {
        xE = xMin;
        yE = p1_down * xMin + p0_down;
      }

      // Point F
      Double_t xF, yF = 0;
      xF = ( yMin - p0_down ) / p1_down;
      if ( xF > xMax) {
        xF = xMax;
        yF = p1_down * xMax + p0_down;
      } else if ( xF > xMin && xF <= xMax ) {
        yF = yMin;
      } else if ( xF < xMin ) {
        xF = xMin;
        yF = p1_down * xMin + p0_down;
      }

      // cout << "E = (" << xE << ", " << yE << ") and F = (" << xF << ", " << yF << ")" << endl;

      string cutName = "xCorrDiagonal_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TCutG * diagonal = new TCutG(cutName.c_str(),5);
      diagonal->SetPoint(0, xC, yC );
      diagonal->SetPoint(1, xD, yD );
      diagonal->SetPoint(2, xE, yE );
      diagonal->SetPoint(3, xF, yF );
      diagonal->SetPoint(4, xC, yC );
      diagonal->SetLineColor( kRed );
      diagonal->SetLineStyle( 2 );
      diagonal->Draw("LP");


      string correlationName = "xCorrClu_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      string option = "[" + cutName + "]";
      TProfile * correlation = histo->ProfileX(correlationName.c_str(), 1, -1, option.c_str());
      newTitle = "Fit (Det "  + toString( iDetector ) + " vs Det" + toString( iDetector + 1 ) + ")";
      correlation->SetTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      correlation->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      correlation->Fit("pol1");
      TF1 * fitFunc = (TF1*) correlation->GetFunction( "pol1" );
      fitFunc->SetLineColor( kBlue ) ;
      fitFunc->SetLineStyle( 2 );
      padVec[iPad]->Update();
      padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      TPaveStats * stats = (TPaveStats*) correlation->GetListOfFunctions()->FindObject("stats");
      stats->SetX1NDC( 0.62 );
      stats->SetY1NDC( 0.65 );

      padVec[iPad]->Modified( true );

      ++iPad;
    }
  }


  TDirectoryFile * clusterYFolder = (TDirectoryFile*) correlatorFolder->Get("ClusterY");
  if ( clusterYFolder != 0x0 ) {

    TString canvasBaseName    = "ClusterYCanvas";
    closeCanvases( canvasBaseName );

    // ok there is a cluster x folder, let's do it!
    padVec.clear();

    UInt_t nHisto = clusterYFolder->GetListOfKeys()->GetSize();
    set< int > sensorIDSet;
    string separator = "-_";

    for ( size_t i = 0 ; i < nHisto; ++i ) {
      string name  = clusterYFolder->GetListOfKeys()->At( i )->GetName();
      size_t last  = name.find_last_not_of( separator );
      size_t first = name.find_last_not_of( separator, last );
      sensorIDSet.insert( atoi( name.substr( first, last ).c_str() ) );
    }

    nDetector = sensorIDSet.size();
    nCanvas = nDetector / nDetPerCanvas;
    if ( nDetector % nDetPerCanvas != 0 ) {
      ++nCanvas;
    }

    Double_t titleHeight = 0.10;
    Int_t canvasWidth  = 800;
    Int_t canvasHeight = 800;
    for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

      string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
      string canvasTitle = string(runName) + " - Cluster correlation (y) " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

      TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
      c->Range(0,0,1,1);
      c->SetBorderSize(0);
      c->SetFrameFillColor(0);
      c->SetBorderMode(0);
      canvasVec.push_back( c );

      // title pad
      TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
      titlePad->Draw();
      titlePad->SetBorderMode(0);
      titlePad->SetBorderSize(0);
      titlePad->SetFrameFillColor(0);
      titlePad->cd();
      TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
      title->SetBorderSize(1);
      title->SetLabel( canvasTitle.c_str() );
      title->Draw();
      c->cd();

      // big pad for the rest
      TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
      bigPad->Draw();
      bigPad->cd();
      bigPad->SetBorderMode(0);
      bigPad->SetBorderSize(0);
      bigPad->SetFrameFillColor(0);

      // divide the bigPad in 2 x 3 TPad and add them to the subpad list
      Int_t nX = 2, nY = 3;
      bigPad->Divide(nX, nY);


      for ( Int_t i = 0; i < nX * nY; i++ ) {
        TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
        smallPad->SetBorderMode(0);
        smallPad->SetBorderSize(0);
        smallPad->SetFrameFillColor(0);
        if ( padVec.size() < 2 * (nDetector - 1) ) {
          padVec.push_back( smallPad);
        }
      }
    }

    gStyle->SetPalette(1,0);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(11);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.35);
    gStyle->SetTitleFontSize( 0.05 );
    gStyle->SetTitleFillColor( kCyan - 9 );

    iPad = 0;
    string newTitle;
    set< int >::iterator iter = sensorIDSet.begin();
    for ( size_t i = 0 ; i < sensorIDSet.size() - 1 ; ++i, ++iter ) {
      int iDetector = *iter;
      string histoName = "ClusterYCorrelationHisto_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TH2D * histo = (TH2D*) clusterYFolder->Get(histoName.c_str());

      newTitle = "Cluster Y: " + toString( iDetector ) + " vs " + toString( iDetector + 1 );
      histo->SetTitle( newTitle.c_str() );

      newTitle = "Detector " + toString( iDetector ) + " [pixel] ";
      histo->SetXTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      histo->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      histo->Draw("colz");
      padVec[iPad]->Update();
      TPaveText * padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      padVec[iPad]->Modified( true );

      ++iPad;

      // find xMin, xMax, yMin, yMax
      Double_t xMin = histo->GetXaxis()->GetXmin();
      Double_t xMax = histo->GetXaxis()->GetXmax();
      Double_t yMin = histo->GetYaxis()->GetXmin();
      Double_t yMax = histo->GetYaxis()->GetXmax();

      // now make a slice of the histo in a 1 / 3 of the height
      Int_t yBin = histo->GetYaxis()->FindBin( ( yMax - yMin ) / 3. ) ;
      string profName = string(histo->GetName()) + "_p1";
      TH1D * profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      Int_t xBin = profile->GetMaximumBin();

      Double_t xA = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yA = (yMax - yMin) / 3;

      // now at 2 / 3
      yBin = histo->GetYaxis()->FindBin( 2 * ( yMax - yMin ) / 3. ) ;
      profName = string(histo->GetName()) + "_p2";
      profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      xBin = profile->GetMaximumBin();

      Double_t xB = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yB = 2 * ( yMax - yMin ) / 3.;

      // cout << "A = (" << xA << ", " << yA << ") and B = (" << xB << ", " << yB << ")" << endl;

      // this is the candidate for the correlation line
      Double_t p1 = (yA - yB) / (xA - xB);
      Double_t p0 = yA - p1 * xA;

      // shift some pixels above
      Double_t shiftUp = 15;
      Double_t p1_up = p1;
      Double_t p0_up = p0 + shiftUp;

      // Point C:
      Double_t xC, yC = 0;
      xC = ( yMin - p0_up) / p1_up;
      if ( xC < xMin ) {
        xC = xMin;
        yC = p1_up * xMin + p0_up;
      } else if ( xC > xMin && xC <= xMax ) {
        yC = yMin;
      } else if ( xC > xMax ) {
        xC = xMax;
        yC = p1_up * xMax + p0_up;
      }

      // Point D
      Double_t xD, yD = 0;
      xD = (yMax - p0_up) / p1_up;
      if ( xD > xMax ) {
        xD = xMax;
        yD = p1_up * xMax + p0_up;
      } else if ( xD > xMin && xD <= xMax ) {
        yD = yMax;
      } else if ( xD < xMin ) {
        xD = xMin;
        yD = p1_up * xMin + p0_up;
      }

      // cout << "C = (" << xC << ", " << yC << ") and D = (" << xD << ", " << yD << ")" << endl;

      // shift some pixels above
      Double_t shiftDown = 15;
      Double_t p1_down = p1;
      Double_t p0_down = p0 - shiftDown;

      // Point E
      Double_t xE, yE = 0;
      xE = ( yMax - p0_down ) / p1_down ;
      if ( xE > xMax ) {
        xE = xMax;
        yE = p1_down * xMax + p0_down;
      } else if ( xE > xMin && xE < xMax ) {
        yE = yMax;
      } else if ( xE < xMin ) {
        xE = xMin;
        yE = p1_down * xMin + p0_down;
      }

      // Point F
      Double_t xF, yF = 0;
      xF = ( yMin - p0_down ) / p1_down;
      if ( xF > xMax) {
        xF = xMax;
        yF = p1_down * xMax + p0_down;
      } else if ( xF > xMin && xF <= xMax ) {
        yF = yMin;
      } else if ( xF < xMin ) {
        xF = xMin;
        yF = p1_down * xMin + p0_down;
      }

      // cout << "E = (" << xE << ", " << yE << ") and F = (" << xF << ", " << yF << ")" << endl;

      string cutName = "yCorrDiagonal_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TCutG * diagonal = new TCutG(cutName.c_str(),5);
      diagonal->SetPoint(0, xC, yC );
      diagonal->SetPoint(1, xD, yD );
      diagonal->SetPoint(2, xE, yE );
      diagonal->SetPoint(3, xF, yF );
      diagonal->SetPoint(4, xC, yC );
      diagonal->SetLineColor( kRed );
      diagonal->SetLineStyle( 2 );
      diagonal->Draw("LP");


      string correlationName = "yCorrClu_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      string option = "[" + cutName + "]";
      TProfile * correlation = histo->ProfileX(correlationName.c_str(), 1, -1, option.c_str());
      newTitle = "Fit (Det "  + toString( iDetector ) + " vs Det" + toString( iDetector + 1 ) + ")";
      correlation->SetTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      correlation->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      correlation->Fit("pol1");
      TF1 * fitFunc = (TF1*) correlation->GetFunction( "pol1" );
      fitFunc->SetLineColor( kBlue ) ;
      fitFunc->SetLineStyle( 2 );
      padVec[iPad]->Update();
      padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      TPaveStats * stats = (TPaveStats*) correlation->GetListOfFunctions()->FindObject("stats");
      stats->SetX1NDC( 0.62 );
      stats->SetY1NDC( 0.65 );

      padVec[iPad]->Modified( true );

      ++iPad;
    }
  }

  TDirectoryFile * hitXFolder = (TDirectoryFile*) correlatorFolder->Get("HitX");
  if ( hitXFolder != 0x0 ) {

    TString canvasBaseName    = "HitXCanvas";
    closeCanvases( canvasBaseName );

    // ok there is a cluster x folder, let's do it!
    padVec.clear();

    UInt_t nHisto = hitXFolder->GetListOfKeys()->GetSize();
    set< int > sensorIDSet;
    string separator = "-_";

    for ( size_t i = 0 ; i < nHisto; ++i ) {
      string name  = hitXFolder->GetListOfKeys()->At( i )->GetName();
      size_t last  = name.find_last_not_of( separator );
      size_t first = name.find_last_not_of( separator, last );
      sensorIDSet.insert( atoi( name.substr( first, last ).c_str() ) );
    }

    nDetector = sensorIDSet.size();
    nCanvas = nDetector / nDetPerCanvas;
    if ( nDetector % nDetPerCanvas != 0 ) {
      ++nCanvas;
    }


    Double_t titleHeight = 0.10;
    Int_t canvasWidth  = 800;
    Int_t canvasHeight = 800;
    for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

      string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
      string canvasTitle = string(runName) + " - Hit correlation (x) " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

      TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
      c->Range(0,0,1,1);
      c->SetBorderSize(0);
      c->SetFrameFillColor(0);
      c->SetBorderMode(0);
      canvasVec.push_back( c );

      // title pad
      TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
      titlePad->Draw();
      titlePad->SetBorderMode(0);
      titlePad->SetBorderSize(0);
      titlePad->SetFrameFillColor(0);
      titlePad->cd();
      TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
      title->SetBorderSize(1);
      title->SetLabel( canvasTitle.c_str() );
      title->Draw();
      c->cd();

      // big pad for the rest
      TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
      bigPad->Draw();
      bigPad->cd();
      bigPad->SetBorderMode(0);
      bigPad->SetBorderSize(0);
      bigPad->SetFrameFillColor(0);

      // divide the bigPad in 2 x 3 TPad and add them to the subpad list
      Int_t nX = 2, nY = 3;
      bigPad->Divide(nX, nY);


      for ( Int_t i = 0; i < nX * nY; i++ ) {
        TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
        smallPad->SetBorderMode(0);
        smallPad->SetBorderSize(0);
        smallPad->SetFrameFillColor(0);
        if ( padVec.size() < 2 * (nDetector - 1) ) {
          padVec.push_back( smallPad);
        }
      }
    }

    gStyle->SetPalette(1,0);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(11);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.35);
    gStyle->SetTitleFontSize( 0.05 );
    gStyle->SetTitleFillColor( kCyan - 9 );

    iPad = 0;
    string newTitle;
    set< int >::iterator iter = sensorIDSet.begin();
    for ( size_t i = 0 ; i < sensorIDSet.size() - 1 ; ++i, ++iter ) {
      int iDetector = *iter;
      string histoName = "HitXCorrelatioHisto_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TH2D * histo = (TH2D*) hitXFolder->Get(histoName.c_str());

      newTitle = "Hit X: " + toString( iDetector ) + " vs " + toString( iDetector + 1 );
      histo->SetTitle( newTitle.c_str() );

      newTitle = "Detector " + toString( iDetector ) + " [pixel] ";
      histo->SetXTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      histo->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      histo->Draw("colz");
      padVec[iPad]->Update();
      TPaveText * padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      padVec[iPad]->Modified( true );

      ++iPad;

      // find xMin, xMax, yMin, yMax
      Double_t xMin = histo->GetXaxis()->GetXmin();
      Double_t xMax = histo->GetXaxis()->GetXmax();
      Double_t yMin = histo->GetYaxis()->GetXmin();
      Double_t yMax = histo->GetYaxis()->GetXmax();

      // now make a slice of the histo in a 1 / 3 of the height
      Int_t yBin = histo->GetYaxis()->FindBin( ( yMax - yMin ) / 3.  + yMin ) ;
      cout << "yBin =" << yBin  << endl;
      string profName = string(histo->GetName()) + "_p1";
      TH1D * profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      Int_t xBin = profile->GetMaximumBin();
      cout << "xBin = " << xBin << endl;
      Double_t xA = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yA = (yMax - yMin) / 3 + yMin;

      // now at 2 / 3
      yBin = histo->GetYaxis()->FindBin( 2 * ( yMax - yMin ) / 3. + yMin ) ;
      cout << "yBin =" << yBin  << endl;
      profName = string(histo->GetName()) + "_p2";
      profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      xBin = profile->GetMaximumBin();
      cout << "xBin = " << xBin << endl;
      Double_t xB = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yB = 2 * ( yMax - yMin ) / 3. + yMin;

      cout << "A = (" << xA << ", " << yA << ") and B = (" << xB << ", " << yB << ")" << endl;

      // this is the candidate for the correlation line
      Double_t p1 = (yA - yB) / (xA - xB);
      Double_t p0 = yA - p1 * xA;

      // shift some pixels above
      Double_t shiftUp = 0.25; // mm
      Double_t p1_up = p1;
      Double_t p0_up = p0 + shiftUp;

      // Point C:
      Double_t xC, yC = 0;
      xC = ( yMin - p0_up) / p1_up;
      if ( xC < xMin ) {
        xC = xMin;
        yC = p1_up * xMin + p0_up;
      } else if ( xC > xMin && xC <= xMax ) {
        yC = yMin;
      } else if ( xC > xMax ) {
        xC = xMax;
        yC = p1_up * xMax + p0_up;
      }

      // Point D
      Double_t xD, yD = 0;
      xD = (yMax - p0_up) / p1_up;
      if ( xD > xMax ) {
        xD = xMax;
        yD = p1_up * xMax + p0_up;
      } else if ( xD > xMin && xD <= xMax ) {
        yD = yMax;
      } else if ( xD < xMin ) {
        xD = xMin;
        yD = p1_up * xMin + p0_up;
      }

      cout << "C = (" << xC << ", " << yC << ") and D = (" << xD << ", " << yD << ")" << endl;

      // shift some pixels above
      Double_t shiftDown = 0.25; // mm
      Double_t p1_down = p1;
      Double_t p0_down = p0 - shiftDown;

      // Point E
      Double_t xE, yE = 0;
      xE = ( yMax - p0_down ) / p1_down ;
      if ( xE > xMax ) {
        xE = xMax;
        yE = p1_down * xMax + p0_down;
      } else if ( xE > xMin && xE < xMax ) {
        yE = yMax;
      } else if ( xE < xMin ) {
        xE = xMin;
        yE = p1_down * xMin + p0_down;
      }

      // Point F
      Double_t xF, yF = 0;
      xF = ( yMin - p0_down ) / p1_down;
      if ( xF > xMax) {
        xF = xMax;
        yF = p1_down * xMax + p0_down;
      } else if ( xF > xMin && xF <= xMax ) {
        yF = yMin;
      } else if ( xF < xMin ) {
        xF = xMin;
        yF = p1_down * xMin + p0_down;
      }

      cout << "E = (" << xE << ", " << yE << ") and F = (" << xF << ", " << yF << ")" << endl;

      string cutName = "xHitCorrDiagonal_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TCutG * diagonal = new TCutG(cutName.c_str(),5);
      diagonal->SetPoint(0, xC, yC );
      diagonal->SetPoint(1, xD, yD );
      diagonal->SetPoint(2, xE, yE );
      diagonal->SetPoint(3, xF, yF );
      diagonal->SetPoint(4, xC, yC );
      diagonal->SetLineColor( kRed );
      diagonal->SetLineStyle( 2 );
      diagonal->Draw("LP");


      string correlationName = "xCorrHit_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      string option = "[" + cutName + "]";
      TProfile * correlation = histo->ProfileX(correlationName.c_str(), 1, -1, option.c_str());
      newTitle = "Fit (Det "  + toString( iDetector ) + " vs Det" + toString( iDetector + 1 ) + ")";
      correlation->SetTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      correlation->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      correlation->Fit("pol1");
      TF1 * fitFunc = (TF1*) correlation->GetFunction( "pol1" );
      fitFunc->SetLineColor( kBlue ) ;
      fitFunc->SetLineStyle( 2 );
      padVec[iPad]->Update();
      padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      TPaveStats * stats = (TPaveStats*) correlation->GetListOfFunctions()->FindObject("stats");
      stats->SetX1NDC( 0.62 );
      stats->SetY1NDC( 0.65 );

      padVec[iPad]->Modified( true );

      ++iPad;
    }
  }


  TDirectoryFile * hitYFolder = (TDirectoryFile*) correlatorFolder->Get("HitY");
  if ( hitYFolder != 0x0 ) {

    TString canvasBaseName    = "HitYCanvas";
    closeCanvases( canvasBaseName );

    // ok there is a cluster x folder, let's do it!
    padVec.clear();

    UInt_t nHisto = hitYFolder->GetListOfKeys()->GetSize();
    set< int > sensorIDSet;
    string separator = "-_";

    for ( size_t i = 0 ; i < nHisto; ++i ) {
      string name  = hitYFolder->GetListOfKeys()->At( i )->GetName();
      size_t last  = name.find_last_not_of( separator );
      size_t first = name.find_last_not_of( separator, last );
      sensorIDSet.insert( atoi( name.substr( first, last ).c_str() ) );
    }

    nDetector = sensorIDSet.size();
    nCanvas = nDetector / nDetPerCanvas;
    if ( nDetector % nDetPerCanvas != 0 ) {
      ++nCanvas;
    }


    Double_t titleHeight = 0.10;
    Int_t canvasWidth  = 800;
    Int_t canvasHeight = 800;
    for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

      string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
      string canvasTitle = string(runName) + " - Hit correlation (y) " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

      TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
      c->Range(0,0,1,1);
      c->SetBorderSize(0);
      c->SetFrameFillColor(0);
      c->SetBorderMode(0);
      canvasVec.push_back( c );

      // title pad
      TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
      titlePad->Draw();
      titlePad->SetBorderMode(0);
      titlePad->SetBorderSize(0);
      titlePad->SetFrameFillColor(0);
      titlePad->cd();
      TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
      title->SetBorderSize(1);
      title->SetLabel( canvasTitle.c_str() );
      title->Draw();
      c->cd();

      // big pad for the rest
      TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
      bigPad->Draw();
      bigPad->cd();
      bigPad->SetBorderMode(0);
      bigPad->SetBorderSize(0);
      bigPad->SetFrameFillColor(0);

      // divide the bigPad in 2 x 3 TPad and add them to the subpad list
      Int_t nX = 2, nY = 3;
      bigPad->Divide(nX, nY);


      for ( Int_t i = 0; i < nX * nY; i++ ) {
        TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
        smallPad->SetBorderMode(0);
        smallPad->SetBorderSize(0);
        smallPad->SetFrameFillColor(0);
        if ( padVec.size() < 2 * (nDetector - 1) ) {
          padVec.push_back( smallPad);
        }
      }
    }

    gStyle->SetPalette(1,0);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(11);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.35);
    gStyle->SetTitleFontSize( 0.05 );
    gStyle->SetTitleFillColor( kCyan - 9 );

    iPad = 0;
    string newTitle;
    set< int >::iterator iter = sensorIDSet.begin();
    for ( size_t i = 0 ; i < sensorIDSet.size() - 1 ; ++i, ++iter ) {
      int iDetector = *iter;
      string histoName = "HitYCorrelationHisto_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TH2D * histo = (TH2D*) hitYFolder->Get(histoName.c_str());

      newTitle = "Hit Y: " + toString( iDetector ) + " vs " + toString( iDetector + 1 );
      histo->SetTitle( newTitle.c_str() );

      newTitle = "Detector " + toString( iDetector ) + " [pixel] ";
      histo->SetXTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      histo->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      histo->Draw("colz");
      padVec[iPad]->Update();
      TPaveText * padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      padVec[iPad]->Modified( true );

      ++iPad;

      // find xMin, xMax, yMin, yMax
      Double_t xMin = histo->GetXaxis()->GetXmin();
      Double_t xMax = histo->GetXaxis()->GetXmax();
      Double_t yMin = histo->GetYaxis()->GetXmin();
      Double_t yMax = histo->GetYaxis()->GetXmax();

      // now make a slice of the histo in a 1 / 3 of the height
      Int_t yBin = histo->GetYaxis()->FindBin( ( yMax - yMin ) / 3.  + yMin ) ;
      cout << "yBin =" << yBin  << endl;
      string profName = string(histo->GetName()) + "_p1";
      TH1D * profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      Int_t xBin = profile->GetMaximumBin();
      cout << "xBin = " << xBin << endl;
      Double_t xA = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yA = (yMax - yMin) / 3 + yMin;

      // now at 2 / 3
      yBin = histo->GetYaxis()->FindBin( 2 * ( yMax - yMin ) / 3. + yMin ) ;
      cout << "yBin =" << yBin  << endl;
      profName = string(histo->GetName()) + "_p2";
      profile = histo->ProjectionX( profName.c_str(), yBin, 1 + yBin );
      xBin = profile->GetMaximumBin();
      cout << "xBin = " << xBin << endl;
      Double_t xB = profile->GetXaxis()->GetBinCenter( xBin );
      Double_t yB = 2 * ( yMax - yMin ) / 3. + yMin;

      cout << "A = (" << xA << ", " << yA << ") and B = (" << xB << ", " << yB << ")" << endl;

      // this is the candidate for the correlation line
      Double_t p1 = (yA - yB) / (xA - xB);
      Double_t p0 = yA - p1 * xA;

      // shift some pixels above
      Double_t shiftUp = 0.25; // mm
      Double_t p1_up = p1;
      Double_t p0_up = p0 + shiftUp;

      // Point C:
      Double_t xC, yC = 0;
      xC = ( yMin - p0_up) / p1_up;
      if ( xC < xMin ) {
        xC = xMin;
        yC = p1_up * xMin + p0_up;
      } else if ( xC > xMin && xC <= xMax ) {
        yC = yMin;
      } else if ( xC > xMax ) {
        xC = xMax;
        yC = p1_up * xMax + p0_up;
      }

      // Point D
      Double_t xD, yD = 0;
      xD = (yMax - p0_up) / p1_up;
      if ( xD > xMax ) {
        xD = xMax;
        yD = p1_up * xMax + p0_up;
      } else if ( xD > xMin && xD <= xMax ) {
        yD = yMax;
      } else if ( xD < xMin ) {
        xD = xMin;
        yD = p1_up * xMin + p0_up;
      }

      cout << "C = (" << xC << ", " << yC << ") and D = (" << xD << ", " << yD << ")" << endl;

      // shift some pixels above
      Double_t shiftDown = 0.25; // mm
      Double_t p1_down = p1;
      Double_t p0_down = p0 - shiftDown;

      // Point E
      Double_t xE, yE = 0;
      xE = ( yMax - p0_down ) / p1_down ;
      if ( xE > xMax ) {
        xE = xMax;
        yE = p1_down * xMax + p0_down;
      } else if ( xE > xMin && xE < xMax ) {
        yE = yMax;
      } else if ( xE < xMin ) {
        xE = xMin;
        yE = p1_down * xMin + p0_down;
      }

      // Point F
      Double_t xF, yF = 0;
      xF = ( yMin - p0_down ) / p1_down;
      if ( xF > xMax) {
        xF = xMax;
        yF = p1_down * xMax + p0_down;
      } else if ( xF > xMin && xF <= xMax ) {
        yF = yMin;
      } else if ( xF < xMin ) {
        xF = xMin;
        yF = p1_down * xMin + p0_down;
      }

      cout << "E = (" << xE << ", " << yE << ") and F = (" << xF << ", " << yF << ")" << endl;

      string cutName = "yHitCorrDiagonal_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      TCutG * diagonal = new TCutG(cutName.c_str(),5);
      diagonal->SetPoint(0, xC, yC );
      diagonal->SetPoint(1, xD, yD );
      diagonal->SetPoint(2, xE, yE );
      diagonal->SetPoint(3, xF, yF );
      diagonal->SetPoint(4, xC, yC );
      diagonal->SetLineColor( kRed );
      diagonal->SetLineStyle( 2 );
      diagonal->Draw("LP");


      string correlationName = "yCorrHit_d" + toString( iDetector ) + "_d" + toString( iDetector + 1 );
      string option = "[" + cutName + "]";
      TProfile * correlation = histo->ProfileX(correlationName.c_str(), 1, -1, option.c_str());
      newTitle = "Fit (Det "  + toString( iDetector ) + " vs Det" + toString( iDetector + 1 ) + ")";
      correlation->SetTitle( newTitle.c_str() );
      newTitle = "Detector " + toString( iDetector + 1) + " [pixel] ";
      correlation->SetYTitle( newTitle.c_str() );
      padVec[iPad]->cd();
      correlation->Fit("pol1");
      TF1 * fitFunc = (TF1*) correlation->GetFunction( "pol1" );
      fitFunc->SetLineColor( kBlue ) ;
      fitFunc->SetLineStyle( 2 );
      padVec[iPad]->Update();
      padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
      padTitle->SetX2NDC( 0.60 );
      padTitle->SetY1NDC( 0.90 );
      TPaveStats * stats = (TPaveStats*) correlation->GetListOfFunctions()->FindObject("stats");
      stats->SetX1NDC( 0.62 );
      stats->SetY1NDC( 0.65 );

      padVec[iPad]->Modified( true );

      ++iPad;
    }
  }

  // save every canvases
  string path( prepareOutputFolder( "Correlation" ) );
  for ( UInt_t iCanvas = 0; iCanvas < canvasVec.size() ; ++iCanvas ) {
    string figName = path + canvasVec[iCanvas]->GetName() + pictureOutputFormat.Data();
    canvasVec[iCanvas]->SaveAs( figName.c_str() );
  }

}

void showClusterPlot( const char * filename ) {

  // check if this file is already open
  TFile * inputFile = closeAndReopenFile( filename );

  if ( inputFile == 0x0 ) {
    cerr << "Problems opening file " << filename << endl;
    return;
  }

  // typical names for the processor folder
  vector< string > folderNames;
  folderNames.push_back( clusterHistoFolderName.Data() );
  folderNames.push_back( "Clustering"   );
  folderNames.push_back( "FilterHisto"  );
  folderNames.push_back( "ClusterHisto" );

  string dutFolderName =  clusterHistoFolderName.Data()  + toString( "_dut" );
  folderNames.push_back( dutFolderName );
  dutFolderName = clusterHistoFolderName.Data()  + toString( "_dep" );
  folderNames.push_back( dutFolderName );
  dutFolderName = clusterHistoFolderName.Data()  + toString( "_taki" );
  folderNames.push_back( dutFolderName );
  folderNames.push_back( toString("FilterHisto") + toString("_dut") );
  folderNames.push_back( toString("FilterHisto") + toString("_dep") );
  folderNames.push_back( toString("FilterHisto") + toString("_tel") );
  folderNames.push_back( toString("Clustering") +  toString("_dut") );
  folderNames.push_back( toString("Clustering") +  toString("_tel") );
  folderNames.push_back( toString("Clustering") +  toString("_dep") );


  // this function will create the following canvases:
  //
  // --> 1 canvas with SN plot (one every 3 detectors)
  // --> 1 canvas with hitmaps (one every 6 detectors)
  // --> 1 canvas with cluster noise (one every 6 detectors)
  // --> 1 canvas with the cluster multiplicity (one every 6 detectors)
  UInt_t nDetPerCanvas = 3;

  TString snrCanvasBaseName    = "SNRCanvas";

  // close all canvases with this name
  closeCanvases( snrCanvasBaseName );

  // look into the input file for a folder named
  TDirectoryFile * clusterHistoFolder = checkFolder( folderNames, inputFile );
  if ( clusterHistoFolder == 0x0 ) {
    return ;
  }

  // this folder should contain one subfolder for each detector
  // so guess the number of detectors
  UInt_t nDetector = clusterHistoFolder->GetListOfKeys()->GetSize();
  vector< string > detectorFolderNames;
  vector< int >    sensorIDs;
  string separator = "-_";

  for ( size_t i = 0 ; i < nDetector; ++i ) {
    string name( clusterHistoFolder->GetListOfKeys()->At( i )->GetName() );
    detectorFolderNames.push_back( name ) ;
    sensorIDs.push_back( atoi(name.substr( name.find_last_of( separator ) + 1,  name.length() ).c_str()) );
  }
  UInt_t nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  vector<TCanvas * > canvasVec;
  vector<TPad * >    padVec;
  Double_t titleHeight = 0.10;
  Int_t canvasWidth  = 800;
  Int_t canvasHeight = 800;
  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(snrCanvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - SNR histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 2 * nDetector ) {
        padVec.push_back( smallPad);
      }
    }
  }

  // set the style for these canvases
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.35);
  gStyle->SetTitleFontSize( 0.05 );
  gStyle->SetTitleFillColor( kCyan - 9 );

  // now plot the histograms in the right pad and make the fits
  UInt_t iPad = 0;
  string newTitle;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) clusterHistoFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > firstHistoNames;
    firstHistoNames.push_back( "seedSNR-d" + toString( sensorIDs.at(iDetector) ) );
    firstHistoNames.push_back( "seedSNR_d" + toString( sensorIDs.at(iDetector) ) );

    TH1D * firstHisto = NULL;
    for ( size_t iName = 0 ; iName < firstHistoNames.size(); ++iName ) {
      firstHisto =  (TH1D*) detectorFolder->Get( firstHistoNames.at( iName ).c_str() );
      if ( firstHisto != NULL ) break;
    }
    if ( firstHisto == NULL ) {
      cerr << "None of the histo name possibilities: " << endl;
      for ( size_t iName = 0; iName < firstHistoNames.size() ; ++iName ) {
        cerr << "\t" << firstHistoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return ;
    }

    newTitle = string(firstHisto->GetTitle()) + " - Detector " + toString( sensorIDs.at( iDetector ) );
    firstHisto->SetTitle( newTitle.c_str() );
    firstHisto->SetXTitle("SNR");
    firstHisto->SetFillColor( kCyan - 5 );
    padVec[iPad]->cd();
    firstHisto->Fit("landau");
    firstHisto->GetFunction("landau")->SetLineColor( kRed );
    ++iPad;


    vector< string > secondHistoNames;
    secondHistoNames.push_back( "clusterSNR3x3-d" + toString( sensorIDs.at(iDetector) ) );
    secondHistoNames.push_back( "clusterSNR3x3_d" + toString( sensorIDs.at(iDetector) ) );

    TH1D * secondHisto     = NULL;
    for ( size_t iName = 0; iName < secondHistoNames.size(); ++iName ) {
      secondHisto =   (TH1D*) detectorFolder->Get( secondHistoNames.at( iName ).c_str() );
      if ( secondHisto != NULL ) break;
    }
    if ( secondHisto == NULL ) {
      cerr << "None of the histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < secondHistoNames.size() ; ++iName ) {
        cerr << "\t" << secondHistoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string(secondHisto->GetTitle()) + " - Detector " + toString(  sensorIDs.at(iDetector) );
    secondHisto->SetTitle( newTitle.c_str() );
    secondHisto->SetXTitle("SNR");
    secondHisto->SetFillColor( kCyan - 5 );
    padVec[iPad]->cd();
    secondHisto->Fit("landau");
    secondHisto->GetFunction("landau")->SetLineColor( kRed );
    ++iPad;

  }

  // now display the hit maps
  nDetPerCanvas = 6;
  TString hitMapCanvasBaseName = "HitMapCanvas";

  // close all canvases with this name
  closeCanvases( hitMapCanvasBaseName );

  nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {
    string canvasName  = string(hitMapCanvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - HitMap histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c =  new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() <  nDetector ) {
        padVec.push_back( smallPad );
      }
    }
  }

  gStyle->SetPalette( 1, 0 );

  // now plot the histograms in the right pad
  iPad = 0;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) clusterHistoFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > histoNames;
    histoNames.push_back( "hitMap-d" + toString( sensorIDs.at( iDetector )   ) );
    histoNames.push_back( "hitMap_d" + toString( sensorIDs.at( iDetector )   ) );

    TH2D * histo = NULL;
    for ( size_t iName = 0; iName < histoNames.size(); ++iName ) {
      histo =  (TH2D*) detectorFolder->Get( histoNames.at( iName ).c_str() );
      if ( histo != NULL ) break;
    }
    if ( histo == NULL ) {
      cerr << "None of the map histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < histoNames.size() ; ++iName ) {
        cerr << "\t" << histoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string( histo->GetTitle() ) + " - Detector " + toString(  sensorIDs.at( iDetector ) );
    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle("x [pixel]");
    histo->SetYTitle("y [pixel]");
    histo->SetStats( false );
    padVec[iPad]->cd();
    histo->Draw("colz");
    ++iPad;
  }

  // now display the noise histograms
  nDetPerCanvas = 6;

  TString noiseCanvasBaseName = "NoiseCanvas";

  // close all canvases with this name
  closeCanvases( noiseCanvasBaseName );

  nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(noiseCanvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - Cluster noise histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < nDetector ) {
        padVec.push_back( smallPad);
      }
    }
  }


  iPad = 0;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {
    TDirectoryFile * detectorFolder = (TDirectoryFile*) clusterHistoFolder->Get( detectorFolderNames.at( iDetector ).c_str() );

    vector< string > histoNames;
    histoNames.push_back( "clusterNoise-d" + toString( sensorIDs.at( iDetector )   ) );
    histoNames.push_back( "clusterNoise_d" + toString( sensorIDs.at( iDetector )   ) );

    TH1D * histo = NULL;
    for ( size_t iName = 0; iName < histoNames.size(); ++iName ) {
      histo =  (TH1D*) detectorFolder->Get( histoNames.at( iName ).c_str() );
      if ( histo != NULL ) break;
    }
    if ( histo == NULL ) {
      cerr << "None of the map histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < histoNames.size() ; ++iName ) {
        cerr << "\t" << histoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string(histo->GetTitle()) + " - Detector " + toString( sensorIDs.at( iDetector ) );
    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle("Noise [ADC]");
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad]->cd();
    histo->Draw();
    padVec[iPad]->Update();
    TPaveStats * st = (TPaveStats*)  histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptStat(1110);
    padVec[iPad]->Modified( true );
    ++iPad;

  }

  // now display the event multiplicity canvas
  nDetPerCanvas = 6;

  TString multiplicityCanvasBaseName = "MultiplicityCanvas";

  // close all canvases with this name
  closeCanvases( multiplicityCanvasBaseName );

  nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(multiplicityCanvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - Cluster multiplicity histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < nDetector ) {
        padVec.push_back( smallPad);
      }
    }
  }


  iPad = 0;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) clusterHistoFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > histoNames;
    histoNames.push_back( "eventMultiplicity-d" + toString( sensorIDs.at( iDetector )   ) );
    histoNames.push_back( "eventMultiplicity_d" + toString( sensorIDs.at( iDetector )   ) );

    TH1D * histo = NULL;
    for ( size_t iName = 0; iName < histoNames.size(); ++iName ) {
      histo =  (TH1D*) detectorFolder->Get( histoNames.at( iName ).c_str() );
      if ( histo != NULL ) break;
    }
    if ( histo == NULL ) {
      cerr << "None of the map histo name possibilities:" << endl;
      for ( size_t iName = 0; iName < histoNames.size() ; ++iName ) {
        cerr << "\t" << histoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return;
    }

    newTitle = string(histo->GetTitle()) + " - Detector " + toString( sensorIDs.at( iDetector )  );
    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle("Multiplicity [#]");
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad]->cd();
    histo->Draw();
    padVec[iPad]->Update();
    TPaveStats * st = (TPaveStats*)  histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptStat(1110);
    padVec[iPad]->Modified( true );
    ++iPad;

  }

  // save every canvases
  string path( prepareOutputFolder( "Cluster" ) );
  for ( UInt_t iCanvas = 0; iCanvas < canvasVec.size() ; ++iCanvas ) {
    string figName = path + canvasVec[iCanvas]->GetName() + pictureOutputFormat.Data();
    canvasVec[iCanvas]->SaveAs( figName.c_str() );
  }

  return;
}

void showEtaPlot( const char * filename ) {

  // check if this file is already open
  TFile * inputFile = closeAndReopenFile( filename );

  if ( inputFile == 0x0 ) {
    cerr << "Problems opening file " << filename << endl;
    return;
  }

  vector< string > folderNames;
  folderNames.push_back( etaFolderName.Data() );

  string dutFolderName = etaFolderName.Data() + toString( "_dut" );
  folderNames.push_back( dutFolderName );
  dutFolderName = etaFolderName.Data() + toString( "_dep" );
  folderNames.push_back( dutFolderName );
  dutFolderName = etaFolderName.Data() + toString( "_taki" );
  folderNames.push_back( dutFolderName );

  // this function will create the following canvases:
  //
  // --> 1 canvas with CoG and ETA along x (one every 3 detectors)
  // --> 1 canvas with CoG and ETA along y (one every 3 detectors)
  // --> 1 canvas with 2D CoG map (one every 6 detectors);

  UInt_t nDetPerCanvas   = 3;
  TString canvasBaseName = "EtaXCanvas";

  // close all canvases with this name
  closeCanvases( canvasBaseName );

  // look into the input file for a folder named ETA
  TDirectoryFile * etaFolder = checkFolder( folderNames, inputFile );
  if ( etaFolder == 0x0 ) {
    return ;
  }

  // this folder should contain one subfolder for each detector
  // so guess the number of detectors
  UInt_t nDetector = etaFolder->GetListOfKeys()->GetSize();
  vector< string > detectorFolderNames;
  vector< int >    sensorIDs;
  string separator = "-_";

  for ( size_t i = 0; i < nDetector; ++i ) {
    string name( etaFolder->GetListOfKeys()->At( i )->GetName() );
    detectorFolderNames.push_back( name );
    sensorIDs.push_back( atoi(name.substr( name.find_last_of( separator ) + 1,  name.length() ).c_str()) );
  }

  UInt_t nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  vector<TCanvas * > canvasVec;
  vector<TPad * >    padVec;
  Double_t titleHeight = 0.10;
  Int_t canvasWidth  = 800;
  Int_t canvasHeight = 800;

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - ETA X histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 2 * nDetector ) {
        padVec.push_back( smallPad);
      }
    }
  }

  // set the style for these canvases
  gStyle->SetOptStat("emr");
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.35);
  gStyle->SetTitleFontSize( 0.05 );
  gStyle->SetTitleFillColor( kCyan - 9 );

  UInt_t iPad = 0;
  string newTitle;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) etaFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > cogHistoNames;
    cogHistoNames.push_back( "CoG-X-" + toString( sensorIDs.at(iDetector)) );
    cogHistoNames.push_back( "CoG_X_" + toString( sensorIDs.at(iDetector)) );

    TH1D * cogHisto = NULL;

    for ( size_t iName = 0 ; iName < cogHistoNames.size(); ++iName ) {
      cogHisto =  (TH1D*) detectorFolder->Get( cogHistoNames.at( iName ).c_str() );
      if ( cogHisto != NULL ) break;
    }
    if ( cogHisto == NULL ) {
      cerr << "None of the histo name possibilities: " << endl;
      for ( size_t iName = 0; iName < cogHistoNames.size() ; ++iName ) {
        cerr << "\t" << cogHistoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return ;
    }

    newTitle = string( cogHisto->GetTitle()) + " - Detector " + toString( iDetector );
    cogHisto->SetTitle( newTitle.c_str() );
    cogHisto->SetXTitle(" x [pitch unit] ");
    cogHisto->SetFillColor( kCyan - 5 );
    padVec[iPad++]->cd();
    cogHisto->Draw();

    vector< string > etaHistoNames;
    etaHistoNames.push_back( "EtaProfile-X-" + toString( sensorIDs.at(iDetector) ));
    etaHistoNames.push_back( "EtaProfile_X_" + toString( sensorIDs.at(iDetector) ));

    TProfile * etaHisto = NULL;
    for ( size_t iName = 0 ; iName < etaHistoNames.size(); ++iName ) {
      etaHisto =  (TProfile*) detectorFolder->Get( etaHistoNames.at( iName ).c_str() );
      if ( etaHisto != NULL ) break;
    }
    if ( etaHisto == NULL ) {
      cerr << "None of the histo name possibilities: " << endl;
      for ( size_t iName = 0; iName < etaHistoNames.size() ; ++iName ) {
        cerr << "\t" << etaHistoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return ;
    }

    newTitle = string( etaHisto->GetTitle()) + " - Detector " + toString( iDetector );
    etaHisto->SetTitle( newTitle.c_str() );
    etaHisto->SetXTitle(" x [picth unit] " );
    etaHisto->SetLineColor( kRed );
    padVec[iPad++]->cd();
    etaHisto->Draw();
  }

  nDetPerCanvas   = 3;
  canvasBaseName = "EtaYCanvas";

  // close all canvases with this name
  closeCanvases( canvasBaseName );

  nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - ETA Y histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 2 * nDetector ) {
        padVec.push_back( smallPad);
      }
    }
  }

  // set the style for these canvases
  gStyle->SetOptStat("emr");
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.35);
  gStyle->SetTitleFontSize( 0.05 );
  gStyle->SetTitleFillColor( kCyan - 9 );

  iPad = 0;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) etaFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > cogHistoNames;
    cogHistoNames.push_back( "CoG-Y-" + toString( sensorIDs.at(iDetector)) );
    cogHistoNames.push_back( "CoG_Y_" + toString( sensorIDs.at(iDetector)) );

    TH1D * cogHisto = NULL;

    for ( size_t iName = 0 ; iName < cogHistoNames.size(); ++iName ) {
      cogHisto =  (TH1D*) detectorFolder->Get( cogHistoNames.at( iName ).c_str() );
      if ( cogHisto != NULL ) break;
    }
    if ( cogHisto == NULL ) {
      cerr << "None of the histo name possibilities: " << endl;
      for ( size_t iName = 0; iName < cogHistoNames.size() ; ++iName ) {
        cerr << "\t" << cogHistoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return ;
    }

    newTitle = string( cogHisto->GetTitle()) + " - Detector " + toString( iDetector );
    cogHisto->SetTitle( newTitle.c_str() );
    cogHisto->SetXTitle(" y [pitch unit] ");
    cogHisto->SetFillColor( kCyan - 5 );
    padVec[iPad++]->cd();
    cogHisto->Draw();

 vector< string > etaHistoNames;
    etaHistoNames.push_back( "EtaProfile-Y-" + toString( sensorIDs.at(iDetector) ));
    etaHistoNames.push_back( "EtaProfile_Y_" + toString( sensorIDs.at(iDetector) ));

    TProfile * etaHisto = NULL;
    for ( size_t iName = 0 ; iName < etaHistoNames.size(); ++iName ) {
      etaHisto =  (TProfile*) detectorFolder->Get( etaHistoNames.at( iName ).c_str() );
      if ( etaHisto != NULL ) break;
    }
    if ( etaHisto == NULL ) {
      cerr << "None of the histo name possibilities: " << endl;
      for ( size_t iName = 0; iName < etaHistoNames.size() ; ++iName ) {
        cerr << "\t" << etaHistoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return ;
    }

    newTitle = string( etaHisto->GetTitle()) + " - Detector " + toString( iDetector );
    etaHisto->SetTitle( newTitle.c_str() );
    etaHisto->SetXTitle(" y [picth unit] " );
    etaHisto->SetLineColor( kRed ) ;
    padVec[iPad++]->cd();
    etaHisto->Draw();
  }

  nDetPerCanvas = 6;
  canvasBaseName = "CoGCanvas";

  // close all canvases with this name
  closeCanvases( canvasBaseName );

  nCanvas   = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - CoG histograms " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 1 * nDetector ) {
        padVec.push_back( smallPad);
      }
    }
  }

  gStyle->SetPalette(1,0);

  iPad = 0;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    TDirectoryFile * detectorFolder = (TDirectoryFile*) etaFolder->Get( detectorFolderNames.at(iDetector).c_str() );

    vector< string > cogHistoNames;
    cogHistoNames.push_back( "CoG-Histo2D-" + toString( sensorIDs.at(iDetector)) );
    cogHistoNames.push_back( "CoG_Histo2D_" + toString( sensorIDs.at(iDetector)) );

    TH2D * cogHisto = NULL;

    for ( size_t iName = 0 ; iName < cogHistoNames.size(); ++iName ) {
      cogHisto =  (TH2D*) detectorFolder->Get( cogHistoNames.at( iName ).c_str() );
      if ( cogHisto != NULL ) break;
    }
    if ( cogHisto == NULL ) {
      cerr << "None of the histo name possibilities: " << endl;
      for ( size_t iName = 0; iName < cogHistoNames.size() ; ++iName ) {
        cerr << "\t" << cogHistoNames.at( iName ) << endl;
      }
      cerr << "was found in " << filename << endl;
      return ;
    }

    newTitle = string( cogHisto->GetTitle()) + " - Detector " + toString( iDetector );
    cogHisto->SetTitle( newTitle.c_str() );
    cogHisto->SetXTitle(" x [pitch unit] ");
    cogHisto->SetYTitle(" y [pitch unit] ");
    cogHisto->SetStats( false );
    padVec[iPad++]->cd();
    cogHisto->Draw("colz");

  }

  string path( prepareOutputFolder( "ETA" ));
  for ( UInt_t iCanvas = 0 ; iCanvas < canvasVec.size(); ++iCanvas ) {
    string figName = path + canvasVec[iCanvas]->GetName() + pictureOutputFormat.Data();
    canvasVec[iCanvas]->SaveAs( figName.c_str() );
  }

}

void showTrackerPlot( const char * filename ) {

  // check if this file is already open
  TFile * inputFile = closeAndReopenFile( filename );

  if ( inputFile == 0x0 ) {
    cerr << "Problems opening file " << filename << endl;
    return;
  }
  vector< string > folderNames;
  folderNames.push_back( trackerFolderName.Data() );

  string dutFolderName = trackerFolderName.Data() + toString( "_dut" );
  folderNames.push_back( dutFolderName );
  dutFolderName = trackerFolderName.Data() + toString( "_dep" );
  folderNames.push_back( dutFolderName );
  dutFolderName = trackerFolderName.Data() + toString( "_taki" );
  folderNames.push_back( dutFolderName );

  // --> 1 canvas with reconstructed hit position (one every 6 det.)
  // --> 1 canvas with rotation check plots (one every 3 det.)
  // --> 1 canvas with shift check plot (one every 3 det. )

  UInt_t nDetPerCanvas = 6;

  // start from the hit reconstructed
  TString canvasBaseName = "RecoHitCanvas";

  // close all canvases with this name
  closeCanvases( canvasBaseName );

  // look into the input file for a folder named
  TDirectoryFile * trackerFolder = checkFolder( folderNames, inputFile );

  if ( trackerFolder == 0x0 ) {
    return;
  }

  // guess the number of sensors
  UInt_t nDetector = 0;
  while ( true ) {
    string name = "fittedXY_" + toString(nDetector);
    TH2D * histo = (TH2D*) trackerFolder->Get(name.c_str());
    if ( histo != 0x0 ) ++nDetector;
    else break;
  }
  if ( nDetector == 0 ) {
    cerr << "Something wrong with the number of detectors" << endl;
    return;
  }

  UInt_t nCanvas = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  vector<TCanvas * > canvasVec;
  vector<TPad * >    padVec;
  Double_t titleHeight = 0.10;
  Int_t canvasWidth  = 800;
  Int_t canvasHeight = 800;
  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - Reconstructed hit position " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 1 * nDetector ) {
        padVec.push_back( smallPad);
      }
    }
  }

  // set the style for these canvases
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.35);
  gStyle->SetTitleFontSize( 0.05 );
  gStyle->SetTitleFillColor( kCyan - 9 );
  gStyle->SetPalette(1,0);
  gStyle->SetLabelOffset( 0.01 );
  gStyle->SetLabelSize( 0.05 );
  gStyle->SetTitleOffset( 0.96 );
  gStyle->SetTitleSize( 0.05 );

  // now plot the histograms in the right pad and make the fits
  UInt_t iPad = 0;
  string newTitle;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector; ++iDetector ) {

    string histoName = "fittedXY_" + toString( iDetector );
    TH2D * histo = (TH2D*) trackerFolder->Get( histoName.c_str() );
    newTitle = "Reco hit map on detector " + toString( iDetector );
    histo->SetTitle( newTitle.c_str());
    histo->SetXTitle( "x [mm]" );
    histo->SetYTitle( "y [mm]" );
    padVec[iPad]->cd();
    histo->Draw("colz");
    padVec[iPad]->Update();
    setDefaultAxis( histo->GetXaxis() );
    setDefaultAxis( histo->GetYaxis() );
    TPaveText * padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle->SetX1NDC( 0.053 );
    padTitle->SetY1NDC( 0.867 );
    padTitle->SetX2NDC( 0.601 );
    padTitle->SetY2NDC( 0.998 );
    padVec[iPad]->Modified( true );
    ++iPad;
  }

  nDetPerCanvas = 3;

  canvasBaseName = "RotationCanvas";
  closeCanvases( canvasBaseName );

  // the number of histograms to plot is different from the number of
  // but I can get the name from the input list
  vector< string > rotX2DVec;
  vector< string > rotY2DVec;
  {
    TIter next( trackerFolder->GetListOfKeys() );
    while ( TObject * obj = next() ){
      TString objName = obj->GetName();

      if ( objName.BeginsWith( "relRotX2D_" )  ) {
        if ( find( rotX2DVec.begin(), rotX2DVec.end(), objName.Data()) == rotX2DVec.end() ) {
          rotX2DVec.push_back( objName.Data() );
        }
      }
      if ( objName.BeginsWith( "relRotY2D_" ) ) {
        if (  find( rotY2DVec.begin(), rotY2DVec.end(), objName.Data()) == rotY2DVec.end() ) {
          rotY2DVec.push_back( objName.Data() );
        }
      }
    }
    sort( rotX2DVec.begin(), rotX2DVec.end() );
    sort( rotY2DVec.begin(), rotY2DVec.end() );
  }

  nCanvas = rotX2DVec.size() / nDetPerCanvas;
  if ( rotX2DVec.size() % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - Rotation " + toString(iCanvas + 1) + " / " +  toString( nCanvas );


    Int_t xShift = 50 * canvasVec.size();
    if ( xShift > 600 ) xShift = 0;
    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), xShift, 0, canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 2 * rotX2DVec.size() ) {
        padVec.push_back( smallPad);
      }
    }
  }

  iPad = 0;
  TString oldTitle;
  for ( UInt_t iDetector = 0 ; iDetector < rotX2DVec.size(); ++iDetector ) {

    TH2D * histo = (TH2D*) trackerFolder->Get( rotX2DVec[iDetector].c_str() );
    setDefaultAxis( histo->GetXaxis() );
    setDefaultAxis( histo->GetYaxis() );

    // dirty game to retype the title...
    oldTitle = histo->GetTitle();
    cout << "rotX " << oldTitle << endl;
    if ( !oldTitle.BeginsWith( "#Delta" ) ) {
      TObjArray * objArray = oldTitle.Tokenize(" ");
      for ( Int_t iToken = 0 ; iToken < objArray->GetSize(); ++iToken ) {
        string objString = ((TObjString*) objArray->At(iToken))->GetString().Data();
        if ( objString == "plane" ) {
          newTitle = "#Deltax vs y on detector " + toString(  ((TObjString*) objArray->At(iToken + 1))->GetString().Data() );
          break;
        }
      }
      objArray->Delete();
    } else {
      newTitle = oldTitle.Data();
    }

    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle( "y [mm]" );
    histo->SetYTitle( "#Deltax [mm]");
    padVec[iPad]->cd();
    histo->Draw("colz");
    padVec[iPad]->Update();
    TPaveText * padTitle1 = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle1->SetX1NDC( 0.053 );
    padTitle1->SetY1NDC( 0.867 );
    padTitle1->SetX2NDC( 0.601 );
    padTitle1->SetY2NDC( 0.998 );
    padVec[iPad]->Modified( true );
    ++iPad;

    histo = (TH2D*) trackerFolder->Get( rotY2DVec[iDetector].c_str() );
    setDefaultAxis( histo->GetXaxis() );
    setDefaultAxis( histo->GetYaxis() );

    oldTitle = histo->GetTitle();
    if ( !oldTitle.BeginsWith( "#Delta" ) ) {
      TObjArray * objArray = oldTitle.Tokenize(" ");
      for ( Int_t iToken = 0 ; iToken < objArray->GetSize(); ++iToken ) {
        string objString = ((TObjString*) objArray->At(iToken))->GetString().Data();
        if ( objString == "plane" ) {
          newTitle = "#Deltay vs x on detector " + toString(  ((TObjString*) objArray->At(iToken + 1))->GetString().Data() );
          break;
        }
      }
      objArray->Delete();
    } else {
      newTitle = oldTitle.Data();
    }

    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle( "x [mm]" );
    histo->SetYTitle( "#Deltay [mm]");
    padVec[iPad]->cd();
    histo->Draw("colz");
    padVec[iPad]->Update();
    TPaveText * padTitle2 = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle2->SetX1NDC( 0.053 );
    padTitle2->SetY1NDC( 0.867 );
    padTitle2->SetX2NDC( 0.601 );
    padTitle2->SetY2NDC( 0.998 );
    padVec[iPad]->Modified( true );

    ++iPad;
  }


  canvasBaseName = "ShiftCanvas";
  closeCanvases( canvasBaseName );

  // the number of histograms to plot is different from the number of
  // but I can get the name from the input list

  vector< string > shiftXVec;
  vector< string > shiftYVec;

  {
    TIter next( trackerFolder->GetListOfKeys() );
    while ( TObject * obj = next() ){
      TString objName = obj->GetName();
      if ( objName.BeginsWith( "relShiftX_" ) ) {
        if ( find( shiftXVec.begin(), shiftXVec.end(), objName.Data() ) == shiftXVec.end() ) {
          shiftXVec.push_back( objName.Data() );
        }
      }
      if ( objName.BeginsWith( "relShiftY_" ) ) {
        if ( find( shiftYVec.begin(), shiftYVec.end(), objName.Data() ) == shiftYVec.end() ) {
          shiftYVec.push_back( objName.Data() );
        }
      }
    }
    sort( shiftXVec.begin(), shiftXVec.end() );
    sort( shiftYVec.begin(), shiftYVec.end() );
  }


  nCanvas = shiftXVec.size() /  nDetPerCanvas;
  if ( rotX2DVec.size() % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - Shift " + toString(iCanvas + 1) + " / " +  toString( nCanvas );


    Int_t xShift = 50 * canvasVec.size();
    if ( xShift > 600 ) xShift = 0;
    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), xShift, 0, canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 2 * rotX2DVec.size() ) {
        padVec.push_back( smallPad);
      }
    }
  }

  iPad = 0;
  for ( UInt_t iDetector = 0 ; iDetector < shiftXVec.size(); ++iDetector ) {

    TH1D * histo = (TH1D*) trackerFolder->Get( shiftXVec[iDetector].c_str() );
    setDefaultAxis( histo->GetXaxis() );

    // dirty game to retype the title...
    oldTitle = histo->GetTitle();
    if ( !oldTitle.BeginsWith("#Delta") ) {
      TObjArray * objArray = oldTitle.Tokenize(" ");
      for ( Int_t iToken = 0 ; iToken < objArray->GetSize(); ++iToken ) {
        string objString = ((TObjString*) objArray->At(iToken))->GetString().Data();
        if ( objString == "plane" ) {
          newTitle = "#Deltax on detector " + toString(  ((TObjString*) objArray->At(iToken + 1))->GetString().Data() );
          break;
        }
      }
      objArray->Delete();
    } else {
      newTitle = oldTitle.Data();
    }

    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle( "#Deltax [mm]" );
    padVec[iPad]->cd();
    histo->Fit("gaus");
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad]->Update();
    TPaveText * padTitle1 = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle1->SetX1NDC( 0.050 );
    padTitle1->SetY1NDC( 0.862 );
    padTitle1->SetX2NDC( 0.509 );
    padTitle1->SetY2NDC( 0.998 );

    TPaveStats * st = (TPaveStats*) histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptFit(111);
    st->SetX1NDC( 0.5476 );
    st->SetY1NDC( 0.4790 );
    st->SetX2NDC( 0.9792 );
    st->SetY2NDC( 0.9983 );

    TF1 * fitFunc = (TF1*) histo->GetListOfFunctions()->FindObject("gaus");
    fitFunc->SetLineColor( kBlue );
    Double_t mean = fitFunc->GetParameter(1);
    Double_t rms  = fitFunc->GetParameter(2);

    histo->GetXaxis()->SetRangeUser( mean - 5 * rms, mean + 5 * rms );
    padVec[iPad]->Modified( true );
    ++iPad;

    histo = (TH1D*) trackerFolder->Get( shiftYVec[iDetector].c_str() );
    setDefaultAxis( histo->GetXaxis() );

    // dirty game to retype the title...
    oldTitle = histo->GetTitle();
    if ( !oldTitle.BeginsWith( "#Delta" ) ) {
      TObjArray * objArray = oldTitle.Tokenize(" ");
      for ( Int_t iToken = 0 ; iToken < objArray->GetSize(); ++iToken ) {
        string objString = ((TObjString*) objArray->At(iToken))->GetString().Data();
        if ( objString == "plane" ) {
          newTitle = "#Deltay on detector " + toString(  ((TObjString*) objArray->At(iToken + 1))->GetString().Data() );
          break;
        }
      }
      objArray->Delete();
    } else {
      newTitle = oldTitle.Data();
    }

    histo->SetTitle( newTitle.c_str() );
    histo->SetXTitle( "#Deltay [mm]" );
    padVec[iPad]->cd();
    histo->Fit("gaus");
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad]->Update();
    padTitle1 = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle1->SetX1NDC( 0.050 );
    padTitle1->SetY1NDC( 0.862 );
    padTitle1->SetX2NDC( 0.509 );
    padTitle1->SetY2NDC( 0.998 );

    st = (TPaveStats*) histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptFit(111);
    st->SetX1NDC( 0.5476 );
    st->SetY1NDC( 0.4790 );
    st->SetX2NDC( 0.9792 );
    st->SetY2NDC( 0.9983 );

    fitFunc   = (TF1*) histo->GetListOfFunctions()->FindObject("gaus");
    fitFunc->SetLineColor( kBlue );
    mean      = fitFunc->GetParameter(1);
    rms       = fitFunc->GetParameter(2);

    histo->GetXaxis()->SetRangeUser( mean - 5 * rms, mean + 5 * rms );
    padVec[iPad]->Modified( true );
    ++iPad;
  }

  nDetPerCanvas = 3;
  canvasBaseName = "ResidualCanvas";
  closeCanvases( canvasBaseName );

  nCanvas = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ){
    ++nCanvas;
  }

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - Residual " + toString(iCanvas + 1) + " / " +  toString( nCanvas );

    Int_t xShift = 50 * canvasVec.size();
    if ( xShift > 600 ) xShift = 0;
    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), xShift, 0, canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 2 * nDetector ) {
        cout << smallPad << endl;
        padVec.push_back( smallPad);
      }
    }
  }

  iPad = 0;

  for ( UInt_t iDetector = 0; iDetector < nDetector ; ++iDetector ) {

    // residual x
    string histoName   = "residualX_" + toString( iDetector );
    string histoTitle  = "Residual X on plane " + toString( iDetector );
    TH1D * histo = (TH1D*)  trackerFolder->Get( histoName.c_str() );
    setDefaultAxis( histo->GetXaxis() );
    histo->SetTitle( histoTitle.c_str() );
    histo->SetXTitle( "#Deltax [mm]" );
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad]->cd();
    histo->Fit("gaus");
    padVec[iPad]->Update();

    TPaveText * padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle->SetX1NDC( 0.050 );
    padTitle->SetY1NDC( 0.862 );
    padTitle->SetX2NDC( 0.509 );
    padTitle->SetY2NDC( 0.998 );

    TPaveStats * st = (TPaveStats*) histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptFit(111);
    st->SetOptStat(10);
    st->SetOptFit(111);
    st->SetFitFormat("3.2e");
    st->SetX1NDC( 0.5476 );
    st->SetY1NDC( 0.4790 );
    st->SetX2NDC( 0.9792 );
    st->SetY2NDC( 0.9983 );

    TF1 * fitFunc = (TF1*) histo->GetListOfFunctions()->FindObject("gaus");
    fitFunc->SetLineColor( kBlue );
    Double_t mean      = fitFunc->GetParameter(1);
    Double_t rms       = fitFunc->GetParameter(2);

    histo->GetXaxis()->SetRangeUser( mean - 5 * rms, mean + 5 * rms );
    padVec[iPad]->Modified( true );
    ++iPad;

    // residual Y
    histoName   = "residualY_" + toString( iDetector );
    histoTitle  = "Residual Y on plane " + toString( iDetector );
    histo = (TH1D*)  trackerFolder->Get( histoName.c_str() );
    setDefaultAxis( histo->GetXaxis() );
    histo->SetTitle( histoTitle.c_str() );
    histo->SetXTitle( "#Deltay [mm]" );
    histo->SetFillColor( kCyan - 5 );
    padVec[iPad]->cd();
    histo->Fit("gaus");
    padVec[iPad]->Update();

    padTitle = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle->SetX1NDC( 0.050 );
    padTitle->SetY1NDC( 0.862 );
    padTitle->SetX2NDC( 0.509 );
    padTitle->SetY2NDC( 0.998 );


    st = (TPaveStats*) histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptFit(111);
    st->SetOptStat(10);
    st->SetFitFormat("3.2e");
    st->SetX1NDC( 0.5476 );
    st->SetY1NDC( 0.4790 );
    st->SetX2NDC( 0.9792 );
    st->SetY2NDC( 0.9983 );

    fitFunc = (TF1*) histo->GetListOfFunctions()->FindObject("gaus");
    fitFunc->SetLineColor( kBlue );
    mean      = fitFunc->GetParameter(1);
    rms       = fitFunc->GetParameter(2);

    histo->GetXaxis()->SetRangeUser( mean - 5 * rms, mean + 5 * rms );
    padVec[iPad]->Modified( true );
    ++iPad;

  }

  string path( prepareOutputFolder( trackerFolderName.Data() ));
  for ( UInt_t iCanvas = 0 ; iCanvas < canvasVec.size(); ++iCanvas ) {
    string figName = path + canvasVec[iCanvas]->GetName() + pictureOutputFormat.Data();
    canvasVec[iCanvas]->SaveAs( figName.c_str() );
  }

}


void showMillePlot( const char * filename ) {

  // check if this file is already open
  TFile * inputFile = closeAndReopenFile( filename );

  if ( inputFile == 0x0 ) {
    cerr << "Problems opening file " << filename << endl;
    return;
  }

  // typical names for the processor folder
  vector< string > folderNames;
  folderNames.push_back( milleFolderName.Data() );
  folderNames.push_back( "Mille" );

  string dutFolderName = milleFolderName.Data() + toString( "_dut" );
  folderNames.push_back( dutFolderName );
  dutFolderName = milleFolderName.Data() + toString( "_dep" );
  folderNames.push_back( dutFolderName );
  dutFolderName = milleFolderName.Data() + toString( "_taki" );
  folderNames.push_back( dutFolderName );

  // --> 1 canvas with residuals (one every 3 det. )

  UInt_t nDetPerCanvas = 3;

  // canvas name
  TString canvasBaseName = "MilleResiduals";

  // close all canvases with this name
  closeCanvases( canvasBaseName );

  // look into the input file for a folder named
  TDirectoryFile * milleFolder = checkFolder( folderNames, inputFile );
  if ( milleFolder == 0x0 ) {
    return ;
  }

  // guess the number of sensors
  UInt_t nDetector = 0;
  vector< int > sensorIDs;
  string separator = "_-" ;


  for ( int i = 0 ; i < milleFolder->GetListOfKeys()->GetSize(); ++i ) {
    string name( milleFolder->GetListOfKeys()->At( i )->GetName() );
    string compName = "ResidualX_d";

    if ( name.compare( 0, compName.length(), compName ) == 0 ) {
      sensorIDs.push_back( atoi( name.substr( name.find_last_not_of( separator ), name.length()).c_str()) );
    }
  }
  nDetector = sensorIDs.size();

  if ( nDetector == 0 ) {
    cerr << "Something wrong with the number of detectors" << endl;
    return;
  }

  UInt_t nCanvas = nDetector / nDetPerCanvas;
  if ( nDetector % nDetPerCanvas != 0 ) {
    ++nCanvas;
  }

  vector<TCanvas * > canvasVec;
  vector<TPad * >    padVec;
  Double_t titleHeight = 0.10;
  Int_t canvasWidth  = 800;
  Int_t canvasHeight = 800;

  padVec.clear();

  for ( UInt_t iCanvas = 0; iCanvas < nCanvas; iCanvas++ ) {

    string canvasName  = string(canvasBaseName.Data()) + "_" + toString(iCanvas);
    string canvasTitle = string(runName) + " - Mille residuals " + toString(iCanvas + 1) + " / " +  toString( nCanvas );


    Int_t xShift = 50 * canvasVec.size();
    if ( xShift > 600 ) xShift = 0;
    TCanvas * c = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), xShift, 0, canvasWidth, canvasHeight);
    c->Range(0,0,1,1);
    c->SetBorderSize(0);
    c->SetFrameFillColor(0);
    c->SetBorderMode(0);
    canvasVec.push_back( c );

    // title pad
    TPad * titlePad = new TPad("title","title",0, 1 - titleHeight,1,1);
    titlePad->Draw();
    titlePad->SetBorderMode(0);
    titlePad->SetBorderSize(0);
    titlePad->SetFrameFillColor(0);
    titlePad->cd();
    TPaveLabel * title = new TPaveLabel(0.10,0.10,0.90,0.90,"arc");
    title->SetBorderSize(1);
    title->SetLabel( canvasTitle.c_str() );
    title->Draw();
    c->cd();

    // big pad for the rest
    TPad * bigPad = new TPad("bigPad","bigPad", 0, 0, 1, 1 - titleHeight );
    bigPad->Draw();
    bigPad->cd();
    bigPad->SetBorderMode(0);
    bigPad->SetBorderSize(0);
    bigPad->SetFrameFillColor(0);

    // divide the bigPad in 2 x 3 TPad and add them to the subpad list
    Int_t nX = 2, nY = 3;
    bigPad->Divide(nX, nY);

    for ( Int_t i = 0; i < nX * nY; i++ ) {
      TPad * smallPad =  dynamic_cast<TPad*> (bigPad->cd( 1 + i ));
      smallPad->SetBorderMode(0);
      smallPad->SetBorderSize(0);
      smallPad->SetFrameFillColor(0);
      if ( padVec.size() < 2 * nDetector ) {
        padVec.push_back( smallPad );
      }
    }
  }

  UInt_t iPad = 0;
  vector< Double_t > xMeanVec, yMeanVec;
  for ( UInt_t iDetector = 0 ; iDetector < nDetector ; ++iDetector ) {

    string histoName  = "ResidualX_d" + toString( sensorIDs.at( iDetector ) );
    string histoTitle = "Residual along X for detector " + toString( sensorIDs.at( iDetector ) );
    TH1D * histo = (TH1D*) milleFolder->Get( histoName.c_str() );
    histo->SetTitle( histoTitle.c_str() );
    setDefaultAxis( histo->GetXaxis() );
    histo->SetXTitle( "#Deltax [#mum]" );
    padVec[iPad]->cd();
    histo->SetFillColor( kCyan - 5 );
    histo->Draw();
    padVec[iPad]->Update();
    Double_t xMax = histo->GetBinCenter( histo->GetMaximumBin() );
    histo->GetXaxis()->SetRangeUser( xMax - 300, xMax + 300 );
    xMeanVec.push_back( histo->GetMean() );

    TPaveText * padTitle1 = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle1->SetX1NDC( 0.050 );
    padTitle1->SetY1NDC( 0.862 );
    padTitle1->SetX2NDC( 0.509 );
    padTitle1->SetY2NDC( 0.998 );

    TPaveStats * st = (TPaveStats*) histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptFit(111);
    st->SetX1NDC( 0.5476 );
    st->SetY1NDC( 0.4790 );
    st->SetX2NDC( 0.9792 );
    st->SetY2NDC( 0.9983 );
    padVec[iPad]->Modified( true );

    ++iPad;

    histoName = "ResidualY_d" + toString( sensorIDs.at( iDetector ) );
    histoTitle = "Residual along Y for detector " + toString( sensorIDs.at( iDetector ) );
    histo = (TH1D*) milleFolder->Get(  histoName.c_str() );
    histo->SetTitle( histoTitle.c_str() );
    setDefaultAxis( histo->GetXaxis() );
    histo->SetXTitle( "#Deltay [#mum]" );
    padVec[iPad]->cd();
    histo->SetFillColor( kCyan - 5 );
    histo->Draw();
    padVec[iPad]->Update();
    xMax = histo->GetBinCenter( histo->GetMaximumBin() );
    histo->GetXaxis()->SetRangeUser( xMax - 300, xMax + 300 );
    yMeanVec.push_back( histo->GetMean() );


    padTitle1 = (TPaveText*) padVec[iPad]->GetListOfPrimitives()->FindObject("title");
    padTitle1->SetX1NDC( 0.050 );
    padTitle1->SetY1NDC( 0.862 );
    padTitle1->SetX2NDC( 0.509 );
    padTitle1->SetY2NDC( 0.998 );

    st = (TPaveStats*) histo->GetListOfFunctions()->FindObject("stats");
    st->SetOptFit(111);
    st->SetX1NDC( 0.5476 );
    st->SetY1NDC( 0.4790 );
    st->SetX2NDC( 0.9792 );
    st->SetY2NDC( 0.9983 );
    padVec[iPad]->Modified( true );

    ++iPad;
  }

  // printing the expected starting value.
  vector<Double_t >::iterator xMeanBegin = xMeanVec.begin();
  vector<Double_t >::iterator xMeanEnd   = xMeanVec.end();
  vector<Double_t >::iterator yMeanBegin = yMeanVec.begin();
  vector<Double_t >::iterator yMeanEnd   = yMeanVec.end();

  vector<Double_t >::iterator iter = xMeanVec.begin();

  Double_t xMean = (*xMeanBegin + *(xMeanEnd - 1) ) / 2;
  Double_t yMean = (*yMeanBegin + *(yMeanEnd - 1) ) / 2;

  cout << setw(20) << "Detector" << setw(20) << "x" << setw(20) << "y" << endl;
  for ( size_t iPos = 0 ; iPos < xMeanVec.size(); ++iPos ) {
    if ( ( iPos == 0 ) || ( iPos == xMeanVec.size() - 1 ) ) {
      cout << setw(20) << sensorIDs.at( iPos ) << setw(20)<< "0" << setw(20) << "0" << endl;
    } else {
      cout << setw(20) << sensorIDs.at( iPos ) << setw(20) << xMean - xMeanVec.at( iPos ) << setw(20) << yMean - yMeanVec.at( iPos ) << endl;
    }
  }



  string path( prepareOutputFolder( milleFolderName.Data() ));
  for ( UInt_t iCanvas = 0 ; iCanvas < canvasVec.size(); ++iCanvas ) {
    string figName = path + canvasVec[iCanvas]->GetName() + pictureOutputFormat.Data();
    canvasVec[iCanvas]->SaveAs( figName.c_str() );
  }

}

TFile * closeAndReopenFile( const char * filename ) {

  // check if the file is already opened
  TList * listOfFiles = (TList*) gROOT->GetListOfFiles();

  // loop over all files
  for ( Int_t iFile = 0 ; iFile < listOfFiles->GetSize()  ; ++iFile ) {
    TFile * currentFile = (TFile*) listOfFiles->At( iFile ) ;
    TString currentFileName = currentFile->GetName() ;
    TString filenameString( filename );
    if ( currentFileName.EndsWith( filename ) ){
      if ( currentFile != 0x0 ) {
        return currentFile;
      }
    }
  }

  return TFile::Open( filename );

}

void closeCanvases( const char * canvasBaseName ) {

  TList * listOfOpenedCanvas = (TList*) gROOT->GetListOfCanvases();
  for ( int i = 0 ; i < listOfOpenedCanvas->GetSize() ; ++i ) {
    TCanvas * canvas = (TCanvas*) listOfOpenedCanvas->At( i );
    TString canvasName2 = canvas->GetName();
    if ( canvasName2.Contains( canvasBaseName ) ) {
      canvas->Close();
    }
  }

}

void setRunName( const char * name ) {
  runName = name ;
}

const char * prepareOutputFolder( const char * type) {

  // check if a folder "pics" exists
  string pwd( gSystem->pwd() );


  if ( !gSystem->cd("pics") ) {
    // it doesn't exist, so create it
    gSystem->mkdir( "pics" );
    gSystem->cd("pics");
  }

  // replace blanks with "_" in the runname;
  runName.ReplaceAll(" ", 1, "_", 1);

  // check if a folder with this name already exists
  if ( !gSystem->cd( runName.Data() ) ) {
    // it doesn't exist, so create it
    gSystem->mkdir( runName.Data() );
    gSystem->cd( runName.Data() );
  }

  // check if a folder with the type specified already exists
  if ( !gSystem->cd( type ) ) {
    // it doesn't, create it!
    gSystem->mkdir( type );
    gSystem->cd( type );
  }


  string newFolder( gSystem->pwd() + string("/") );
  gSystem->cd( pwd.c_str() );

  return newFolder.c_str();
}

void usage() {

  vector< string > listOfFunction;
  listOfFunction.push_back( "void showPedeNoisePlot( const char * filename, const char * detector = \"mimotel\")" );
  listOfFunction.push_back( "void showClusterPlot( const char * filename )" );
  listOfFunction.push_back( "void showEtaPlot( const char * filename )" );
  listOfFunction.push_back( "void showCorrelationPlot( const char * filename )");
  listOfFunction.push_back( "void showMillePlot( const char * filename )" );
  listOfFunction.push_back( "void showTrackerPlot( const char * filename )" );


  cout << endl;
  cout << "First set the overall run name using " << endl;
  cout << "--->   setRunName( const char * ) " << endl;
  cout << "Then you can use the following function to display commonly used histograms:" << endl;

  vector< string >::iterator iter = listOfFunction.begin();
  while ( iter!= listOfFunction.end() ) {
    cout << "--->   " << (*iter++) << endl;
  }

}


void setDefaultAxis(TAxis * axis ) {

  axis->SetLabelOffset( 0.01 );
  axis->SetLabelSize( 0.05 );
  axis->SetTitleOffset( 0.96 );
  axis->SetTitleSize( 0.05 );

}

TDirectoryFile * checkFolder( vector< string >& candidateNames, TFile * inputFile ) {

  vector< string >::iterator iter = candidateNames.begin();
  vector< string > goodFolders;
  while ( iter != candidateNames.end() ) {
    if ( dynamic_cast< TDirectoryFile * > ( inputFile->Get( (*iter).c_str() ) ) ) {
      goodFolders.push_back( *iter );
    }
    ++iter;
  }

  if ( goodFolders.size() == 0 ) {
    cerr << "None of the following folders were found in the input file: " << endl;
    iter = candidateNames.begin();
    while ( iter != candidateNames.end() ) {
      cerr << " --> " << (*iter) << endl;
      ++iter;
    }

    // print out the directory contained into this file:
    TIter iterator( inputFile->GetListOfKeys() );
    cerr << endl << "This file contains the following folders: " << endl;
    size_t counter = 0;
    vector< string > availableFolders;
    while ( TKey * key = (TKey*) iterator() ) {
      string keyClassName = key->GetClassName();
      if ( keyClassName == "TDirectoryFile" ) {
        cerr << " " << counter << " --> " << key->GetName() << endl;
        availableFolders.push_back( key->GetName() );
        ++counter;
      }
    }
    cerr << " q --> to exit " << endl;

    bool badChoice = true;
    string selectedName;
    string ans;
    while ( badChoice ) {
      cout << "Choose one by the corresponding number [ 0 - " << counter - 1 << " ] " << endl;
      cin >> ans;

      if ( ans == "q" ) {
        cout << "Exiting " << endl;
        return NULL;
      } else {

        TString choice( ans );
        choice.ReplaceAll( " ", 1 , "", 0 );
        if ( choice.IsDigit() ) {
          try {
            selectedName = availableFolders.at( choice.Atoi() );
            badChoice = false;
            return  dynamic_cast< TDirectoryFile * > ( inputFile->Get( selectedName.c_str() ) );
          } catch ( std::exception  output_of_range ) {
            cout << "Invalid answer " << endl;
            badChoice = true;
          }

        }
      }
    }

  } else if ( goodFolders.size() == 1 ) {

    return dynamic_cast< TDirectoryFile * > ( inputFile->Get( goodFolders.at( 0 ).c_str() ) );

  } else if ( goodFolders.size() > 1 ) {

    cout << "Found more than one folder candidate. Plz choose one of the following:" << endl;
    iter = goodFolders.begin();
    size_t counter = 0;
    while ( iter != goodFolders.end() ) {
      cout << " " << counter << ". --> " << (*iter ) << endl;
      ++iter;
      ++counter;
    }
    cerr << " q --> to exit " << endl;

    bool badChoice = true;
    string selectedName;
    string ans;
    while ( badChoice ) {
      cout << "Choose one by the corresponding number [ 0 - " << counter - 1 << " ] " << endl;
      cin >> ans;

      if ( ans == "q" ) {
        cout << "Exiting " << endl;
        return NULL;
      } else {

        TString choice( ans );
        choice.ReplaceAll( " ", 1 , "", 0 );
        if ( choice.IsDigit() ) {
          try {
            selectedName = goodFolders.at( choice.Atoi() );
            badChoice = false;
            return  dynamic_cast< TDirectoryFile * > ( inputFile->Get( selectedName.c_str() ) );
          } catch ( std::exception  output_of_range ) {
            cout << "Invalid answer " << endl;
            badChoice = true;
          }

        }
      }
    }
  } else {
    return NULL;
  }


  return NULL;
}

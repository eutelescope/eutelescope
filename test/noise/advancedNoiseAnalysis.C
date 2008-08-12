#include <iostream>
#include <sstream>
#include <string>


#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
#include <TGraphErrors.h>
#include <TColor.h>
#include <TColorWheel.h>
#include <Rtypes.h>
#include <TLegend.h>
#include <TPaveLabel.h>

using namespace std;

// a few global constants
const size_t kXPixel = 256;
const size_t kYPixel = 256;

const size_t kNChan  = 4;
size_t xLimit[ kNChan + 1  ] = { 0, 64, 128 , 192, 256 };
int kColor[10] = { 
  kRed,
  kBlue,
  kBlack,
  kGreen,
  kTeal,
  kOrange,
  kPink
};



template<typename T>
string toString(T data ) {
  stringstream ss;
  ss << data ;
  return ss.str();

}

double * noiseMean = NULL, * noiseRMS = NULL, * channel = NULL;


void advancedNoiseAnalysis( unsigned int runNumber, unsigned int loop = 1) {
  
  string inputFileName = "./histo/run00" + toString( runNumber ) + "-ped-histo.root";
  string outputFileName = "./histo/run00" + toString( runNumber ) + "-adv-noise.root";
  

  // before opening the input and the output files, try to see if they
  // are not opened yet and in case close them before continue   
  TList * listOfOpenedFile = (TList*) gROOT->GetListOfFiles();
  for ( int i = 0; i < listOfOpenedFile->GetSize() ; ++i ) {
    TFile * file = (TFile*) listOfOpenedFile->At( i ) ;
    TString fileName(file->GetName());
    TString inputFileName1( inputFileName.c_str() );
    TString outputFileName1( outputFileName.c_str() );

    if (  ( fileName.Contains( inputFileName1 ) ) ||
	  ( inputFileName1.Contains( fileName ) ) ||
	  ( fileName.Contains( outputFileName1 ) ) ||
	  ( outputFileName1.Contains( fileName ) ) ) {
      cout << "Closing " << fileName << " before reopen " << endl;
      file->Close();
    }
  }


  // close also all the previously opened canvas
  TList * listOfOpenedCanvas = (TList*) gROOT->GetListOfCanvases();
  for ( int i = 0 ; i < listOfOpenedCanvas->GetSize() ; ++i ) {
    TCanvas * canvas = (TCanvas*) listOfOpenedCanvas->At( i );
    TString canvasName2 = canvas->GetName();
    if ( canvasName2.Contains( "det" ) ) {
      canvas->Close();
    }
  }

	 
  // now safely open the file
  TFile * inputFile = TFile::Open( inputFileName.c_str() ) ;
  TFile * outputFile = TFile::Open( outputFileName.c_str(), "RECREATE") ;
  TList * outputHistoList = new TList;

  // look into the inputFile for a folder named
  string pedeProcessorFolderName = "PedestalAndNoiseCalculator";
  TDirectoryFile * pedeProcessorFolder = (TDirectoryFile*) inputFile->Get( pedeProcessorFolderName.c_str() );
  
  if ( pedeProcessorFolder == 0 ) { 
    cerr << "No pedestal processor folder found in file " << inputFileName << endl;
    return ;
  }

  // this folder should contain one folder for each loop.
  string loopFolderName = "loop-" + toString( loop );
  TDirectoryFile * loopFolder = (TDirectoryFile *) pedeProcessorFolder->Get( loopFolderName.c_str() );
  
  if ( loopFolder == 0 ) {
    cerr << "No " << loopFolderName << " found in file " << inputFileName << endl;
    return ;
  }

  // guess the number of sensors from the number of subfolder in the loopfolder
  size_t nDetector = loopFolder->GetListOfKeys()->GetSize();
  cout << "This file contains " << nDetector << " detectors" << endl;

  // prepare arrays to store the mean and the rms of the noise distribution
  if ( noiseMean == NULL ) {
    delete [] noiseMean;
    noiseMean = NULL;
  }
  if ( noiseRMS == NULL ) {
    delete [] noiseRMS;
    noiseRMS = NULL;
  }
  if ( channel == NULL ) {
    delete [] channel;
    channel = NULL;
  }

  noiseMean = new double[ nDetector * kNChan ];
  noiseRMS  = new double[ nDetector * kNChan ];
  channel   = new double[ kNChan ];

  string canvasName = "comparison";
  string canvasTitle = "Noise comparison";

  TCanvas * comparisonCanvas = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), 1000, 500 );
  comparisonCanvas->Divide(1,2);
  
  TPad * topPad = (TPad*) comparisonCanvas->cd(1);
  topPad->Divide( nDetector );
  
  TPad * middlePad = (TPad *) comparisonCanvas->cd(2);
  middlePad->Divide( kNChan );


  // for each detector we have to get the noise map and to prepare 4
  // separe histos and maps
  for ( unsigned int iDetector = 0; iDetector < nDetector; iDetector++ ) {

    // get the noise map.
    string noiseMapName = "detector-" + toString( iDetector ) ;
    noiseMapName += "/NoiseMap-d" + toString( iDetector )  ;
    noiseMapName += "-l" + toString( loop ) ;
 
    TH2D * noiseMap = ( TH2D* ) loopFolder->Get( noiseMapName.c_str() ); 
    

    // create a folder in the output file
    TDirectory * subfolder = outputFile->mkdir( string( "detector_" + toString( iDetector ) ).c_str(),
						string( "detector_" + toString( iDetector ) ).c_str()
						);
    subfolder->cd();

    
    string canvasName = "det" + toString( iDetector );
    string canvasTitle = "Detector " + toString( iDetector );
    
    TCanvas * canvas = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), 1000, 500 );
    canvas->Divide( kNChan, 2 );

    // ok now start the loop on channels
    for ( size_t iChan = 0 ; iChan < kNChan ; ++iChan ) { 
	
      if ( iDetector == 0 ) channel[iChan] = iChan - 0.5;

      string tempName = "NoiseMap_d" + toString( iDetector ) + "_l" + toString( loop )	+ "_ch" + toString( iChan ) ;
      string tempTitle = "NoiseMap Det. " + toString( iDetector ) + " - Ch. " + toString( iChan ) ;

      TH2D * noiseMapCh = new TH2D ( tempName.c_str() , tempTitle.c_str(), 
				     kXPixel / kNChan , -0.5 + xLimit[ iChan ] , -0.5 + xLimit[ iChan + 1 ],
				     kYPixel, -0.5, -0.5 + kYPixel );
      noiseMapCh->SetXTitle("X [pixel]");
      noiseMapCh->SetYTitle("Y [pixel]");
      noiseMapCh->SetZTitle("Noise [ADC]");
      noiseMapCh->SetStats( false );
      outputHistoList->Add( noiseMapCh ) ;

      tempName = "NoiseDist_d" + toString( iDetector ) + "_l" + toString( loop )	+ "_ch" + toString( iChan ) ;
      
      tempTitle = "NoiseDist Det. " + toString( iDetector ) + " - Ch. " + toString( iChan ) ; 

      TH1D * noiseDistCh = new TH1D( tempName.c_str(), tempTitle.c_str(), 50, 0., 10. );
      noiseDistCh->SetXTitle("Noise [ADC]");
      noiseDistCh->SetLineColor( kColor[iDetector]  );
      noiseDistCh->SetLineStyle( iChan + 2 );
      noiseDistCh->SetLineWidth( 2 );
      outputHistoList->Add( noiseDistCh );

      // let's start looping on pixels now
      for ( size_t yPixel = 1 ; yPixel <= kYPixel ; ++yPixel ) {
	for ( size_t xPixel = xLimit[ iChan ] + 1; xPixel <= xLimit[ iChan +1 ] ; ++xPixel ) {
	  double noise = noiseMap->GetBinContent( xPixel , yPixel );
	  noiseMapCh->Fill( xPixel - 1 , yPixel - 1, noise );
	  noiseDistCh->Fill( noise );
	  
	}
      }

      canvas->cd( iChan + 1 ) ;
      noiseMapCh->Draw("colz");
      canvas->cd( iChan + kNChan + 1  );
      noiseDistCh->Draw();
      
      topPad->cd( iDetector + 1 );
      if ( iChan == 0 ) {
	noiseDistCh->Draw();
      } else {
	noiseDistCh->Draw("same");
      }

      middlePad->cd( iChan + 1 );
      if ( iDetector == 0 ) {
	noiseDistCh->Draw();
      } else {
	noiseDistCh->Draw("same");
      }


      noiseMean[ kNChan * iDetector + iChan ] = noiseDistCh->GetMean();
      noiseRMS[ kNChan * iDetector  + iChan ] = noiseDistCh->GetRMS();
    }
    canvas->Write();

  }

  canvasName = "summary";
  canvasTitle = "Noise summary";

  TCanvas * summaryCanvas = new TCanvas( canvasName.c_str(), canvasTitle.c_str(), 1000, 500 );
  summaryCanvas->SetGridx(1);
  TLegend * legend = new TLegend(0.5, 4.8, 1.5, 4.3,"","br");;
  

  for ( size_t iDetector = 0 ; iDetector < nDetector ; ++iDetector ) {
    
    TGraphErrors * gr = new TGraphErrors( kNChan, channel, &noiseMean[ iDetector * kNChan ], NULL, &noiseRMS[ iDetector * kNChan ] );
    gr->SetName( string( "NoisePerChannel_d" + toString( iDetector )).c_str());
    gr->SetTitle(string("Detector " + toString( iDetector )).c_str());
    gr->GetXaxis()->SetTitle("Channel #");
    gr->GetYaxis()->SetTitle("Noise [ADC]");
    gr->GetXaxis()->SetNdivisions( 5 );
    gr->GetXaxis()->SetLabelSize( 0 );
    gr->SetMarkerStyle( iDetector + 1 );
    gr->SetMarkerColor( kColor[iDetector] );
    gr->SetLineColor( kColor[iDetector] );
    gr->SetLineWidth( 2 );

    
    legend->AddEntry( gr, string("Detector " + toString( iDetector )).c_str(), "LP");

    if ( iDetector == 0 ) {
      gr->Draw("ALP");
    } else {
      gr->Draw("LP");
    }
    

  }

  
  
  legend->Draw();

  for ( size_t iChan = 0 ; iChan < kNChan ; ++iChan ) {
    
    TPaveLabel * label = new TPaveLabel( iChan - 0.75 , 3.2 , iChan -0.25 , 3, string("Ch " + toString( iChan ) ).c_str());
    label->Draw();
  }


  summaryCanvas->Write();
  comparisonCanvas->Write();

  outputHistoList->Write();

  

 
} 

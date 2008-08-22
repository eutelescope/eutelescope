// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

#include <TString.h>
#include <TFile.h>
#include <TROOT.h>
#include <TDirectoryFile.h>
#include <TH2.h>
#include <TCanvas.h>

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

//! Converts to string
/*! This template function converts to string any type of data that
 *  has a proper ostream streamer defined.
 */
template<typename T>
string toString(T data ) {
  stringstream ss;
  ss << data ;
  return ss.str();

}


//! This is the real function
/*! This function is plotting the results of the EUTelCorrelator
 *  processor.
 *
 *  It produces also an output ROOT file with a canvas having all the
 *  correlation plots nicely displayed.
 *
 *  @param inputFileName The ROOT file containing the original
 *  histograms
 */

void correlationAnalysis(const char * filename) {

  // the inputFileName should already have the proper path, so that
  // the new output file can have exactly the same name as the input
  // with appended something to make it different.
  // 
  // I don't like the idea of opening the inputFileName in append
  // mode! 
  
  
  string inputFileName( filename ); // converted to string ! ! !
  const string extension( ".root" );
  string outputFileName = inputFileName.substr( 0, inputFileName.rfind( extension, inputFileName.length() )) + "-new" + extension;

  // before opening the input and the output files, try to see if they
  // are not opened yet and in case close them before continue   
  TList * listOfOpenedFile = (TList*) gROOT->GetListOfFiles();
  for ( int i = 0; i < listOfOpenedFile->GetSize() ; ++i ) {
    TFile * file = (TFile*) listOfOpenedFile->At( i ) ;
    TString fileName(file->GetName());
    TString inputFileName1( inputFileName.c_str() );
    TString outputFileName1( outputFileName.c_str() );

    if ( ( fileName.Contains( inputFileName1 ) ) ||
         ( fileName.Contains( outputFileName1 ) ) ) {
      cout << "Closing " << fileName << " before reopen " << endl;
      file->Close();
    }
  }



  // now we are finally ready to reopen the input file in read only and
  // the output one in r/w
  TFile * inputFile = TFile::Open( inputFileName.c_str() ) ;
  TFile * outputFile = TFile::Open( outputFileName.c_str(), "RECREATE") ;

  const string correlatorFolderName = "Correlator";

  TDirectoryFile * correlatorFolder =  (TDirectoryFile*) inputFile->Get( correlatorFolderName.c_str() );

  if ( correlatorFolder == 0x0 ) {
    cerr << "No correlation folder found." << endl;
    return;
  }

  // within this folder there can be up two 4 subfolders. Anyway it's
  // a even number.
  TList * listOfFolder = (TList*) correlatorFolder->GetListOfKeys();
  
  // loop over all the entries in this subfolder
  for ( short iSubFolder = 0; iSubFolder < listOfFolder->GetSize() ; ++iSubFolder ) {
    TDirectoryFile * subFolder = (TDirectoryFile*) correlatorFolder->Get( listOfFolder->At( iSubFolder )->GetName() );
    cout << "Processing " << subFolder->GetName() << endl;

    // within the subfolder we should have the histograms with the
    // correlation plots. 
    unsigned short noOfHistos  = subFolder->GetListOfKeys()->GetSize();
    unsigned short noOfSensors = 0.5* ( 1 + sqrt( 1 + 4 * noOfHistos ) );
    cout << "The numebr of sensors is " << noOfSensors << endl;
    TCanvas * c = new TCanvas(subFolder->GetName(), subFolder->GetName(), 1000, 500 );
    c->Divide( noOfSensors, noOfSensors );

    TList * listOfHisto = (TList*) subFolder->GetListOfKeys();

    int iPos = 0;
    int iPos2 = 0;
    for ( int x = 0 ; x < noOfSensors; ++x ) {
      for ( int y = 0; y < noOfSensors; ++y ) {
        c->cd( iPos2 + 1 );

        if ( x != y ) {
          TH2D * histo = (TH2D*) subFolder->Get( listOfHisto->At( iPos )->GetName() );
          histo->Draw("colz");
          
          ++iPos;
        }
        ++iPos2;
      }
    }
  }


  return ;
}

#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <cstdlib>
#include "TROOT.h"



#include "langaufit.h"

using namespace std;

int main(int argc, char *argv[]){
	// argv[1] is the root file
	// argv[2] is the optinal name
	
	char tmpchar[100];

	TFile *f = TFile::Open(argv[1]);
	TH1D *h0 = (TH1D*)f->Get("hHitAmp_run2006_id6");
	TH1D *h1 = (TH1D*)f->Get("hHitAmp_run2006_id7");
	
	string aname = "./";
	
	
	Double_t chisqr;
	Int_t    ndf;
	Double_t LangauPeak, LangauFWHM;
	Double_t fitrange[2];
   Double_t startvalues[4], parlimitslow[4], parlimitshigh[4];
	Double_t fitparameters[4], fitparameters_error[4];
	
	TCanvas *cc0 =  new TCanvas("cc0","cc0",800,600);
	TCanvas *cc1 =  new TCanvas("cc1","cc1",800,600);
	
	//gStyle->SetOptStat(11111111);
	//gStyle->SetOptFit(1111);

	TF1 *final_fit0;
	TF1 *final_fit1;

	//////////  ---- chip0 -----

	//h0->Rebin(2);
	h0->SetAxisRange(10,100,"X");
	fitrange[0]=12.0;
	fitrange[1]=35.0;
	
	//   par[0]=Width (scale) parameter of Landau density
	startvalues[0]=1.0; parlimitslow[0]=0.05; parlimitshigh[0]=20.0;
   //   par[1]=Most Probable (MP, location) parameter of Landau density
	startvalues[1]=17.0; parlimitslow[1]=10.0; parlimitshigh[1]=30.0;
   //   par[2]=Total area (integral -inf to inf, normalization constant)
	startvalues[2]=2000.0; parlimitslow[2]=1000.; parlimitshigh[2]=6000.0;
   //   par[3]=Width (sigma) of convoluted Gaussian function
	startvalues[3]=4.; parlimitslow[3]=0.001; parlimitshigh[3]=20.0;

	cc0->cd();
	
	final_fit0 = langaufit(h0,fitrange,startvalues,parlimitslow,parlimitshigh,fitparameters,fitparameters_error,&chisqr,&ndf);
	langaupro(fitparameters,LangauPeak,LangauFWHM);
	
	std::cout<<"Fit_results_for chip0 "<<argv[2]<<std::endl;
	// final_fit->Print();
	
	cout<<"Width "<<fitparameters[0]<<" Error "<<fitparameters_error[0]<<endl;
	cout<<"MP "<<fitparameters[1]<<" Error "<<fitparameters_error[1]<<endl;
	cout<<"Area "<<fitparameters[2]<<" Error "<<fitparameters_error[2]<<endl;
	cout<<"GSigma "<<fitparameters[3]<<" Error "<<fitparameters_error[3]<<endl;
	
	cout<<"chisqr "<<chisqr<<" ndf "<<ndf<<" chisqr/ndf "<<chisqr/ndf<<endl;

	sprintf(tmpchar,"%s/LangausFit_chip0_%s.C",aname.c_str(),argv[2]);
	cc0->SaveAs(tmpchar);
	
	sprintf(tmpchar,"%s/LangausFit_chip0_%s.pdf",aname.c_str(),argv[2]);
	cc0->SaveAs(tmpchar);
	
	
	//////////  ---- chip1 -----
	
	h1->Rebin(2);
	h1->SetAxisRange(15,100,"X");
	fitrange[0]=20.0;
	fitrange[1]=80.0;
	cc1->cd();
	
	//   par[0]=Width (scale) parameter of Landau density
	startvalues[0]=1.0; parlimitslow[0]=0.5; parlimitshigh[0]=20.0;
   //   par[1]=Most Probable (MP, location) parameter of Landau density
	startvalues[1]=23.0; parlimitslow[1]=10.0; parlimitshigh[1]=30.0;
   //   par[2]=Total area (integral -inf to inf, normalization constant)
	startvalues[2]=20000.0; parlimitslow[2]=10000.; parlimitshigh[2]=100000.0;
   //   par[3]=Width (sigma) of convoluted Gaussian function
	startvalues[3]=4.0; parlimitslow[3]=0.1; parlimitshigh[3]=10.0;
	

	
	final_fit1 = langaufit(h1,fitrange,startvalues,parlimitslow,parlimitshigh,fitparameters,fitparameters_error,&chisqr,&ndf);
	langaupro(fitparameters,LangauPeak,LangauFWHM);
	
	std::cout<<"Fit_results_for chip1 "<<argv[2]<<std::endl;
	// final_fit->Print();
	
	cout<<"Width "<<fitparameters[0]<<" Error "<<fitparameters_error[0]<<endl;
	cout<<"MP "<<fitparameters[1]<<" Error "<<fitparameters_error[1]<<endl;
	cout<<"Area "<<fitparameters[2]<<" Error "<<fitparameters_error[2]<<endl;
	cout<<"GSigma "<<fitparameters[3]<<" Error "<<fitparameters_error[3]<<endl;
	
	cout<<"chisqr "<<chisqr<<" ndf "<<ndf<<" chisqr/ndf "<<chisqr/ndf<<endl;
	
	
	sprintf(tmpchar,"%s/LangausFit_chip1_%s.C",aname.c_str(),argv[2]);
	cc1->SaveAs(tmpchar);
	
	sprintf(tmpchar,"%s/LangausFit_chip1_%s.pdf",aname.c_str(),argv[2]);
	cc1->SaveAs(tmpchar);
	
	

	return 0;
	
}


   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf
	

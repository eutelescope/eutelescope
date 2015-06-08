#ifndef CROSSSECTION_HPP
#define CROSSSECTION_HPP

#include <vector>
#include <iostream>


#include "TObject.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TColor.h"
#include "TString.h"
#include "TMath.h"
#include "TDirectory.h"

class CrossSection: public TObject {

    private:
        TH2F *fH2;               // input 2D histogram
        Int_t fNbinsH2;        // num bins of x-axis 2D hist

        TH1F *fHcsXY;           // cross section along central x or y axis of 2D histogram
        TH1F *fHcsBorder;       // cross section along border of 2D histogram
        TH1F *fHcsDiag;         // cross section along central diagonal of 2D histogram

        std::vector< TH1F* > fVcs;
    
        TLegend *fLeg;            // legend

        Int_t fNbins;           // num bins of x-axis cs hists
        Float_t fBinWidth;      // bin width of h2 
        Float_t *fBins;          // array containing bin boundaries of cs hists (variable bin width)

        Bool_t fNbinsEven;     // n_bins even or odd
       
        Bool_t fCScreated;     // n_bins even or odd
        Bool_t fDrawLegend;    // flag for drawing legend

    public:
        CrossSection (TH2F *h2, TString name_outhist_prefix="");
        ~CrossSection ();

        void CreateCrossSections();                     // get vector with histograms for the three sections
        void RecreateCrossSections();                     // get vector with histograms for the three sections

        std::vector< TH1F* > GetCrossSections();                     // get vector with histograms for the three sections
        TLegend *GetLegend();                                 

        void DrawCrossSection(TString option="", Bool_t create_canvas=kFALSE);     // draw the histrograms 
       
        void SetCSlineColor(Int_t index_cs=0, Color_t color=kBlack);
        void SetDrawLegend(Bool_t draw_legend=kTRUE);

//        ClassDef(CrossSection,1);
};

#endif


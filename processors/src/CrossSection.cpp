#include <vector>
#include <iostream>

#include "CrossSection.hpp"

using namespace std;

//ClassImp(CrossSection)

// constructor
//_______________________________________________________________________________
CrossSection::CrossSection(TH2F *h2, TString name_outhist_prefix)
    :TObject(),
     fH2(h2),
     fNbinsH2(h2->GetNbinsX()),
     fHcsXY(0), 
     fHcsBorder(0),
     fHcsDiag(0),
     fVcs(0),
     fLeg(0),
     fNbins(-1),
     fBinWidth(h2->GetXaxis()->GetBinWidth(1)),
     fBins(0),
     fNbinsEven(0),
     fCScreated(0),
     fDrawLegend(0)
{
    // check if bin number even or odd
    fNbinsEven = kFALSE;
    if (fNbinsH2%2 == 0) {
        fNbinsEven = kTRUE;
    }
     
    // number of bins for 1D hists
    if (fNbinsEven) {
        fNbins = fNbinsH2/2;
    }
    else {
        fNbins = fNbinsH2/2 + 1;
    }

    // different bin width for the diagonal cross section
    fBins = new Float_t[fNbins*3+1];
    for (Int_t i_bin=0; i_bin<=fNbins; i_bin++) {
        fBins[i_bin] = i_bin*fBinWidth;
        fBins[i_bin+fNbins] = (i_bin+fNbins)*fBinWidth;
        fBins[i_bin+fNbins*2] = fNbins*2*fBinWidth + (i_bin)*fBinWidth*TMath::Sqrt(2);            
    }
    
    // get name of provided 2D hist and make it base for name of 1D hist
    if (name_outhist_prefix.CompareTo("")==0) {  
        TString name_h2 = h2->GetName();
        Int_t i_hist = 0; 
        while(gDirectory->Get(Form("%s_border_%i",name_h2.Data(), i_hist))!=NULL) {
            i_hist++;
        }
        fHcsBorder = new TH1F(Form("%s_border_%i", name_h2.Data(), i_hist),"section along border of 2D histogram", fNbins*3, fBins);
        fHcsXY = new TH1F(Form("%s_center_%i", name_h2.Data(), i_hist),"section along central x or y axis of 2D histogram", fNbins*3, fBins);
        fHcsDiag = new TH1F(Form("%s_diag_%i", name_h2.Data(), i_hist),"section along central diagonal of 2D histogram", fNbins*3, fBins);
    }
    else {
        fHcsBorder = new TH1F(Form("%s_border", name_outhist_prefix.Data()),"section along border of 2D histogram", fNbins*3, fBins);
        fHcsXY = new TH1F(Form("%s_center", name_outhist_prefix.Data()),"section along central x or y axis of 2D histogram", fNbins*3, fBins);
        fHcsDiag = new TH1F(Form("%s_diag", name_outhist_prefix.Data()),"section along central diagonal of 2D histogram", fNbins*3, fBins);
    }

    fCScreated = kFALSE;
    fDrawLegend = kTRUE;
}

// destructor
//_______________________________________________________________________________
CrossSection::~CrossSection() {
    fVcs.clear(); 
    delete fLeg;        fLeg=0;
    delete fHcsXY;      fHcsXY=0;
    delete fHcsDiag;    fHcsDiag=0;
    delete fHcsBorder;  fHcsBorder=0;
    delete fBins;       fBins=0;
}

// create histograms for cross sections
//_______________________________________________________________________________
void  CrossSection::CreateCrossSections() {
    // assuming symmetry, average over similar pixels
    for (Int_t i_col=1; i_col<=fNbins; i_col++) {
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(i_col, 1));
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(1, i_col));
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(fNbinsH2-i_col+1, 1));
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(1, fNbinsH2-i_col+1));
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(i_col, fNbinsH2));
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(fNbinsH2, i_col));
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(fNbinsH2-i_col+1, fNbinsH2));
        fHcsBorder->AddBinContent(i_col, fH2->GetBinContent(fNbinsH2, fNbinsH2-i_col+1));

        fHcsDiag->AddBinContent(fNbins*3 - i_col + 1, fH2->GetBinContent(i_col, i_col));
        fHcsDiag->AddBinContent(fNbins*3 - i_col + 1, fH2->GetBinContent(fNbinsH2-i_col+1, fNbinsH2-i_col+1));
        fHcsDiag->AddBinContent(fNbins*3 - i_col + 1, fH2->GetBinContent(i_col, fNbinsH2-i_col+1));
        fHcsDiag->AddBinContent(fNbins*3 - i_col + 1, fH2->GetBinContent(fNbinsH2-i_col+1, i_col));

        fHcsXY->AddBinContent(i_col + fNbins, fH2->GetBinContent(i_col, fNbins));
        fHcsXY->AddBinContent(i_col + fNbins, fH2->GetBinContent(fNbinsH2-i_col+1, fNbins));
        fHcsXY->AddBinContent(i_col + fNbins, fH2->GetBinContent(fNbins, i_col));
        fHcsXY->AddBinContent(i_col + fNbins, fH2->GetBinContent(fNbins, fNbinsH2-i_col+1));
    }
    // normalize
    fHcsXY->Scale(1./4);
    fHcsDiag->Scale(1./4);
    fHcsBorder->Scale(1./8);

    fHcsXY->GetXaxis()->SetTitle("Distance along cross section");
    fHcsDiag->GetXaxis()->SetTitle("Distance along cross section");
    fHcsBorder->GetXaxis()->SetTitle("Distance along cross section");

    for (UInt_t i=0; i<3;i++) {
        SetCSlineColor(i,i+2);
    }

    fCScreated = kTRUE;
}

// recreate histograms for cross sections
//_______________________________________________________________________________
void  CrossSection::RecreateCrossSections() {
    fHcsXY->Reset();
    fHcsDiag->Reset();
    fHcsBorder->Reset();
   
    CreateCrossSections();
}

// get vector with histograms for cross sections
//_______________________________________________________________________________
vector< TH1F* >  CrossSection::GetCrossSections() {
    if (!fCScreated) {
        CreateCrossSections();
    }

    fVcs.push_back(fHcsBorder);
    fVcs.push_back(fHcsXY);
    fVcs.push_back(fHcsDiag);

    return fVcs;
}

// get legend
//_______________________________________________________________________________
TLegend *CrossSection::GetLegend() {

    fLeg = new TLegend(0.15, 0.15, 0.5, 0.35); 
    fLeg->AddEntry(fHcsBorder, "border","lp");
    fLeg->AddEntry(fHcsXY, "central axis","lp");
    fLeg->AddEntry(fHcsDiag, "central diagonal","lp");
    
    fLeg->SetFillColor(kWhite);
    fLeg->SetBorderSize(0);
    fLeg->SetTextSize(fHcsBorder->GetXaxis()->GetTitleSize());
    
    return fLeg;
}

// draw the cross sections in summary plot
//_______________________________________________________________________________
void CrossSection::DrawCrossSection(TString option, Bool_t create_canvas) { 
    if (!fCScreated) {
        CreateCrossSections();
    }

    if (create_canvas) {
        TCanvas *c=new TCanvas("c", "c", 1200, 600);
         c->cd();
    } 

    fHcsBorder->SetTitle("");
    //fHcsBorder->GetYaxis()->SetRangeUser(0.0,1.05);
    fHcsBorder->Draw(option.Data());
    fHcsXY->Draw("same");
    fHcsDiag->Draw("same");
    fHcsBorder->Draw("axis same");
    
    if (fDrawLegend) {
        GetLegend()->Draw();
    }
}

// change the line colors of the hists
//_______________________________________________________________________________
void CrossSection::SetCSlineColor(Int_t index_cs, Color_t color) { 

    if (index_cs==0) {
        fHcsBorder->SetLineColor(color);
        //return kTRUE;
    } 
    else if (index_cs==1)  {
        fHcsXY->SetLineColor(color);
        //return kTRUE;
    }
    else if (index_cs==2)  {
        fHcsDiag->SetLineColor(color);
        //return kTRUE;
    }
    else {
        cout<<"index not within range, please check!!"<<endl;
        //return kFALSE;
    }

}

// set flag to draw legend
//_______________________________________________________________________________
void CrossSection::SetDrawLegend(Bool_t draw_legend) { 
    fDrawLegend = draw_legend;
}




// Mainframe macro generated from application: /opt/root/bin/root.exe
// By ROOT version 5.20/00 on 2008-08-10 17:39:30

#if !defined( __CINT__) || defined (__MAKECINT__)

#ifndef ROOT_TStyle
#include "TStyle.h"
#endif

#endif

void Style_AdvNoise()
{
   // Add the saved style to the current ROOT session.

   delete gROOT->GetStyle("AdvNoise");

   TStyle *tmpStyle = new TStyle("AdvNoise", "Style for the advanced noise study");
   tmpStyle->SetNdivisions(510, "x");
   tmpStyle->SetNdivisions(510, "y");
   tmpStyle->SetNdivisions(510, "z");
   tmpStyle->SetAxisColor(1, "x");
   tmpStyle->SetAxisColor(1, "y");
   tmpStyle->SetAxisColor(1, "z");
   tmpStyle->SetLabelColor(1, "x");
   tmpStyle->SetLabelColor(1, "y");
   tmpStyle->SetLabelColor(1, "z");
   tmpStyle->SetLabelFont(62, "x");
   tmpStyle->SetLabelFont(62, "y");
   tmpStyle->SetLabelFont(62, "z");
   tmpStyle->SetLabelOffset(0.005, "x");
   tmpStyle->SetLabelOffset(0.005, "y");
   tmpStyle->SetLabelOffset(0.005, "z");
   tmpStyle->SetLabelSize(0.04, "x");
   tmpStyle->SetLabelSize(0.04, "y");
   tmpStyle->SetLabelSize(0.04, "z");
   tmpStyle->SetTickLength(0.03, "x");
   tmpStyle->SetTickLength(0.03, "y");
   tmpStyle->SetTickLength(0.03, "z");
   tmpStyle->SetTitleOffset(1, "x");
   tmpStyle->SetTitleOffset(1, "y");
   tmpStyle->SetTitleOffset(1, "z");
   tmpStyle->SetTitleSize(0.04, "x");
   tmpStyle->SetTitleSize(0.04, "y");
   tmpStyle->SetTitleSize(0.04, "z");
   tmpStyle->SetTitleColor(1, "x");
   tmpStyle->SetTitleColor(1, "y");
   tmpStyle->SetTitleColor(1, "z");
   tmpStyle->SetTitleFont(62, "x");
   tmpStyle->SetTitleFont(62, "y");
   tmpStyle->SetTitleFont(62, "z");
   tmpStyle->SetBarWidth(1000);
   tmpStyle->SetBarOffset(0);
   tmpStyle->SetDrawBorder(0);
   tmpStyle->SetOptLogx(0);
   tmpStyle->SetOptLogy(0);
   tmpStyle->SetOptLogz(0);
   tmpStyle->SetOptDate(0);
   tmpStyle->SetOptStat(1110);
   tmpStyle->SetOptTitle(kTRUE);
   tmpStyle->SetOptFit(0);
   tmpStyle->SetNumberContours(20);
   tmpStyle->GetAttDate()->SetTextFont(62);
   tmpStyle->GetAttDate()->SetTextSize(0.025);
   tmpStyle->GetAttDate()->SetTextAngle(0);
   tmpStyle->GetAttDate()->SetTextAlign(11);
   tmpStyle->GetAttDate()->SetTextColor(1);
   tmpStyle->SetDateX(0.01);
   tmpStyle->SetDateY(0.01);
   tmpStyle->SetEndErrorSize(2);
   tmpStyle->SetErrorX(0.5);
   tmpStyle->SetFuncColor(1);
   tmpStyle->SetFuncStyle(1);
   tmpStyle->SetFuncWidth(3);
   tmpStyle->SetGridColor(0);
   tmpStyle->SetGridStyle(3);
   tmpStyle->SetGridWidth(1);
   tmpStyle->SetLegendBorderSize(4);
   tmpStyle->SetHatchesLineWidth(1);
   tmpStyle->SetHatchesSpacing(1);
   tmpStyle->SetFrameFillColor(19);
   tmpStyle->SetFrameLineColor(1);
   tmpStyle->SetFrameFillStyle(1001);
   tmpStyle->SetFrameLineStyle(1);
   tmpStyle->SetFrameLineWidth(1);
   tmpStyle->SetFrameBorderSize(1);
   tmpStyle->SetFrameBorderMode(1);
   tmpStyle->SetHistFillColor(0);
   tmpStyle->SetHistLineColor(1);
   tmpStyle->SetHistFillStyle(1001);
   tmpStyle->SetHistLineStyle(1);
   tmpStyle->SetHistLineWidth(1);
   tmpStyle->SetHistMinimumZero(kFALSE);
   tmpStyle->SetCanvasPreferGL(kFALSE);
   tmpStyle->SetCanvasColor(19);
   tmpStyle->SetCanvasBorderSize(2);
   tmpStyle->SetCanvasBorderMode(1);
   tmpStyle->SetCanvasDefH(500);
   tmpStyle->SetCanvasDefW(700);
   tmpStyle->SetCanvasDefX(10);
   tmpStyle->SetCanvasDefY(10);
   tmpStyle->SetPadColor(19);
   tmpStyle->SetPadBorderSize(2);
   tmpStyle->SetPadBorderMode(1);
   tmpStyle->SetPadBottomMargin(0.1);
   tmpStyle->SetPadTopMargin(0.1);
   tmpStyle->SetPadLeftMargin(0.1);
   tmpStyle->SetPadRightMargin(0.1);
   tmpStyle->SetPadGridX(kFALSE);
   tmpStyle->SetPadGridY(kFALSE);
   tmpStyle->SetPadTickX(0);
   tmpStyle->SetPadTickY(0);
   tmpStyle->SetPaperSize(20, 26);
   tmpStyle->SetScreenFactor(0.8558);
   tmpStyle->SetStatColor(19);
   tmpStyle->SetStatTextColor(1);
   tmpStyle->SetStatBorderSize(2);
   tmpStyle->SetStatFont(62);
   tmpStyle->SetStatFontSize(0);
   tmpStyle->SetStatStyle(1001);
   tmpStyle->SetStatFormat("6.4g");
   tmpStyle->SetStatX(0.98);
   tmpStyle->SetStatY(0.995);
   tmpStyle->SetStatW(0.4);
   tmpStyle->SetStatH(0.36);
   tmpStyle->SetStripDecimals(kTRUE);
   tmpStyle->SetTitleAlign(13);
   tmpStyle->SetTitleFillColor(19);
   tmpStyle->SetTitleTextColor(1);
   tmpStyle->SetTitleBorderSize(2);
   tmpStyle->SetTitleFont(62);
   tmpStyle->SetTitleFontSize(0);
   tmpStyle->SetTitleStyle(1001);
   tmpStyle->SetTitleX(0.01);
   tmpStyle->SetTitleY(0.995);
   tmpStyle->SetTitleW(0);
   tmpStyle->SetTitleH(0);
   tmpStyle->SetLegoInnerR(0.5);

   Int_t fPaletteColor[50] = {51, 52, 53, 54, 55, 56, 57, 58, 59, 
                             60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 
                             70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
                             80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 
                             90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
   tmpStyle->SetPalette(50, fPaletteColor);

   TString fLineStyleArrayTmp[30] = {"", "    ", "   12 12", "   4 8", 
                             "   12 16 4 16", "   20 12 4 12", "   20 12 4 12 4 12 4 12", "   20 20", "   20 12 4 12 4 12", 
                             "   80 20", "   80 40 4 40", "    ", "    ", "    ", 
                             "    ", "    ", "    ", "    ", "    ", 
                             "    ", "    ", "    ", "    ", "    ", 
                             "    ", "    ", "    ", "    ", "    ", "    "};
   for (Int_t i=0; i<30; i++)
      tmpStyle->SetLineStyleString(i, fLineStyleArrayTmp[i]);

   tmpStyle->SetHeaderPS("");
   tmpStyle->SetTitlePS("");
   tmpStyle->SetFitFormat("5.4g");
   tmpStyle->SetPaintTextFormat("g");
   tmpStyle->SetLineScalePS(3);
   tmpStyle->SetColorModelPS(0);
   tmpStyle->SetTimeOffset(788918400);

   tmpStyle->SetLineColor(1);
   tmpStyle->SetLineStyle(1);
   tmpStyle->SetLineWidth(1);
   tmpStyle->SetFillColor(19);
   tmpStyle->SetFillStyle(1001);
   tmpStyle->SetMarkerColor(1);
   tmpStyle->SetMarkerSize(1);
   tmpStyle->SetMarkerStyle(1);
   tmpStyle->SetTextAlign(11);
   tmpStyle->SetTextAngle(0);
   tmpStyle->SetTextColor(1);
   tmpStyle->SetTextFont(62);
   tmpStyle->SetTextSize(0);



   gStyle = tmpStyle;
}

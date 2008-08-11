{

  gStyle->SetPalette(1 , 0);
  
  gROOT->ProcessLine(".x Style_AdvNoise.C");
  gROOT->ProcessLine(".L advancedNoiseAnalysis.C+");
  cout << "Loading  advancedNoiseAnalysis.C " << endl;

}

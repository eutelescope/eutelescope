{

  gStyle->SetPalette(1 , 0);
  
  gROOT->ProcessLine(".L correlationAnalysis.C+");
  cout << "Loading  correlationAnalysis.C " << endl;

}

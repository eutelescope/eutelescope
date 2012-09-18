/*
 Macro to create the files used in example.
 The file contains two histograms and a TTree
 usage from ROOT:
 .L createExampleFile.C
 createFile( "testfile.root" , 123456 )
 createFile( "referencefile.root" , 7890 )
 
 */

void createTree(const char* fn , Int_t seed) {
  TFile* f = TFile::Open(fn,"update");
  f->cd();
  Double_t aD;
  Float_t aF;
  Int_t aI;
  Double_t aVD[5];
  vector<Double_t> aSD(10);
  vector<Double_t>* pSD = &aSD;
  TTree *t = new TTree("TestTree","TestTree");
  t->Branch("aDouble",&aD,"aDouble/D");
  t->Branch("aFloat",&aF,"aFloat/F");
  t->Branch("aInt",&aI,"aInt/I");
  t->Branch("aArrayD",&(aVD[0]),"aArrayD[5]/D");
  t->Branch("aStdVecD","vector<Double_t>",&pSD);

  TRandom rn(seed);
  for ( Int_t evt = 0 ; evt < 1000 ; ++ evt ) {
    aD = rn.Gaus();
    aF = rn.Uniform();
    aI = rn.Integer(100);
    for ( Int_t i = 0 ; i<10 ; ++i ){
      if ( i < 5 ) aVD[i] = rn.Gaus( i , 0.1 );
      aSD[i] = rn.Gaus( -i , 0.1 );
    }
    t->Fill();
  }
  t->Print();
  t->Write();
  f->Close();
}

void createHistos(const char* fn) {
    TFile* f = TFile::Open(fn,"update");
    TH1F* h2 = new TH1F("h2","Random Guas -2",100,-4,4);
    h2->FillRandom("gaus");
    h2->Write();
    f->mkdir("ADirectory");
    f->cd("ADirectory");
    TH1F* h1= new TH1F("h1","Random Gaus",100,-10,10);
    h1->FillRandom("gaus");
    h1->Write();
    TH1F* h3 = new TH1F("h3","Random Gaus -3",100,0,4);
    h3->FillRandom("gaus");
    h2->Write();
    f->Close();
}

void createFile( const char* fn = "test.root" , Int_t seed = 123456) {
    createTree(fn,seed);
    createHistos(fn);
}

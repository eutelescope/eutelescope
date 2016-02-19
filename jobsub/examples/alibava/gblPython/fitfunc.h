Double_t fit_func(Double_t *x, Double_t *par){

  const Int_t    n     = 50; // Number of integration step
  const Double_t xmin = x[0] - 5.0*par[3];
  const Double_t xmax = x[0] + 5.0*par[3];
  const Double_t h    = (xmax - xmin)/(n - 1);

  Double_t sum = 0;
  Double_t x_  = xmin;

  Double_t Lari_sum = 0;
  Double_t Gaus_sum = 0;

  for (Int_t i=1; i<n; i++)
    {
      Lari_sum =  par[1]*TMath::Abs(TMath::Tan(x_*pi/180) - TMath::Tan(par[0]*pi/180)) + par[2];
      Gaus_sum =  1/TMath::Sqrt(2*pi)/(par[3])*TMath::Exp(-(x[0] - x_ )*(x[0] - x_ )/(2*par[3]*par[3]))*h;
      sum += Lari_sum*Gaus_sum;
      x_ += h;
    }
  return   sum;
}


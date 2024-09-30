#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
 
void fft()
{
  // Histograms
  // =========
  //prepare the canvas for drawing
  Int_t k=2000;
  Int_t n=4000;
  TCanvas *canvas1 = new TCanvas("canvas1", "Fast Fourier Transform", 1500, 1000);
  //! Function
  TF1* fLogis = new TF1("fLogis", "1./(1.+exp(1.*(x-20.0)))", 0, n);
  fLogis->SetNpx(10000);
  fLogis->Draw();

  TCanvas *canvas = new TCanvas("canvas", "Fast Fourier Transform", 1500, 1000);
  canvas->Divide(2, 3);

  //! Example function and noise
  canvas->cd(1);

  //A function to sample
  TF1 *fcos = new TF1("fcos", "1.0*cos(1.0*x)+0.1*cos(10*x)+0.01*cos(100*x)", -4*TMath::Pi(), 4*TMath::Pi());
  // fcos->Draw();

  TH1D *hcos = new TH1D("hcos", "hcos", k, 0., 4*TMath::Pi());
   
  Double_t x;
  //Fill the histogram with function values
  for (Int_t i=1; i<=k; i++){
    x = (Double_t(i)/k)*(4*TMath::Pi());
    hcos->SetBinContent(i, fcos->Eval(x));
  }

  TH1D* hcoseven = new TH1D("hcoseven", "hcoseven", n, -4*TMath::Pi(), 4*TMath::Pi());

  double binContent;
  for (int i=1; i<=n/2; i++)
  {
    binContent = hcos->GetBinContent(i);
    hcoseven->SetBinContent(i, binContent);
    hcoseven->SetBinContent(n - i, binContent);
  }
  hcoseven->Draw();

  hcos->SetLineColor(kRed);
  hcos->SetMarkerStyle(4);
  hcos->Draw("same");

  //! Magnitude
  canvas->cd(2);
  //Compute the transform and look at the magnitude of the output
  TH1 *hm =nullptr;
  TVirtualFFT::SetTransform(nullptr);
  hm = hcoseven->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the transform");
  hm->Draw();
  
  //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
  //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
  TH1D* newhm = new TH1D("newhm", "newhm", n, -4*TMath::Pi(), 4*TMath::Pi());
  for(int i = 1; i <= n; i++)
  {
    newhm->SetBinContent(i, (hm->GetBinContent(i))/n);
  }
  // newhm->Draw();

  TH1D* hLogis = new TH1D("hLogis", "Logis", n, 0, n);
  for(int i = 1; i <= n/2; i++)
  {
    hLogis->SetBinContent(i, fLogis->Eval(i));
    hLogis->SetBinContent(n-i, fLogis->Eval(i));
  }

  hLogis->SetLineColor(kRed);
  hLogis->Draw("same");

  //! Real part
  canvas->cd(3);
  TH1 *hr =nullptr;
  hr = hcoseven->FFT(hr, "RE");
  hr->SetTitle("Real part of the transform");
  hr->Draw();

  TH1D* newhr = new TH1D("newhr", "newhr", n, -4*TMath::Pi(), 4*TMath::Pi());
  for(int i = 1; i <= n; i++)
  {
    newhr->SetBinContent(i, hr->GetBinContent(i)/n);
  }
  // newhr->Draw();

  //! Imaginary part
  canvas->cd(4);
  TH1 *hi =nullptr;
  hi = hcoseven->FFT(hi, "IM");
  hi->SetTitle("Imaginary part of the transform");
  hi->Draw();

  TH1D* newhi = new TH1D("newhi", "newhi", n, -4*TMath::Pi(), 4*TMath::Pi());
  for(int i = 1; i <= n; i++)
  {
    newhi->SetBinContent(i, hi->GetBinContent(i)/n);
  }
  // newhi->Draw();

  //! Phase
  canvas->cd(5);
  //Look at the phase of the output
  TH1 *hp = nullptr;
  hp = hcoseven->FFT(hp, "PH");
  hp->SetTitle("Phase of the transform");
  hp->Draw();
  TH1D* newhp = new TH1D("newhp", "newhp", n, -4*TMath::Pi(), 4*TMath::Pi());
  for(int i = 1; i <= n; i++)
  {
    newhp->SetBinContent(i, hp->GetBinContent(i)/n);
  }
  // newhp->Draw();

  //! Apply threshold
  canvas->cd(6);
  Double_t re, im;
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();  
  //Use the following method to get the full output:
  Double_t *re_full = new Double_t[n];
  Double_t *im_full = new Double_t[n];
 
  fft->GetPointsComplex(re_full,im_full);

  for(int i = 1; i <= n; i++)
  {
    std::cout << "old real part [" << i << "] = " << re_full[i] << "\n";
    std::cout << "old im part [" << i << "] = " << im_full[i] << "\n";
    re_full[i] = re_full[i]*hLogis->GetBinContent(i);
    im_full[i] = im_full[i]*hLogis->GetBinContent(i);
    std::cout << "new real part [" << i << "] = " << re_full[i] << "\n";
    std::cout << "new im part [" << i << "] = " << im_full[i] << "\n";
  }

  //Now let's make a backward transform:
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = nullptr;
  //Let's look at the output
  hb = TH1::TransformHisto(fft_back,hb,"Re");

  TH1D* newhb = new TH1D("newhb", "newhb", n, -4*TMath::Pi(), 4*TMath::Pi());
  for(int i = 1; i <= n; i++)
  {
    newhb->SetBinContent(i, (hb->GetBinContent(i))/n);
  }
  newhb->SetTitle("The backward transform result");
  newhb->Draw();

  // hb->Draw();

  //NOTE: here you get at the x-axes number of bins and not real values
  //in this case 25 bins has to be rescaled to a range between 0 and 4*Pi;
  //also here the y-axes has to be rescaled (factor 1/bins)
  canvas->cd(6);

  delete fft_back;
  fft_back=nullptr;

  delete [] re_full;
  delete [] im_full;
}
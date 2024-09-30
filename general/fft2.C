#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"
#include "TF1.h"

void fft2()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Canvas 1
  TCanvas *canvas1 = new TCanvas("canvas1", "Exp data show", 1500, 1000);
  canvas1->cd();
  canvas1->Divide(4, 2);

  //! Pad1: Measurement
  canvas1->cd(1);
  TFile* fileOpen = new TFile("save_new_passive_divider_Na22_no_coinc_ch0HV1925_ch2HV1550_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");

  TH1D* hChannel = new TH1D();
  hChannel = fileOpen->Get<TH1D>("hCh0");
  hChannel->SetLineColor(kRed);
  hChannel->Draw();

  //! Pad2: Derivative
  canvas1->cd(2);
  int binExp = 700;
  int binExpEven = 2*binExp ;
  double xminExp = 100;
  double xmaxExp = 800;

  TH1D* hDiffChannel = new TH1D("hDiffChannel", "Exp data Diff", binExp, xminExp, xmaxExp);

  double binContent, binContentPlus, binLength, newBinContent;
  binLength = (xmaxExp - xminExp)/ binExp;

  // for(int i = 1; i <= binExp; i++)
  // {
  //   binContent = hChannel->GetBinContent(i);
  //   binContentPlus = hChannel->GetBinContent(i+1);
  //   newBinContent = (binContentPlus - binContent)/binLength;
  //   hDiffChannel->SetBinContent(i, newBinContent);
  // }

  for(int i = 3; i <= binExp - 2; i++)
  {
    newBinContent = (-hChannel->GetBinContent(i+2) 
    + 8*hChannel->GetBinContent(i+1) 
    - 8*hChannel->GetBinContent(i-1) 
    + hChannel->GetBinContent(i-2))/12*binLength;
    hDiffChannel->SetBinContent(i, newBinContent);
  }

  hDiffChannel->Draw();

  //! Pad3: Even function
  canvas1->cd(3);
  TH1D* hEvenChannel = new TH1D("hEvenChannel", "Transformed to even function", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);

  for(int i = 1; i <= binExpEven/2; i++)
  {
    // binContent = hChannel->GetBinContent(i);
    binContent = hDiffChannel->GetBinContent(i);
    hEvenChannel->SetBinContent(i, binContent);
    hEvenChannel->SetBinContent(binExpEven - i, binContent);
  }
  hEvenChannel->SetBinContent(binExp, hEvenChannel->GetBinContent(binExp + 1));
  hEvenChannel->Draw();

  //! Pad4: Magnitude
  canvas1->cd(4);
  double k = 0.1;
  double c = 30.0;

  TF1* fLogis = new TF1("fLogis", "1./(1.+exp([0]*(x-[1])))", 0, binExpEven);
  fLogis->SetParameter(0, k);
  fLogis->SetParameter(1, c);
  fLogis->SetNpx(10000);
  // fLogis->Draw();
  //Compute the transform and look at the magnitude of the output
  TH1 *hm =nullptr;
  TVirtualFFT::SetTransform(nullptr);
  hm = hEvenChannel->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the transform");
  hm->Draw();
  
  //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
  //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
  TH1D* newhm = new TH1D("newhm", "newhm", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    newhm->SetBinContent(i, hm->GetBinContent(i)/binExpEven);
  }
  // newhm->Draw();

  TH1D* hLogis = new TH1D("hLogis", "Logis", binExpEven, 0, binExpEven);
  for(int i = 1; i <= binExpEven/2; i++)
  {
    hLogis->SetBinContent(i, fLogis->Eval(i));
    hLogis->SetBinContent(binExpEven-i, fLogis->Eval(i));
  }

  hLogis->SetLineColor(kRed);
  // hLogis->Scale(1000000, "noSW2");
  hLogis->Draw("same");

  //! Pad5: Real part
  canvas1->cd(5);
  TH1 *hr =nullptr;
  hr = hEvenChannel->FFT(hr, "RE");
  hr->SetTitle("Real part of the transform");
  hr->Draw();

  TH1D* newhr = new TH1D("newhr", "newhr", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    newhr->SetBinContent(i, hr->GetBinContent(i)/binExpEven);
  }
  // newhr->Draw();

  //! Pad6: Imaginary part
  canvas1->cd(6);
  TH1 *hi =nullptr;
  hi = hEvenChannel->FFT(hi, "IM");
  hi->SetTitle("Imaginary part of the transform");
  hi->Draw();

  TH1D* newhi = new TH1D("newhi", "newhi", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    newhi->SetBinContent(i, hi->GetBinContent(i)/binExpEven);
  }
  // newhi->Draw();

  //! Pad7: Phase
  canvas1->cd(7);
  //Look at the phase of the output
  TH1 *hp = nullptr;
  hp = hEvenChannel->FFT(hp, "PH");
  hp->SetTitle("Phase of the transform");
  hp->Draw();

  TH1D* newhp = new TH1D("newhp", "newhp", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    newhp->SetBinContent(i, hp->GetBinContent(i)/binExpEven);
  }
  // newhp->Draw();

  //! Pad8: Apply threshold
  canvas1->cd(8);

  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();  
  //Use the following method to get the full output:
  Double_t *re_full = new Double_t[binExpEven];
  Double_t *im_full = new Double_t[binExpEven];
 
  fft->GetPointsComplex(re_full, im_full);

  for(int i = 1; i <= binExpEven; i++)
  {
    re_full[i] = re_full[i]*hLogis->GetBinContent(i);
    im_full[i] = im_full[i]*hLogis->GetBinContent(i);
  }

  //Now let's make a backward transform:
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &binExpEven, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = nullptr;
  //Let's look at the output
  hb = TH1::TransformHisto(fft_back,hb,"RE");

  TH1D* newhb = new TH1D("newhb", "newhb", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    // std::cout << "value = " << hb->GetBinContent(i)/binExpEven << "\n";
    newhb->SetBinContent(i, hb->GetBinContent(i)/binExpEven);
  }
  newhb->SetTitle("The backward transform result");
  newhb->Draw();

  //! Canvas 2
  TCanvas* canvas2 = new TCanvas("canvas2", "Compare", 800, 600);
  canvas2->Divide(2, 1);

  //! Pad 1: Measurement
  canvas2->cd(1);
  hChannel->Draw();

  //! Pad 2: Filtered
  canvas2->cd(2);
  TH1D* hFiltered = new TH1D("hFiltered", "Exp data Diff filtered", binExp, xminExp, xmaxExp);
  for(int i = 1; i <= binExp; i++)
  {
    hFiltered->SetBinContent(i, newhb->GetBinContent(i));
  }
  // hDiffChannel->Draw();

  hFiltered->SetLineColor(kBlue);
  hFiltered->Draw();

  TF1* fit = new TF1("fit", "[0]*x+[1]", 280, 380);
  hFiltered->Fit("fit", "", "", 280, 380);

  // double sum = 0.;
  // int first = 280;
  // int last = 360;
  // for (int i = first; i <= last; i++)
  // {
  //   sum += atan(hFiltered->GetBinContent(i));
  // }
  // std::cout << "average angle = " << sum/(last - first + 1) << "\n";
  // std::cout << "last to first angle = " << (hFiltered->GetBinContent(last) - hFiltered->GetBinContent(first))/(last - first + 1) << "\n";

  std::cout << "time: " << timer->RealTime() << " seconds \n";
}
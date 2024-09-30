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
#include "TGraph.h"
#include "TMarker.h"

void fitBG(){
  TFile* file1 = new TFile("./expfiles/old/save_old_divider_bg_no_coinc_ch2HV1746_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
  TH1D* hBg = new TH1D();
  hBg = file1->Get<TH1D>("hCh1");
  
  int binExp = 700;
  int binExpEven = 2*binExp ;
  double xminExp = 100;
  double xmaxExp = 800;

  double binContent, binContentPlus, binLength, newBinContent;
  binLength = (xmaxExp - xminExp)/ binExp;

  TH1D* hEvenChannel = new TH1D("hEvenChannel", "Transformed to even function", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);

  for(int i = 1; i <= binExpEven/2; i++)
  {
    binContent = hBg->GetBinContent(i);
    // binContent = hDiffChannel->GetBinContent(i);
    hEvenChannel->SetBinContent(i, binContent);
    hEvenChannel->SetBinContent(binExpEven - i, binContent);
  }
  hEvenChannel->SetBinContent(binExp, hEvenChannel->GetBinContent(binExp + 1));

  double para_k = 0.1;
  double para_c = 500.0;

  TF1* fLogis = new TF1("fLogis", "1./(1.+exp([0]*(x-[1])))", 0, binExpEven);
  fLogis->SetParameter(0, para_k);
  fLogis->SetParameter(1, para_c);
  fLogis->SetNpx(10000);

  TH1D* hLogis = new TH1D("hLogis", "Logis", binExpEven, 0, binExpEven);
  for(int i = 1; i <= binExpEven/2; i++)
  {
    hLogis->SetBinContent(i, fLogis->Eval(i));
    hLogis->SetBinContent(binExpEven-i, fLogis->Eval(i));
  }

  TH1 *hm =nullptr;
  TVirtualFFT::SetTransform(nullptr);
  hm = hEvenChannel->FFT(hm, "MAG");
  
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

  TH1D* hFiltered = new TH1D("hFiltered", "Exp data Diff filtered", binExp, xminExp, xmaxExp);
  for(int i = 1; i <= binExp; i++)
  {
    hFiltered->SetBinContent(i, newhb->GetBinContent(i));
  }

  double xminCal = -16.1112;
  double xmaxCal = 2203.38;
  
  TH1D* hCal = new TH1D("hCal", "Calibrated", binExp, xminCal, xmaxCal);
  for(int i = 1; i <= binExp; i++){
    hCal->SetBinContent(i, hFiltered->GetBinContent(i));
  }
  hCal->SetXTitle("Energy(keV)");
  hCal->SetYTitle("Count");
  
  TCanvas* canvas1 = new TCanvas("canvas1", "Fit", 800, 600);
  // TF1* fit = new TF1("fit", "[0]*pow(x,3) + [1]*x*x + [2]*x + [3]");
  TF1* fit = new TF1("fit", "[3]*x*x*x + [0]*x*x + [1]*x + [2]");
  hCal->Fit("fit", "", "", 250., 700.);
}
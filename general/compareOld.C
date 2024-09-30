#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TVirtualFFT.h"

//! histogram options
double xminSim = 0.;
double xmaxSim = 1.5;
int xminExp = 100;
int xmaxExp = 800;
int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;

//! ENERGY CALIBRATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ECs137 = 0.477;
double ENa22 = 1.061;
double ECo60 = 1.117;

//? Start
int ChCs137[3] = {268, 221, 285};
int ChNa22[3] = {468, 368, 520};
int ChCo60[3] = {474, 375, 514};

//? End
// int ChCs137[3] = {275, 218, 286};
// int ChNa22[3] = {477, 360, 515};
// int ChCo60[3] = {489, 349, 506};

//! RESOLUTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//? Start
double deltaSigCs137[3] = {8.40, 806, 7.82};
double deltaSigNa22[3] = {5.53, 5.48, 5.32};
double deltaSigCo60[3] = {4.63, 4.60, 4.66};

//? End
// double deltaSigCs137[3] = {8.10, 8.83, 7.82};
// double deltaSigNa22[3] = {5.53, 6.53, 5.17};
// double deltaSigCo60[3] = {4.74, 5.12, 4.32};

double ScaleCs137[3] = {5.5, 4.0, 2.9};
double ScaleNa22[3] = {2.0, 1.40, 1.65};
double ScaleCo60[3] = {2.0, 1.45, 1.65};
//! CHANNEL
int channel = 0;
// int channel = 1;
// int channel = 2;

const char* TH1namechar;
const char* TH1titlechar;
TH1D* fft(TH1D* hChannel, double para_k, double para_c, std::string TH1name, std::string TH1title, int bin, double xmin, double xmax)
{
  int binEven = 2*bin;
  TH1namechar = TH1name.c_str();
  TH1titlechar = TH1title.c_str();

  TH1D* hFiltered = new TH1D(TH1namechar, TH1titlechar, bin, xmin, xmax);
  TH1D* hEvenChannel = new TH1D("hEvenChannel", "Transformed to even function", binEven, -(xmax - xmin), xmax - xmin);

  for(int i = 1; i <= binEven/2; i++)
  {
    hEvenChannel->SetBinContent(i, hChannel->GetBinContent(i));
    hEvenChannel->SetBinContent(binEven - i, hChannel->GetBinContent(i));
  }
  hEvenChannel->SetBinContent(binExp, hEvenChannel->GetBinContent(binExp + 1));

  TF1* fLogis = new TF1("fLogis", "1./(1.+exp([0]*(x-[1])))", 0, binEven);
  fLogis->SetParameter(0, para_k);
  fLogis->SetParameter(1, para_c);
  fLogis->SetNpx(10000);

  //! Magnitude
  TH1 *hm = nullptr;
  TVirtualFFT::SetTransform(nullptr);
  hm = hEvenChannel->FFT(hm, "MAG");
  
  TH1D* newhm = new TH1D("newhm", "newhm", binEven, -(xmax - xmin), xmax - xmin);
  for(int i = 1; i <= binEven; i++)
  {
    newhm->SetBinContent(i, hm->GetBinContent(i)/binEven);
  }

  TH1D* hLogis = new TH1D("hLogis", "Logis", binEven, 0, binEven);
  for(int i = 1; i <= binEven/2; i++)
  {
    hLogis->SetBinContent(i, fLogis->Eval(i));
    hLogis->SetBinContent(binEven-i, fLogis->Eval(i));
  }

  //! Apply threshold
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  Double_t *re_full = new Double_t[binEven];
  Double_t *im_full = new Double_t[binEven];
 
  fft->GetPointsComplex(re_full, im_full);

  for(int i = 1; i <= binEven; i++)
  {
    re_full[i] = re_full[i]*hLogis->GetBinContent(i);
    im_full[i] = im_full[i]*hLogis->GetBinContent(i);
  }

  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &binEven, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = nullptr;
  hb = TH1::TransformHisto(fft_back,hb,"RE");

  TH1D* newhb = new TH1D("newhb", "newhb", binEven, -(xmax - xmin), xmax - xmin);
  for(int i = 1; i <= binEven; i++)
  {
    newhb->SetBinContent(i, hb->GetBinContent(i)/binEven);
  }

  for(int i = 1; i <= bin; i++)
  {
    hFiltered->SetBinContent(i, newhb->GetBinContent(i));
  }

  delete hEvenChannel;
  delete hm;
  delete newhm;
  delete fLogis;
  delete hLogis;
  delete[] re_full;
  delete[] im_full;
  delete fft_back;
  delete hb;
  delete newhb;

  return hFiltered;
}

TH1D* Diff(TH1D* hChannel, std::string TH1name, std::string TH1title, int bin, double xmin, double xmax)
{
  TH1namechar = TH1name.c_str();
  TH1titlechar = TH1title.c_str();
  double newBinContent;
  double binLength = (xmax-xmin)/bin;

  TH1D* hDiff = new TH1D(TH1namechar, TH1titlechar, bin, xmin, xmax);
  for(int i = 3; i <= binExp - 2; i++)
  {
    newBinContent = (-hChannel->GetBinContent(i+2) 
    + 8*hChannel->GetBinContent(i+1) 
    - 8*hChannel->GetBinContent(i-1) 
    + hChannel->GetBinContent(i-2))/12*binLength;
    hDiff->SetBinContent(i, newBinContent);
  }
  return hDiff;
}

void compareOld()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* fsimCs137 = new TFile("./simfiles/withAl_simCs137_1.root", "read");
  // TFile* fexpCs137 = new TFile("./expfiles/time/plastic_detectors_Cs137_start_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");
  TFile* fexpCs137 = new TFile("./expfiles/time/plastic_detectors_Cs137_end_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");

  TTree* tsimCs137 =  (TTree*) fsimCs137->Get("dEEtree");
  TTree* texpCs137 =  (TTree*) fexpCs137->Get("AnalysisxTree");

  TFile* fsimNa22 = new TFile("./simfiles/withAl_simNa22_2.root", "read");
  // TFile* fexpNa22 = new TFile("./expfiles/time/plastic_detectors_Na22_start_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");
  TFile* fexpNa22 = new TFile("./expfiles/time/plastic_detectors_Na22_end_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");

  TTree* tsimNa22 =  (TTree*) fsimNa22->Get("dEEtree");
  TTree* texpNa22 =  (TTree*) fexpNa22->Get("AnalysisxTree");

  TFile* fsimCo60 = new TFile("./simfiles/withAl_simCo60_1.root", "read");
  // TFile* fexpCo60 = new TFile("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");
  TFile* fexpCo60 = new TFile("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0008.root", "read");

  TTree* tsimCo60 =  (TTree*) fsimCo60->Get("dEEtree");
  TTree* texpCo60 =  (TTree*) fexpCo60->Get("AnalysisxTree");

  TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
  hCali->SetBinContent(ChNa22[channel] - xminExp, ENa22);
  hCali->SetBinContent(ChCo60[channel] - xminExp, ECo60);
  hCali->SetBinContent(ChCs137[channel] - xminExp, ECs137);

  TF1* fLinear = new TF1("fLinear", "[0]*x + [1]", xminExp, xmaxExp);

  TCanvas* cFit = new TCanvas("cFit", "Linear calibration", 800, 600);
  cFit->Divide(2,1);
  cFit->cd(1);
  std::cout << "\nENERGY CALIBRATION\n";
  hCali->Fit("fLinear");

  double a = fLinear->GetParameter(0);
  double b = fLinear->GetParameter(1);

  double xminCal = xminExp*a+b;
  double xmaxCal = xmaxExp*a+b;

  TH1D* hsimCs137res = new TH1D("hsimCs137res", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
  TH1D* hcalCs137 = new TH1D("hcalCs137", "Cs137 Experiment Calibrated", binExp, xminCal, xmaxCal);

  hsimCs137res->GetXaxis()->SetTitle("Energy(MeV)");
  hsimCs137res->GetYaxis()->SetTitle("Count");

  hcalCs137->GetXaxis()->SetTitle("Energy(MeV)");
  hcalCs137->GetYaxis()->SetTitle("Count");

  TH1D* hsimNa22res = new TH1D("hsimNa22res", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);
  TH1D* hcalNa22 = new TH1D("hcalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);

  hsimNa22res->GetXaxis()->SetTitle("Energy(MeV)");
  hsimNa22res->GetYaxis()->SetTitle("Count");

  hcalNa22->GetXaxis()->SetTitle("Energy(MeV)");
  hcalNa22->GetYaxis()->SetTitle("Count");

  TH1D* hsimCo60res = new TH1D("hsimCo60res", "Co60 Simulation with resolution", binExp, xminCal, xmaxCal);
  TH1D* hcalCo60 = new TH1D("hcalCo60", "Co60 Experiment Calibrated", binExp, xminCal, xmaxCal);

  hsimCo60res->GetXaxis()->SetTitle("Energy(MeV)");
  hsimCo60res->GetYaxis()->SetTitle("Count");

  hcalCo60->GetXaxis()->SetTitle("Energy(MeV)");
  hcalCo60->GetYaxis()->SetTitle("Count");
  //! Fill experiment histograms
  UShort_t chX[48];

  texpCs137->SetBranchAddress("NeEvent.neutAmp[48]", chX);
  
  TH1D* hexpCs137 = new TH1D("hexpCs137", "Cs137 Measurement", xmaxExp-xminExp, xminExp, xmaxExp);
  Long64_t entries = texpCs137->GetEntries();

  for(int i = 0; i < entries; i++){
    texpCs137->GetEntry(i);
    hcalCs137->Fill(a*(chX[channel]+0.5) + b);
    hexpCs137->Fill(chX[channel]);
  }

  UShort_t chY[48];

  texpNa22->SetBranchAddress("NeEvent.neutAmp[48]", chY);
  
  entries = texpNa22->GetEntries();

  for(int i = 0; i < entries; i++){
    texpNa22->GetEntry(i);
    hcalNa22->Fill(a*(chY[channel]+0.5) + b);
  }

  UShort_t chZ[48];

  texpCo60->SetBranchAddress("NeEvent.neutAmp[48]", chZ);
  
  entries = texpCo60->GetEntries();

  for(int i = 0; i < entries; i++){
    texpCo60->GetEntry(i);
    hcalCo60->Fill(a*(chZ[channel]+0.5) + b);
  }
  //! Resolution for simulation histograms
  std::cout << "\nENERGY RESOLUTION\n";
  // std::cout << "delta1 = " << deltaSig1*100 << "%; delta2 = " << deltaSig2*100 << "%\n";

  TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal, xmaxCal);
  hRes->SetBinContent(hcalCs137->FindBin(ECs137), deltaSigCs137[channel]/100.);
  hRes->SetBinContent(hcalNa22->FindBin(ENa22), deltaSigNa22[channel]/100.);
  hRes->SetBinContent(hcalCo60->FindBin(ECo60), deltaSigCo60[channel]/100.);

  //! Choose fit function
  TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");

  cFit->cd(2);
  hRes->Fit("fRes", "", "", ECs137-0.2, ECo60+0.2);

  //! Correspond to chosen function
  double coA = fRes->GetParameter(0);
  double coB = 0.;
  double coC = fRes->GetParameter(1);

  // coA = 0.;
  // coB = 0.;
  // coC = 0.;

  double x_sim;

  tsimCs137->SetBranchAddress("Scintillator", &x_sim);

  TRandom3* ranGen = new TRandom3();

  entries = tsimCs137->GetEntries();

  // for (int i = 0; i < entries; i++)
  for (int i = 0; i < 1000000; i++)
  {
    tsimCs137->GetEntry(i);
    double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
    hsimCs137res->Fill(ranGen->Gaus(x_sim,sigma));
  }
  
  double y_sim;

  tsimNa22->SetBranchAddress("Scintillator", &y_sim);

  entries = tsimNa22->GetEntries();

  for (int i = 0; i < entries; i++)
  {
    tsimNa22->GetEntry(i);
    double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
    hsimNa22res->Fill(ranGen->Gaus(y_sim,sigma));
  }

  double z_sim;

  tsimCo60->SetBranchAddress("Scintillator", &z_sim);

  entries = tsimCo60->GetEntries();

  for (int i = 0; i < entries; i++)
  {
    tsimCo60->GetEntry(i);
    double sigma = z_sim*sqrt( pow(coA,2) + pow(coB/sqrt(z_sim),2) + pow(coC/z_sim,2) );
    hsimCo60res->Fill(ranGen->Gaus(z_sim,sigma));
  }
  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  c1->cd();
  c1->Divide(1, 3);

  c1->cd(1);
  hcalCs137->SetLineColor(kRed);
  hcalCs137->Draw();
  hsimCs137res->Scale(ScaleCs137[channel], "noSW2");
  hsimCs137res->Draw("same");

  TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg1->SetHeader("Cs137", "C");
  leg1->SetBorderSize(2);
  leg1->AddEntry(hsimCs137res, "simulation", "l");
  leg1->AddEntry(hcalCs137, "experiment", "l");
  leg1->Draw();

  c1->cd(2);
  hcalNa22->SetLineColor(kRed);
  hcalNa22->Draw();
  hsimNa22res->Scale(ScaleNa22[channel], "noSW2");
  hsimNa22res->Draw("same");

  TLegend *leg2 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg2->SetHeader("Na22", "C");
  leg2->SetBorderSize(2);
  leg2->AddEntry(hsimNa22res, "simulation", "l");
  leg2->AddEntry(hcalNa22, "experiment", "l");
  leg2->Draw();

  c1->cd(3);
  hcalCo60->SetLineColor(kRed);
  hcalCo60->Draw();
  hsimCo60res->Scale(ScaleCo60[channel], "noSW2");
  hsimCo60res->Draw("same");

  TLegend *leg3 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg3->SetHeader("Co60", "C");
  leg3->SetBorderSize(2);
  leg3->AddEntry(hsimCo60res, "simulation", "l");
  leg3->AddEntry(hcalCo60, "experiment", "l");
  leg3->Draw();

  TCanvas* CanvasSimAndExp = new TCanvas("CanvasSimAndExp", "CanvasSimAndExp", 800, 600);
  CanvasSimAndExp->Divide(1,2);
  CanvasSimAndExp->cd(1)->SetGridx();
  CanvasSimAndExp->cd(1)->SetGridy();
  CanvasSimAndExp->cd(1);
  hexpCs137->SetLineColor(kRed);
  hexpCs137->Draw();
  CanvasSimAndExp->cd(2)->SetGridx();
  CanvasSimAndExp->cd(2)->SetGridy();
  CanvasSimAndExp->cd(2);
  hsimCs137res->Draw();

  TCanvas* CanvasFFT = new TCanvas("CanvasFFT", "CanvasFFT", 800, 600);
  TH1D* hexpCs137Filtered = fft(hexpCs137, 0.1, 100., "hexpCs137Filtered", "Cs137 Measurement Filtered", hexpCs137->GetNbinsX(), xminExp, xmaxExp);
  CanvasFFT->Divide(1,2);
  CanvasFFT->cd(1);
  hexpCs137->Draw();
  CanvasFFT->cd(2);
  hexpCs137Filtered->Draw();

  TCanvas* CanvasDiff = new TCanvas("CanvasDiff", "CanvasDiff", 800, 600);
  TH1D* hexpCs137Diff = Diff(hexpCs137Filtered, "hexpCs137Diff", "Differentiated Cs137 Spectrum", hexpCs137->GetNbinsX(), xminExp, xmaxExp);
  CanvasDiff->cd();
  hexpCs137Diff->Draw();

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
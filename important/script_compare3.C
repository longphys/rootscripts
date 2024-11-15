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
#include "TGraph.h"
#include "TVirtualFFT.h"

//! histogram options
double xminSim = 0.;
double xmaxSim = 1.3;
double xminExp = 100.;
double xmaxExp = 1500.;
int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;

double E1 = 0.477;
double E2 = 1.061;

//! Zoom for GRADIENT DESCENT
double minTestCs137 = 0.35; //MeV
double maxTestCs137 = 0.7;
double minTestNa22 = 0.8; //MeV
double maxTestNa22 = 1.3;

//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//! ZOOM FOR FFT
int minRangeCs137 = 200; int maxRangeCs137 = 1500;
int minRangeNa22 = 400; int maxRangeNa22 = 1500;
// int minRangeCo60 = 300; int maxRangeCo60 = 1500;

//! FFT Threshold FOR INITIAL VALUES (30-1500)
//? To increase FFT effect, lower the threshold
int inithreshCs137 = 100.;
int inithreshNa22 = 100.; 

//! SIMULATION AND MEASUREMENT FILES
TFile* fsimCs137 = new TFile("../../simfiles/withAl_simCs137_1.root", "read");
// TFile* fexpCs137 = new TFile("./expfiles/time/4/plastic_detectors_Cs137_ch0HV1787_ch1HV1821_ch2HV1850_ch0_9849_ch1_9759_ch2_9825_run_0_0001.root", "read");
// TFile* fexpCs137 = new TFile("./expfiles/time/5/plastic_detectors_Cs137_ch0HV1671_ch1HV1885_ch2HV1800_ch0_9837_ch1_9729_ch2_9419_run_0_0001.root", "read");

TFile* fsimNa22 = new TFile("../../simfiles/withAl_simNa22_1.root", "read");
// TFile* fexpNa22 = new TFile("./expfiles/time/4/plastic_detectors_Na22_ch0HV1787_ch1HV1821_ch2HV1850_ch0_9849_ch1_9759_ch2_9825_run_0_0001.root", "read");
// TFile* fexpNa22 = new TFile("./expfiles/time/5/plastic_detectors_Na22_ch0HV1671_ch1HV1885_ch2HV1800_ch0_9837_ch1_9729_ch2_9419_run_0_0001.root", "read");

//!Neutron measurement calibration
TFile* fexpCs137 = new TFile("../../EfficiencyPlasticFLNR/declib/output/test12_2_Plast_BC404_Cs137.root","read");
TFile* fexpNa22 = new TFile("../../EfficiencyPlasticFLNR/declib/output/test13_2_Plast_BC404_Na22.root", "read");
std::vector <double> *deCs137 = nullptr;
std::vector <double> *deNa22 = nullptr;

//! CHANNEL CHANGE
int channel = 0;
// int channel = 1;
// int channel = 2;

//! NUMBER OF ENTRIES FOR GRADIENT DESCENT
int entriesExp = 500000;
int entriesSim = 500000;

// int entriesExp = 200000;
// int entriesSim = 200000;
//! FFT Threshold FOR GRADIENT DESCENT (30-1500)
double threshExp = 800.;
double threshSim = 800.;

//! TUNING RATE FOR GRADIENT DESCENT
//? To increase the speed of gradient descent, increase the tuning rate
//? Ch rates are for energy calibration changes, Sig rates are for energy resolution changes

double tuningRateCh1 = 2.;
double tuningRateCh2 = 20.;
double tuningRateSig1 = 0.05;
double tuningRateSig2 = 0.05;
//!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! RESOLUTION
double dSig1 = 8.;
double dSig2 = 5.5;
// double dSig3 = 5.;
//? In percentage

//! SCALING
//? ch0_9803_ch1_9842_ch2_9854
double ScaleCs137[3] = {8.5, 8.5, 10.7};
double ScaleNa22[3] = {3.4, 3.6, 3.4};

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

  hLogis->SetLineColor(kRed);
  hLogis->Draw("same");

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
  for(int i = 3; i <= bin - 2; i++)
  {
    newBinContent = (-hChannel->GetBinContent(i+2) 
    + 8*hChannel->GetBinContent(i+1) 
    - 8*hChannel->GetBinContent(i-1) 
    + hChannel->GetBinContent(i-2))/12*binLength;
    hDiff->SetBinContent(i, newBinContent);
  }
  return hDiff;
}

double chi2(TH1D* hChannel1, TH1D* hChannel2){
  double chi2Result = hChannel1->Chi2Test(hChannel2, "UU CHI2/NDF");
  return chi2Result;
}

//! Linear Energy Calibration
double fitLinear(int channelCs137, int channelNa22, TCanvas* canvasFit, int i){
  TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
  hCali->SetBinContent(channelCs137 - xminExp, E1);
  hCali->SetBinContent(channelNa22 - xminExp, E2);

  TF1* fLinear = new TF1("fLinear", "[0]*x + [1]", xminExp, xmaxExp);
  canvasFit->cd(1);
  hCali->Fit("fLinear", "Q");

  double EnergyA = fLinear->GetParameter(0);
  double EnergyB = fLinear->GetParameter(1);

  canvasFit->Modified();
  canvasFit->Update();

  delete hCali;
  delete fLinear;

  if(i == 0){return EnergyA;}
  else{return EnergyB;}
}

//! Resolution for simulation histograms
double fitRes(double Sig1, double Sig2, TH1D* hCs137, TH1D* hNa22, double xmin, double xmax, TCanvas* canvasFit, int i){
  TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xmin, xmax);
  hRes->SetBinContent(hCs137->FindBin(E1), Sig1/100.);
  hRes->SetBinContent(hNa22->FindBin(E2), Sig2/100.);

  //! Choose fit function
  TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");
  canvasFit->cd(2);
  hRes->Fit("fRes", "Q", "", E1-0.2, E2+0.3);

  //! Correspond to chosen function
  double EnergyCoA = std::abs(fRes->GetParameter(0));
  double EnergyCoB = 0.;
  double EnergyCoC = std::abs(fRes->GetParameter(1));

  canvasFit->Modified();
  canvasFit->Update();

  delete hRes;
  delete fRes;

  if(i == 0){return EnergyCoA;}
  else if(i == 1){return EnergyCoB;}
  else{return EnergyCoC;}
}

void script_compare3()
{
  TCanvas* cFit = new TCanvas("cFit", "cFit", 1000, 500);
  cFit->Divide(2,1);

  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TTree* tsimCs137 = (TTree*) fsimCs137->Get("dEEtree");
  // TTree* texpCs137 = (TTree*) fexpCs137->Get("AnalysisxTree");

  //! Neutron Calibration
  TTree* texpCs137 = (TTree*)fexpCs137->Get("ETree");
  texpCs137->SetBranchAddress("EDep", &deCs137);

  UShort_t chX[48];
  // texpCs137->SetBranchAddress("NeEvent.neutAmp[48]", chX);
  // Long64_t entriesExpCs137 = texpCs137->GetEntries();
  Long64_t entriesExpCs137 = entriesExp;

  double x_sim;
  tsimCs137->SetBranchAddress("Scintillator", &x_sim);
  // Long64_t entriesSimCs137 = tsimCs137->GetEntries();
  Long64_t entriesSimCs137 = entriesSim;

  TTree* tsimNa22 =  (TTree*) fsimNa22->Get("dEEtree");
  // TTree* texpNa22 =  (TTree*) fexpNa22->Get("AnalysisxTree");

  //! Neutron Calibration
  TTree* texpNa22 = (TTree*)fexpNa22->Get("ETree");
  texpNa22->SetBranchAddress("EDep", &deNa22);

  UShort_t chY[48];
  // texpNa22->SetBranchAddress("NeEvent.neutAmp[48]", chY);
  // Long64_t entriesExpNa22 = texpNa22->GetEntries();
  Long64_t entriesExpNa22 = entriesExp;

  double y_sim;
  tsimNa22->SetBranchAddress("Scintillator", &y_sim);
  // Long64_t entriesSimNa22 = tsimNa22->GetEntries();
  Long64_t entriesSimNa22 = entriesSim;

  //! FFT to get Initial Compton Edge position
  TH1D* hExpCs137 = new TH1D("hExpCs137", "Cs137 Spectrum", binExp, xminExp, xmaxExp);
  TH1D* hExpNa22 = new TH1D("hExpNa22", "Na22 Spectrum", binExp, xminExp, xmaxExp);

  for(int i = 0; i < 100000; i++){
    texpCs137->GetEntry(i);
    // hExpCs137->Fill(chX[channel]);

    int de_sizeCsFFT = (int)deCs137->size(); //! Neutron Calibration
    double *de_dataCsFFT = deCs137->data();
    // if(de_dataCsFFT[channel]>5.){
    hExpCs137->Fill(de_dataCsFFT[channel]); 
    // }
  }
  
  TH1D* hExpCs137Filtered = fft(hExpCs137, 0.1, inithreshCs137, "hExpCs137Filtered", "Cs137 Spectrum Filtered", binExp, xminExp, xmaxExp);
  hExpCs137Filtered->GetXaxis()->SetTitle("Channel");
  hExpCs137Filtered->GetYaxis()->SetTitle("Count");
  TH1D* hExpCs137Diff = Diff(hExpCs137Filtered, "hExpCs137Diff", "Cs137 Spectrum Derivative", binExp, xminExp, xmaxExp);

  // std::cout << "~~~~~~~~~TEST~~~~~~~~~~~\n";
  for(int i = 0; i < 100000; i++){
    texpNa22->GetEntry(i);
    // hExpNa22->Fill(chY[channel]);

    int de_sizeNaFFT = (int)deNa22->size(); //! Neutron Calibration
    double *de_dataNaFFT = deNa22->data();
    // if(de_dataNaFFT[channel]>5.){
      hExpNa22->Fill(de_dataNaFFT[channel]); 
    // }
  }

  for(int i = 0; i<5; i++){
    hExpCs137->SetBinContent(i,0.);
    hExpNa22->SetBinContent(i,0.);
  }

  TH1D* hExpNa22Filtered = fft(hExpNa22, 0.1, inithreshNa22, "hExpNa22Filtered", "Na22 Spectrum Filtered", binExp, xminExp, xmaxExp);
  hExpNa22Filtered->GetXaxis()->SetTitle("Channel");
  hExpNa22Filtered->GetYaxis()->SetTitle("Count");
  TH1D* hExpNa22Diff = Diff(hExpNa22Filtered, "hExpNa22Diff", "Na22 Spectrum Derivative", binExp, xminExp, xmaxExp);

  TCanvas* cShow2 = new TCanvas("cShow2", "Derivative check", 1000, 500);
  cShow2->Divide(2,2);
  cShow2->cd(1);
  hExpCs137Filtered->Draw();
  cShow2->cd(3);
  hExpCs137Diff->Draw();
  hExpCs137Diff->GetXaxis()->SetRangeUser(minRangeCs137, maxRangeCs137);
  int Ch1 = hExpCs137Diff->GetMinimumBin() + xminExp;
  // std::cout << "Ch1 = " << Ch1 << "\n";

  cShow2->cd(2);
  hExpNa22Filtered->Draw();
  cShow2->cd(4);
  hExpNa22Diff->Draw();
  hExpNa22Diff->GetXaxis()->SetRangeUser(minRangeNa22, maxRangeNa22);
  int Ch2 = hExpNa22Diff->GetMinimumBin() + xminExp;
  // std::cout << "Ch2 = " << Ch2 << "\n";

  cShow2->Modified();
  cShow2->Update();
  
  double Ch1step = 1;
  double Ch2step = 1;

  double dSig1step = dSig1*0.1;
  double dSig2step = dSig2*0.1;
  
  TRandom3* ranGen = new TRandom3();


  TCanvas* cShow_SAVE = new TCanvas("cShow_SAVE", "cShow_SAVE", 1200, 1000);
  cShow_SAVE->SetLeftMargin(0.15);
  cShow_SAVE->SetTicks();
  cShow_SAVE->SetGrid();

  TCanvas* cShow = new TCanvas("cShow", "Compare Spectrum", 1800, 1000);
  cShow->Divide(2,2);
  cShow->cd(1)->SetLeftMargin(0.15);
  cShow->cd(3)->SetLeftMargin(0.15);
  cShow->cd(2)->SetLeftMargin(0.15);
  cShow->cd(4)->SetLeftMargin(0.15);

  TCanvas* cShow1 = new TCanvas("cShow1", "Gradient Check", 1500, 600);
  cShow1->Divide(4,2);

  double chi2_Cs137, chi2_Na22;
  TGraph* grChi2_Cs137 = new TGraph();
  TGraph* grChi2_Na22 = new TGraph();

  double deltaChi2 = 1000000.;
  double deltaChi2_1 = 0., deltaChi2_2 = 0., deltaChi2_3 = 0.;

  double thresh = 0.3;
  double chi2limit = 0.5;
  int count = 1;

  TCanvas* results = new TCanvas("Preliminary results", "Preliminary results", 1500, 500);
  results->Divide(3,1);

  while (deltaChi2>thresh){
    //! Energy Calibration
    std::cout << "\nIteration " << count << ":\n";

    double a = fitLinear(Ch1, Ch2, cFit, 0);
    double b = fitLinear(Ch1, Ch2, cFit, 1);

    double aUpCh1 = fitLinear(Ch1+Ch1step, Ch2, cFit, 0);
    double bUpCh1 = fitLinear(Ch1+Ch1step, Ch2, cFit, 1);

    double aUpCh2 = fitLinear(Ch1, Ch2+Ch2step, cFit, 0);
    double bUpCh2 = fitLinear(Ch1, Ch2+Ch2step, cFit, 1);

    double xminCal = xminExp*a+b;
    double xmaxCal = xmaxExp*a+b;

    double xminCalUpCh1 = xminExp*aUpCh1+bUpCh1;
    double xmaxCalUpCh1 = xmaxExp*aUpCh1+bUpCh1;

    double xminCalUpCh2 = xminExp*aUpCh2+bUpCh2;
    double xmaxCalUpCh2 = xmaxExp*aUpCh2+bUpCh2;
    
    //! Histograms (Cs137)
    TH1D* hsimCs137res = new TH1D("hsimCs137res", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalCs137 = new TH1D("hcalCs137", "Cs137 Experiment Calibrated", binExp, xminCal, xmaxCal);
    
    TH1D* hsimCs137resUpCh1 = new TH1D("hsimCs137resUpCh1", "Cs137 Simulation with resolution", binExp, xminCalUpCh1, xmaxCalUpCh1);
    TH1D* hcalCs137UpCh1 = new TH1D("hcalCs137UpCh1", "Cs137 Experiment Calibrated", binExp, xminCalUpCh1, xmaxCalUpCh1);

    TH1D* hsimCs137resUpCh2 = new TH1D("hsimCs137resUpCh2", "Cs137 Simulation with resolution", binExp, xminCalUpCh2, xmaxCalUpCh2);
    TH1D* hcalCs137UpCh2 = new TH1D("hcalCs137UpCh2", "Cs137 Experiment Calibrated", binExp, xminCalUpCh2, xmaxCalUpCh2);

    TH1D* hsimCs137resUpSig1 = new TH1D("hsimCs137resUpSig1", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);

    TH1D* hsimCs137resUpSig2 = new TH1D("hsimCs137resUpSig2", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);

    //! Histograms (Na22)
    TH1D* hsimNa22res = new TH1D("hsimNa22res", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalNa22 = new TH1D("hcalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);
    
    TH1D* hsimNa22resUpCh1 = new TH1D("hsimNa22resUpCh1", "Na22 Simulation with resolution", binExp, xminCalUpCh1, xmaxCalUpCh1);
    TH1D* hcalNa22UpCh1 = new TH1D("hcalNa22UpCh1", "Na22 Experiment Calibrated", binExp, xminCalUpCh1, xmaxCalUpCh1);

    TH1D* hsimNa22resUpCh2 = new TH1D("hsimNa22resUpCh2", "Na22 Simulation with resolution", binExp, xminCalUpCh2, xmaxCalUpCh2);
    TH1D* hcalNa22UpCh2 = new TH1D("hcalNa22UpCh2", "Na22 Experiment Calibrated", binExp, xminCalUpCh2, xmaxCalUpCh2);

    TH1D* hsimNa22resUpSig1 = new TH1D("hsimNa22resUpSig1", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    TH1D* hsimNa22resUpSig2 = new TH1D("hsimNa22resUpSig2", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    //! Fill experiment histograms (Cs137)
    for(int i = 0; i < entriesExpCs137; i++){
      texpCs137->GetEntry(i);
      // hcalCs137->Fill(a*(chX[channel]+0.5) + b);
      // hcalCs137UpCh1->Fill(aUpCh1*(chX[channel]+0.5) + bUpCh1);
      // hcalCs137UpCh2->Fill(aUpCh2*(chX[channel]+0.5) + bUpCh2);

      int de_sizeCs137 = (int)deCs137->size(); //! Neutron Calibration
      double *de_dataCs137 = deCs137->data();
      hcalCs137->Fill(a*(de_dataCs137[channel]+0.5) + b);
      hcalCs137UpCh1->Fill(aUpCh1*(de_dataCs137[channel]+0.5) + bUpCh1);
      hcalCs137UpCh2->Fill(aUpCh2*(de_dataCs137[channel]+0.5) + bUpCh2); 
    }
    TH1D* hcalCs137Filtered = fft(hcalCs137, 0.1, threshExp, "hcalCs137Filtered", "hcalCs137Filtered", binExp, xminCal, xmaxCal);
    TH1D* hcalCs137UpCh1Filtered = fft(hcalCs137UpCh1, 0.1, threshExp, "hcalCs137UpCh1Filtered", "hcalCs137UpCh1Filtered", binExp, xminCalUpCh1, xmaxCalUpCh1);
    TH1D* hcalCs137UpCh2Filtered = fft(hcalCs137UpCh2, 0.1, threshExp, "hcalCs137UpCh2Filtered", "hcalCs137UpCh2Filtered", binExp, xminCalUpCh2, xmaxCalUpCh2);

    //! Fill experiment histograms (Na22)
    for(int i = 0; i < entriesExpNa22; i++){
      texpNa22->GetEntry(i);
      hcalNa22->Fill(a*(chY[channel]+0.5) + b);
      hcalNa22UpCh1->Fill(aUpCh1*(chY[channel]+0.5) + bUpCh1);
      hcalNa22UpCh2->Fill(aUpCh2*(chY[channel]+0.5) + bUpCh2);

      int de_sizeNa22 = (int)deNa22->size(); //! Neutron Calibration
      double *de_dataNa22 = deNa22->data();
      hcalNa22->Fill(a*(de_dataNa22[channel]+0.5) + b);
      hcalNa22UpCh1->Fill(aUpCh1*(de_dataNa22[channel]+0.5) + bUpCh1);
      hcalNa22UpCh2->Fill(aUpCh2*(de_dataNa22[channel]+0.5) + bUpCh2); 
    }
    TH1D* hcalNa22Filtered = fft(hcalNa22, 0.1, threshExp, "hcalNa22Filtered", "hcalNa22Filtered", binExp, xminCal, xmaxCal);
    TH1D* hcalNa22UpCh1Filtered = fft(hcalNa22UpCh1, 0.1, threshExp, "hcalNa22UpCh1Filtered", "hcalNa22UpCh1Filtered", binExp, xminCalUpCh1, xmaxCalUpCh1);
    TH1D* hcalNa22UpCh2Filtered = fft(hcalNa22UpCh2, 0.1, threshExp, "hcalNa22UpCh2Filtered", "hcalNa22UpCh2Filtered", binExp, xminCalUpCh2, xmaxCalUpCh2);

    //! Prepare Energy Resolution coefficients
    double coA = fitRes(dSig1, dSig2, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 0);
    double coB = fitRes(dSig1, dSig2, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 1);
    double coC = fitRes(dSig1, dSig2, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 2);

    double coAUpSig1 = fitRes(dSig1+dSig1step, dSig2, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 0);
    double coBUpSig1 = fitRes(dSig1+dSig1step, dSig2, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 1);
    double coCUpSig1 = fitRes(dSig1+dSig1step, dSig2, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 2);

    double coAUpSig2 = fitRes(dSig1, dSig2+dSig2step, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 0);
    double coBUpSig2 = fitRes(dSig1, dSig2+dSig2step, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 1);
    double coCUpSig2 = fitRes(dSig1, dSig2+dSig2step, hcalCs137Filtered, hcalNa22Filtered, xminCal, xmaxCal, cFit, 2);

    //! Fill Simulation histograms (Cs137)
    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
      double energy = ranGen->Gaus(x_sim,sigma);
      hsimCs137res->Fill(energy);
      hsimCs137resUpCh1->Fill(energy);
      hsimCs137resUpCh2->Fill(energy);
    }
    TH1D* hsimCs137resFiltered = fft(hsimCs137res, 0.1, threshSim, "hsimCs137resFiltered", "hsimCs137resFiltered", binExp, xminCal, xmaxCal);
    TH1D* hsimCs137resUpCh1Filtered = fft(hsimCs137resUpCh1, 0.1, threshSim, "hsimCs137resUpCh1Filtered", "hsimCs137resUpCh1Filtered", binExp, xminCalUpCh1, xmaxCalUpCh1);
    TH1D* hsimCs137resUpCh2Filtered = fft(hsimCs137resUpCh2, 0.1, threshSim, "hsimCs137resUpCh2Filtered", "hsimCs137resUpCh2Filtered", binExp, xminCalUpCh2, xmaxCalUpCh2);    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coAUpSig1,2) + pow(coBUpSig1/sqrt(x_sim),2) + pow(coCUpSig1/x_sim,2) );
      hsimCs137resUpSig1->Fill(ranGen->Gaus(x_sim,sigma));
    }
    TH1D* hsimCs137resUpSig1Filtered = fft(hsimCs137resUpSig1, 0.1, threshSim, "hsimCs137resUpSig1Filtered", "hsimCs137resUpSig1Filtered", binExp, xminCal, xmaxCal);
  
    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coAUpSig2,2) + pow(coBUpSig2/sqrt(x_sim),2) + pow((coCUpSig2)/x_sim,2) );
      hsimCs137resUpSig2->Fill(ranGen->Gaus(x_sim,sigma));
    }
    TH1D* hsimCs137resUpSig2Filtered = fft(hsimCs137resUpSig2, 0.1, threshSim, "hsimCs137resUpSig2Filtered", "hsimCs137resUpSig2Filtered", binExp, xminCal, xmaxCal);

    hcalCs137Filtered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137resFiltered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);

    //! Fill Simulation histograms (Na22)
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      double energy = ranGen->Gaus(y_sim,sigma);
      hsimNa22res->Fill(energy);
      hsimNa22resUpCh1->Fill(energy);
      hsimNa22resUpCh2->Fill(energy);
    }
    TH1D* hsimNa22resFiltered = fft(hsimNa22res, 0.1, threshSim, "hsimNa22resFiltered", "hsimNa22resFiltered", binExp, xminCal, xmaxCal);
    TH1D* hsimNa22resUpCh1Filtered = fft(hsimNa22resUpCh1, 0.1, threshSim, "hsimNa22resUpCh1Filtered", "hsimNa22resUpCh1Filtered", binExp, xminCalUpCh1, xmaxCalUpCh1);
    TH1D* hsimNa22resUpCh2Filtered = fft(hsimNa22resUpCh2, 0.1, threshSim, "hsimNa22resUpCh2Filtered", "hsimNa22resUpCh2Filtered", binExp, xminCalUpCh2, xmaxCalUpCh2);

    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coAUpSig1,2) + pow(coBUpSig1/sqrt(y_sim),2) + pow(coCUpSig1/y_sim,2) );
      hsimNa22resUpSig1->Fill(ranGen->Gaus(y_sim,sigma));
    }
    TH1D* hsimNa22resUpSig1Filtered = fft(hsimNa22resUpSig1, 0.1, threshSim, "hsimNa22resUpSig1Filtered", "hsimNa22resUpSig1Filtered", binExp, xminCal, xmaxCal);
  
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coAUpSig2,2) + pow(coBUpSig2/sqrt(y_sim),2) + pow((coCUpSig2)/y_sim,2) );
      hsimNa22resUpSig2->Fill(ranGen->Gaus(y_sim,sigma));
    }
    TH1D* hsimNa22resUpSig2Filtered = fft(hsimNa22resUpSig2, 0.1, threshSim, "hsimNa22resUpSig2Filtered", "hsimNa22resUpSig2Filtered", binExp, xminCal, xmaxCal);

    hcalNa22Filtered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22resFiltered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    
    //! Delta Chi2
    if(count == 2){deltaChi2_1 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_Cs137) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_Na22);}
    else if(count == 3){deltaChi2_2 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_Cs137) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_Na22);}
    else if(count == 4){deltaChi2_3 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_Cs137) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_Na22);
    deltaChi2 = deltaChi2_1 + deltaChi2_2 + deltaChi2_3;}
    else if (count > 4){
    deltaChi2_1 = deltaChi2_2;
    deltaChi2_2 = deltaChi2_3;
    deltaChi2_3 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_Cs137) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_Na22);
    deltaChi2 = deltaChi2_1 + deltaChi2_2 + deltaChi2_3;
    }

    std::cout << "\ndeltaChi2_1 = " << deltaChi2_1 
    << "\ndeltaChi2_2 = " << deltaChi2_2
    << "\ndeltaChi2_3 = " << deltaChi2_3
    << "\nTotal deltaChi2 = " << deltaChi2 << "\n";

    chi2_Cs137 = chi2(hcalCs137Filtered, hsimCs137resFiltered);
    chi2_Na22 = chi2(hcalNa22Filtered, hsimNa22resFiltered);
    std::cout << "~~~~~~~~~~~CHECK~~~~~~~~~~~~\n";

    grChi2_Cs137->SetTitle("#chi^{2} distance by iteration for {}^{137}Cs");
    cShow->cd(1)->SetTicks();
    cShow->cd(1)->SetGrid();
    grChi2_Cs137->GetXaxis()->SetTitle("Iteration");
    grChi2_Cs137->GetYaxis()->SetTitle("#chi^{2} distance");

    grChi2_Cs137->SetStats(0);
    grChi2_Cs137->SetLineColor(kBlack);
    grChi2_Cs137->SetLineWidth(3);

    grChi2_Cs137->GetXaxis()->SetLabelFont(42);
    grChi2_Cs137->GetXaxis()->SetTitleFont(52);
    grChi2_Cs137->GetXaxis()->SetTitleSize(0.04);
    grChi2_Cs137->GetXaxis()->CenterTitle(true);

    grChi2_Cs137->GetYaxis()->SetLabelFont(42);
    grChi2_Cs137->GetYaxis()->SetTitleFont(52);
    grChi2_Cs137->GetYaxis()->SetTitleSize(0.04);
    grChi2_Cs137->GetYaxis()->CenterTitle(true);

    grChi2_Cs137->AddPoint(count, chi2_Cs137);

    grChi2_Na22->SetTitle("#chi^{2} distance by iteration for {}^{22}Na");
    cShow->cd(2)->SetTicks();
    cShow->cd(2)->SetGrid();
    grChi2_Na22->GetXaxis()->SetTitle("Iteration");
    grChi2_Na22->GetYaxis()->SetTitle("#chi^{2} distance");

    grChi2_Na22->SetStats(0);
    grChi2_Na22->SetLineColor(kBlack);
    grChi2_Na22->SetLineWidth(3);

    grChi2_Na22->GetXaxis()->SetLabelFont(42);
    grChi2_Na22->GetXaxis()->SetTitleFont(52);
    grChi2_Na22->GetXaxis()->SetTitleSize(0.04);
    grChi2_Na22->GetXaxis()->CenterTitle(true);

    grChi2_Na22->GetYaxis()->SetLabelFont(42);
    grChi2_Na22->GetYaxis()->SetTitleFont(52);
    grChi2_Na22->GetYaxis()->SetTitleSize(0.04);
    grChi2_Na22->GetYaxis()->CenterTitle(true);

    grChi2_Na22->AddPoint(count, chi2_Na22);

    cShow->cd(1);
    grChi2_Cs137->Draw();
    cShow->cd(1)->Modified();
    cShow->cd(1)->Update();
    
    cShow->cd(2);
    grChi2_Na22->Draw();
    cShow->cd(2)->Modified();
    cShow->cd(2)->Update();

    //! Gradient(Cs137)
    hcalCs137UpCh1Filtered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137resUpCh1Filtered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    double chi2_Cs137_UpCh1 = chi2(hcalCs137UpCh1Filtered, hsimCs137resUpCh1Filtered);

    hcalCs137UpCh2Filtered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137resUpCh2Filtered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    double chi2_Cs137_UpCh2 = chi2(hcalCs137UpCh2Filtered, hsimCs137resUpCh2Filtered);

    hsimCs137resUpSig1Filtered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    double chi2_Cs137_UpSig1 = chi2(hcalCs137Filtered, hsimCs137resUpSig1Filtered);
    
    hsimCs137resUpSig2Filtered->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    double chi2_Cs137_UpSig2 = chi2(hcalCs137Filtered, hsimCs137resUpSig2Filtered);

    //! Gradient(Na22)
    hcalNa22UpCh1Filtered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22resUpCh1Filtered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    double chi2_Na22_UpCh1 = chi2(hcalNa22UpCh1Filtered, hsimNa22resUpCh1Filtered);

    hcalNa22UpCh2Filtered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22resUpCh2Filtered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    double chi2_Na22_UpCh2 = chi2(hcalNa22UpCh2Filtered, hsimNa22resUpCh2Filtered);

    hsimNa22resUpSig1Filtered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    double chi2_Na22_UpSig1 = chi2(hcalNa22Filtered, hsimNa22resUpSig1Filtered);
    
    hsimNa22resUpSig2Filtered->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    double chi2_Na22_UpSig2 = chi2(hcalNa22Filtered, hsimNa22resUpSig2Filtered);

    //! Derivative (Cs137&Na220)
    double devChi2UpCh1 = ((chi2_Cs137_UpCh1+chi2_Na22_UpCh1)-(chi2_Cs137+chi2_Na22))/Ch1step;
    double devChi2UpCh2 = ((chi2_Cs137_UpCh2+chi2_Na22_UpCh2)-(chi2_Cs137+chi2_Na22))/Ch2step;
    double devChi2UpSig1 = ((chi2_Cs137_UpSig1+chi2_Na22_UpSig1)-(chi2_Cs137+chi2_Na22))/dSig1step;
    double devChi2UpSig2 = ((chi2_Cs137_UpSig2+chi2_Na22_UpSig2)-(chi2_Cs137+chi2_Na22))/dSig2step;

    std::cout << "devChi2UpCh1 = " << devChi2UpCh1 << "; devChi2UpCh2 = " << devChi2UpCh2 
    << "\ndevChi2UpSig1 = " << devChi2UpSig1 << "; devChi2UpSig2 = " << devChi2UpSig2 << "\n";

    Ch1 = Ch1 - tuningRateCh1*devChi2UpCh1;
    Ch2 = Ch2 - tuningRateCh2*devChi2UpCh2;

    dSig1 = dSig1 - tuningRateSig1*devChi2UpSig1;
    dSig2 = dSig2 - tuningRateSig2*devChi2UpSig2;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "Ch1 = " << Ch1 << "; Ch2 = " << Ch2 
    << "\nSig1 = " << dSig1 << "; Sig2 = " << dSig2 << "\n";

    double energy_resolution_at_1MeV = sqrt( pow(coA,2.) + std::pow((coC/1.),2.) ) / 1;
    std::cout << "E_resolution at 1MeV = " << energy_resolution_at_1MeV*100. << " (%)\n";

    cShow->cd(3);
    cShow->cd(3)->SetTicks();
    hcalCs137Filtered->SetTitle("");
    hcalCs137Filtered->GetXaxis()->SetTitle("E(MeVee)");
    hcalCs137Filtered->GetYaxis()->SetTitle("Events");

    hcalCs137Filtered->SetStats(0);
    hcalCs137Filtered->SetLineWidth(2);
    hsimCs137resFiltered->SetLineWidth(2);

    hcalCs137Filtered->GetXaxis()->SetLabelFont(42);
    hcalCs137Filtered->GetXaxis()->SetTitleFont(52);
    hcalCs137Filtered->GetXaxis()->SetTitleSize(0.04);
    hcalCs137Filtered->GetXaxis()->CenterTitle(true);

    hcalCs137Filtered->GetYaxis()->SetLabelFont(42);
    hcalCs137Filtered->GetYaxis()->SetTitleFont(52);
    hcalCs137Filtered->GetYaxis()->SetTitleSize(0.04);
    hcalCs137Filtered->GetYaxis()->CenterTitle(true);

    hcalCs137Filtered->GetXaxis()->UnZoom();
    hsimCs137resFiltered->GetXaxis()->UnZoom();
    hcalCs137Filtered->SetLineColor(kRed);
    hcalCs137Filtered->Draw();
    hsimCs137resFiltered->Scale(ScaleCs137[channel], "noSW2");
    hsimCs137resFiltered->Draw("same");

    TLegend *legend = new TLegend(0.4, 0.55, 0.8, 0.85);
    legend->SetBorderSize(0);
    legend->SetLineWidth(2);
    legend->SetTextSize(0.06);
    legend->AddEntry(hcalCs137Filtered, "{}^{137}Cs Measurement", "l");
    legend->AddEntry(hsimCs137resFiltered, "{}^{137}Cs Simulation", "l");

    legend->Draw();

    cShow->cd(4);
    cShow->cd(4)->SetTicks();
    hcalNa22Filtered->SetTitle("");
    hcalNa22Filtered->GetXaxis()->SetTitle("E(MeVee)");
    hcalNa22Filtered->GetYaxis()->SetTitle("Events");

    hcalNa22Filtered->SetStats(0);
    hcalNa22Filtered->SetLineWidth(2);
    hsimNa22resFiltered->SetLineWidth(2);

    hcalNa22Filtered->GetXaxis()->SetLabelFont(42);
    hcalNa22Filtered->GetXaxis()->SetTitleFont(52);
    hcalNa22Filtered->GetXaxis()->SetTitleSize(0.04);
    hcalNa22Filtered->GetXaxis()->CenterTitle(true);

    hcalNa22Filtered->GetYaxis()->SetLabelFont(42);
    hcalNa22Filtered->GetYaxis()->SetTitleFont(52);
    hcalNa22Filtered->GetYaxis()->SetTitleSize(0.04);
    hcalNa22Filtered->GetYaxis()->CenterTitle(true);

    hcalNa22Filtered->GetXaxis()->UnZoom();
    hsimNa22resFiltered->GetXaxis()->UnZoom();
    hcalNa22Filtered->SetLineColor(kRed);
    hcalNa22Filtered->Draw();
    hsimNa22resFiltered->Scale(ScaleNa22[channel], "noSW2");
    hsimNa22resFiltered->Draw("same");

    TLegend *legend1 = new TLegend(0.4, 0.55, 0.8, 0.85);
    legend1->SetBorderSize(0);
    legend1->SetLineWidth(10);
    legend1->SetTextSize(0.06);
    legend1->AddEntry(hcalNa22Filtered, "{}^{22}Na Measurement", "l");
    legend1->AddEntry(hsimNa22resFiltered, "{}^{22}Na Simulation", "l");

    legend1->Draw();

    cShow->cd(3)->Modified();
    cShow->cd(3)->Update();
    cShow->cd(4)->Modified();
    cShow->cd(4)->Update();

    cShow1->cd(1);    
    hcalCs137UpCh1Filtered->SetLineColor(kRed);
    hcalCs137UpCh1Filtered->Draw();
    hsimCs137resUpCh1Filtered->Scale(ScaleCs137[channel], "noSW2");
    hsimCs137resUpCh1Filtered->Draw("same");
    cShow1->cd(1)->Modified();
    cShow1->cd(1)->Update();

    cShow1->cd(2);    
    hcalCs137UpCh2Filtered->SetLineColor(kRed);
    hcalCs137UpCh2Filtered->Draw();
    hsimCs137resUpCh2Filtered->Scale(ScaleCs137[channel], "noSW2");
    hsimCs137resUpCh2Filtered->Draw("same");
    cShow1->cd(2)->Modified();
    cShow1->cd(2)->Update();

    cShow1->cd(3);
    hcalCs137Filtered->Draw();
    hsimCs137resUpSig1Filtered->Scale(ScaleCs137[channel], "noSW2");
    hsimCs137resUpSig1Filtered->Draw("same");
    cShow1->cd(3)->Modified();
    cShow1->cd(3)->Update();

    cShow1->cd(4);
    hcalCs137Filtered->Draw();
    hsimCs137resUpSig2Filtered->Scale(ScaleCs137[channel], "noSW2");
    hsimCs137resUpSig2Filtered->Draw("same");
    cShow1->cd(4)->Modified();
    cShow1->cd(4)->Update();

    cShow1->cd(5);    
    hcalNa22UpCh1Filtered->SetLineColor(kRed);
    hcalNa22UpCh1Filtered->Draw();
    hsimNa22resUpCh1Filtered->Scale(ScaleNa22[channel], "noSW2");
    hsimNa22resUpCh1Filtered->Draw("same");
    cShow1->cd(5)->Modified();
    cShow1->cd(5)->Update();

    cShow1->cd(6);    
    hcalNa22UpCh2Filtered->SetLineColor(kRed);
    hcalNa22UpCh2Filtered->Draw();
    hsimNa22resUpCh2Filtered->Scale(ScaleNa22[channel], "noSW2");
    hsimNa22resUpCh2Filtered->Draw("same");
    cShow1->cd(6)->Modified();
    cShow1->cd(6)->Update();

    cShow1->cd(7);
    hcalNa22Filtered->Draw();
    hsimNa22resUpSig1Filtered->Scale(ScaleNa22[channel], "noSW2");
    hsimNa22resUpSig1Filtered->Draw("same");
    cShow1->cd(7)->Modified();
    cShow1->cd(7)->Update();

    cShow1->cd(8);
    hcalNa22Filtered->Draw();
    hsimNa22resUpSig2Filtered->Scale(ScaleNa22[channel], "noSW2");
    hsimNa22resUpSig2Filtered->Draw("same");
    cShow1->cd(8)->Modified();
    cShow1->cd(8)->Update();

    cShow->cd();
    // cShow->Print("script_compare3_final3.gif+50");
    
    if (count == 1){
      // cShow_SAVE->cd();
      // hcalNa22Filtered->Draw();
      // hsimNa22resFiltered->Draw("same");
      // cShow_SAVE->Modified();
      // cShow_SAVE->Update();

      // TLegend *legend_SAVE = new TLegend(0.45, 0.55, 0.75, 0.85);
      // legend_SAVE->SetBorderSize(0);
      // legend_SAVE->SetLineWidth(2);
      // legend_SAVE->SetTextSize(0.05);
      // legend_SAVE->AddEntry(hcalNa22Filtered, "{}^{22}Na Measurement", "l");
      // legend_SAVE->AddEntry(hsimNa22resFiltered, "{}^{22}Na Simulation", "l");

      // legend_SAVE->Draw();

      // cShow_SAVE->Print("Na22_firstfit.png");
    }

    if(count > 1){
      if (chi2_Cs137 < chi2limit){
        tuningRateCh1 = 0.;
        tuningRateSig1 = 0.;
      }

      if (chi2_Na22 < chi2limit){
        tuningRateCh2 = 0.;
        tuningRateSig2 = 0.;
      }
    }

    if(deltaChi2 > thresh){
      delete hsimCs137res;
      delete hcalCs137;
      delete hsimCs137resUpCh1;
      delete hcalCs137UpCh1;
      delete hsimCs137resUpCh2;
      delete hcalCs137UpCh2;
      delete hsimCs137resUpSig1;
      delete hsimCs137resUpSig2;

      delete hsimCs137resFiltered;
      delete hcalCs137Filtered;
      delete hsimCs137resUpCh1Filtered;
      delete hcalCs137UpCh1Filtered;
      delete hsimCs137resUpCh2Filtered;
      delete hcalCs137UpCh2Filtered;
      delete hsimCs137resUpSig1Filtered;
      delete hsimCs137resUpSig2Filtered;

      delete hsimNa22res;
      delete hcalNa22;
      delete hsimNa22resUpCh1;
      delete hcalNa22UpCh1;
      delete hsimNa22resUpCh2;
      delete hcalNa22UpCh2;
      delete hsimNa22resUpSig1;
      delete hsimNa22resUpSig2;

      delete hsimNa22resFiltered;
      delete hcalNa22Filtered;
      delete hsimNa22resUpCh1Filtered;
      delete hcalNa22UpCh1Filtered;
      delete hsimNa22resUpCh2Filtered;
      delete hcalNa22UpCh2Filtered;
      delete hsimNa22resUpSig1Filtered;
      delete hsimNa22resUpSig2Filtered;

    }
    else{
      cShow_SAVE->cd();
      // hcalNa22Filtered->Draw();
      // hsimNa22resFiltered->Draw("same");
      // cShow_SAVE->Modified();
      // cShow_SAVE->Update();

      // cShow_SAVE->Print("Na22_lastfit.png");
      
      grChi2_Na22->Draw();

      cShow_SAVE->Modified();
      cShow_SAVE->Update();
      cShow_SAVE->Print("chi2_descent_Na22.png");
    }

    count++;
  } 

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
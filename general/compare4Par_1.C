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
double xminExp = 100;
double xmaxExp = 1500;
int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;

// double E1 = 0.341;
double E1 = 0.477;
double E2 = 1.061;
double E3 = 1.117;

//! Zoomed in Histograms options
double minTest662 = 0.3; //MeV
double maxTest662 = 0.7;
double minTest1274 = 0.6; //MeV
double maxTest1274 = 1.3;
double minTest1332 = 0.7; //MeV
double maxTest1332 = 1.4;

//! ENERGY CALIBRATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// double Ch1 = 315;
// double Ch1 = 410;
// double Ch2 = 800;
//? Now determined by FFT
int minRangeNa22 = 300;
int minRangeCo60 = 300;

//! RESOLUTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double deltaSig1 = 8.;
double deltaSig2 = 5.5;
double deltaSig3 = 5.;
//? In percentage

//! SCALING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ScaleCs137 = 8.;
double ScaleNa22 = 2.2;
double ScaleCo60 = 3.;

//! CHANNELS
// int channel = 0;
int channel = 1;
// int channel = 2;

//! Entries for checking
int entriesExp = 500000;
int entriesSim = 500000;

//! FFT Threshold
double threshExp = 100.;
double threshSim = 50.;

//! Tuning Rate
double tuningRateA = 0.00000000007;
double tuningRateB = 0.00000007;
double tuningRateCoA = 0.0000002;
double tuningRateCoC = 0.0000002;

//! CHANNEL
TFile* fsimCs137 = new TFile("./simfiles/withAl_simCs137_1.root", "read");
// TFile* fexpCs137 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Cs137_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
// TFile* fexpCs137 = new TFile("./expfiles/usable/stilbene_divider_Cs137_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
// TFile* fexpCs137 = new TFile("./expfiles/time/plastic_detectors_Cs137_start_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");
TFile* fexpCs137 = new TFile("./expfiles/time/plastic_detectors_Cs137_end_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");

TFile* fsimNa22 = new TFile("./simfiles/withAl_simNa22_2.root", "read");
// TFile* fexpNa22 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Na22_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
// TFile* fexpNa22 = new TFile("./expfiles/usable/stilbene_divider_Na22_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
// TFile* fexpNa22 = new TFile("./expfiles/time/plastic_detectors_Na22_start_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");
TFile* fexpNa22 = new TFile("./expfiles/time/plastic_detectors_Na22_end_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");

TFile* fsimCo60 = new TFile("./simfiles/withAl_simCo60_1.root", "read");
// TFile* fexpCo60 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
// TFile* fexpCo60 = new TFile("./expfiles/usable/stilbene_divider_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
// TFile* fexpCo60 = new TFile("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");
TFile* fexpCo60 = new TFile("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0008.root", "read");



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

double chi2(TH1D* hChannel1, TH1D* hChannel2){
  double chi2Result = abs(hChannel1->Chi2Test(hChannel2, "UU CHI2/NDF") - 1);
  return chi2Result;
}

void compare4Par_1()
{
  auto timer = new TStopwatch();
  timer->Start();

  TCanvas* cFit = new TCanvas("cFit", "cFit", 1000, 500);
  cFit->Divide(2,1);

  //! Files and trees
  TTree* tsimCs137 = (TTree*) fsimCs137->Get("dEEtree");
  TTree* texpCs137 = (TTree*) fexpCs137->Get("AnalysisxTree");

  UShort_t chX[48];
  texpCs137->SetBranchAddress("NeEvent.neutAmp[48]", chX);
  // Long64_t entriesExpCs137 = texpCs137->GetEntries();
  Long64_t entriesExpCs137 = entriesExp;

  double x_sim;
  tsimCs137->SetBranchAddress("Scintillator", &x_sim);
  // Long64_t entriesSimCs137 = tsimCs137->GetEntries();
  Long64_t entriesSimCs137 = entriesSim;

  TTree* tsimNa22 =  (TTree*) fsimNa22->Get("dEEtree");
  TTree* texpNa22 =  (TTree*) fexpNa22->Get("AnalysisxTree");

  UShort_t chY[48];
  texpNa22->SetBranchAddress("NeEvent.neutAmp[48]", chY);
  // Long64_t entriesExpNa22 = texpNa22->GetEntries();
  Long64_t entriesExpNa22 = entriesExp;

  double y_sim;
  tsimNa22->SetBranchAddress("Scintillator", &y_sim);
  // Long64_t entriesSimNa22 = tsimNa22->GetEntries();
  Long64_t entriesSimNa22 = entriesSim;

  TTree* tsimCo60 =  (TTree*) fsimCo60->Get("dEEtree");
  TTree* texpCo60 =  (TTree*) fexpCo60->Get("AnalysisxTree");

  UShort_t chZ[48];
  texpCo60->SetBranchAddress("NeEvent.neutAmp[48]", chZ);
  // Long64_t entriesExpCo60 = texpCo60->GetEntries();
  Long64_t entriesExpCo60 = entriesExp;

  double z_sim;
  tsimCo60->SetBranchAddress("Scintillator", &z_sim);
  // Long64_t entriesSimCo60 = tsimCo60->GetEntries();
  Long64_t entriesSimCo60 = entriesSim;

  //! FFT to get Compton Edge position
  TH1D* hExpCs137 = new TH1D("hExpCs137", "Cs137 Spectrum", binExp, xminExp, xmaxExp);
  TH1D* hExpNa22 = new TH1D("hExpNa22", "Na22 Spectrum", binExp, xminExp, xmaxExp);
  TH1D* hExpCo60 = new TH1D("hExpCo60", "Co60 Spectrum", binExp, xminExp, xmaxExp);

  for(int i = 0; i < 1000000; i++){
    texpCs137->GetEntry(i);
    hExpCs137->Fill(chX[channel]);
  }
  TH1D* hExpCs137Filtered = fft(hExpCs137, 0.1, 100, "hExpCs137Filtered", "Cs137 Spectrum Filtered", binExp, xminExp, xmaxExp);
  TH1D* hExpCs137Diff = Diff(hExpCs137Filtered, "hExpCs137Diff", "Cs137 Spectrum Derivative", binExp, xminExp, xmaxExp);
  TH1D* hExpCs137DiffFiltered = fft(hExpCs137Diff, 0.1, 50, "hExpCs137DiffFiltered", "Cs137 Spectrum Derivative Filtered", binExp, xminExp, xmaxExp);

  for(int i = 0; i < 1000000; i++){
    texpNa22->GetEntry(i);
    hExpNa22->Fill(chY[channel]);
  }
  TH1D* hExpNa22Filtered = fft(hExpNa22, 0.1, 100, "hExpNa22Filtered", "Na22 Spectrum Filtered", binExp, xminExp, xmaxExp);
  TH1D* hExpNa22Diff = Diff(hExpNa22Filtered, "hExpNa22Diff", "Na22 Spectrum Derivative", binExp, xminExp, xmaxExp);
  TH1D* hExpNa22DiffFiltered = fft(hExpNa22Diff, 0.1, 50, "hExpNa22DiffFiltered", "Na22 Spectrum Derivative Filtered", binExp, xminExp, xmaxExp);

  for(int i = 0; i < 1000000; i++){
    texpCo60->GetEntry(i);
    hExpCo60->Fill(chZ[channel]);
  }
  TH1D* hExpCo60Filtered = fft(hExpCo60, 0.1, 100, "hExpCo60Filtered", "Co60 Spectrum Filtered", binExp, xminExp, xmaxExp);
  TH1D* hExpCo60Diff = Diff(hExpCo60Filtered, "hExpCo60Diff", "Co60 Spectrum Derivative", binExp, xminExp, xmaxExp);
  TH1D* hExpCo60DiffFiltered = fft(hExpCo60Diff, 0.1, 50, "hExpCo60DiffFiltered", "Co60 Spectrum Derivative Filtered", binExp, xminExp, xmaxExp);

  TCanvas* cShow2 = new TCanvas("cShow2", "Derivative check", 1000, 500);
  cShow2->Divide(3,2);
  cShow2->cd(1);
  hExpCs137Filtered->Draw();
  cShow2->cd(4);
  hExpCs137DiffFiltered->Draw();
  // hExpCs137Diff->Draw();
  int Ch1 = hExpCs137DiffFiltered->GetMinimumBin() + xminExp;
  std::cout << "Ch1 = " << Ch1 << "\n";

  cShow2->cd(2);
  hExpNa22Filtered->Draw();
  cShow2->cd(5);
  hExpNa22DiffFiltered->Draw();
  hExpNa22DiffFiltered->GetXaxis()->SetRangeUser(minRangeNa22, xmaxExp);
  int Ch2 = hExpNa22DiffFiltered->GetMinimumBin() + xminExp;
  // hExpNa22DiffFiltered->GetXaxis()->UnZoom();
  std::cout << "Ch2 = " << Ch2 << "\n";
  // hExpNa22Diff->Draw();

  cShow2->cd(3);
  hExpCo60Filtered->Draw();
  cShow2->cd(6);
  hExpCo60DiffFiltered->Draw();
  hExpCo60DiffFiltered->GetXaxis()->SetRangeUser(minRangeCo60, xmaxExp);
  int Ch3 = hExpCo60DiffFiltered->GetMinimumBin() + xminExp;
  // hExpCo60DiffFiltered->GetXaxis()->UnZoom();
  std::cout << "Ch3 = " << Ch3 << "\n";
  // hExpCo60Diff->Draw();

  cShow2->Modified();
  cShow2->Update();
  
  TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
  hCali->SetBinContent(Ch1 - xminExp, E1);
  hCali->SetBinContent(Ch2 - xminExp, E2);
  hCali->SetBinContent(Ch3 - xminExp, E3);

  TF1* fLinear = new TF1("fLinear", "[0]*x + [1]", xminExp, xmaxExp);
  cFit->cd(1);
  hCali->Fit("fLinear");

  double a = fLinear->GetParameter(0);
  double b = fLinear->GetParameter(1);

  double aStep = 0.01*a;
  double bStep = 0.01*b;

  //! Resolution for simulation histograms

  std::cout << "delta1 = " << deltaSig1 << "%; delta2 = " << deltaSig2 << "%; delta3 = " << deltaSig3 << "%\n";

  double xminCal = xminExp*a+b;
  double xmaxCal = xmaxExp*a+b;

  TH1D* hFirstCalCs137 = new TH1D("hFirstCalCs137", "Relative Resolution Calibrated", binExp, xminCal, xmaxCal);
  for(int i = 0; i < entriesExpCs137; i++){
    texpCs137->GetEntry(i);
    hFirstCalCs137->Fill(a*(chX[channel]+0.5) + b);
  }

  TH1D* hFirstCalNa22 = new TH1D("hFirstCalNa22", "Relative Resolution Calibrated", binExp, xminCal, xmaxCal);
  for(int i = 0; i < entriesExpNa22; i++){
    texpNa22->GetEntry(i);
    hFirstCalNa22->Fill(a*(chY[channel]+0.5) + b);
  }

  TH1D* hFirstCalCo60 = new TH1D("hFirstCalCo60", "Relative Resolution Calibrated", binExp, xminCal, xmaxCal);
  for(int i = 0; i < entriesExpCo60; i++){
    texpCo60->GetEntry(i);
    hFirstCalCo60->Fill(a*(chZ[channel]+0.5) + b);
  }

  TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal, xmaxCal);
  hRes->SetBinContent(hFirstCalCs137->FindBin(E1), deltaSig1/100.);
  hRes->SetBinContent(hFirstCalNa22->FindBin(E2), deltaSig2/100.);
  hRes->SetBinContent(hFirstCalCo60->FindBin(E3), deltaSig3/100.);

  //! Choose fit function
  TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");
  cFit->cd(2);
  hRes->Fit("fRes", "", "", E1-0.2, E2+0.3);
  cFit->Modified();
  cFit->Update();

  //! Correspond to chosen function
  double coA = std::abs(fRes->GetParameter(0));
  double coB = 0.;
  double coC = std::abs(fRes->GetParameter(1));

  double coAStep = coA*0.01;
  // double coBStep = coB*0.2;
  double coCStep = coC*0.01;
  
  TRandom3* ranGen = new TRandom3();

  TCanvas* cShow = new TCanvas("cShow", "Compare Spectrum", 1800, 600);
  cShow->Divide(3,2);

  TCanvas* cShow1 = new TCanvas("cShow1", "Gradient Check", 1500, 600);
  cShow1->Divide(4,3);

  double chi2_662, chi2_1274, chi2_1332;
  TGraph* grChi2_662 = new TGraph();
  TGraph* grChi2_1274 = new TGraph();
  TGraph* grChi2_1332 = new TGraph();

  double deltaChi2 = 1000000.;
  double deltaChi2_1 = 0., deltaChi2_2 = 0., deltaChi2_3 = 0.;

  double thresh = 1.0;
  int count = 1;
  while (deltaChi2>thresh){
    std::cout << "\nIteration " << count << ":\n";

    xminCal = xminExp*a+b;
    xmaxCal = xmaxExp*a+b;

    double xminCalUpA = xminExp*(a + aStep)+b;
    double xmaxCalUpA = xmaxExp*(a + aStep)+b;

    double xminCalUpB = xminExp*a+(b + bStep);
    double xmaxCalUpB = xmaxExp*a+(b + bStep);

    //! Histograms (Na22)
    TH1D* hsimNa22res = new TH1D("hsimNa22res", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalNa22 = new TH1D("hcalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);
    
    TH1D* hsimNa22resUpA = new TH1D("hsimNa22resUpA", "Na22 Simulation with resolution", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalNa22UpA = new TH1D("hcalNa22UpA", "Na22 Experiment Calibrated", binExp, xminCalUpA, xmaxCalUpA);
   
    TH1D* hsimNa22resUpB = new TH1D("hsimNa22resUpB", "Na22 Simulation with resolution", binExp, xminCalUpB, xmaxCalUpB);
    TH1D* hcalNa22UpB = new TH1D("hcalNa22UpB", "Na22 Experiment Calibrated", binExp, xminCalUpB, xmaxCalUpB);

    TH1D* hsimNa22resUpCoA = new TH1D("hsimNa22resUpCoA", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    TH1D* hsimNa22resUpCoC = new TH1D("hsimNa22resUpCoC", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    //! Histograms (Cs137)
    TH1D* hsimCs137res = new TH1D("hsimCs137res", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalCs137 = new TH1D("hcalCs137", "Cs137 Experiment Calibrated", binExp, xminCal, xmaxCal);
    
    TH1D* hsimCs137resUpA = new TH1D("hsimCs137resUpA", "Cs137 Simulation with resolution", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalCs137UpA = new TH1D("hcalCs137UpA", "Cs137 Experiment Calibrated", binExp, xminCalUpA, xmaxCalUpA);
   
    TH1D* hsimCs137resUpB = new TH1D("hsimCs137resUpB", "Cs137 Simulation with resolution", binExp, xminCalUpB, xmaxCalUpB);
    TH1D* hcalCs137UpB = new TH1D("hcalCs137UpB", "Cs137 Experiment Calibrated", binExp, xminCalUpB, xmaxCalUpB);

    TH1D* hsimCs137resUpCoA = new TH1D("hsimCs137resUpCoA", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);

    TH1D* hsimCs137resUpCoC = new TH1D("hsimCs137resUpCoC", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);

    //! Histograms (Co60)
    TH1D* hsimCo60res = new TH1D("hsimCo60res", "Co60 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalCo60 = new TH1D("hcalCo60", "Co60 Experiment Calibrated", binExp, xminCal, xmaxCal);
    
    TH1D* hsimCo60resUpA = new TH1D("hsimCo60resUpA", "Co60 Simulation with resolution", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalCo60UpA = new TH1D("hcalCo60UpA", "Co60 Experiment Calibrated", binExp, xminCalUpA, xmaxCalUpA);
   
    TH1D* hsimCo60resUpB = new TH1D("hsimCo60resUpB", "Co60 Simulation with resolution", binExp, xminCalUpB, xmaxCalUpB);
    TH1D* hcalCo60UpB = new TH1D("hcalCo60UpB", "Co60 Experiment Calibrated", binExp, xminCalUpB, xmaxCalUpB);

    TH1D* hsimCo60resUpCoA = new TH1D("hsimCo60resUpCoA", "Co60 Simulation with resolution", binExp, xminCal, xmaxCal);

    TH1D* hsimCo60resUpCoC = new TH1D("hsimCo60resUpCoC", "Co60 Simulation with resolution", binExp, xminCal, xmaxCal);

    //! Fill experiment histograms (Cs137)
    for(int i = 0; i < entriesExpCs137; i++){
      texpCs137->GetEntry(i);
      hcalCs137->Fill(a*(chX[channel]+0.5) + b);
      hcalCs137UpA->Fill((a+aStep)*(chX[channel]+0.5) + b);
      hcalCs137UpB->Fill(a*(chX[channel]+0.5) + (b+bStep));
    }
    TH1D* hcalCs137Filtered = fft(hcalCs137, 0.1, threshExp, "hcalCs137Filtered", "hcalCs137Filtered", binExp, xminCal, xmaxCal);
    TH1D* hcalCs137UpAFiltered = fft(hcalCs137UpA, 0.1, threshExp, "hcalCs137UpAFiltered", "hcalCs137UpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalCs137UpBFiltered = fft(hcalCs137UpB, 0.1, threshExp, "hcalCs137UpBFiltered", "hcalCs137UpBFiltered", binExp, xminCalUpB, xmaxCalUpB);

    //! Fill experiment histograms (Na22)
    for(int i = 0; i < entriesExpNa22; i++){
      texpNa22->GetEntry(i);
      hcalNa22->Fill(a*(chY[channel]+0.5) + b);
      hcalNa22UpA->Fill((a+aStep)*(chY[channel]+0.5) + b);
      hcalNa22UpB->Fill(a*(chY[channel]+0.5) + (b+bStep));
    }
    TH1D* hcalNa22Filtered = fft(hcalNa22, 0.1, threshExp, "hcalNa22Filtered", "hcalNa22Filtered", binExp, xminCal, xmaxCal);
    TH1D* hcalNa22UpAFiltered = fft(hcalNa22UpA, 0.1, threshExp, "hcalNa22UpAFiltered", "hcalNa22UpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalNa22UpBFiltered = fft(hcalNa22UpB, 0.1, threshExp, "hcalNa22UpBFiltered", "hcalNa22UpBFiltered", binExp, xminCalUpB, xmaxCalUpB);

    //! Fill experiment histograms (Co60)
    for(int i = 0; i < entriesExpCo60; i++){
      texpCo60->GetEntry(i);
      hcalCo60->Fill(a*(chZ[channel]+0.5) + b);
      hcalCo60UpA->Fill((a+aStep)*(chZ[channel]+0.5) + b);
      hcalCo60UpB->Fill(a*(chZ[channel]+0.5) + (b+bStep));
    }
    TH1D* hcalCo60Filtered = fft(hcalCo60, 0.1, threshExp, "hcalCo60Filtered", "hcalCo60Filtered", binExp, xminCal, xmaxCal);
    TH1D* hcalCo60UpAFiltered = fft(hcalCo60UpA, 0.1, threshExp, "hcalCo60UpAFiltered", "hcalCo60UpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalCo60UpBFiltered = fft(hcalCo60UpB, 0.1, threshExp, "hcalCo60UpBFiltered", "hcalCo60UpBFiltered", binExp, xminCalUpB, xmaxCalUpB);

    //! Fill Simulation histograms (Cs137)
    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
      double energy = ranGen->Gaus(x_sim,sigma);
      hsimCs137res->Fill(energy);
      hsimCs137resUpA->Fill(energy);
      hsimCs137resUpB->Fill(energy);
    }
    TH1D* hsimCs137resFiltered = fft(hsimCs137res, 0.1, threshSim, "hsimCs137resFiltered", "hsimCs137resFiltered", binExp, xminCal, xmaxCal);
    TH1D* hsimCs137resUpAFiltered = fft(hsimCs137resUpA, 0.1, threshSim, "hsimCs137resUpAFiltered", "hsimCs137resUpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hsimCs137resUpBFiltered = fft(hsimCs137resUpB, 0.1, threshSim, "hsimCs137resUpBFiltered", "hsimCs137resUpBFiltered", binExp, xminCalUpB, xmaxCalUpB);

    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA+coAStep,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
      hsimCs137resUpCoA->Fill(ranGen->Gaus(x_sim,sigma));
    }
    TH1D* hsimCs137resUpCoAFiltered = fft(hsimCs137resUpCoA, 0.1, threshSim, "hsimCs137resUpCoAFiltered", "hsimCs137resUpCoAFiltered", binExp, xminCal, xmaxCal);
  
    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow((coC+coCStep)/x_sim,2) );
      hsimCs137resUpCoC->Fill(ranGen->Gaus(x_sim,sigma));
    }
    TH1D* hsimCs137resUpCoCFiltered = fft(hsimCs137resUpCoA, 0.1, threshSim, "hsimCs137resUpCoCFiltered", "hsimCs137resUpCoCFiltered", binExp, xminCal, xmaxCal);

    hcalCs137Filtered->GetXaxis()->SetRangeUser(minTest662, maxTest662);
    hsimCs137resFiltered->GetXaxis()->SetRangeUser(minTest662, maxTest662);

    //! Fill Simulation histograms (Na22)
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      double energy = ranGen->Gaus(y_sim,sigma);
      hsimNa22res->Fill(energy);
      hsimNa22resUpA->Fill(energy);
      hsimNa22resUpB->Fill(energy);
    }
    TH1D* hsimNa22resFiltered = fft(hsimNa22res, 0.1, threshSim, "hsimNa22resFiltered", "hsimNa22resFiltered", binExp, xminCal, xmaxCal);
    TH1D* hsimNa22resUpAFiltered = fft(hsimNa22resUpA, 0.1, threshSim, "hsimNa22resUpAFiltered", "hsimNa22resUpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hsimNa22resUpBFiltered = fft(hsimNa22resUpB, 0.1, threshSim, "hsimNa22resUpBFiltered", "hsimNa22resUpBFiltered", binExp, xminCalUpB, xmaxCalUpB);    
    
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA+coAStep,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      hsimNa22resUpCoA->Fill(ranGen->Gaus(y_sim,sigma));
    }
    TH1D* hsimNa22resUpCoAFiltered = fft(hsimNa22resUpCoA, 0.1, threshSim, "hsimNa22resUpCoAFiltered", "hsimNa22resUpCoAFiltered", binExp, xminCal, xmaxCal);
  
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow((coC+coCStep)/y_sim,2) );
      hsimNa22resUpCoC->Fill(ranGen->Gaus(y_sim,sigma));
    }
    TH1D* hsimNa22resUpCoCFiltered = fft(hsimNa22resUpCoA, 0.1, threshSim, "hsimNa22resUpCoCFiltered", "hsimNa22resUpCoCFiltered", binExp, xminCal, xmaxCal);

    hcalNa22Filtered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    hsimNa22resFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    
    //! Fill Simulation histograms (Co60)
    for (int i = 0; i < entriesSimCo60; i++)
    {
      tsimCo60->GetEntry(i);
      double sigma = z_sim*sqrt( pow(coA,2) + pow(coB/sqrt(z_sim),2) + pow(coC/z_sim,2) );
      double energy = ranGen->Gaus(z_sim,sigma);
      hsimCo60res->Fill(energy);
      hsimCo60resUpA->Fill(energy);
      hsimCo60resUpB->Fill(energy);
    }
    TH1D* hsimCo60resFiltered = fft(hsimCo60res, 0.1, threshSim, "hsimCo60resFiltered", "hsimCo60resFiltered", binExp, xminCal, xmaxCal);
    TH1D* hsimCo60resUpAFiltered = fft(hsimCo60resUpA, 0.1, threshSim, "hsimCo60resUpAFiltered", "hsimCo60resUpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hsimCo60resUpBFiltered = fft(hsimCo60resUpB, 0.1, threshSim, "hsimCo60resUpBFiltered", "hsimCo60resUpBFiltered", binExp, xminCalUpB, xmaxCalUpB);    
    
    for (int i = 0; i < entriesSimCo60; i++)
    {
      tsimCo60->GetEntry(i);
      double sigma = z_sim*sqrt( pow(coA+coAStep,2) + pow(coB/sqrt(z_sim),2) + pow(coC/z_sim,2) );
      hsimCo60resUpCoA->Fill(ranGen->Gaus(z_sim,sigma));
    }
    TH1D* hsimCo60resUpCoAFiltered = fft(hsimCo60resUpCoA, 0.1, threshSim, "hsimCo60resUpCoAFiltered", "hsimCo60resUpCoAFiltered", binExp, xminCal, xmaxCal);
  
    for (int i = 0; i < entriesSimCo60; i++)
    {
      tsimCo60->GetEntry(i);
      double sigma = z_sim*sqrt( pow(coA,2) + pow(coB/sqrt(z_sim),2) + pow((coC+coCStep)/z_sim,2) );
      hsimCo60resUpCoC->Fill(ranGen->Gaus(z_sim,sigma));
    }
    TH1D* hsimCo60resUpCoCFiltered = fft(hsimCo60resUpCoA, 0.1, threshSim, "hsimCo60resUpCoCFiltered", "hsimCo60resUpCoCFiltered", binExp, xminCal, xmaxCal);

    hcalCo60Filtered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    hsimCo60resFiltered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    
    //! Delta Chi2
    if(count == 2){deltaChi2_1 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_662) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_1274)
    + abs(chi2(hcalCo60Filtered, hsimCo60resFiltered) -chi2_1332);}
    else if(count == 3){deltaChi2_2 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_662) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_1274)
    + abs(chi2(hcalCo60Filtered, hsimCo60resFiltered) -chi2_1332);}
    else if(count == 4){deltaChi2_3 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_662) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_1274)
    + abs(chi2(hcalCo60Filtered, hsimCo60resFiltered) -chi2_1332);
    deltaChi2 = deltaChi2_1 + deltaChi2_2 + deltaChi2_3;}
    else if (count > 4){
    deltaChi2_1 = deltaChi2_2;
    deltaChi2_2 = deltaChi2_3;
    deltaChi2_3 = abs(chi2(hcalCs137Filtered, hsimCs137resFiltered) - chi2_662) 
    + abs(chi2(hcalNa22Filtered, hsimNa22resFiltered) -chi2_1274)
    + abs(chi2(hcalCo60Filtered, hsimCo60resFiltered) -chi2_1332);
    deltaChi2 = deltaChi2_1 + deltaChi2_2 + deltaChi2_3;
    }

    std::cout << "\ndeltaChi2_1 = " << deltaChi2_1 
    << "\ndeltaChi2_2 = " << deltaChi2_2
    << "\ndeltaChi2_3 = " << deltaChi2_3
    << "\nTotal deltaChi2 = " << deltaChi2 << "\n";

    chi2_662 = chi2(hcalCs137Filtered, hsimCs137resFiltered);
    chi2_1274 = chi2(hcalNa22Filtered, hsimNa22resFiltered);
    chi2_1332 = chi2(hcalCo60Filtered, hsimCo60resFiltered);

    grChi2_662->SetTitle("|r-Chi2 - 1| by iteration (Cs137)");
    grChi2_662->GetXaxis()->SetTitle("Iteration");
    grChi2_662->GetYaxis()->SetTitle("|r-Chi2 - 1|");
    grChi2_662->AddPoint(count, chi2_662);

    grChi2_1274->SetTitle("|r-Chi2 - 1| by iteration (Na22)");
    grChi2_1274->GetXaxis()->SetTitle("Iteration");
    grChi2_1274->GetYaxis()->SetTitle("|r-Chi2 - 1|");
    grChi2_1274->AddPoint(count, chi2_1274);

    grChi2_1332->SetTitle("|r-Chi2 - 1| by iteration (Co60)");
    grChi2_1332->GetXaxis()->SetTitle("Iteration");
    grChi2_1332->GetYaxis()->SetTitle("|r-Chi2 - 1|");
    grChi2_1332->AddPoint(count, chi2_1332);

    cShow->cd(1);
    grChi2_662->Draw();
    cShow->cd(1)->Modified();
    cShow->cd(1)->Update();
    
    cShow->cd(2);
    grChi2_1274->Draw();
    cShow->cd(2)->Modified();
    cShow->cd(2)->Update();

    cShow->cd(3);
    grChi2_1332->Draw();
    cShow->cd(3)->Modified();
    cShow->cd(3)->Update();

    //! Gradient(Cs137)
    hcalCs137UpAFiltered->GetXaxis()->SetRangeUser(minTest662, maxTest662);
    hsimCs137resUpAFiltered->GetXaxis()->SetRangeUser(minTest662, maxTest662);
    double chi2_662_UpA = chi2(hcalCs137UpAFiltered, hsimCs137resUpAFiltered);

    hcalCs137UpBFiltered->GetXaxis()->SetRangeUser(minTest662, maxTest662);
    hsimCs137resUpBFiltered->GetXaxis()->SetRangeUser(minTest662, maxTest662);
    double chi2_662_UpB = chi2(hcalCs137UpBFiltered, hsimCs137resUpBFiltered);

    hsimCs137resUpCoAFiltered->GetXaxis()->SetRangeUser(minTest662, maxTest662);
    double chi2_662_UpCoA = chi2(hcalCs137Filtered, hsimCs137resUpCoAFiltered);
    
    hsimCs137resUpCoCFiltered->GetXaxis()->SetRangeUser(minTest662, maxTest662);
    double chi2_662_UpCoC = chi2(hcalCs137Filtered, hsimCs137resUpCoCFiltered);

    //! Gradient(Na22)
    hcalNa22UpAFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    hsimNa22resUpAFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpA = chi2(hcalNa22UpAFiltered, hsimNa22resUpAFiltered);

    hcalNa22UpBFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    hsimNa22resUpBFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpB = chi2(hcalNa22UpBFiltered, hsimNa22resUpBFiltered);

    hsimNa22resUpCoAFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpCoA = chi2(hcalNa22Filtered, hsimNa22resUpCoAFiltered);
    
    hsimNa22resUpCoCFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpCoC = chi2(hcalNa22Filtered, hsimNa22resUpCoCFiltered);

    //! Gradient(Co60)
    hcalCo60UpAFiltered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    hsimCo60resUpAFiltered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    double chi2_1332_UpA = chi2(hcalCo60UpAFiltered, hsimCo60resUpAFiltered);

    hcalCo60UpBFiltered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    hsimCo60resUpBFiltered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    double chi2_1332_UpB = chi2(hcalCo60UpBFiltered, hsimCo60resUpBFiltered);

    hsimCo60resUpCoAFiltered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    double chi2_1332_UpCoA = chi2(hcalCo60Filtered, hsimCo60resUpCoAFiltered);
    
    hsimCo60resUpCoCFiltered->GetXaxis()->SetRangeUser(minTest1332, maxTest1332);
    double chi2_1332_UpCoC = chi2(hcalCo60Filtered, hsimCo60resUpCoCFiltered);

    //! Derivative (Cs137&Na22&Co60)
    double devChi2UpA = ((chi2_1274_UpA+chi2_662_UpA+chi2_1332_UpA)-(chi2_1274+chi2_662+chi2_1332))/aStep;
    double devChi2UpB = ((chi2_1274_UpB+chi2_662_UpB+chi2_1332_UpB)-(chi2_1274+chi2_662+chi2_1332))/bStep;
    double devChi2UpCoA = ((chi2_1274_UpCoA+chi2_662_UpCoA+chi2_1332_UpCoA)-(chi2_1274+chi2_662+chi2_1332))/coAStep;
    double devChi2UpCoC = ((chi2_1274_UpCoC+chi2_662_UpCoC+chi2_1332_UpCoC)-(chi2_1274+chi2_662+chi2_1332))/coCStep;

    std::cout << "devUpA = " << devChi2UpA << "; devUpB = " << devChi2UpB
    << "\ndevUpCoA = " << devChi2UpCoA << "; devUpCoC = " << devChi2UpCoC << "\n";

    a = a - tuningRateA*devChi2UpA;
    b = b - tuningRateB*devChi2UpB;
    coA = coA - tuningRateCoA*devChi2UpCoA;
    coC = coC - tuningRateCoC*devChi2UpCoC;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "a = " << a << "; b = " << b << "\n";
    std::cout << "coA = " << coA << "; coC = " << coC << "\n"; 

    double finalChannel1 = (E1-b)/a;
    double finalChannel2 = (E2-b)/a;
    double finalChannel3 = (E3-b)/a;

    double finalSig1 = 100.*sqrt(coA*coA + coC*coC/(E1*E1));
    double finalSig2 = 100.*sqrt(coA*coA + coC*coC/(E2*E2));
    double finalSig3 = 100.*sqrt(coA*coA + coC*coC/(E3*E3));

    std::cout << "Ch1 = " << finalChannel1 << "; Ch2 = " << finalChannel2 << "; Ch3 = " << finalChannel3
    << "\nSig1 = " << finalSig1 << "; Sig2 = " << finalSig2 << "; Sig3 = " << finalSig3 << "\n";

    cShow->cd(4);
    hcalCs137Filtered->GetXaxis()->UnZoom();
    hsimCs137resFiltered->GetXaxis()->UnZoom();
    hcalCs137Filtered->SetLineColor(kRed);
    hcalCs137Filtered->Draw();
    hsimCs137resFiltered->Scale(ScaleCs137, "noSW2");
    hsimCs137resFiltered->Draw("same");
    cShow->cd(4)->Modified();
    cShow->cd(4)->Update();

    cShow->cd(5);
    hcalNa22Filtered->GetXaxis()->UnZoom();
    hsimNa22resFiltered->GetXaxis()->UnZoom();
    hcalNa22Filtered->SetLineColor(kRed);
    hcalNa22Filtered->Draw();
    hsimNa22resFiltered->Scale(ScaleNa22, "noSW2");
    hsimNa22resFiltered->Draw("same");
    cShow->cd(5)->Modified();
    cShow->cd(5)->Update();

    cShow->cd(6);
    hcalCo60Filtered->GetXaxis()->UnZoom();
    hsimCo60resFiltered->GetXaxis()->UnZoom();
    hcalCo60Filtered->SetLineColor(kRed);
    hcalCo60Filtered->Draw();
    hsimCo60resFiltered->Scale(ScaleCo60, "noSW2");
    hsimCo60resFiltered->Draw("same");
    cShow->cd(6)->Modified();
    cShow->cd(6)->Update();

    cShow1->cd(1);    
    hcalNa22UpAFiltered->SetLineColor(kRed);
    hcalNa22UpAFiltered->Draw();
    hsimNa22resUpAFiltered->Scale(ScaleNa22, "noSW2");
    hsimNa22resUpAFiltered->Draw("same");
    cShow1->cd(1)->Modified();
    cShow1->cd(1)->Update();

    cShow1->cd(2);    
    hcalNa22UpBFiltered->SetLineColor(kRed);
    hcalNa22UpBFiltered->Draw();
    hsimNa22resUpBFiltered->Scale(ScaleNa22, "noSW2");
    hsimNa22resUpBFiltered->Draw("same");
    cShow1->cd(2)->Modified();
    cShow1->cd(2)->Update();

    cShow1->cd(3);
    hcalNa22Filtered->Draw();
    hsimNa22resUpCoAFiltered->Scale(ScaleNa22, "noSW2");
    hsimNa22resUpCoAFiltered->Draw("same");
    cShow1->cd(3)->Modified();
    cShow1->cd(3)->Update();

    cShow1->cd(4);
    hcalNa22Filtered->Draw();
    hsimNa22resUpCoCFiltered->Scale(ScaleNa22, "noSW2");
    hsimNa22resUpCoCFiltered->Draw("same");
    cShow1->cd(4)->Modified();
    cShow1->cd(4)->Update();

    cShow1->cd(5);    
    hcalCs137UpAFiltered->SetLineColor(kRed);
    hcalCs137UpAFiltered->Draw();
    hsimCs137resUpAFiltered->Scale(ScaleCs137, "noSW2");
    hsimCs137resUpAFiltered->Draw("same");
    cShow1->cd(5)->Modified();
    cShow1->cd(5)->Update();

    cShow1->cd(6);    
    hcalCs137UpBFiltered->SetLineColor(kRed);
    hcalCs137UpBFiltered->Draw();
    hsimCs137resUpBFiltered->Scale(ScaleCs137, "noSW2");
    hsimCs137resUpBFiltered->Draw("same");
    cShow1->cd(6)->Modified();
    cShow1->cd(6)->Update();

    cShow1->cd(7);
    hcalCs137Filtered->Draw();
    hsimCs137resUpCoAFiltered->Scale(ScaleCs137, "noSW2");
    hsimCs137resUpCoAFiltered->Draw("same");
    cShow1->cd(7)->Modified();
    cShow1->cd(7)->Update();

    cShow1->cd(8);
    hcalCs137Filtered->Draw();
    hsimCs137resUpCoCFiltered->Scale(ScaleCs137, "noSW2");
    hsimCs137resUpCoCFiltered->Draw("same");
    cShow1->cd(8)->Modified();
    cShow1->cd(8)->Update();

    cShow1->cd(9);    
    hcalCo60UpAFiltered->SetLineColor(kRed);
    hcalCo60UpAFiltered->Draw();
    hsimCo60resUpAFiltered->Scale(ScaleCo60, "noSW2");
    hsimCo60resUpAFiltered->Draw("same");
    cShow1->cd(9)->Modified();
    cShow1->cd(9)->Update();

    cShow1->cd(10);    
    hcalCo60UpBFiltered->SetLineColor(kRed);
    hcalCo60UpBFiltered->Draw();
    hsimCo60resUpBFiltered->Scale(ScaleCo60, "noSW2");
    hsimCo60resUpBFiltered->Draw("same");
    cShow1->cd(10)->Modified();
    cShow1->cd(10)->Update();

    cShow1->cd(11);
    hcalCo60Filtered->Draw();
    hsimCo60resUpCoAFiltered->Scale(ScaleCo60, "noSW2");
    hsimCo60resUpCoAFiltered->Draw("same");
    cShow1->cd(11)->Modified();
    cShow1->cd(11)->Update();

    cShow1->cd(12);
    hcalCo60Filtered->Draw();
    hsimCo60resUpCoCFiltered->Scale(ScaleCo60, "noSW2");
    hsimCo60resUpCoCFiltered->Draw("same");
    cShow1->cd(12)->Modified();
    cShow1->cd(12)->Update();

    if(deltaChi2 > thresh){
      delete hsimNa22res;
      delete hcalNa22;
      delete hsimNa22resUpA;
      delete hcalNa22UpA;
      delete hsimNa22resUpB;
      delete hcalNa22UpB;
      delete hsimNa22resUpCoA;
      delete hsimNa22resUpCoC;

      delete hsimNa22resFiltered;
      delete hcalNa22Filtered;
      delete hsimNa22resUpAFiltered;
      delete hcalNa22UpAFiltered;
      delete hsimNa22resUpBFiltered;
      delete hcalNa22UpBFiltered;
      delete hsimNa22resUpCoAFiltered;
      delete hsimNa22resUpCoCFiltered;

      delete hsimCs137res;
      delete hcalCs137;
      delete hsimCs137resUpA;
      delete hcalCs137UpA;
      delete hsimCs137resUpB;
      delete hcalCs137UpB;
      delete hsimCs137resUpCoA;
      delete hsimCs137resUpCoC;

      delete hsimCs137resFiltered;
      delete hcalCs137Filtered;
      delete hsimCs137resUpAFiltered;
      delete hcalCs137UpAFiltered;
      delete hsimCs137resUpBFiltered;
      delete hcalCs137UpBFiltered;
      delete hsimCs137resUpCoAFiltered;
      delete hsimCs137resUpCoCFiltered;

      delete hsimCo60res;
      delete hcalCo60;
      delete hsimCo60resUpA;
      delete hcalCo60UpA;
      delete hsimCo60resUpB;
      delete hcalCo60UpB;
      delete hsimCo60resUpCoA;
      delete hsimCo60resUpCoC;

      delete hsimCo60resFiltered;
      delete hcalCo60Filtered;
      delete hsimCo60resUpAFiltered;
      delete hcalCo60UpAFiltered;
      delete hsimCo60resUpBFiltered;
      delete hcalCo60UpBFiltered;
      delete hsimCo60resUpCoAFiltered;
      delete hsimCo60resUpCoCFiltered;
    }
    else{}
    
    count++;
  }

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
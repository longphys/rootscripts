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
#include "TChain.h"

//! histogram options
double xminSim = 0.;
double xmaxSim = 1.5;
double xminExp = 100;
double xmaxExp = 1500;
int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;

double ECs137 = 0.477;
double ENa22 = 1.061;
double ECo60 = 1.117;

//! ch0_9751_ch1_9800_ch2_9421
// int ChCs137[3] = {249, 239, 404};
// int ChNa22[3] = {398, 400, 787};
// int ChCo60[3] = {386, 425, 815};

// double deltaSigCs137[3] = {9.8, 8.2, 7.6};
// double deltaSigNa22[3] = {5.6, 5.4, 5.2};
// double deltaSigCo60[3] = {4.2, 4.3, 4.5};

// double ScaleCs137[3] = {5.5, 4.7, 2.9};
// double ScaleNa22[3] = {2.0, 1.45, 1.65};
// double SimScaleCo60[3] = {1.9, 1.45, 5.0};

// double CalScaleCo60[3] = {1., 1., 1.};

//! ch0_9841_ch1_9748_ch2_9702
// int ChCs137[3] = {280, 226, 281};
// int ChNa22[3] = {473, 360, 502};
// int ChCo60[3] = {486, 373, 501};

// double deltaSigCs137[3] = {7.4, 7.6, 7.6};
// double deltaSigNa22[3] = {4.9, 5.0, 5.0};
// double deltaSigCo60[3] = {4.8, 4.9, 4.9};

// double ScaleCs137[3] = {5.5, 4.7, 2.9};
// double ScaleNa22[3] = {2.0, 1.45, 1.65};
// double SimScaleCo60[3] = {2.0, 1.45, 1.65};

// double CalScaleCo60[3] = {1., 1., 1.};

//! ch0_9843_ch1_9805_ch2_9421
int ChCs137[3] = {337, 269, 298};
int ChNa22[3] = {634, 484, 597};
int ChCo60[3] = {667, 481, 608};

double deltaSigCs137[3] = {7.5, 7.7, 7.0};
double deltaSigNa22[3] = {5.2, 5.2, 4.7};
double deltaSigCo60[3] = {4.4, 4.5, 3.7};

double ScaleCs137[3] = {5.5, 4.7, 2.9};
double ScaleNa22[3] = {2.0, 1.45, 1.65};
double SimScaleCo60[3] = {2.0, 1.45, 1.65};

double CalScaleCo60[3] = {1., 1., 1.};

//! ch0_9803_ch1_9842_ch2_9854
// int ChCs137[3] = {329, 342, 362};
// int ChNa22[3] = {610, 583, 615};
// int ChCo60[3] = {613, 616, 728};

// double deltaSigCs137[3] = {7.4, 8.4, 7.6};
// double deltaSigNa22[3] = {4.8, 5.8, 5.5};
// double deltaSigCo60[3] = {4.0, 5.3, 4.8};

// double ScaleCs137[3] = {5.5, 4.7, 2.9};
// double ScaleNa22[3] = {2.0, 1.45, 1.65};
// double SimScaleCo60[3] = {2.0, 1.45, 1.65};

// double CalScaleCo60[3] = {1.05, 1., 2.4};

//! energy range to analyze timing resolution (MeV)
double min_time = 0.35;
double max_time = 0.4;

//! Channels
int channel[3] = {0, 1, 2};

void time()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* fsimCo60 = new TFile("~/data/simfiles/withAl_simCo60_1.root", "read");
  // TFile* fexpCo60 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
  // TFile* fexpCo60 = new TFile("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");
  TFile* fexpCo60 = new TFile("~/data//expfiles/time/2/plastic_detectors_Co60_ch0HV1900_ch1HV2250_ch2HV1850_ch0_9843_ch1_9805_ch2_9421_run_0_0001.root", "read");
  // TFile* fexpCo60 = new TFile("./expfiles/time/plastic_detectors_Co60_ch0HV2050_ch1HV2250_ch2HV1850_ch0_9803_ch1_9842_ch2_9854_run_0_0001.root", "read");

  TChain* chExpCo60 = new TChain("chexpCo60");
  // chExpCo60->Add("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0002.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0003.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0004.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0005.root?#AnalysisxTree");

  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0002.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0003.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0004.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0005.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0006.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0007.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0008.root?#AnalysisxTree");

  chExpCo60->Add("~/data/expfiles/time/2/plastic_detectors_Co60_ch0HV1900_ch1HV2250_ch2HV1850_ch0_9843_ch1_9805_ch2_9421_run_0_0001.root?#AnalysisxTree");
  chExpCo60->Add("~/data/expfiles/time/2/plastic_detectors_Co60_ch0HV1900_ch1HV2250_ch2HV1850_ch0_9843_ch1_9805_ch2_9421_run_0_0002.root?#AnalysisxTree");
  chExpCo60->Add("~/data/expfiles/time/2/plastic_detectors_Co60_ch0HV1900_ch1HV2250_ch2HV1850_ch0_9843_ch1_9805_ch2_9421_run_0_0003.root?#AnalysisxTree");
  chExpCo60->Add("~/data/expfiles/time/2/plastic_detectors_Co60_ch0HV1900_ch1HV2250_ch2HV1850_ch0_9843_ch1_9805_ch2_9421_run_0_0004.root?#AnalysisxTree");
  chExpCo60->Add("~/data/expfiles/time/2/plastic_detectors_Co60_ch0HV1900_ch1HV2250_ch2HV1850_ch0_9843_ch1_9805_ch2_9421_run_0_0005.root?#AnalysisxTree");
  chExpCo60->Add("~/data/expfiles/time/2/plastic_detectors_Co60_ch0HV1900_ch1HV2250_ch2HV1850_ch0_9843_ch1_9805_ch2_9421_run_0_0006.root?#AnalysisxTree");

  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV2050_ch1HV2250_ch2HV1850_ch0_9803_ch1_9842_ch2_9854_run_0_0001.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV2050_ch1HV2250_ch2HV1850_ch0_9803_ch1_9842_ch2_9854_run_0_0002.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV2050_ch1HV2250_ch2HV1850_ch0_9803_ch1_9842_ch2_9854_run_0_0003.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV2050_ch1HV2250_ch2HV1850_ch0_9803_ch1_9842_ch2_9854_run_0_0004.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV2050_ch1HV2250_ch2HV1850_ch0_9803_ch1_9842_ch2_9854_run_0_0005.root?#AnalysisxTree");
  // chExpCo60->Add("./expfiles/time/plastic_detectors_Co60_ch0HV2050_ch1HV2250_ch2HV1850_ch0_9803_ch1_9842_ch2_9854_run_0_0006.root?#AnalysisxTree");

  TTree* tsimCo60 =  (TTree*) fsimCo60->Get("dEEtree");
  TTree* texpCo60 =  (TTree*) fexpCo60->Get("AnalysisxTree");

  double a[3], b[3], xminCal[3], xmaxCal[3], coA[3], coB[3], coC[3];

  TH1D* hcalCo60[3];
  TH1D* hsimCo60[3];

  UShort_t chX[48];
  texpCo60->SetBranchAddress("NeEvent.neutAmp[48]", chX);
  Long64_t entriesExp = texpCo60->GetEntries();

  char* nameCal = new char[20];
  char* nameSim = new char[20];

  double x_sim;
  tsimCo60->SetBranchAddress("Scintillator", &x_sim);
  TRandom3* ranGen = new TRandom3();
  Long64_t entriesSim = tsimCo60->GetEntries();

  TCanvas* cFit = new TCanvas("cFit", "Linear calibration", 800, 600);
  for(int i = 0; i < 3; i++){
    TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
    hCali->SetBinContent(ChCs137[i] - xminExp, ECs137);
    hCali->SetBinContent(ChNa22[i] - xminExp, ENa22);
    hCali->SetBinContent(ChCo60[i] - xminExp, ECo60);

    TF1* fLinear = new TF1("fLinear", "[0]*x + [1]", xminExp, xmaxExp);

    cFit->Divide(2,1);
    cFit->cd(1);
    std::cout << "\nENERGY CALIBRATION\n";
    hCali->Fit("fLinear");

    a[i] = fLinear->GetParameter(0);
    b[i] = fLinear->GetParameter(1);

    xminCal[i] = xminExp*a[i]+b[i];
    xmaxCal[i] = xmaxExp*a[i]+b[i];

    sprintf(nameCal,"hcalCo60_Ch%d",i); 

    hcalCo60[i] = new TH1D(nameCal, nameCal, binExp, xminCal[i], xmaxCal[i]);
    hcalCo60[i]->GetXaxis()->SetTitle("Energy(MeV)");
    hcalCo60[i]->GetYaxis()->SetTitle("Count");

    for(int k = 0; k < entriesExp; k++){
    // for(int k = 0; k < 500000; k++){
      texpCo60->GetEntry(k);
      hcalCo60[i]->Fill(a[i]*(chX[channel[i]]+0.5) + b[i]);
    }

    TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal[i], xmaxCal[i]);
    hRes->SetBinContent(hcalCo60[i]->FindBin(ECs137), deltaSigCs137[i]/100.);
    hRes->SetBinContent(hcalCo60[i]->FindBin(ENa22), deltaSigNa22[i]/100.);
    hRes->SetBinContent(hcalCo60[i]->FindBin(ECo60), deltaSigCo60[i]/100.);

    TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");

    cFit->cd(2);
    hRes->Fit("fRes", "", "", ECs137-0.2, ECo60+0.2);

    coA[i] = fRes->GetParameter(0);
    coB[i] = 0.;
    coC[i] = fRes->GetParameter(1);

    sprintf(nameSim,"hsimCo60_Ch%d",i); 

    hsimCo60[i] = new TH1D(nameSim, nameSim, binExp, xminCal[i], xmaxCal[i]);
    hsimCo60[i]->GetXaxis()->SetTitle("Energy(MeV)");
    hsimCo60[i]->GetYaxis()->SetTitle("Count");

    for (int k = 0; k < entriesSim; k++)
    {
      tsimCo60->GetEntry(k);
      double sigma = x_sim*sqrt( pow(coA[i],2) + pow(coB[i]/sqrt(x_sim),2) + pow(coC[i]/x_sim,2) );
      hsimCo60[i]->Fill(ranGen->Gaus(x_sim,sigma));
    }
  
    delete hCali;
    delete fLinear;
    delete hRes;
    delete fRes;
  }

  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  c1->cd();
  c1->Divide(1, 1);

  c1->cd(1);
  hcalCo60[0]->SetLineColor(kBlack);
  hcalCo60[1]->SetLineColor(kRed);
  hcalCo60[2]->SetLineColor(kBlue);
  
  hcalCo60[0]->Scale(CalScaleCo60[0], "noSW2");
  hcalCo60[1]->Scale(CalScaleCo60[1], "noSW2");
  hcalCo60[2]->Scale(CalScaleCo60[2], "noSW2");

  TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg1->SetHeader("Co60", "C");
  leg1->SetBorderSize(2);

  for(int i = 0; i<3; i++){
    hcalCo60[i]->Draw("same");
    sprintf(nameCal,"Channel %d",i);
    leg1->AddEntry(hcalCo60[i], nameCal, "l");
  }
  leg1->Draw();

  TCanvas* c3 = new TCanvas("c3", "c3", 1500, 600);
  c3->cd();
  c3->Divide(3,1);

  hsimCo60[0]->Scale(SimScaleCo60[0], "no SW2");
  hsimCo60[1]->Scale(SimScaleCo60[1], "no SW2");
  hsimCo60[2]->Scale(SimScaleCo60[2], "no SW2");
  for(int i = 1; i <= 3; i++){
    c3->cd(i);
    hsimCo60[i-1]->Draw();
    hcalCo60[i-1]->Draw("same");
  }

  //! Time
  UShort_t chY[48];
  texpCo60->SetBranchAddress("NeEvent.neutTDC[48]", chY);

  TH1D* hTimeCo60[3];
  TH1D* hTimeSubCo60[3];

  std::cout << "a[0] = " << a[0] << "; b[0] = " << b[0] <<"\n";
  std::cout << "a[1] = " << a[1] << "; b[1] = " << b[1] <<"\n";
  std::cout << "a[2] = " << a[2] << "; b[2] = " << b[2] <<"\n";

  TCanvas* c2 = new TCanvas("c2", "c2", 1500, 600);
  c2->Divide(3,1);

  TH1D* hTime1 = new TH1D("hTime1", "hTime1", 10000, -5000, 5000);
  TH1D* hTime2 = new TH1D("hTime2", "hTime2", 10000, -5000, 5000);
  TH1D* hTime3 = new TH1D("hTime3", "hTime3", 10000, -5000, 5000);
  hTime3->SetLineWidth(3);
  // hTime3->SetStats(0);

  hTime3->GetYaxis()->SetTitle("Events");
  hTime3->GetYaxis()->SetLabelFont(42);
  hTime3->GetYaxis()->SetTitleFont(52);
  hTime3->GetYaxis()->SetTitleSize(0.04);
  hTime3->GetYaxis()->CenterTitle(true);

  hTime3->GetXaxis()->SetTitle("Time difference (a. units)");
  hTime3->GetXaxis()->SetLabelFont(42);
  hTime3->GetXaxis()->SetTitleFont(52);
  hTime3->GetXaxis()->SetTitleSize(0.04);
  hTime3->GetXaxis()->CenterTitle(true);
  TF1* fGaus = new TF1("fGaus","gaus");

  auto a0 = std::to_string(a[0]);
  auto b0 = std::to_string(-b[0]);
  auto a1 = std::to_string(a[1]);
  auto b1 = std::to_string(-b[1]);
  auto a2 = std::to_string(a[2]);
  auto b2 = std::to_string(-b[2]);

  auto min_time_string = std::to_string(min_time);
  auto max_time_string = std::to_string(max_time);

  std::string prompt0_full_str = 
  "NeEvent.neutAmp[0] > (" + min_time_string + "+" + b0 + ")/" + a0 + " && " + 
  "NeEvent.neutAmp[0] < (" + max_time_string + "+" + b0 + ")/" + a0 + " && " + 
  "NeEvent.neutAmp[1] > (" + min_time_string + "+" + b1 + ")/" + a1 + " && " + 
  "NeEvent.neutAmp[1] < (" + max_time_string + "+" + b1 + ")/" + a1;
  const char* prompt0_full = prompt0_full_str.c_str();

  std::string prompt1_full_str = 
  "NeEvent.neutAmp[1] > (" + min_time_string + "+" + b1 + ")/" + a1 + " && " + 
  "NeEvent.neutAmp[1] < (" + max_time_string + "+" + b1 + ")/" + a1 + " && " + 
  "NeEvent.neutAmp[2] > (" + min_time_string + "+" + b2 + ")/" + a2 + " && " + 
  "NeEvent.neutAmp[2] < (" + max_time_string + "+" + b2 + ")/" + a2;
  const char* prompt1_full = prompt1_full_str.c_str();

  std::string prompt2_full_str = 
  "NeEvent.neutAmp[2] > (" + min_time_string + "+" + b2 + ")/" + a2 + " && " + 
  "NeEvent.neutAmp[2] < (" + max_time_string + "+" + b2 + ")/" + a2 + " && " + 
  "NeEvent.neutAmp[0] > (" + min_time_string + "+" + b0 + ")/" + a0 + " && " + 
  "NeEvent.neutAmp[0] < (" + max_time_string + "+" + b0 + ")/" + a0;
  const char* prompt2_full = prompt2_full_str.c_str();

  TCanvas* CanvasTime = new TCanvas("CanvasTime", "CanvasTime", 1100, 1000);
  CanvasTime->SetLeftMargin(0.15);
  CanvasTime->cd();
  chExpCo60->Draw("NeEvent.neutTDC[2]-NeEvent.neutTDC[0]>>hTime3", prompt2_full, "");
  hTime3->Fit("fGaus", "", "", -200., 200.);

  c2->cd(1);
  chExpCo60->Draw("NeEvent.neutTDC[0]-NeEvent.neutTDC[1]>>hTime1", prompt0_full, "");
  hTime1->Fit("fGaus", "", "", -200., 200.);
  double timeResCh[3];
  timeResCh[0] = fGaus->GetParameter(2);

  c2->cd(2);
  chExpCo60->Draw("NeEvent.neutTDC[1]-NeEvent.neutTDC[2]>>hTime2", prompt1_full, "");
  hTime2->Fit("fGaus", "", "", -200., 200.);
  timeResCh[1] = fGaus->GetParameter(2);

  c2->cd(3);
  chExpCo60->Draw("NeEvent.neutTDC[2]-NeEvent.neutTDC[0]>>hTime3", prompt2_full, "");
  hTime3->Fit("fGaus", "", "", -200., 200.);
  timeResCh[2] = fGaus->GetParameter(2);

  std::cout << "Time resolution 0 - 1 = " << timeResCh[0]*31.25 << "(ps)\n";
  std::cout << "Time resolution 1 - 2 = " << timeResCh[1]*31.25 << "(ps)\n";
  std::cout << "Time resolution 2 - 0 = " << timeResCh[2]*31.25 << "(ps)\n";

  double A = pow(timeResCh[0]*31.25,2);
  double B = pow(timeResCh[1]*31.25,2);
  double C = pow(timeResCh[2]*31.25,2);

  double x = (A-B+C)/2;
  double z = C-x;
  double y = B-z;

  double timeRes0 = sqrt(x);
  double timeRes1 = sqrt(y);
  double timeRes2 = sqrt(z);

  std::cout << "Time resolution Channel 0 = " << timeRes0 << "(ps)\n";
  std::cout << "Time resolution Channel 1 = " << timeRes1 << "(ps)\n";
  std::cout << "Time resolution Channel 2 = " << timeRes2 << "(ps)\n";

  std::cout << "Energy resolution at 1MeV Channel 0 = " << 
  sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(1.0),2) + pow(coC[0]/1.0,2) )*100 << "(%)\n";
  std::cout << "Energy resolution at 1MeV Channel 1 = " << 
  sqrt( pow(coA[1],2) + pow(coB[1]/sqrt(1.0),2) + pow(coC[1]/1.0,2) )*100 << "(%)\n";
  std::cout << "Energy resolution at 1MeV Channel 2 = " << 
  sqrt( pow(coA[2],2) + pow(coB[2]/sqrt(1.0),2) + pow(coC[2]/1.0,2) )*100 << "(%)\n";

  int threshChannel = 180;
  std::cout << "Energy Threshold = " << a[2]*threshChannel + b[2] << "(keV)\n";

  std::cout << "CURRENT ENERGY RANGE FOR TIMING RESOLUTION: " << min_time << "-" << max_time << " (MeV)\n";

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
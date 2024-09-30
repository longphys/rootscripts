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

int ChCs137[3] = {277, 226, 285};
int ChNa22[3] = {467, 360, 491};
int ChCo60[3] = {486, 373, 511};
double deltaSigCs137[3] = {7.4, 7.6, 7.6};
double deltaSigNa22[3] = {4.9, 5.0, 5.0};
double deltaSigCo60[3] = {4.8, 4.9, 4.9};
double ScaleCs137[3] = {5.5, 4.7, 2.9};
double ScaleNa22[3] = {2.0, 1.45, 1.65};
double ScaleCo60[3] = {2.0, 1.45, 1.65};

int channel[3] = {0, 1, 2};

void compareCo60()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* fsimCo60 = new TFile("./simfiles/withAl_simCo60_1.root", "read");
  // TFile* fexpCo60 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
  TFile* fexpCo60 = new TFile("./expfiles/time/plastic_detectors_Co60_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");


  TTree* tsimCo60 =  (TTree*) fsimCo60->Get("dEEtree");
  TTree* texpCo60 =  (TTree*) fexpCo60->Get("AnalysisxTree");

  double a[3], b[3], xminCal[3], xmaxCal[3], coA[3], coB[3], coC[3];

  TH1D* hcalCo60[3];
  TH1D* hsimCo60res = new TH1D("hsimCo60res", "Co60 Simulation with resolution", binExp, xminCal[0], xmaxCal[0]);
  hsimCo60res->GetXaxis()->SetTitle("Energy(MeV)");
  hsimCo60res->GetYaxis()->SetTitle("Count");

  UShort_t chX[48];
  texpCo60->SetBranchAddress("NeEvent.neutAmp[48]", chX);
  Long64_t entries = texpCo60->GetEntries();

  char* name = new char[20];

  TCanvas* cFit = new TCanvas("cFit", "Linear calibration", 800, 600);
  for(int i = 0; i < 3; i++){
    TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
    hCali->SetBinContent(ChNa22[i] - xminExp, ENa22);
    hCali->SetBinContent(ChCs137[i] - xminExp, ECs137);
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

    sprintf(name,"hcalCo60_Ch%d",i); 

    hcalCo60[i] = new TH1D(name, "Co60 Experiment Calibrated Ch0", binExp, xminCal[i], xmaxCal[i]);
    hcalCo60[i]->GetXaxis()->SetTitle("Energy(MeV)");
    hcalCo60[i]->GetYaxis()->SetTitle("Count");

    for(int k = 0; k < entries; k++){
      texpCo60->GetEntry(k);
      hcalCo60[i]->Fill(a[i]*(chX[channel[i]]+0.5) + b[i]);
    }

    TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal[i], xmaxCal[i]);
    hRes->SetBinContent(hcalCo60[i]->FindBin(ENa22), deltaSigNa22[i]/100.);
    hRes->SetBinContent(hcalCo60[i]->FindBin(ECs137), deltaSigCs137[i]/100.);
    hRes->SetBinContent(hcalCo60[i]->FindBin(ECo60), deltaSigCo60[i]/100.);

    TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");

    cFit->cd(2);
    hRes->Fit("fRes", "", "", ENa22-0.2, ECo60+0.2);

    coA[i] = fRes->GetParameter(0);
    coB[i] = 0.;
    coC[i] = fRes->GetParameter(1);
 
    delete hCali;
    delete fLinear;
    delete hRes;
    delete fRes;
  }

  double x_sim;

  tsimCo60->SetBranchAddress("Scintillator", &x_sim);

  TRandom3* ranGen = new TRandom3();

  entries = tsimCo60->GetEntries();

  for (int i = 0; i < entries; i++)
  {
    tsimCo60->GetEntry(i);
    double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
    hsimCo60res->Fill(ranGen->Gaus(x_sim,sigma));
  }

  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  c1->cd();
  c1->Divide(1, 1);

  c1->cd(1);
  hcalCo60[0]->SetLineColor(kBlack);
  hcalCo60[1]->SetLineColor(kRed);
  hcalCo60[2]->SetLineColor(kBlue);
  
  hcalCo60[0]->Scale(1., "noSW2");
  hcalCo60[1]->Scale(1.5, "noSW2");
  hcalCo60[2]->Scale(4.0, "noSW2");

  TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg1->SetHeader("Co60", "C");
  leg1->SetBorderSize(2);

  for(int i = 0; i<3; i++){
    hcalCo60[i]->Draw("same");
    sprintf(name,"Channel %d",i);
    leg1->AddEntry(hcalCo60[i], name, "l");
  }
  leg1->Draw();

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
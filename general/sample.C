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

void sample()
{
  auto timer = new TStopwatch();
  timer->Start();

  TFile* fexpCs137 = new TFile("./expfiles/time/plastic_detectors_Cs137_end_ch0HV1800_ch1HV2350_ch2HV2150_ch0_9841_ch1_9748_ch2_9702_run_0_0001.root", "read");

  TTree* texpCs137 =  (TTree*) fexpCs137->Get("AnalysisxTree");

  TH1D* hsimCs137res = new TH1D("hsimCs137res", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
  TH1D* hcalCs137 = new TH1D("hcalCs137", "Cs137 Experiment Calibrated", binExp, xminCal, xmaxCal);

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

  //! Resolution
  // for (int i = 0; i < entries; i++)
  for (int i = 0; i < 1000000; i++)
  {
    tsimCs137->GetEntry(i);
    double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
    hsimCs137res->Fill(ranGen->Gaus(x_sim,sigma));
  }
  
  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  c1->cd();

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

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
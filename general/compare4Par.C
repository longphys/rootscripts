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
double xmaxSim = 1.3;
double xminExp = 100;
double xmaxExp = 1500;
int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;

double E1 = 0.477;
double E2 = 1.061;

//! Zoomed in Histograms options
double minTestCs137 = 0.22; //MeV
double maxTestCs137 = 0.7;

double minTestNa22 = 0.6; //MeV
double maxTestNa22 = 1.4;

//! ENERGY CALIBRATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//? new active dividers channel 1
// double Ch1 = 240;
// double Ch2 = 418;

//? new active dividers channel 2
// double Ch1 = 316;
// double Ch2 = 563;

//? stilbene dividers channel 2
// double Ch1 = 315;
// double Ch2 = 563;

//? NEW stilbene dividers channel 0
// double Ch1 = 246;
// double Ch2 = 393;

//? NEW stilbene dividers channel 1
// double Ch1 = 248;
// double Ch2 = 420;

double Ch1 = 260;
double Ch2 = 440;

// double Ch1 = 148;
// double Ch2 = 520;

//? NEW stilbene dividers channel 2
// double Ch1 = 402;
// double Ch2 = 785;

//! RESOLUTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//? new active dividers channel 1
// double sig1 = 0.037;
// double sig2 = 0.040;

//? new active dividers channel 2
// double sig1 = 0.037;
// double sig2 = 0.057;

//? stilbene dividers channel 2
// double sig1 = 0.040;
// double sig2 = 0.058;

//? NEW stilbene dividers channel 0
// double sig1 = 0.052;
// double sig2 = 0.078;

//? NEW stilbene dividers channel 1
// double deltaSig1 = 0.0734;
// double deltaSig2 = 0.0518;

double deltaSig1 = 0.08;
double deltaSig2 = 0.04;

//? NEW stilbene dividers channel 2
// double sig1 = 0.025;
// double sig2 = 0.052;

//! SCALING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//? new active dividers channel 1
// double ScaleCs137 = 2.6;
// double ScaleNa22 = 1.6;

//? new active dividers channel 2
// double ScaleCs137 = 4.5;
// double ScaleNa22 = 2.3;

//? stilbene dividers channel 2
// double ScaleCs137 = 4.5;
// double ScaleNa22 = 2.3;

// //? NEW stilbene dividers channel 0
// double ScaleCs137 = 5.7;
// double ScaleNa22 = 2.7;

// //? NEW stilbene dividers channel 1
double ScaleCs137 = 4.7;
double ScaleNa22 = 1.9;

//? NEW stilbene dividers channel 2
// double ScaleCs137 = 2.9;
// double ScaleNa22 = 2.2;

//! CHANNEL
// int channel = 0;
int channel = 1;
// int channel = 2;

void compare4Par()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* fsimCs137 = new TFile("./simfiles/withAl_simCs137_1.root", "read");
  TFile* fexpCs137 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Cs137_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");

  TTree* tsimCs137 =  (TTree*) fsimCs137->Get("dEEtree");
  TTree* texpCs137 =  (TTree*) fexpCs137->Get("AnalysisxTree");

  TFile* fsimNa22 = new TFile("./simfiles/withAl_simNa22_1.root", "read");
  TFile* fexpNa22 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Na22_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");

  TTree* tsimNa22 =  (TTree*) fsimNa22->Get("dEEtree");
  TTree* texpNa22 =  (TTree*) fexpNa22->Get("AnalysisxTree");

  UShort_t chX[48];
  texpCs137->SetBranchAddress("NeEvent.neutAmp[48]", chX);
  Long64_t entriesExpCs137 = texpCs137->GetEntries();

  UShort_t chY[48];
  texpNa22->SetBranchAddress("NeEvent.neutAmp[48]", chY);
  Long64_t entriesExpNa22 = texpNa22->GetEntries();

  double x_sim;
  tsimCs137->SetBranchAddress("Scintillator", &x_sim);
  Long64_t entriesSimCs137 = tsimCs137->GetEntries();

  double y_sim;
  tsimNa22->SetBranchAddress("Scintillator", &y_sim);
  Long64_t entriesSimNa22 = tsimNa22->GetEntries();

  TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
  hCali->SetBinContent(Ch1 - xminExp, E1);
  hCali->SetBinContent(Ch2 - xminExp, E2);

  TF1* fLinear = new TF1("fLinear", "[0]*x + [1]", xminExp, xmaxExp);
  hCali->Fit("fLinear");

  double a = fLinear->GetParameter(0);
  double b = fLinear->GetParameter(1);

  double aStep = 0.1*a;
  double bStep = 0.1*b;

  //! Resolution for simulation histograms

  std::cout << "delta1 = " << deltaSig1*100 << "%; delta2 = " << deltaSig2*100 << "%\n";

  double xminCal = xminExp*a+b;
  double xmaxCal = xmaxExp*a+b;

  TH1D* hFirstCalCs137 = new TH1D("hFirstCalCs137", "Cs137 Experiment Calibrated", binExp, xminCal, xmaxCal);
  for(int i = 0; i < entriesExpCs137; i++){
    texpCs137->GetEntry(i);
    hFirstCalCs137->Fill(a*(chX[channel]+0.5) + b);
  }

  TH1D* hFirstCalNa22 = new TH1D("hFirstCalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);
  for(int i = 0; i < entriesExpNa22; i++){
    texpNa22->GetEntry(i);
    hFirstCalNa22->Fill(a*(chY[channel]+0.5) + b);
  }

  TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal, xmaxCal);
  hRes->SetBinContent(hFirstCalCs137->FindBin(E1), deltaSig1);
  hRes->SetBinContent(hFirstCalNa22->FindBin(E2), deltaSig2);

  //! Choose fit function
  // TF1* fRes = new TF1("fRes", "sqrt([0]*[0]/x + [1]*[1]/(x*x))");
  TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");

  hRes->Fit("fRes", "", "", E1-0.2, E2+0.2);

  //! Correspond to chosen function
  double coA = std::abs(fRes->GetParameter(0));
  double coB = 0.;
  double coC = std::abs(fRes->GetParameter(1));

  double coAStep = coA*0.2;
  double coBStep = coB*0.2;
  double coCStep = coC*0.2;

  int tuningTimes = 10;
  double tuningRate = 0.00001;
  
  TRandom3* ranGen = new TRandom3();

  for (int i = 1; i<=tuningTimes; i++){

    xminCal = xminExp*a+b;
    xmaxCal = xmaxExp*a+b;

    double xminCalUpA = xminExp*(a + aStep)+b;
    double xmaxCalUpA = xmaxExp*(a + aStep)+b;

    double xminCalUpB = xminExp*a+(b + bStep);
    double xmaxCalUpB = xmaxExp*a+(b + bStep);

    //! Histograms
    TH1D* hsimCs137res = new TH1D("hsimCs137res", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalCs137 = new TH1D("hcalCs137", "Cs137 Experiment Calibrated", binExp, xminCal, xmaxCal);
    
    TH1D* hsimNa22res = new TH1D("hsimNa22res", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalNa22 = new TH1D("hcalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);

    TH1D* hsimCs137resUpA = new TH1D("hsimCs137resUpA", "Cs137 Simulation with resolution", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalCs137UpA = new TH1D("hcalCs137UpA", "Cs137 Experiment Calibrated", binExp, xminCalUpA, xmaxCalUpA);
    
    TH1D* hsimNa22resUpA = new TH1D("hsimNa22resUpA", "Na22 Simulation with resolution", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalNa22UpA = new TH1D("hcalNa22UpA", "Na22 Experiment Calibrated", binExp, xminCalUpA, xmaxCalUpA);

    TH1D* hsimCs137resUpB = new TH1D("hsimCs137resUpB", "Cs137 Simulation with resolution", binExp, xminCalUpB, xmaxCalUpB);
    TH1D* hcalCs137UpB = new TH1D("hcalCs137UpB", "Cs137 Experiment Calibrated", binExp, xminCalUpB, xmaxCalUpB);
    
    TH1D* hsimNa22resUpB = new TH1D("hsimNa22resUpB", "Na22 Simulation with resolution", binExp, xminCalUpB, xmaxCalUpB);
    TH1D* hcalNa22UpB = new TH1D("hcalNa22UpB", "Na22 Experiment Calibrated", binExp, xminCalUpB, xmaxCalUpB);

    TH1D* hsimCs137resUpCoA = new TH1D("hsimCs137resUpCoA", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hsimNa22resUpCoA = new TH1D("hsimNa22resUpCoA", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    TH1D* hsimCs137resUpCoC = new TH1D("hsimCs137resUpCoC", "Cs137 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hsimNa22resUpCoC = new TH1D("hsimNa22resUpCoC", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    //! Fill experiment histograms
    for(int i = 0; i < entriesExpCs137; i++){
      texpCs137->GetEntry(i);
      hcalCs137->Fill(a*(chX[channel]+0.5) + b);
    }

    for(int i = 0; i < entriesExpNa22; i++){
      texpNa22->GetEntry(i);
      hcalNa22->Fill(a*(chY[channel]+0.5) + b);
    }

    for(int i = 0; i < entriesExpCs137; i++){
      texpCs137->GetEntry(i);
      hcalCs137UpA->Fill((a+aStep)*(chX[channel]+0.5) + b);
    }

    for(int i = 0; i < entriesExpNa22; i++){
      texpNa22->GetEntry(i);
      hcalNa22UpA->Fill((a+aStep)*(chY[channel]+0.5) + b);
    }

    for(int i = 0; i < entriesExpCs137; i++){
      texpCs137->GetEntry(i);
      hcalCs137UpB->Fill(a*(chX[channel]+0.5) + (b+bStep));
    }

    for(int i = 0; i < entriesExpNa22; i++){
      texpNa22->GetEntry(i);
      hcalNa22UpB->Fill(a*(chY[channel]+0.5) + (b+bStep));
    }

    //! Fill Simulation histograms
    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
      hsimCs137res->Fill(ranGen->Gaus(x_sim,sigma));
    }
    
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      hsimNa22res->Fill(ranGen->Gaus(y_sim,sigma));
    }

    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
      hsimCs137resUpA->Fill(ranGen->Gaus(x_sim,sigma));
    }
    
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      hsimNa22resUpA->Fill(ranGen->Gaus(y_sim,sigma));
    }

    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
      hsimCs137resUpB->Fill(ranGen->Gaus(x_sim,sigma));
    }
    
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      hsimNa22resUpB->Fill(ranGen->Gaus(y_sim,sigma));
    }

    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA+coAStep,2) + pow(coB/sqrt(x_sim),2) + pow(coC/x_sim,2) );
      hsimCs137resUpCoA->Fill(ranGen->Gaus(x_sim,sigma));
    }
    
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA+coAStep,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      hsimNa22resUpCoA->Fill(ranGen->Gaus(y_sim,sigma));
    }

    for (int i = 0; i < entriesSimCs137; i++)
    {
      tsimCs137->GetEntry(i);
      double sigma = x_sim*sqrt( pow(coA,2) + pow(coB/sqrt(x_sim),2) + pow((coC+coCStep)/x_sim,2) );
      hsimCs137resUpCoC->Fill(ranGen->Gaus(x_sim,sigma));
    }
    
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow((coC+coCStep)/y_sim,2) );
      hsimNa22resUpCoC->Fill(ranGen->Gaus(y_sim,sigma));
    }

    hcalCs137->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137res->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hcalNa22->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22res->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);

    double chi2Cs137 = hcalCs137->Chi2Test(hsimCs137res, "UU CHI2/NDF");
    double chi2Na22 = hcalNa22->Chi2Test(hsimNa22res, "UU CHI2/NDF");

    hcalCs137UpA->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137resUpA->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hcalNa22UpA->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22resUpA->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);

    double chi2Cs137UpA = hcalCs137UpA->Chi2Test(hsimCs137resUpA, "UU CHI2/NDF");
    double chi2Na22UpA = hcalNa22UpA->Chi2Test(hsimNa22resUpA, "UU CHI2/NDF");

    hcalCs137UpB->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137resUpB->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hcalNa22UpB->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22resUpB->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);

    double chi2Cs137UpB = hcalCs137UpB->Chi2Test(hsimCs137resUpB, "UU CHI2/NDF");
    double chi2Na22UpB = hcalNa22UpB->Chi2Test(hsimNa22resUpB, "UU CHI2/NDF");

    hcalCs137->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137resUpCoA->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hcalNa22->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22resUpCoA->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);

    double chi2Cs137UpCoA = hcalCs137->Chi2Test(hsimCs137resUpCoA, "UU CHI2/NDF");
    double chi2Na22UpCoA = hcalNa22->Chi2Test(hsimNa22resUpCoA, "UU CHI2/NDF");
    
    hcalCs137->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hsimCs137resUpCoC->GetXaxis()->SetRangeUser(minTestCs137, maxTestCs137);
    hcalNa22->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);
    hsimNa22resUpCoC->GetXaxis()->SetRangeUser(minTestNa22, maxTestNa22);

    double chi2Cs137UpCoC = hcalCs137->Chi2Test(hsimCs137resUpCoC, "UU CHI2/NDF");
    double chi2Na22UpCoC = hcalNa22->Chi2Test(hsimNa22resUpCoC, "UU CHI2/NDF");

    double devChi2Cs137UpA = (chi2Cs137UpA-chi2Cs137)/aStep;
    double devChi2Na22UpA = (chi2Na22UpA-chi2Na22)/aStep;
    double devChi2Cs137UpB = (chi2Cs137UpB-chi2Cs137)/bStep;
    double devChi2Na22UpB = (chi2Na22UpB-chi2Na22)/bStep;
    double devChi2Cs137UpCoA = (chi2Cs137UpCoA-chi2Cs137)/coAStep;
    double devChi2Na22UpCoA = (chi2Na22UpCoA-chi2Na22)/coAStep;
    double devChi2Cs137UpCoC = (chi2Cs137UpCoC-chi2Cs137)/coCStep;
    double devChi2Na22UpCoC = (chi2Na22UpCoC-chi2Na22)/coCStep;

    // a = a - tuningRate*devChi2Cs137UpA;
    // b = b - tuningRate*devChi2Cs137UpB;
    // coA = coA - tuningRate*devChi2Cs137UpCoA;
    // coB = coB - tuningRate*devChi2Cs137UpCoC;

    a = a - tuningRate*devChi2Na22UpA;
    b = b - tuningRate*devChi2Na22UpB;
    coA = coA - tuningRate*devChi2Na22UpCoA;
    coB = coB - tuningRate*devChi2Na22UpCoC;

    double finalChannel1 = (E1-b)/a;
    double finalChannel2 = (E2-b)/a;

    double finalSig1 = 100.*sqrt(coA*coA + coC*coC/(E1*E1));
    double finalSig2 = 100.*sqrt(coA*coA + coC*coC/(E2*E2));

    std::cout << "Ch1 = " << finalChannel1 << "; Ch2 = " << finalChannel2
    << "\nSig1 = " << finalSig1 << "; Sig2 = " << finalSig1 << "\n";
    // delete hsimCs137res;
    // delete hcalCs137;
    // delete hsimNa22res;
    // delete hcalNa22;
  }

  //! Canvas and draw
  // TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  // c1->cd();
  // c1->Divide(1, 2);

  // c1->cd(1);
  // hcalCs137->SetLineColor(kRed);
  // hcalCs137->Draw();
  // hsimCs137res->Scale(ScaleCs137, "noSW2");
  // hsimCs137res->Draw("same");

  // TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  // leg1->SetHeader("Cs137", "C");
  // leg1->SetBorderSize(2);
  // leg1->AddEntry(hsimCs137res, "simulation", "l");
  // leg1->AddEntry(hcalCs137, "experiment", "l");
  // leg1->Draw();

  // c1->cd(2);
  // hcalNa22->SetLineColor(kRed);
  // // hcalNa22_511->Scale(0.34, "noSW2");
  // hcalNa22->Draw();
  // hsimNa22res->Scale(ScaleNa22, "noSW2");
  // hsimNa22res->Draw("same");

  // TLegend *leg2 = new TLegend(0.75, 0.6, 0.98, 0.75);
  // leg2->SetHeader("Na22", "C");
  // leg2->SetBorderSize(2);
  // leg2->AddEntry(hsimNa22res, "simulation", "l");
  // leg2->AddEntry(hcalNa22, "experiment", "l");
  // leg2->Draw();

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
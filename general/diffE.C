#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"

void diffE()
{
  //! Files and trees
  TFile* fsimCs137 = new TFile("withAl_simCs137_new2.root", "read");
  TFile* fsimNa22 = new TFile("withAl_simNa22_new2.root", "read");
  TFile* fsimCo60 = new TFile("withAl_simCo60_new2.root", "read");

  TFile* fexpCs137 = new TFile("old_divider_Cs137_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
  TFile* fexpNa22 = new TFile("old_divider_Na22_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");
  TFile* fexpCo60 = new TFile("old_divider_Co60_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");

  TTree* tsimCs137 =  (TTree*) fsimCs137->Get("dEEtree");
  TTree* tsimNa22 =  (TTree*) fsimNa22->Get("dEEtree");
  TTree* tsimCo60 =  (TTree*) fsimCo60->Get("dEEtree");
  
  TTree* texpCs137 =  (TTree*) fexpCs137->Get("AnalysisxTree");
  TTree* texpNa22 =  (TTree*) fexpNa22->Get("AnalysisxTree");
  TTree* texpCo60 =  (TTree*) fexpCo60->Get("AnalysisxTree");

  //! Calibration

  double a = 0.00274; //fourth
  double b = -0.25;

  //! Resolution coefficients
  
  double coA = 0.0249800419143;
  double coB = 0.;
  double coC = 0.018076222210897;

  //! histogram options
  int binSim = 1300;
  int binExp = 700;
  double xminSim = 0.;
  double xmaxSim = 1.3;
  double xminExp = 100;
  double xmaxExp = 800;
  double xminCal = xminExp*a+b;
  double xmaxCal = xmaxExp*a+b;

  TH1D* hsimCs137 = new TH1D("hsimCs137", "Cs137 Simulation", binExp, xminCal, xmaxCal);
  TH1D* hsimNa22 = new TH1D("hsimNa22", "Na22 Simulation", binExp, xminCal, xmaxCal);
  TH1D* hsimCo60 = new TH1D("hsimCo60", "Co60 Simulation", binExp, xminCal, xmaxCal);

  TH1D* hexpCs137 = new TH1D("hexpCs137", "Cs137 Experiment", binExp, xminExp, xmaxExp);
  TH1D* hexpNa22 = new TH1D("hexpNa22", "Na22 Experiment", binExp, xminExp, xmaxExp);
  TH1D* hexpCo60 = new TH1D("hexpCo60", "Co60 Experiment", binExp, xminExp, xmaxExp);

  TH1D* hcalCs137 = new TH1D("hCalexpCs137", "Cs137 Experiment Calibrated", binExp, xminCal, xmaxCal);
  TH1D* hcalNa22 = new TH1D("hcalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);
  TH1D* hcalCo60 = new TH1D("hcalCo60", "Co60 Experiment Calibrated", binExp, xminCal, xmaxCal);

  hsimCs137->GetXaxis()->SetTitle("Energy(MeV)");
  hsimCs137->GetYaxis()->SetTitle("Count");

  hsimNa22->GetXaxis()->SetTitle("Energy(MeV)");
  hsimNa22->GetYaxis()->SetTitle("Count");

  hsimNa22->GetXaxis()->SetTitle("Energy(MeV)");
  hsimNa22->GetYaxis()->SetTitle("Count");

  //! Fill experiment histograms
  UShort_t x_exp[48], y_exp[48], z_exp[48];

  texpCs137->SetBranchAddress("NeEvent.neutAmp[48]", x_exp);
  texpNa22->SetBranchAddress("NeEvent.neutAmp[48]", y_exp);
  texpCo60->SetBranchAddress("NeEvent.neutAmp[48]", z_exp);

  Long64_t entries = texpCs137->GetEntries();

  for(int i = 0; i < entries; i++){
    texpCs137->GetEntry(i);
    hcalCs137->Fill(a*(x_exp[0]+0.5) + b);
  }
  
  entries = texpNa22->GetEntries();

  for(int i = 0; i < entries; i++){
    texpNa22->GetEntry(i);
    hcalNa22->Fill(a*(y_exp[0]+0.5) + b);
  }

  entries = texpCo60->GetEntries();

  for(int i = 0; i < entries; i++){
    texpCo60->GetEntry(i);
    hcalCo60->Fill(a*(z_exp[0]+0.5) + b);
  }

  //! Resolution for simulation histograms
  double x_sim, y_sim, z_sim;

  tsimCs137->SetBranchAddress("Scintillator", &x_sim);
  tsimNa22->SetBranchAddress("Scintillator", &y_sim);
  tsimCo60->SetBranchAddress("Scintillator", &z_sim);

  TRandom3* ranGen = new TRandom3();

  entries = tsimCs137->GetEntries();

  for (int i = 0; i < entries; i++)
  {
    tsimCs137->GetEntry(i);
    double sigma = x_sim*sqrt( pow(coA, 2) + pow(coB/sqrt(x_sim), 2) + pow(coC/x_sim, 2) );
    if(ranGen->Gaus(x_sim, sigma) <= 0.00001)
    {
      continue;
    }
    hsimCs137->Fill(ranGen->Gaus(x_sim, sigma));
  }

  entries = tsimNa22->GetEntries();

  for (int i = 0; i < entries; i++)
  {
    tsimNa22->GetEntry(i);
    double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
    // sigma = 0.;
    hsimNa22->Fill(ranGen->Gaus(y_sim,sigma));
  }

  entries =  tsimCo60->GetEntries();

  for (int i = 0; i < entries; i++)
  {
    tsimCo60->GetEntry(i);
    double sigma = z_sim*sqrt( pow(coA,2) + pow(coB/sqrt(z_sim),2) + pow(coC/z_sim,2) );
    // sigma = 0.;
    hsimCo60->Fill(ranGen->Gaus(z_sim,sigma));
  }


  //! Canvas and draw
  // TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1000);
  // c1->cd();
  // c1->Divide(3, 2);

  // c1->cd(1);
  // hcalCs137->SetLineColor(kRed);
  // hcalCs137->Draw();
  // hsimCs137->Scale(0.68, "noSW2");
  // hsimCs137->Draw("same");

  // TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  // leg1->SetHeader("Cs137", "C");
  // leg1->SetBorderSize(2);
  // leg1->AddEntry(hsimCs137, "simulation", "l");
  // leg1->AddEntry(hcalCs137, "experiment", "l");
  // leg1->Draw();

  // c1->cd(2);
  // hsimNa22->Draw();
  // hcalNa22->SetLineColor(kRed);
  // // hcalNa22->Scale(0.175, "noSW2");
  // hcalNa22->Scale(0.175*2.8, "noSW2");
  // hcalNa22->Draw("same");

  // TLegend *leg2 = new TLegend(0.75, 0.6, 0.98, 0.75);
  // leg2->SetHeader("Na22", "C");
  // leg2->SetBorderSize(2);
  // leg2->AddEntry(hsimNa22, "simulation", "l");
  // leg2->AddEntry(hcalNa22, "experiment", "l");
  // leg2->Draw();

  // c1->cd(3);
  // hsimCo60->Draw();
  // hcalCo60->SetLineColor(kRed);
  // // hcalCo60->Scale(0.138, "noSW2");
  // hcalCo60->Scale(0.138/1.07, "noSW2");
  // hcalCo60->Draw("same");

  // TLegend *leg3 = new TLegend(0.75, 0.6, 0.98, 0.75);
  // leg3->SetHeader("Co60", "C");
  // leg3->SetBorderSize(2);
  // leg3->AddEntry(hsimCo60, "simulation", "l");
  // leg3->AddEntry(hcalCo60, "experiment", "l");
  // leg3->Draw();

  // c1->cd(4);
  // c1->cd(4)->SetLogy();
  // hcalCs137->SetLineColor(kRed);
  // hcalCs137->Draw();
  // hsimCs137->Draw("same");

  // c1->cd(5);
  // c1->cd(5)->SetLogy();
  // hcalNa22->SetLineColor(kRed);
  // hsimNa22->Draw();
  // hcalNa22->Draw("same");

  // c1->cd(6);
  // c1->cd(6)->SetLogy();
  // hcalCo60->SetLineColor(kRed);
  // hcalCo60->Draw();
  // hsimCo60->Draw("same");

  //! Differentiate
  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
  c2->cd();
  c2->Divide(3, 2);

  double binContent, binContentPlus, binLength, newBinContent;

  // hcalNa22->Rebin(4);
  // hcalCs137->Rebin(4);
  // binExp = binExp/4;

  hcalNa22->Smooth(2);
  hcalCs137->Smooth(2);
  hcalCo60->Smooth(2);

  binLength = (xmaxCal - xminCal)/ binExp;
  // std::cout << "binLength = " << binLength << "\n";

  TH1D* hDifcalNa22 = new TH1D("hDifcalNa22", "Na22 Exp Differentiated", binExp, xminCal, xmaxCal);
  TH1D* hDifcalCs137 = new TH1D("hDifcalCs137", "Cs137 Exp Differentiated", binExp, xminCal, xmaxCal);
  TH1D* hDifcalCo60 = new TH1D("hDifcalCo60", "Co60 Exp Differentiated", binExp, xminCal, xmaxCal);

  for(int i = 1; i < binExp; i++)
  {
    binContent = hcalNa22->GetBinContent(i);
    binContentPlus = hcalNa22->GetBinContent(i+1);
    newBinContent = (binContentPlus - binContent)/binLength;
    hDifcalNa22->SetBinContent(i, newBinContent);

    binContent = hcalCs137->GetBinContent(i);
    binContentPlus = hcalCs137->GetBinContent(i+1);
    newBinContent = (binContentPlus - binContent)/binLength;
    hDifcalCs137->SetBinContent(i, newBinContent);

    binContent = hcalCo60->GetBinContent(i);
    binContentPlus = hcalCo60->GetBinContent(i+1);
    newBinContent = (binContentPlus - binContent)/binLength;
    hDifcalCo60->SetBinContent(i, newBinContent);
  }

  // hDifcalCs137->Rebin(4);
  // hDifcalNa22->Rebin(4);

  c2->cd(1);
  hcalNa22->Draw();

  c2->cd(4);
  hDifcalNa22->SetLineColor(kRed);
  hDifcalNa22->Draw();

  c2->cd(2);
  hcalCs137->Draw();

  c2->cd(5);
  hDifcalCs137->SetLineColor(kRed);
  hDifcalCs137->Draw();

  c2->cd(3);
  hcalCo60->Draw();

  c2->cd(6);
  hDifcalCo60->SetLineColor(kRed);
  hDifcalCo60->Draw();
}
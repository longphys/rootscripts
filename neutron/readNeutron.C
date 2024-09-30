#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStopwatch.h"

//! histogram options
double xminSim = 0.;
double xmaxSim = 1.5;
double xminExp = 0;
double xmaxExp = 20000;
int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;

double ECs137 = 0.477;
double ENa22 = 1.061;
double ECo60 = 1.117;

//! ch0_9751_ch1_9800_ch2_9421
int ChCs137[2] = {388, 283};
int ChNa22[2] = {885, 642};

double deltaSigCs137[2] = {7.6, 7.6};
double deltaSigNa22[2] = {5.5, 5.5};

//! Channels
int channel[2] = {0, 1};

void readNeutron()
{
  auto timer = new TStopwatch();
  timer->Start();

  std::vector <double> *EDep = nullptr;
  std::vector <double> *Time = nullptr;

  //! Files and trees
  // TFile* f = new TFile("../../EfficiencyPlasticFLNR/declib/output/test08_2_Plast_BC404_14x_14y_Pos_x0_y0_out.root", "read");
  TFile* f = new TFile("../../EfficiencyPlasticFLNR/declib/output/test09_2_Plast_BC404_14x_14y_Pos_x20_y-20_out.root", "read");
  // TFile* f = new TFile("../../EfficiencyPlasticFLNR/declib/output/test10_2_Plast_BC404_14x_14y_Pos_x0_y-92_out.root", "read");
  // TFile* f = new TFile("../../EfficiencyPlasticFLNR/declib/output/test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1_out.root", "read");
  
  TTree* t =  (TTree*) f->Get("ETree");

  double a[2], b[2], xminCal[2], xmaxCal[2], coA[2], coB[2], coC[2];

  // TH1D* hcal = new TH1D("hcal", "Calibrated spectrum", binExp, xminExp, xmaxExp);
  // TH1D* hexp = new TH1D("hexp", "hexp", );
  TH1D* hcal[2];
  char* nameCal = new char[20];

  t->SetBranchAddress("EDep", &EDep);
  t->SetBranchAddress("Time", &Time);
  Long64_t entries = t->GetEntries();
  std::cout << "Number of entries: " << entries << "\n";

  TCanvas* cFit = new TCanvas("cFit", "Linear calibration", 800, 600);
  for(int i = 0; i < 2; i++){
    TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
    hCali->SetBinContent(ChCs137[i] - xminExp, ECs137);
    hCali->SetBinContent(ChNa22[i] - xminExp, ENa22);

    TF1* fCali = new TF1("fCali", "[0]*x + [1]", xminExp, xmaxExp);

    cFit->Divide(2,1);
    cFit->cd(1);
    std::cout << "\nENERGY CALIBRATION\n";
    hCali->Fit("fCali");

    a[i] = fCali->GetParameter(0);
    b[i] = fCali->GetParameter(1);

    xminCal[i] = xminExp*a[i]+b[i];
    xmaxCal[i] = xmaxExp*a[i]+b[i];

    sprintf(nameCal,"channel%d",i); 

    hcal[i] = new TH1D(nameCal, nameCal, binExp, xminCal[i], xmaxCal[i]);
    hcal[i]->GetXaxis()->SetTitle("Energy(MeV)");
    hcal[i]->GetYaxis()->SetTitle("Count");

    // for(int k = 0; k < entries; k++){
    for(int k = 0; k < 5000000; k++){
      if(k%1000000==0){std::cout << "entry: " << k << "\n";}
      t->GetEntry(k);
      int EDep_size = EDep->size();
      double* EDep_data = EDep->data();
      hcal[i]->Fill(a[i]*(EDep_data[channel[i]]+0.5) + b[i]);
    }

    TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal[i], xmaxCal[i]);
    hRes->SetBinContent(hcal[i]->FindBin(ECs137), deltaSigCs137[i]/100.);
    hRes->SetBinContent(hcal[i]->FindBin(ENa22), deltaSigNa22[i]/100.);

    TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");

    cFit->cd(2);
    hRes->Fit("fRes", "", "", ECs137-0.2, ECo60+0.2);

    coA[i] = fRes->GetParameter(0);
    coB[i] = 0.;
    coC[i] = fRes->GetParameter(1);

    delete hCali;
    delete fCali;
    delete hRes;
    delete fRes;
  }

  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
  c1->cd();
  c1->Divide(2, 1);

  hcal[0]->SetLineColor(kBlue);
  hcal[1]->SetLineColor(kRed);
  
  // hcal[0]->Scale(CalScaleCo60[0], "noSW2");
  // hcal[1]->Scale(CalScaleCo60[1], "noSW2");

  TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg1->SetHeader("test", "C");
  leg1->SetBorderSize(2);

  for(int i = 0; i<2; i++){
    c1->cd(i+1);
    hcal[i]->Draw();
    sprintf(nameCal,"Channel %d",i);
    leg1->AddEntry(hcal[i], nameCal, "l");
  }
  leg1->Draw();

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStopwatch.h"

//! histogram options
double xminSim = 0.;
double xmaxSim = 20.;
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
  TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test08_2_Plast_BC404_14x_14y_Pos_x0_y0_out.root", "read");
  // TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test09_2_Plast_BC404_14x_14y_Pos_x20_y-20_out.root", "read");
  // TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test10_2_Plast_BC404_14x_14y_Pos_x0_y-92_out.root", "read");
  // TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1_out.root", "read");
  TTree* texp =  (TTree*) fexp->Get("ETree");
  
  TFile* fsim = new TFile("../../simfiles/neutron/prem_088.root", "read");
  TTree* tsim =  (TTree*) fsim->Get("dEEtree");

  double a[2], b[2], xminCal[2], xmaxCal[2], coA[2], coB[2], coC[2];

  //! Measurement histogram
  TH1D* hcal[2];
  char* nameCal = new char[20];

  texp->SetBranchAddress("EDep", &EDep);
  texp->SetBranchAddress("Time", &Time);
  
  Long64_t entries_exp = texp->GetEntries();
  // std::cout << "Number of entries (exp): " << entries_exp << "\n";

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

    hcal[i] = new TH1D(nameCal, nameCal, binExp/10, xminCal[i], xmaxCal[i]);
    hcal[i]->GetXaxis()->SetTitle("Energy(MeV)");
    hcal[i]->GetYaxis()->SetTitle("Count");

    // for(int k = 0; k < entries_exp; k++){
    for(int k = 0; k < 5000000; k++){
      // if(k%1000000==0){std::cout << "entry: " << k << "\n";}
      texp->GetEntry(k);
      int EDep_size = EDep->size();
      double* EDep_data = EDep->data();
      hcal[i]->Fill(a[i]*(EDep_data[channel[i]]+0.5) + b[i]);
    }

    std::cout << "bin Cs137: " << hcal[i]->FindBin(ECs137) << "\n";
    std::cout << "bin Na22: " << hcal[i]->FindBin(ENa22) << "\n";

    TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal[i], xmaxCal[i]);
    hRes->SetBinContent(hcal[i]->FindBin(ECs137), deltaSigCs137[i]/100.);
    hRes->SetBinContent(hcal[i]->FindBin(ENa22), deltaSigNa22[i]/100.);

    TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");

    cFit->cd(2);
    hRes->Fit("fRes", "", "", ECs137-0.5, ECo60+0.5);

    cFit->Modified();
    cFit->Update();

    coA[i] = fRes->GetParameter(0);
    coB[i] = 0.;
    coC[i] = fRes->GetParameter(1);

    delete hCali;
    delete fCali;
    delete hRes;
    delete fRes;
  }

  //! Simulation histograms
  TH1D* hsim = new TH1D("hsim", "Simulation", binExp/10, xminCal[0], xmaxCal[0]);
  std::vector <double> *neutronE = nullptr;
  std::vector <double> *protonE = nullptr;
  std::vector <double> *gammaE = nullptr;
  std::vector <double> *alphaE = nullptr;
  std::vector <double> *c12E = nullptr;
  std::vector <double> *otherE = nullptr;

  tsim->SetBranchAddress("NeutronEDep", &neutronE);
  tsim->SetBranchAddress("ProtonEDep", &protonE);
  tsim->SetBranchAddress("GammaEDep", &gammaE);
  tsim->SetBranchAddress("AlphaEDep", &alphaE);
  tsim->SetBranchAddress("C12EDep", &c12E);
  tsim->SetBranchAddress("OtherEDep", &otherE);

  TRandom3* ranGen = new TRandom3();

  coA[0] = 0.05;
  coC[0] = 0.002;

  std::cout << "coA = " << coA[0] << "; coB = " << coB[0] << "; coC = " << coC[0] << "\n";


  unsigned int entries_sim = tsim->GetEntriesFast();
  for (int i = 0; i < entries_sim; i++){
    tsim->GetEntry(i);

    int neutronE_size = neutronE->size();
    int protonE_size = protonE->size();
    int gammaE_size = gammaE->size();
    int alphaE_size = alphaE->size();
    int c12E_size = c12E->size();
    int otherE_size = otherE->size();

    double* neutronE_data = neutronE->data();
    double* protonE_data = protonE->data();
    double* gammaE_data = gammaE->data();
    double* alphaE_data = alphaE->data();
    double* c12E_data = c12E->data();
    double* otherE_data = otherE->data();

    double sigma_neutron = neutronE_data[0]*sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(neutronE_data[0]),2) + pow(coC[0]/neutronE_data[0],2) );
    double sigma_proton = protonE_data[0]*sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(protonE_data[0]),2) + pow(coC[0]/protonE_data[0],2) );
    double sigma_gamma = gammaE_data[0]*sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(gammaE_data[0]),2) + pow(coC[0]/gammaE_data[0],2) );
    double sigma_alpha = alphaE_data[0]*sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(alphaE_data[0]),2) + pow(coC[0]/alphaE_data[0],2) );
    double sigma_c12 = c12E_data[0]*sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(c12E_data[0]),2) + pow(coC[0]/c12E_data[0],2) );
    double sigma_other = otherE_data[0]*sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(otherE_data[0]),2) + pow(coC[0]/otherE_data[0],2) );

    double energy_neutron = ranGen->Gaus(neutronE_data[0], sigma_neutron);
    double energy_proton = ranGen->Gaus(protonE_data[0], sigma_proton);
    double energy_gamma = ranGen->Gaus(gammaE_data[0], sigma_gamma);
    double energy_alpha = ranGen->Gaus(alphaE_data[0], sigma_alpha);
    double energy_c12 = ranGen->Gaus(c12E_data[0], sigma_c12);
    double energy_other = ranGen->Gaus(otherE_data[0], sigma_other);

    double energy_abs = neutronE_data[0]+protonE_data[0]+gammaE_data[0]+alphaE_data[0]+c12E_data[0]+otherE_data[0];
    double sigma = energy_abs*sqrt( pow(coA[0],2) + pow(coB[0]/sqrt(energy_abs),2) + pow(coC[0]/energy_abs,2) );
    double energy_res = ranGen->Gaus(energy_abs, sigma);

    // if(i%100000 == 0){
    //   std::cout << "\nEntry: " << i << "\n";
    //   std::cout << "energy_abs = " << energy_abs << "\n";
    //   std::cout << "sigma = " << sigma << "\n";
    //   std::cout << "energy_res = " << energy_res << "\n";
    // }
    // hsim->Fill(energy_abs);
    hsim->Fill(energy_res);
  }
  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
  c1->cd();
  c1->Divide(2, 1);

  hcal[0]->SetLineColor(kRed);
  hcal[1]->SetLineColor(kBlue);
  
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

  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
  c2->cd();
  // c2->Divide(2,1);

  TLegend *leg2 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg2->SetHeader("test", "C");
  leg2->SetBorderSize(2);

  // c2->cd(1);
  hcal[0]->Draw();
  // c2->cd(2);
  hsim->SetLineColor(kBlue);
  hsim->Draw("same");
  leg2->AddEntry(hcal[0], "measurement", "l");
  leg2->AddEntry(hsim, "simulation", "l");

  leg2->Draw();

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
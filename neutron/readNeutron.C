#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"
#include "TCutG.h"

//! histogram options
double xminSim = 0.;
double xmaxSim = 20.;
double xminExp = 0.;
double xmaxExp = 20000;

int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;
int binExpFactor = 10;

double ECs137 = 0.477;
double ENa22 = 1.061;
double ECo60 = 1.117;

int ChCs137[2] = {365, 283};
int ChNa22[2] = {860, 642};
// int ChCs137[2] = {385, 283};
// int ChNa22[2] = {862, 642};

double deltaSigCs137[2] = {7.5, 7.6};
double deltaSigNa22[2] = {5.1, 5.5};

//! Channels
int channel[2] = {0, 1};

void readNeutron()
{
  auto timer = new TStopwatch();
  timer->Start();

  std::vector <double> *EDep = nullptr;
  std::vector <double> *Time = nullptr;
  std::vector <double> *EDepBackground = nullptr;

  //! Files and trees
  TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test08_2_Plast_BC404_14x_14y_Pos_x0_y0_out.root", "read");
  // TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test09_2_Plast_BC404_14x_14y_Pos_x20_y-20_out.root", "read");
  // TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test10_2_Plast_BC404_14x_14y_Pos_x0_y-92_out.root", "read");
  // TFile* fexp = new TFile("../../EfficiencyPlasticFLNR/declib/output/test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1_out.root", "read");

  TFile* fbackground = new TFile("../../EfficiencyPlasticFLNR/declib/output/cal_test04_2_Plast_BC404_Bskgr.root", "read");
  // TFile* fbackground = new TFile("../../EfficiencyPlasticFLNR/declib/output/test14_2_Plast_BC404_Background.root", "read");

  TTree* texp =  (TTree*) fexp->Get("ETree");
  TTree* tbackground =  (TTree*) fbackground->Get("ETree");
  
  // TFile* fsim = new TFile("../../simfiles/neutron/eff_1.45_birk_207.root", "read");
  // TFile* fsim = new TFile("../../simfiles/neutron/eff_0.95_birk_088.root", "read");
  // TFile* fsim = new TFile("../../simfiles/neutron/eff_0.8_birk_0486.root", "read");
  // TFile* fsim = new TFile("../../simfiles/neutron/conference_file.root", "read");
  TFile* fsim = new TFile("../../simfiles/neutron/out.root", "read");
  // TFile* fsim = new TFile("../../simfiles/neutron/test.root", "read");
  // TFile* fsim = new TFile("../../simfiles/neutron/test1.root", "read");
  // TFile* fsim = new TFile("../../../software/geant4/geant4build/sim/build/out.root", "read");
  TTree* tsim =  (TTree*) fsim->Get("dEEtree");

  double a[2], b[2], xminCal[2], xmaxCal[2], coA[2], coB[2], coC[2];

  //! Measurement histogram
  TH1D* hcal_original[2];
  TH1D* hcal[2];
  TH1D* h_sil_time[2];
  TH1D* h_plas_time[2];
  TH2D* htime_diff[2];
  char* nameCal = new char[20];

  texp->SetBranchAddress("EDep", &EDep);
  texp->SetBranchAddress("Time", &Time);

  TH1D* hbackground[2];
  char* nameCalBackground = new char[20];
  tbackground->SetBranchAddress("EDep", &EDepBackground);
  
  Long64_t entries_exp = texp->GetEntries();
  entries_exp = 15460000; //! 10% of data
  // entries_exp = 15460000*3; //!
  // entries_exp = 100000;
  // entries_exp = 5000000;


  // double scale_measurement = 7.0;
  // double scale_measurement = 7.0/3.;
  double scale_measurement = 1.0;
  // double scale_sim = 1.8;
  // double scale_sim = 1.0;
  double scale_sim = 4.0;

  // std::cout << "Number of entries (exp): " << entries_exp << "\n";

  TCutG *cutg = new TCutG("mycut", 0);
  cutg->SetVarX("y");
  cutg->SetVarY("x");

  //! positive peak
  // cutg->SetPoint(0, 0.0, 0.0);
  // cutg->SetPoint(1, 0.5, 0.0);
  // cutg->SetPoint(2, 1.0, 0.14);
  // cutg->SetPoint(3, 1.5, 0.15);
  // cutg->SetPoint(4, 2.0, 0.16);
  // cutg->SetPoint(5, 2.5, 0.165);
  // cutg->SetPoint(6, 3.0, 0.17);
  // cutg->SetPoint(7, 4.0, 0.175);
  // cutg->SetPoint(8, 5.0, 0.18);
  // cutg->SetPoint(9, 6.0, 0.185);
  // cutg->SetPoint(10, 7.0, 0.1875);
  // cutg->SetPoint(11, 8.0, 0.19);
  // cutg->SetPoint(12, 20.0, 0.2);
  // cutg->SetPoint(13, 20.0, 0.2);
  // cutg->SetPoint(14, 8.0, 2.10);
  // cutg->SetPoint(15, 7.0, 2.1025);
  // cutg->SetPoint(16, 6.0, 2.105);
  // cutg->SetPoint(17, 5.0, 2.11);
  // cutg->SetPoint(18, 4.0, 2.115);
  // cutg->SetPoint(19, 3.0, 2.12);
  // cutg->SetPoint(20, 2.5, 2.125);
  // cutg->SetPoint(21, 2.0, 2.13);
  // cutg->SetPoint(22, 1.5, 2.14);
  // cutg->SetPoint(23, 1.0, 2.25);
  // cutg->SetPoint(24, 0.5, 2.35);
  // cutg->SetPoint(25, 0.0, 2.6);
  // cutg->SetPoint(26, 0.0, 0.0);

  //! negative peak
  // cutg->SetPoint(0, 0.0, 0.0);
  // cutg->SetPoint(1, 10.0, 0.0);
  // cutg->SetPoint(2, 10.0, -3.0);
  // cutg->SetPoint(3, 0.0, -3.0);
  // cutg->SetPoint(4, 0.0, 0.0);

  //! 2 peaks
  cutg->SetPoint(0, 0.0, 3.0);
  cutg->SetPoint(1, 10.0, 3.0);
  cutg->SetPoint(2, 10.0, -3.0);
  cutg->SetPoint(3, 0.0, -3.0);
  cutg->SetPoint(4, 0.0, 0.0);

  TCanvas* cFit = new TCanvas("cFit", "Linear calibration", 800, 600);
  for(int i = 0; i < 1; i++){
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

    hcal_original[i] = new TH1D(nameCal, nameCal, binExp/binExpFactor, xminCal[i], xmaxCal[i]);
    hcal_original[i]->GetXaxis()->SetTitle("Energy(MeV)");
    hcal_original[i]->GetYaxis()->SetTitle("Count");

    hcal[i] = new TH1D(nameCal, nameCal, binExp/binExpFactor, xminCal[i], xmaxCal[i]);
    hcal[i]->GetXaxis()->SetTitle("Energy(MeV)");
    hcal[i]->GetYaxis()->SetTitle("Count");

    h_plas_time[i] = new TH1D("h_plas_time", "h_plas_time", 10000, -100., 100.);
    h_sil_time[i] = new TH1D("h_sil_time", "h_sil_time", 10000, -100., 100.);
    htime_diff[i] = new TH2D("htime_diff", "htime_diff", binExp/binExpFactor/10, xminCal[i], xmaxCal[i], 2000, -20, 20);

    int beam_xc = 8;
    int beam_yc = 8;
    int time_diff_strip = 8;

    double threshold_plas_0 = 87.;
    double threshold_plas_1 = 81.;

    for(int k = 0; k < entries_exp; k++){
    // for(int k = 0; k < 10000000; k++){
      if(k%100000==0){std::cout << "entry: " << k << " of " << entries_exp << "\n";}
      texp->GetEntry(k);
      int EDep_size = EDep->size();
      double* EDep_data = EDep->data();
      int Time_size = Time->size();
      double* Time_data = Time->data();

      double energy = a[i]*(EDep_data[channel[i]]+0.05) + b[i];
      hcal_original[i]->Fill(energy);

      std::vector <int> xy_allalphas;
      std::vector <int> xy_plst0;
      std::vector <int> xy_plst1;
      for(int j = 2; j < EDep_size; j++){
        // if (EDep_data[j] > 1 ) //! Signal on silicon
        if (EDep_data[j] > 0 ) //! Signal on silicon
        {
          // xy_allalphas.push_back(j);
          //! Signal on silicon and Plastic0
          if(EDep_data[0] > threshold_plas_0){
            xy_plst0.push_back(j);
          }
        }
        // if(j == time_diff_strip){
        //   double timediff = Time_data[i]-Time_data[j];
        //   htime_diff[i]->Fill(energy, timediff);
        // }
      }

      if (xy_plst0.size() == 2) //! SIGNAL IN SILICONS
      // if (xy_plst0.size() > 1) //!
      {
        int x = xy_plst0.at(0);
        int y = xy_plst0.at(1);

        // double timediff = Time_data[i]-Time_data[time_diff_strip];
        // htime_diff[i]->Fill(energy, timediff);

        if (x < 16 && x > 1 && y > 15 && y < 30)
        {
          x -= 1;
          y -= 15;

          // if ((x) && (y))
          if ((x == beam_xc) && (y == beam_yc))
          // if ((x >= beam_xc-2) && (x <= beam_xc+2) && (y >= beam_yc-2) && (y <= beam_yc+2))
          {
            double timediff = Time_data[i]-Time_data[time_diff_strip];
            h_sil_time[i]->Fill(Time_data[time_diff_strip]);
            h_plas_time[i]->Fill(Time_data[i]);
            htime_diff[i]->Fill(energy, timediff);
            if(Time_data[i] != 0){
              // double timediff = Time_data[i]-Time_data[time_diff_strip];
              // htime_diff[i]->Fill(energy, timediff);
              if (cutg->IsInside(energy, timediff)){ //! if (x = energy, y = timediff is inside cutg) return 1
                hcal[i]->Fill(energy);
              }
            }
          }
        }
      }
    }

    // sprintf(nameCalBackground,"backgroundChannel%d",i); 
    // hbackground[i] = new TH1D(nameCalBackground, nameCalBackground, binExp/binExpFactor, xminCal[i], xmaxCal[i]);
    // hbackground[i]->GetXaxis()->SetTitle("Energy(MeV)");
    // hbackground[i]->GetYaxis()->SetTitle("Count");
    // for(int k = 0; k < entries_exp; k++){
    // // for(int k = 0; k < 10000000; k++){
    //   tbackground->GetEntry(k);
    //   int EDepBackground_size = EDepBackground->size();
    //   double* EDepBackground_data = EDepBackground->data();
    //   hbackground[i]->Fill(a[i]*(EDepBackground_data[channel[i]]+0.5) + b[i]);
    // }

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
  TH1D* hsim[2];
  hsim[0] = new TH1D("hsim", "Simulation", binExp/binExpFactor, xminCal[0], xmaxCal[0]);
  TH1D* hcompare = new TH1D("hcompare", "Simulation", binExp/binExpFactor, xminCal[0], xmaxCal[0]);
  TH1D* htotal = new TH1D("htotal", "Simulation", binExp/binExpFactor, xminCal[0], xmaxCal[0]);

  std::vector <double> *neutronE = nullptr;
  std::vector <double> *protonE = nullptr;
  std::vector <double> *protonERaw = nullptr;
  std::vector <double> *gammaE = nullptr;
  std::vector <double> *alphaE = nullptr;
  std::vector <double> *c12E = nullptr;
  std::vector <double> *otherE = nullptr;

  tsim->SetBranchAddress("NeutronEDep", &neutronE);
  tsim->SetBranchAddress("ProtonEDep", &protonE);
  tsim->SetBranchAddress("ProtonEDepRaw", &protonERaw);
  tsim->SetBranchAddress("GammaEDep", &gammaE);
  tsim->SetBranchAddress("AlphaEDep", &alphaE);
  tsim->SetBranchAddress("C12EDep", &c12E);
  tsim->SetBranchAddress("OtherEDep", &otherE);

  TRandom3* ranGen = new TRandom3();

  int channel_sim = 0;

  // coA[channel_sim] = 0.03;
  // coB[channel_sim] = 0.1;
  // coC[channel_sim] = 0.002;

  // coA[channel_sim] = 0.0001;
  // coB[channel_sim] = 0.04;
  // coC[channel_sim] = 0.004;

  std::cout << "coA = " << coA[channel_sim] << "; coB = " << coB[channel_sim] << "; coC = " << coC[channel_sim] << "\n";

  unsigned int entries_sim = tsim->GetEntriesFast();
  // entries_sim = entries_sim/2;
  for (int i = 0; i < entries_sim; i++){
    tsim->GetEntry(i);

    int neutronE_size = neutronE->size();
    int protonE_size = protonE->size();
    int protonERaw_size = protonERaw->size();
    int gammaE_size = gammaE->size();
    int alphaE_size = alphaE->size();
    int c12E_size = c12E->size();
    int otherE_size = otherE->size();

    double* neutronE_data = neutronE->data();
    double* protonE_data = protonE->data();
    double* protonERaw_data = protonERaw->data();
    double* gammaE_data = gammaE->data();
    double* alphaE_data = alphaE->data();
    double* c12E_data = c12E->data();
    double* otherE_data = otherE->data();

    // double neutron_sigma = neutronE_data[channel_sim]*sqrt( pow(coA[channel_sim],2) + 
    // pow(coB[channel_sim]/sqrt(neutronE_data[channel_sim]),2) + pow(coC[channel_sim]/neutronE_data[channel_sim],2) );
    // double proton_sigma = protonE_data[channel_sim]*sqrt( pow(coA[channel_sim],2) + 
    // pow(coB[channel_sim]/sqrt(protonE_data[channel_sim]),2) + pow(coC[channel_sim]/protonE_data[channel_sim],2) );
    // double gamma_sigma = gammaE_data[channel_sim]*sqrt( pow(coA[channel_sim],2) + 
    // pow(coB[channel_sim]/sqrt(gammaE_data[channel_sim]),2) + pow(coC[channel_sim]/gammaE_data[channel_sim],2) );
    // double alpha_sigma = alphaE_data[channel_sim]*sqrt( pow(coA[channel_sim],2) + 
    // pow(coB[channel_sim]/sqrt(alphaE_data[channel_sim]),2) + pow(coC[channel_sim]/alphaE_data[channel_sim],2) );
    // double c12_sigma = c12E_data[channel_sim]*sqrt( pow(coA[channel_sim],2) + 
    // pow(coB[channel_sim]/sqrt(c12E_data[channel_sim]),2) + pow(coC[channel_sim]/c12E_data[channel_sim],2) );
    // double sigma_other = otherE_data[channel_sim]*sqrt( pow(coA[channel_sim],2) + 
    // pow(coB[channel_sim]/sqrt(otherE_data[channel_sim]),2) + pow(coC[channel_sim]/otherE_data[channel_sim],2) );

    // double neutron_energy = ranGen->Gaus(neutronE_data[channel_sim], neutron_sigma);
    // double proton_energy = ranGen->Gaus(protonE_data[channel_sim], proton_sigma);
    // double gamma_energy = ranGen->Gaus(gammaE_data[channel_sim], gamma_sigma);
    // double alpha_energy = ranGen->Gaus(alphaE_data[channel_sim], alpha_sigma);
    // double c12_energy = ranGen->Gaus(c12E_data[channel_sim], c12_sigma);
    // double other_energy = ranGen->Gaus(otherE_data[channel_sim], sigma_other);

    // double energy_abs = protonERaw_data[channel_sim];
    double energy_abs = protonE_data[channel_sim]+alphaE_data[channel_sim]+c12E_data[channel_sim]+otherE_data[channel_sim];
    // double energy_abs = protonE_data[channel_sim]+alphaE_data[channel_sim]+c12E_data[channel_sim];
    // double energy_abs = protonE_data[channel_sim]+alphaE_data[channel_sim];
    if (energy_abs !=0.){
      double sigma = energy_abs*sqrt( pow(coA[channel_sim],2) + pow(coB[channel_sim]/sqrt(energy_abs),2) + 
      pow(coC[channel_sim]/energy_abs,2) );
      double energy_res = ranGen->Gaus(energy_abs, sigma);

      double energy_abs2 = protonE_data[channel_sim];
      double sigma2 = energy_abs2*sqrt( pow(coA[channel_sim],2) + pow(coB[channel_sim]/sqrt(energy_abs2),2) + 
      pow(coC[channel_sim]/energy_abs2,2) );
      double energy_res2 = ranGen->Gaus(energy_abs2, sigma2);

      double energy_abs3 = alphaE_data[channel_sim];
      double sigma3 = energy_abs3*sqrt( pow(coA[channel_sim],2) + pow(coB[channel_sim]/sqrt(energy_abs3),2) + 
      pow(coC[channel_sim]/energy_abs3,2) );
      double energy_res3 = ranGen->Gaus(energy_abs3, sigma3);
      // if(i == entries_sim){
      //   // std::cout << "\nEntry: " << i << "\n";
      //   std::cout << "energy_abs = " << energy_abs << "\n";
      //   std::cout << "sigma = " << sigma << "\n";
      //   std::cout << "energy_res = " << energy_res << "\n";
      // }
      hsim[channel_sim]->Fill(energy_res2);
      hcompare->Fill(energy_res3);
      htotal->Fill(energy_res);
    }

  }

  //! Canvas and draw
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
  c1->cd();
  c1->Divide(2, 1);

  // hcal[0]->SetLineColor(kRed);
  hcal[0]->SetLineColor(kBlue+2);
  hcal[0]->SetLineWidth(2);
  // hcal[1]->SetLineColor(kRed);
  
  TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg1->SetHeader("test", "C");
  leg1->SetBorderSize(2);

  c1->cd(1);
  hcal[0]->Draw();
  leg1->AddEntry(hcal[0], "plastic 0", "l");

  c1->cd(2);
  // hcal[1]->Draw();
  // leg1->AddEntry(hcal[1], "plastic 1", "l");

  leg1->Draw();

  // c1->cd(3);
  // c1->cd(3)->SetLogy();
  // hbackground[1]->Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
  c2->cd();
  c2->SetLogy();
  c2->SetTicks();
  c2->SetLeftMargin(0.15);
  // c2->Divide(2,1);
  // c2->SetLogy();

  TLegend *leg2 = new TLegend(0.4, 0.6, 0.85, 0.85);
  // leg2->SetHeader("test", "C");
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.06);

  hcal[channel_sim]->SetTitle("");
  hcal[channel_sim]->GetXaxis()->SetTitle("E(MeVee)");
  hcal[channel_sim]->GetXaxis()->SetLabelFont(42);
  hcal[channel_sim]->GetXaxis()->SetTitleFont(52);
  hcal[channel_sim]->GetXaxis()->SetTitleSize(0.04);
  hcal[channel_sim]->GetXaxis()->CenterTitle(true);

  hcal[channel_sim]->GetYaxis()->SetTitle("Events");
  hcal[channel_sim]->GetYaxis()->SetLabelFont(42);
  hcal[channel_sim]->GetYaxis()->SetTitleFont(52);
  hcal[channel_sim]->GetYaxis()->SetTitleSize(0.04);
  hcal[channel_sim]->GetYaxis()->CenterTitle(true);
  
  // c2->cd(1);
  // hcal[channel_sim]->Add(hbackground[channel_sim],-1.);
  hcal[channel_sim]->Scale(scale_measurement, "noSW2");
  hcal[channel_sim]->SetLineColor(kBlue+2);
  hcal[channel_sim]->SetLineWidth(2);
  hcal[channel_sim]->SetStats(0);
  hcal[channel_sim]->Draw();
  // c2->cd(2);
  // hsim[channel_sim]->Scale(1.6, "noSW2");
  hsim[channel_sim]->Scale(scale_sim, "noSW2");
  hsim[channel_sim]->SetLineColor(kBlack);
  hsim[channel_sim]->SetLineWidth(2);
  // hsim[channel_sim]->Draw("same");
  htotal->Scale(scale_sim, "noSW2");
  htotal->SetLineColor(kBlue);
  htotal->SetLineWidth(2);
  htotal->SetStats(0);
  htotal->Draw("same");

  // hcal[channel_sim]->Add(hbackground[channel_sim],-1.);
  hcal_original[channel_sim]->Scale(scale_measurement, "noSW2");
  hcal_original[channel_sim]->SetLineColor(kGreen);
  hcal_original[channel_sim]->SetLineWidth(2);
  // hcal_original[channel_sim]->Draw("same");
  // leg2->AddEntry(hcal_original[channel_sim], "measurement", "l");
  leg2->AddEntry(hcal[channel_sim], "Measurement", "l");
  // leg2->AddEntry(hsim[channel_sim], "simulation only proton", "l");
  // leg2->AddEntry(htotal, "simulation w res", "l");
  leg2->AddEntry(htotal, "Simulation", "l");

  leg2->Draw();

  //! Chi2
  double xmin_chi2 = 1.0;
  double xmax_chi2 = 10.0;
  hcal_original[channel_sim]->GetXaxis()->SetRangeUser(xmin_chi2, xmax_chi2);
  hcal[channel_sim]->GetXaxis()->SetRangeUser(xmin_chi2, xmax_chi2);
  hsim[channel_sim]->GetXaxis()->SetRangeUser(xmin_chi2, xmax_chi2);
  hcal[channel_sim]->GetYaxis()->SetRangeUser(0., 2000.0);
  hsim[channel_sim]->GetYaxis()->SetRangeUser(0., 2000.0);
  double chi2 = hcal[channel_sim]->Chi2Test(hsim[channel_sim], "CHI2/NDF");
  std::cout << "Chi2_proton from " << xmin_chi2 << " to " << xmax_chi2 << " (MeV) = " << chi2 << "\n";

  htotal->GetXaxis()->SetRangeUser(xmin_chi2, xmax_chi2);
  htotal->GetYaxis()->SetRangeUser(0., 2000.0);

  double chi2_total = hcal[channel_sim]->Chi2Test(htotal, "CHI2/NDF");
  std::cout << "Chi2_total from " << xmin_chi2 << " to " << xmax_chi2 << " (MeV) = " << chi2_total << "\n";

  // hcal[channel_sim]->GetXaxis()->UnZoom();
  // hsim[channel_sim]->GetXaxis()->UnZoom();
  // hcal[channel_sim]->GetYaxis()->UnZoom();
  // hsim[channel_sim]->GetYaxis()->UnZoom();

  //! time_data
  TCanvas* ctest = new TCanvas("ctest", "ctest", 800, 600);
  ctest->Divide(2,1);
  ctest->cd(1);
  h_plas_time[0]->Draw();

  ctest->cd(2);
  h_sil_time[0]->Draw();


  double energy0 = 0.15;
  // double energy1 = 0.15;
  double threshold0 = (energy0-b[0])/a[0];
  // double threshold1 = (energy1-b[1])/a[1];

  std::cout << "threshold[0] = " << threshold0 << "\n";
  // std::cout << "threshold[1] = " << threshold1 << "\n";

  TCanvas *c3 = new TCanvas("c3", "c3", 1800, 800);
  c3->Divide(2,1);
  c3->cd(1);
  c3->cd(1)->SetLeftMargin(0.15);
  htime_diff[0]->SetTitle("Time difference over Energy");
  htime_diff[0]->SetStats(0);

  htime_diff[0]->GetXaxis()->SetTitle("E(MeVee)");
  htime_diff[0]->GetXaxis()->SetLabelFont(42);
  htime_diff[0]->GetXaxis()->SetTitleFont(52);
  htime_diff[0]->GetXaxis()->SetTitleSize(0.04);
  htime_diff[0]->GetXaxis()->CenterTitle(true);

  htime_diff[0]->GetYaxis()->SetTitle("Time Difference");
  htime_diff[0]->GetYaxis()->SetLabelFont(42);
  htime_diff[0]->GetYaxis()->SetTitleFont(52);
  htime_diff[0]->GetYaxis()->SetTitleSize(0.04);
  htime_diff[0]->GetYaxis()->CenterTitle(true);
  htime_diff[0]->Draw("col z");
  c3->cd(2);
  c3->cd(2)->SetLeftMargin(0.15);
  TH2D* htime_diff_clone = (TH2D*)htime_diff[0]->Clone();
  htime_diff_clone->SetTitle("Time difference cut-off");
  htime_diff_clone->SetStats(0);

  htime_diff_clone->GetXaxis()->SetTitle("E(MeV)");
  htime_diff_clone->GetXaxis()->SetLabelFont(42);
  htime_diff_clone->GetXaxis()->SetTitleFont(52);
  htime_diff_clone->GetXaxis()->SetTitleSize(0.04);
  htime_diff_clone->GetXaxis()->CenterTitle(true);

  htime_diff_clone->GetYaxis()->SetTitle("Time Difference");
  htime_diff_clone->GetYaxis()->SetLabelFont(42);
  htime_diff_clone->GetYaxis()->SetTitleFont(52);
  htime_diff_clone->GetYaxis()->SetTitleSize(0.04);
  htime_diff_clone->GetYaxis()->CenterTitle(true);
  htime_diff_clone->Draw("col z, [mycut]");

  TCanvas *c4 = new TCanvas("c4", "c4", 1600, 800);
  // c4->Divide(2,1);
  // c4->cd(1);
  // htime_x[0]->Draw();
  // c4->cd(2);
  // c4->cd(2)->SetTitle("htime with cut");
  // htime_y[0]->Draw();

  TH1D* hcut = (TH1D*)hcal[channel_sim]->Clone();
  hcut->Add(hsim[channel_sim], -1.);
  // hcut->Add(hcal[channel_sim], -1.);
  // double diff = hcut->Integral(0,binExp/binExpFactor);
  // std::cout << "diff = " << diff << "\n";
  for (int i = 0; i < hcut->GetNbinsX(); i++){
    if (hcut->GetBinContent(i) < 0){
      hcut->SetBinContent(i, 0);
    }
  }
  c4->cd();
  xmin_chi2 = 0.2;
  xmax_chi2 = 2.0;
  hcompare->Scale(scale_sim, "noSW2");
  hcompare->GetXaxis()->SetRangeUser(xmin_chi2, xmax_chi2);
  hcompare->SetLineColor(kBlack);
  hcompare->Draw();
  hcut->GetXaxis()->SetRangeUser(xmin_chi2, xmax_chi2);
  hcut->Draw("same");
  double chi2_cut = hcut->Chi2Test(hcompare, "CHI2/NDF");
  std::cout << "Chi2_cut from " << xmin_chi2 << " to " << xmax_chi2 << " (MeV) = " << chi2_cut << "\n";

  TLegend *leg3 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg3->SetHeader("test", "C");
  leg3->SetBorderSize(2);
  leg3->AddEntry(hcut, "measurement residue", "l");
  leg3->AddEntry(hcompare, "alpha simulation", "l");

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";

  TCanvas* c5 = new TCanvas("c5", "c5", 1600, 800);
  c5->Divide(2,1);
  c5->cd(1);
  htotal->GetXaxis()->UnZoom();
  htotal->GetYaxis()->UnZoom();
  htotal->Draw();
  c5->cd(2);
  hcal[channel_sim]->GetXaxis()->UnZoom();
  hcal[channel_sim]->GetYaxis()->UnZoom();
  hcal[channel_sim]->Draw();

  double integral_blue = htotal->Integral(htotal->FindBin(0.2), htotal->FindBin(1.2));
  double integral_red = hcal[channel_sim]->Integral(hcal[channel_sim]->FindBin(0.2), hcal[channel_sim]->FindBin(1.2));
  std::cout << "ratio = " << abs(integral_red-integral_blue)/hcal[channel_sim]->Integral(0, hcal[channel_sim]->GetNbinsX()) << "\n";

  TFile* out = new TFile("out.root", "recreate");
  TH1D* hClone = (TH1D*)hcal[0]->Clone();
  out->Write();
  out->Close();

  TCanvas* exportCanvas = new TCanvas();
  exportCanvas->cd();
  hcal[0]->Draw();
}

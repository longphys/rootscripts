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

double E1 = 0.341;
double E2 = 1.061;

//! Zoomed in Histograms options
double minTest1274 = 0.2; //MeV
double maxTest1274 = 1.4;

//! ENERGY CALIBRATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double Ch1 = 315;
double Ch2 = 800;

//! RESOLUTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double deltaSig1 = 0.08;
double deltaSig2 = 0.04;

//! SCALING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ScaleCs137 = 4.7;
double ScaleNa22 = 3.1;

//! CHANNEL
int channel = 2;

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

void compare4Par_0()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Files and trees
  TFile* fsimNa22 = new TFile("./simfiles/withAl_simNa22_2.root", "read");
  TFile* fexpNa22 = new TFile("./expfiles/new/stilbene_dividers_Plastic_ch2_Na22_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");

  TTree* tsimNa22 =  (TTree*) fsimNa22->Get("dEEtree");
  TTree* texpNa22 =  (TTree*) fexpNa22->Get("AnalysisxTree");

  UShort_t chY[48];
  texpNa22->SetBranchAddress("NeEvent.neutAmp[48]", chY);
  // Long64_t entriesExpNa22 = texpNa22->GetEntries();
  Long64_t entriesExpNa22 = 1000000;

  double y_sim;
  tsimNa22->SetBranchAddress("Scintillator", &y_sim);
  // Long64_t entriesSimNa22 = tsimNa22->GetEntries();
  Long64_t entriesSimNa22 = 1000000;

  TH1D* hCali = new TH1D("hCali", "Calibration fit", binExp, xminExp, xmaxExp);
  hCali->SetBinContent(Ch1 - xminExp, E1);
  hCali->SetBinContent(Ch2 - xminExp, E2);

  TF1* fLinear = new TF1("fLinear", "[0]*x + [1]", xminExp, xmaxExp);
  hCali->Fit("fLinear");

  double a = fLinear->GetParameter(0);
  double b = fLinear->GetParameter(1);

  double aStep = 0.0001*a;
  double bStep = 0.0001*b;

  //! Resolution for simulation histograms

  std::cout << "delta1 = " << deltaSig1*100 << "%; delta2 = " << deltaSig2*100 << "%\n";

  double xminCal = xminExp*a+b;
  double xmaxCal = xmaxExp*a+b;

  TH1D* hFirstCalNa22 = new TH1D("hFirstCalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);
  for(int i = 0; i < entriesExpNa22; i++){
    texpNa22->GetEntry(i);
    hFirstCalNa22->Fill(a*(chY[channel]+0.5) + b);
  }

  TH1D* hRes = new TH1D("hRes", "Resolution fit", binExp, xminCal, xmaxCal);
  hRes->SetBinContent(hFirstCalNa22->FindBin(E1), deltaSig1);
  hRes->SetBinContent(hFirstCalNa22->FindBin(E2), deltaSig2);

  //! Choose fit function
  TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");

  hRes->Fit("fRes", "", "", E1-0.2, E2+0.2);

  //! Correspond to chosen function
  double coA = std::abs(fRes->GetParameter(0));
  double coB = 0.;
  double coC = std::abs(fRes->GetParameter(1));

  double coAStep = coA*0.0001;
  // double coBStep = coB*0.2;
  double coCStep = coC*0.0001;

  int tuningTimes = 200;
  double tuningRateA = 0.000000000001;
  double tuningRateB = 0.00000002;
  double tuningRateCoA = 0.00000000001;
  double tuningRateCoC = 0.00000000001;
  
  TRandom3* ranGen = new TRandom3();

  TCanvas* cShow = new TCanvas("cShow", "cShow", 1500, 600);
  cShow->Divide(2,1);

  TCanvas* cShow1 = new TCanvas("cShow1", "cShow1", 1500, 600);
  cShow1->Divide(2,2);

  double chi2_1274;
  TGraph* grChi2 = new TGraph();

  double deltaChi2 = 1000000.;
  double deltaChi2_1 = 0., deltaChi2_2 = 0., deltaChi2_3 = 0.;

  // for (int i = 1; i<=tuningTimes; i++){
  int count = 1;
  while (deltaChi2>50.){
    std::cout << "\nIteration " << count << ":\n";

    xminCal = xminExp*a+b;
    xmaxCal = xmaxExp*a+b;

    double xminCalUpA = xminExp*(a + aStep)+b;
    double xmaxCalUpA = xmaxExp*(a + aStep)+b;

    double xminCalUpB = xminExp*a+(b + bStep);
    double xmaxCalUpB = xmaxExp*a+(b + bStep);

    //! Histograms 
    TH1D* hsimNa22res = new TH1D("hsimNa22res", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);
    TH1D* hcalNa22 = new TH1D("hcalNa22", "Na22 Experiment Calibrated", binExp, xminCal, xmaxCal);
    
    TH1D* hsimNa22resUpA = new TH1D("hsimNa22resUpA", "Na22 Simulation with resolution", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalNa22UpA = new TH1D("hcalNa22UpA", "Na22 Experiment Calibrated", binExp, xminCalUpA, xmaxCalUpA);
   
    TH1D* hsimNa22resUpB = new TH1D("hsimNa22resUpB", "Na22 Simulation with resolution", binExp, xminCalUpB, xmaxCalUpB);
    TH1D* hcalNa22UpB = new TH1D("hcalNa22UpB", "Na22 Experiment Calibrated", binExp, xminCalUpB, xmaxCalUpB);

    TH1D* hsimNa22resUpCoA = new TH1D("hsimNa22resUpCoA", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    TH1D* hsimNa22resUpCoC = new TH1D("hsimNa22resUpCoC", "Na22 Simulation with resolution", binExp, xminCal, xmaxCal);

    //! Fill experiment histograms
    for(int i = 0; i < entriesExpNa22; i++){
      texpNa22->GetEntry(i);
      hcalNa22->Fill(a*(chY[channel]+0.5) + b);
      hcalNa22UpA->Fill((a+aStep)*(chY[channel]+0.5) + b);
      hcalNa22UpB->Fill(a*(chY[channel]+0.5) + (b+bStep));
    }
    TH1D* hcalNa22Filtered = fft(hcalNa22, 0.1, 70., "hcalNa22Filtered", "hcalNa22Filtered", binExp, xminCal, xmaxCal);
    TH1D* hcalNa22UpAFiltered = fft(hcalNa22UpA, 0.1, 70., "hcalNa22UpAFiltered", "hcalNa22UpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hcalNa22UpBFiltered = fft(hcalNa22UpB, 0.1, 70., "hcalNa22UpBFiltered", "hcalNa22UpBFiltered", binExp, xminCalUpB, xmaxCalUpB);

    //! Fill Simulation histograms
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      double energy = ranGen->Gaus(y_sim,sigma);
      hsimNa22res->Fill(energy);
      hsimNa22resUpA->Fill(energy);
      hsimNa22resUpB->Fill(energy);
    }
    TH1D* hsimNa22resFiltered = fft(hsimNa22res, 0.1, 50., "hsimNa22resFiltered", "hsimNa22resFiltered", binExp, xminCal, xmaxCal);
    TH1D* hsimNa22resUpAFiltered = fft(hsimNa22resUpA, 0.1, 50., "hsimNa22resUpAFiltered", "hsimNa22resUpAFiltered", binExp, xminCalUpA, xmaxCalUpA);
    TH1D* hsimNa22resUpBFiltered = fft(hsimNa22resUpB, 0.1, 50., "hsimNa22resUpBFiltered", "hsimNa22resUpBFiltered", binExp, xminCalUpB, xmaxCalUpB);

    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA+coAStep,2) + pow(coB/sqrt(y_sim),2) + pow(coC/y_sim,2) );
      hsimNa22resUpCoA->Fill(ranGen->Gaus(y_sim,sigma));
    }
    TH1D* hsimNa22resUpCoAFiltered = fft(hsimNa22resUpCoA, 0.1, 50., "hsimNa22resUpCoAFiltered", "hsimNa22resUpCoAFiltered", binExp, xminCal, xmaxCal);
  
    for (int i = 0; i < entriesSimNa22; i++)
    {
      tsimNa22->GetEntry(i);
      double sigma = y_sim*sqrt( pow(coA,2) + pow(coB/sqrt(y_sim),2) + pow((coC+coCStep)/y_sim,2) );
      hsimNa22resUpCoC->Fill(ranGen->Gaus(y_sim,sigma));
    }
    TH1D* hsimNa22resUpCoCFiltered = fft(hsimNa22resUpCoA, 0.1, 50., "hsimNa22resUpCoCFiltered", "hsimNa22resUpCoCFiltered", binExp, xminCal, xmaxCal);

    hcalNa22Filtered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    hsimNa22resFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    
    if(count == 2){deltaChi2_1 = abs(hcalNa22Filtered->Chi2Test(hsimNa22resFiltered, "UU CHI2") - chi2_1274);}
    if(count == 3){deltaChi2_2 = abs(hcalNa22Filtered->Chi2Test(hsimNa22resFiltered, "UU CHI2") - chi2_1274);}
    if(count == 4){
      deltaChi2_3 = abs(hcalNa22Filtered->Chi2Test(hsimNa22resFiltered, "UU CHI2") - chi2_1274);
      deltaChi2 = deltaChi2_1 + deltaChi2_2 + deltaChi2_3;
    }
    if(count > 4){
      deltaChi2_1 = deltaChi2_2;
      deltaChi2_2 = deltaChi2_3;
      deltaChi2_3 = abs(hcalNa22Filtered->Chi2Test(hsimNa22resFiltered, "UU CHI2") - chi2_1274);
      deltaChi2 = deltaChi2_1 + deltaChi2_2 + deltaChi2_3;
    }

    std::cout << "\ndeltaChi2_1 = " << deltaChi2_1 
    << "\ndeltaChi2_2 = " << deltaChi2_2
    << "\ndeltaChi2_3 = " << deltaChi2_3
    << "\nTotal deltaChi2 = " << deltaChi2 << "\n";

    // deltaChi2 = abs(hcalNa22Filtered->Chi2Test(hsimNa22resFiltered, "UU CHI2") - chi2_1274);
    // std::cout << "delta Chi2 = " << deltaChi2 << "\n";
    chi2_1274 = hcalNa22Filtered->Chi2Test(hsimNa22resFiltered, "UU CHI2");

    grChi2->AddPoint(count, chi2_1274);
    cShow->cd(1);
    grChi2->Draw();
    cShow->cd(1)->Modified();
    cShow->cd(1)->Update();

    hcalNa22UpAFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    hsimNa22resUpAFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpA = hcalNa22UpAFiltered->Chi2Test(hsimNa22resUpAFiltered, "UU CHI2");

    hcalNa22UpBFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    hsimNa22resUpBFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpB = hcalNa22UpBFiltered->Chi2Test(hsimNa22resUpBFiltered, "UU CHI2");

    hsimNa22resUpCoAFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpCoA = hcalNa22Filtered->Chi2Test(hsimNa22resUpCoAFiltered, "UU CHI2");
    
    hsimNa22resUpCoCFiltered->GetXaxis()->SetRangeUser(minTest1274, maxTest1274);
    double chi2_1274_UpCoC = hcalNa22Filtered->Chi2Test(hsimNa22resUpCoCFiltered, "UU CHI2");

    double devChi2Na22UpA = (chi2_1274_UpA-chi2_1274)/aStep;
    double devChi2Na22UpB = (chi2_1274_UpB-chi2_1274)/bStep;
    double devChi2Na22UpCoA = (chi2_1274_UpCoA-chi2_1274)/coAStep;
    double devChi2Na22UpCoC = (chi2_1274_UpCoC-chi2_1274)/coCStep;

    std::cout << "devUpA = " << devChi2Na22UpA << "; devUpB = " << devChi2Na22UpB
    << "\ndevUpCoA = " << devChi2Na22UpCoA << "; devUpCoC = " << devChi2Na22UpCoC << "\n";

    a = a - tuningRateA*devChi2Na22UpA;
    b = b - tuningRateB*devChi2Na22UpB;
    coA = coA - tuningRateCoA*devChi2Na22UpCoA;
    coC = coC - tuningRateCoC*devChi2Na22UpCoC;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "a = " << a << "; b = " << b << "\n";
    std::cout << "coA = " << coA << "; coC = " << coC << "\n"; 

    double finalChannel1 = (E1-b)/a;
    double finalChannel2 = (E2-b)/a;

    double finalSig1 = 100.*sqrt(coA*coA + coC*coC/(E1*E1));
    double finalSig2 = 100.*sqrt(coA*coA + coC*coC/(E2*E2));

    std::cout << "Ch1 = " << finalChannel1 << "; Ch2 = " << finalChannel2
    << "\nSig1 = " << finalSig1 << "; Sig2 = " << finalSig2 << "\n";

    cShow->cd(2);
    hcalNa22Filtered->SetLineColor(kRed);
    hcalNa22Filtered->Draw();
    hsimNa22resFiltered->Scale(ScaleNa22, "noSW2");
    hsimNa22resFiltered->Draw("same");
    cShow->cd(2)->Modified();
    cShow->cd(2)->Update();

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

    if(deltaChi2 > 30.){
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
    }
    else{}
    
    count++;
  }

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}
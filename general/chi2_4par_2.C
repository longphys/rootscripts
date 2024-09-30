#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStopwatch.h"

void chi2_4par_2(){
  auto timer = new TStopwatch();
  timer->Start();

  TRandom3* rand3 = new TRandom3(0);

  TH1D* hGauss = new TH1D("hGauss", "Gauss Hist", 500, -3., 3.);
  // TF1* fGauss = new TF1("fGauss", "[0]*exp(- ((x-[1])*(x-[1]) / (2*c*c)) )", -5., 5.);

  // fGauss->SetParameter(0, 800.);
  // fGauss->SetParameter(1, 0.);
  // fGauss->SetParameter(2, 2.);
  // fGauss->SetNpx(100000);

  //! Model Gauss
  int fillTimes = 50000;
  for (int i = 1; i<=fillTimes; i++){
    hGauss->Fill(rand3->Gaus(0.82,0.63));
  }

  for (int i = 1; i <= fillTimes/2; i ++){
    hGauss->Fill(rand3->Gaus(-1.37, 0.45));
  }

  double sigmaMin1 = 0.00;
  double sigmaMax1 = 2.00;
  double sigmaStep1 = 0.01;
  int sigmaTimes1 = (sigmaMax1 - sigmaMin1)/sigmaStep1 + 1;
  double sigma1[sigmaTimes1];

  double meanMin1 = 0.00;
  double meanMax1 = 2.00;
  double meanStep1 = 0.01;
  int meanTimes1 = (meanMax1 - meanMin1)/meanStep1 + 1;
  double mean1[meanTimes1];

  double sigmaMin2 = 0.00;
  double sigmaMax2 = 2.00;
  double sigmaStep2 = 0.01;
  int sigmaTimes2 = (sigmaMax2 - sigmaMin2)/sigmaStep2 + 1;
  double sigma2[sigmaTimes2];

  double meanMin2 = -2.00;
  double meanMax2 = 0.00;
  double meanStep2 = 0.01;
  int meanTimes2 = (meanMax2 - meanMin2)/meanStep2 + 1;
  double mean2[meanTimes2];

  double intensity[sigmaTimes1][meanTimes1][sigmaTimes2][meanTimes2];

  int sigmaCurrent1_pos = 3;
  int meanCurrent1_pos = 3;
  int sigmaCurrent2_pos = 3;
  int meanCurrent2_pos = 3;

  double chi2current,
  chi2sigmaUp1, chi2sigmaDown1,
  chi2meanUp1, chi2meanDown1,
  chi2sigmaUp2, chi2sigmaDown2,
  chi2meanUp2, chi2meanDown2;

  double devSigmaUp1, devSigmaDown1,
  devMeanUp1, devMeanDown1,
  devSigmaUp2, devSigmaDown2,
  devMeanUp2, devMeanDown2;

  int rangeSigma1 = 1;
  int rangeMean1 = 1;
  int rangeSigma2 = 1;
  int rangeMean2 = 1;

  int times = 500;

  TH1D* hSave = new TH1D();

  for(int k = 1; k <= times; k++){
    TH1D* hChi2 = new TH1D("hChi2", "hChi2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2SigmaUp1 = new TH1D("hChi2SigmaUp1", "hChi2SigmaUp1", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2SigmaDown1 = new TH1D("hChi2SigmaDown1", "hChi2SigmaDown1", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2MeanUp1 = new TH1D("hChi2MeanUp1", "hChi2MeanUp1", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2MeanDown1 = new TH1D("hChi2MeanDown1", "hChi2MeanDown1", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2SigmaUp2 = new TH1D("hChi2SigmaUp2", "hChi2SigmaUp2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2SigmaDown2 = new TH1D("hChi2SigmaDown2", "hChi2SigmaDown2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2MeanUp2 = new TH1D("hChi2MeanUp2", "hChi2MeanUp2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2MeanDown2 = new TH1D("hChi2MeanDown2", "hChi2MeanDown2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
  
    // TCanvas* checkcanvas = new TCanvas("checkcanvas", "checkcanvas", 800, 600);
    // checkcanvas->Divide(5,2);
    //! CURRENT POSITION
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
      double mean = meanMin1 + meanCurrent1_pos*meanStep1;
      hChi2->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
      double mean = meanMin2 + meanCurrent2_pos*meanStep2;
      hChi2->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(9);
    // hChi2->Draw();

    //! SIGMA 1 UP
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + (sigmaCurrent1_pos+rangeSigma1)*sigmaStep1;
      double mean = meanMin1 + meanCurrent1_pos*meanStep1;
      hChi2SigmaUp1->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
      double mean = meanMin2 + meanCurrent2_pos*meanStep2;
      hChi2SigmaUp1->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(1);
    // hChi2SigmaUp1->Draw();

    //! SIGMA 1 DOWN
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + (sigmaCurrent1_pos-rangeSigma1)*sigmaStep1;
      double mean = meanMin1 + meanCurrent1_pos*meanStep1;
      hChi2SigmaDown1->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
      double mean = meanMin2 + meanCurrent2_pos*meanStep2;
      hChi2SigmaDown1->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(2);
    // hChi2SigmaDown1->Draw();

    //! MEAN 1 UP
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
      double mean = meanMin1 + (meanCurrent1_pos+rangeMean1)*meanStep1;
      hChi2MeanUp1->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
      double mean = meanMin2 + meanCurrent2_pos*meanStep2;
      hChi2MeanUp1->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(3);
    // hChi2MeanUp1->Draw();

    //! MEAN 1 DOWN
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
      double mean = meanMin1 + (meanCurrent1_pos-rangeMean1)*meanStep1;
      hChi2MeanDown1->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
      double mean = meanMin2 + meanCurrent2_pos*meanStep2;
      hChi2MeanDown1->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(4);
    // hChi2MeanDown1->Draw();

    //! SIGMA 2 UP
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
      double mean = meanMin1 + meanCurrent1_pos*meanStep1;
      hChi2SigmaUp2->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + (sigmaCurrent2_pos+rangeSigma2)*sigmaStep2;
      double mean = meanMin2 + meanCurrent2_pos*meanStep2;
      hChi2SigmaUp2->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(5);
    // hChi2SigmaUp2->Draw();

    //! SIGMA 2 DOWN
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
      double mean = meanMin1 + meanCurrent1_pos*meanStep1;
      hChi2SigmaDown2->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + (sigmaCurrent2_pos-rangeSigma2)*sigmaStep2;
      double mean = meanMin2 + meanCurrent2_pos*meanStep2;
      hChi2SigmaDown2->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(6);
    // hChi2SigmaDown2->Draw();

    //! MEAN 2 UP
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
      double mean = meanMin1 + meanCurrent1_pos*meanStep1;
      hChi2MeanUp2->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
      double mean = meanMin2 + (meanCurrent2_pos+rangeMean2)*meanStep2;
      hChi2MeanUp2->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(7);
    // hChi2MeanUp2->Draw();

    //! MEAN 2 DOWN
    for (int i = 1; i <=fillTimes; i++){
      double sigma = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
      double mean = meanMin1 + meanCurrent1_pos*meanStep1;
      hChi2MeanDown2->Fill(rand3->Gaus(mean,sigma));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      double sigma = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
      double mean = meanMin2 + (meanCurrent2_pos-rangeMean2)*meanStep2;
      hChi2MeanDown2->Fill(rand3->Gaus(mean,sigma));
    }

    // checkcanvas->cd(8);
    // hChi2MeanDown2->Draw();

    chi2current = hGauss->Chi2Test(hChi2, "UU CHI2");
    chi2sigmaUp1 = hGauss->Chi2Test(hChi2SigmaUp1, "UU CHI2");
    chi2sigmaDown1 = hGauss->Chi2Test(hChi2SigmaDown1, "UU CHI2");
    chi2meanUp1 = hGauss->Chi2Test(hChi2MeanUp1, "UU CHI2");
    chi2meanDown1 = hGauss->Chi2Test(hChi2MeanDown1, "UU CHI2");
    chi2sigmaUp2 = hGauss->Chi2Test(hChi2SigmaUp2, "UU CHI2");
    chi2sigmaDown2 = hGauss->Chi2Test(hChi2SigmaDown2, "UU CHI2");
    chi2meanUp2 = hGauss->Chi2Test(hChi2MeanUp2, "UU CHI2");
    chi2meanDown2 = hGauss->Chi2Test(hChi2MeanDown2, "UU CHI2");

    devSigmaUp1 = (chi2sigmaUp1-chi2current)/(rangeSigma1*sigmaStep1);
    devSigmaDown1 = (chi2sigmaDown1-chi2current)/(rangeSigma1*sigmaStep1);
    devMeanUp1 = (chi2meanUp1-chi2current)/(rangeMean1*meanStep1);
    devMeanDown1 = (chi2meanDown1-chi2current)/(rangeMean1*meanStep1);
    devSigmaUp2 = (chi2sigmaUp2-chi2current)/(rangeSigma2*sigmaStep2);
    devSigmaDown2 = (chi2sigmaDown2-chi2current)/(rangeSigma2*sigmaStep2);
    devMeanUp2 = (chi2meanUp2-chi2current)/(rangeMean2*meanStep2);
    devMeanDown2 = (chi2meanDown2-chi2current)/(rangeMean2*meanStep2);

    if (devSigmaUp1 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      sigmaCurrent1_pos = sigmaCurrent1_pos + rangeSigma1;
      hSave = (TH1D*)hChi2SigmaUp1->Clone();
    }
    if (devSigmaDown1 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      sigmaCurrent1_pos = sigmaCurrent1_pos - rangeSigma1;
      hSave = (TH1D*)hChi2SigmaDown1->Clone();
    }
    if (devMeanUp1 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      meanCurrent1_pos = meanCurrent1_pos + rangeMean1;
      hSave = (TH1D*)hChi2MeanUp1->Clone();
    }
    if (devMeanDown1 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      meanCurrent1_pos = meanCurrent1_pos - rangeMean1;
      hSave = (TH1D*)hChi2MeanDown1->Clone();
    }
    if (devSigmaUp2 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      sigmaCurrent2_pos = sigmaCurrent2_pos + rangeSigma2;
      hSave = (TH1D*)hChi2SigmaUp2->Clone();
    }
    if (devSigmaDown2 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      sigmaCurrent2_pos = sigmaCurrent2_pos - rangeSigma2;
      hSave = (TH1D*)hChi2SigmaDown2->Clone();
    }
    if (devMeanUp2 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      meanCurrent2_pos = meanCurrent2_pos + rangeMean2;
      hSave = (TH1D*)hChi2MeanUp2->Clone();
    }
    if (devMeanDown2 == std::min(devSigmaUp1, std::min(devSigmaDown1, std::min(devMeanUp1, std::min(devMeanDown1, std::min(devSigmaUp2, std::min(devSigmaDown2, std::min(devMeanUp2, devMeanDown2))))))) ){
      meanCurrent2_pos = meanCurrent2_pos - rangeMean2;
      hSave = (TH1D*)hChi2MeanDown2->Clone();
    }

    delete hChi2;
    delete hChi2SigmaUp1;
    delete hChi2SigmaDown1;
    delete hChi2MeanUp1;
    delete hChi2MeanDown1;
    delete hChi2SigmaUp2;
    delete hChi2SigmaDown2;
    delete hChi2MeanUp2;
    delete hChi2MeanDown2;
  }

  double finalSigma1 = sigmaMin1 + sigmaCurrent1_pos*sigmaStep1;
  double finalMean1 = meanMin1 + meanCurrent1_pos*meanStep1;
  double finalSigma2 = sigmaMin2 + sigmaCurrent2_pos*sigmaStep2;
  double finalMean2 = meanMin2 + meanCurrent2_pos*meanStep2;

  std::cout << "FinalMean1 = " << finalMean1 << "; FinalSigma1 = " << finalSigma1
  << "\nFinalMean2 = " << finalMean2 << "; FinalSigma2 = " << finalSigma2 << "\n";

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
  hGauss->Draw();
  hSave->SetLineColor(kRed);
  hSave->Draw("same");

  std::cout << "time = " << timer->RealTime() << " (s)\n";
}
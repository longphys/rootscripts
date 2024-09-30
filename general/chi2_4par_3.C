#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStopwatch.h"

void chi2_4par_3(){
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
  int fillTimes = 100000;
  for (int i = 1; i<=fillTimes; i++){
    hGauss->Fill(rand3->Gaus(0.82,0.63));
  }

  for (int i = 1; i <= fillTimes/2; i ++){
    hGauss->Fill(rand3->Gaus(-1.37, 0.45));
  }

  int sigmaCurrent1_pos = 3;
  int meanCurrent1_pos = 3;
  int sigmaCurrent2_pos = 3;
  int meanCurrent2_pos = 3;

  double chi2current,
  chi2sigmaUp1, chi2meanUp1, chi2sigmaUp2, chi2meanUp2;

  double devSigmaUp1, devMeanUp1, devSigmaUp2, devMeanUp2;

  double mean1 = 0.5;
  double mean2 = -0.5;
  double sigma1 = 0.5;
  double sigma2 =0.5;

  double sigmaStep1 = 0.01;
  double meanStep1 = 0.01;
  double sigmaStep2 = 0.01;
  double meanStep2 = 0.01;

  int rangeSigma1 = 1;
  int rangeMean1 = 1;
  int rangeSigma2 = 1;
  int rangeMean2 = 1;

  double changeRate = 0.000001;

  int times = 200;

  TH1D* hSave = new TH1D();
  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);

  for(int k = 1; k <= times; k++){
    TH1D* hChi2 = new TH1D("hChi2", "hChi2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2SigmaUp1 = new TH1D("hChi2SigmaUp1", "hChi2SigmaUp1", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2MeanUp1 = new TH1D("hChi2MeanUp1", "hChi2MeanUp1", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2SigmaUp2 = new TH1D("hChi2SigmaUp2", "hChi2SigmaUp2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2MeanUp2 = new TH1D("hChi2MeanUp2", "hChi2MeanUp2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
  
    //! CURRENT POSITION
    for (int i = 1; i <=fillTimes; i++){
      hChi2->Fill(rand3->Gaus(mean1,sigma1));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      hChi2->Fill(rand3->Gaus(mean2,sigma2));
    }

    double sigma1up = sigma1 + sigmaStep1;
    //! SIGMA 1 UP
    for (int i = 1; i <=fillTimes; i++){
      hChi2SigmaUp1->Fill(rand3->Gaus(mean1,sigma1up));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      hChi2SigmaUp1->Fill(rand3->Gaus(mean2,sigma2));
    }

    double mean1up = mean1 + meanStep1;
    //! MEAN 1 UP
    for (int i = 1; i <=fillTimes; i++){
      hChi2MeanUp1->Fill(rand3->Gaus(mean1up,sigma1));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      hChi2MeanUp1->Fill(rand3->Gaus(mean2,sigma2));
    }

    double sigma2up = sigma2 + sigmaStep2;
    //! SIGMA 2 UP
    for (int i = 1; i <=fillTimes; i++){
      hChi2SigmaUp2->Fill(rand3->Gaus(mean1,sigma1));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      hChi2SigmaUp2->Fill(rand3->Gaus(mean2,sigma2up));
    }

    double mean2up = mean2 + meanStep2;
    //! MEAN 2 UP
    for (int i = 1; i <=fillTimes; i++){
      hChi2MeanUp2->Fill(rand3->Gaus(mean1,sigma1));
    }
    for (int i = 1; i <=fillTimes/2; i++){
      hChi2MeanUp2->Fill(rand3->Gaus(mean2up,sigma2));
    }

    chi2current = hGauss->Chi2Test(hChi2, "UU CHI2");
    chi2sigmaUp1 = hGauss->Chi2Test(hChi2SigmaUp1, "UU CHI2");
    chi2meanUp1 = hGauss->Chi2Test(hChi2MeanUp1, "UU CHI2");
    chi2sigmaUp2 = hGauss->Chi2Test(hChi2SigmaUp2, "UU CHI2");
    chi2meanUp2 = hGauss->Chi2Test(hChi2MeanUp2, "UU CHI2");

    devSigmaUp1 = (chi2sigmaUp1-chi2current)/(rangeSigma1*sigmaStep1);
    devMeanUp1 = (chi2meanUp1-chi2current)/(rangeMean1*meanStep1);
    devSigmaUp2 = (chi2sigmaUp2-chi2current)/(rangeSigma2*sigmaStep2);
    devMeanUp2 = (chi2meanUp2-chi2current)/(rangeMean2*meanStep2);

    // std::cout << "devMeanUp1 = " << devMeanUp1 << "; devSigma1 = " << devSigmaUp1
    // << "\ndevMeanUp2 = " << devMeanUp2 << "; devSigmaUp2 = " << devSigmaUp2 << "\n";
   
    // std::cout << "FinalMean1 = " << mean1 << "; FinalSigma1 = " << sigma1
    // << "\nFinalMean2 = " << mean2 << "; FinalSigma2 = " << sigma2 << "\n";

    mean1 = mean1 - changeRate*devMeanUp1;
    mean2 = mean2 - changeRate*devMeanUp2;
    sigma1 = sigma1 - changeRate*devSigmaUp1;
    sigma2 = sigma2 -changeRate*devSigmaUp2;

    hSave = (TH1D*)hChi2->Clone();
        
    canvas2->cd();
    hGauss->Draw();
    hSave->SetLineColor(kRed);
    hSave->Draw("same");
    canvas2->Update();

    // std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    // std::cout << "FinalMean1 = " << mean1 << "; FinalSigma1 = " << sigma1
    // << "\nFinalMean2 = " << mean2 << "; FinalSigma2 = " << sigma2 << "\n";

    delete hChi2;
    delete hChi2SigmaUp1;
    delete hChi2MeanUp1;
    delete hChi2SigmaUp2;
    delete hChi2MeanUp2;
  }



  std::cout << "time = " << timer->RealTime() << " (s)\n";
}
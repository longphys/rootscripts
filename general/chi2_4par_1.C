#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TLegend.h"

void chi2_4par_1(){
  auto timer = new TStopwatch();
  timer->Start();

  TRandom3* rand3 = new TRandom3(0);

  TH1D* hGauss = new TH1D("hGauss", "Gauss Hist", 500, -3., 3.);
  // TF1* fGauss = new TF1("fGauss", "[0]*exp(- ((x-[1])*(x-[1]) / (2*c*c)) )", -5., 5.);

  // fGauss->SetParameter(0, 800.);
  // fGauss->SetParameter(1, 0.);
  // fGauss->SetParameter(2, 2.);
  // fGauss->SetNpx(100000);

  int fillTimes = 50000;
  for (int i = 1; i<=fillTimes; i++){
    hGauss->Fill(rand3->Gaus(0.82,0.63));
  }

  for (int i = 1; i <= fillTimes/2; i ++){
    hGauss->Fill(rand3->Gaus(-1.37, 0.45));
  }

  double currentSigma1;
  double sigmaMin1 = 0.50;
  double sigmaMax1 = 1.00;
  double sigmaStep1 = 0.01;
  int sigmaTimes1 = (sigmaMax1 - sigmaMin1)/sigmaStep1 + 1;
  double sigma1[sigmaTimes1];

  double currentMean1;
  double meanMin1 = 0.50;
  double meanMax1 = 1.00;
  double meanStep1 = 0.01;
  int meanTimes1 = (meanMax1 - meanMin1)/meanStep1 + 1;
  double mean1[meanTimes1];

  double currentSigma2;
  double sigmaMin2 = 0.00;
  double sigmaMax2 = 0.50;
  double sigmaStep2 = 0.01;
  int sigmaTimes2 = (sigmaMax2 - sigmaMin2)/sigmaStep2 + 1;
  double sigma2[sigmaTimes2];

  double currentMean2;
  double meanMin2 = -1.50;
  double meanMax2 = -1.00;
  double meanStep2 = 0.01;
  int meanTimes2 = (meanMax2 - meanMin2)/meanStep2 + 1;
  double mean2[meanTimes2];

  double intensity[sigmaTimes1][meanTimes1][sigmaTimes2][meanTimes2];

  double chi2Save = 10000;

  TH1D* hSave = new TH1D();
  double sigmaSave1, meanSave1, sigmaSave2, meanSave2;

  for (int j = 0; j<sigmaTimes1; j++){
    for (int k = 0; k<meanTimes1; k++){
      for (int l = 0; l<sigmaTimes2; l++){
        for (int m = 0; m<meanTimes2; m++){
          TH1D* hRandom = new TH1D("hRandom", "Random Hist", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
          for (int i = 1; i<=fillTimes; i++){
            currentSigma1 = sigmaMin1 + j*sigmaStep1;
            currentMean1 = meanMin1 + k*meanStep1;
            hRandom->Fill(rand3->Gaus(currentMean1,currentSigma1));
          }
          for (int i = 1; i<=fillTimes/2; i++){
            currentSigma2 = sigmaMin2 + l*sigmaStep2;
            currentMean2 = meanMin2 + m*meanStep2;
            hRandom->Fill(rand3->Gaus(currentMean2,currentSigma2));
          }

          double chi2 = hGauss->Chi2Test(hRandom, "UU CHI2/NDF");
          // double chi2 = hGauss->Chi2Test(hRandom, "UU CHI2/NDF");
          // std::cout << "chi2[" << l << "] = " << chi2 << "\n";

          if (chi2 < chi2Save){
            hSave = (TH1D*)hRandom->Clone();
            chi2Save = chi2;
            sigmaSave1 = currentSigma1;
            meanSave1 = currentMean1;
            sigmaSave2 = currentSigma2;
            meanSave2 = currentMean2;
          }

          sigma1[j] = currentSigma1;
          mean1[k] = currentMean1;
          sigma2[l] = currentSigma2;
          mean2[m] = currentMean2;
          intensity[j][k][l][m] = chi2;

          // for(int i = 1; i<=hRandom->GetNbinsX(); i++){
          //   hRandom->SetBinContent(i, 0.);
          // }
          delete hRandom;
        }
      }
    }
  }

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
  hGauss->Draw();
  hSave->SetLineColor(kRed);
  hSave->Draw("same");

  TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg1->SetHeader("Gaussians", "C");
  leg1->SetBorderSize(2);
  leg1->AddEntry(hGauss, "Model", "l");
  leg1->AddEntry(hSave, "Chi2 map", "l");
  leg1->Draw();

  std::cout << "\nminimum chi2 = " << chi2Save << "\n";
  std::cout << "meanSave1 = " << meanSave1 << "; sigmaSave1 = " << sigmaSave1 << "\n";
  std::cout << "meanSave2 = " << meanSave2 << "; sigmaSave2 = " << sigmaSave2 << "\n";

  std::cout << "time = " << timer->RealTime() << " (s)\n";
}
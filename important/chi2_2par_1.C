#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TLegend.h"

void chi2_2par_1(){
  auto timer = new TStopwatch();
  timer->Start();

  TRandom3* rand3 = new TRandom3(0);

  TH1D* hGauss = new TH1D("hGauss", "Gauss Hist", 500, -3., 4.);

  int fillTimes = 100000;
  for (int i = 1; i<=fillTimes; i++){
    hGauss->Fill(rand3->Gaus(0.82,0.63));
  }

  double sigma;
  double sigmaMin = 0.00;
  double sigmaMax = 2.00;
  double sigmaStep = 0.05;
  int sigmaTimes = (sigmaMax - sigmaMin)/sigmaStep + 1;

  double mean;
  double meanMin = -1.00;
  double meanMax = 1.50;
  double meanStep = 0.05;
  int meanTimes = (meanMax - meanMin)/meanStep + 1;

  double x[sigmaTimes], y[meanTimes], z[sigmaTimes][meanTimes];

  double chi2Save = 10000;

  TH1D* hSave = new TH1D();
  double sigmaSave, meanSave;

  for (int k = 1; k<=meanTimes; k++){
    for (int l = 1; l<=sigmaTimes; l++){
      TH1D* hRandom = new TH1D("hRandom", "Random Hist", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
      for (int i = 1; i<=fillTimes/10; i++){
        mean = meanMin + k*meanStep;
        sigma = sigmaMin + l*sigmaStep;
        hRandom->Fill(rand3->Gaus(mean,sigma));
      }

      double chi2 = hGauss->Chi2Test(hRandom, "UU CHI2");

      if (chi2 < chi2Save){
        hSave = (TH1D*)hRandom->Clone();
        chi2Save = chi2;
        sigmaSave = sigma;
        meanSave = mean;
      }

      x[l-1] = sigma;
      y[k-1] = mean;
      z[l-1][k-1] = chi2;

      delete hRandom;
    }
  }

  TH2D* hChi2 = new TH2D("hChi2", "Chi2 mapping", sigmaTimes, sigmaMin, sigmaMax, meanTimes, meanMin, meanMax);
  for (int k = 1; k<=meanTimes; k++){
    for (int l = 1; l<=sigmaTimes; l++){
      hChi2->SetBinContent(l, k, z[l-1][k-1]);
    }
  }

  double meanNew = -0.7;
  double sigmaNew = 1.8;

  double meanStepNew = 0.01;
  double sigmaStepNew = 0.01;

  double changeRate = 0.000001;

  double chi2current, chi2meanUp, chi2meanDown, chi2sigmaUp, chi2sigmaDown;
  double devmeanUp, devmeanDown, devsigmaUp, devsigmaDown;

  int rangeX = 2;
  int rangeY = 2;

  int times = 100;

  TGraph* grPath = new TGraph();
  grPath->GetXaxis()->SetTitle("Sigma");
  grPath->GetYaxis()->SetTitle("Mean");

  TH1D* hSave1 = new TH1D();
  for (int k = 1; k <= times; k++){
    grPath->AddPoint(sigmaNew, meanNew);

    TH1D* hChi2 = new TH1D("hChi2", "hChi2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2meanUp = new TH1D("hChi2meanUp", "hChi2meanUp", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2meanDown = new TH1D("hChi2meanDown", "hChi2meanDown", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2sigmaUp = new TH1D("hChi2sigmaUp", "hChi2sigmaUp", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2sigmaDown = new TH1D("hChi2sigmaDown", "hChi2sigmaDown", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());

    for (int i = 1; i <=100000; i++){
      hChi2->Fill(rand3->Gaus(meanNew,sigmaNew));
    }

    double meanUp = meanNew + meanStep;
    for (int i = 1; i <=100000; i++){
      hChi2meanUp->Fill(rand3->Gaus(meanUp,sigmaNew));
    }

    double meanDown = meanNew - meanStep;
    for (int i = 1; i <=100000; i++){
      hChi2meanDown->Fill(rand3->Gaus(meanDown,sigmaNew));
    }
    
    double sigmaUp = sigmaNew + sigmaStep;
    for (int i = 1; i <=100000; i++){
      hChi2sigmaUp->Fill(rand3->Gaus(meanNew,sigmaUp));
    }

    double sigmaDown = sigmaNew - sigmaStep;
    for (int i = 1; i <=100000; i++){
      hChi2sigmaDown->Fill(rand3->Gaus(meanNew,sigmaDown));
    }

    chi2current = hGauss->Chi2Test(hChi2, "UU CHI2");
    chi2meanUp = hGauss->Chi2Test(hChi2meanUp, "UU CHI2");
    chi2meanDown = hGauss->Chi2Test(hChi2meanDown, "UU CHI2");
    chi2sigmaUp = hGauss->Chi2Test(hChi2sigmaUp, "UU CHI2");
    chi2sigmaDown = hGauss->Chi2Test(hChi2sigmaDown, "UU CHI2");

    double devMean = (chi2meanUp - chi2meanDown)/(2*meanStep);
    double devSigma = (chi2sigmaUp - chi2sigmaDown)/(2*sigmaStep);

    // std::cout << "dev mean = " << devMean << "; dev sigma = " << devSigma << "\n";

    meanNew = meanNew - changeRate*devMean;
    sigmaNew = sigmaNew - changeRate*devSigma;

    std::cout << "new mean = " << meanNew << "; new sigma = " << sigmaNew << "\n";

    hSave1 = (TH1D*)hChi2->Clone();

    delete hChi2;
    delete hChi2meanUp;
    delete hChi2meanDown;
    delete hChi2sigmaUp;
    delete hChi2sigmaDown;
  }

  TCanvas* canvas1 = new TCanvas("canvas1", "Compare Chi2", 800, 600);
  canvas1->cd();
  hChi2->GetXaxis()->SetTitle("Sigma");
  hChi2->GetYaxis()->SetTitle("Mean");
  hChi2->Draw("colz");
  grPath->SetLineColor(kRed);
  grPath->SetLineWidth(2);
  grPath->SetLineStyle(1);
  grPath->Draw("same");

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);
  hGauss->SetLineColor(kBlack);
  hGauss->Scale(1./hGauss->Integral(1,hGauss->GetNbinsX()), "noSW2");
  hGauss->Draw();
  hSave->Scale(1./hSave->Integral(1,hSave->GetNbinsX()), "noSW2");
  hSave->SetLineColor(kRed);
  hSave->Draw("same");
  hSave1->Scale(1./hSave1->Integral(1,hSave1->GetNbinsX()), "noSW2");
  hSave1->SetLineColor(kBlue);
  hSave1->Draw("same");

  TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  leg1->SetHeader("Gaussians", "C");
  leg1->SetBorderSize(2);
  leg1->AddEntry(hGauss, "Model", "l");
  leg1->AddEntry(hSave1, "Gradient descent", "l");
  leg1->AddEntry(hSave, "Chi2 map", "l");
  leg1->Draw();

  std::cout << "\nminimum chi2 = " << chi2Save << "\n";
  std::cout << "meanSave = " << meanSave << "; sigmaSave = " << sigmaSave << "\n";

  std::cout << "time = " << timer->RealTime() << " (s)\n";
}
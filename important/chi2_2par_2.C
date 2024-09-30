#include <iostream>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TLegend.h"

void chi2_2par_2(){
  auto timer = new TStopwatch();
  timer->Start();

  TRandom3* rand3 = new TRandom3(0);

  TH1D* hGauss = new TH1D("hGauss", "Gauss Hist", 500, -10., 10.);

  int events = 100000;

  for (int i = 1; i<=events; i++){
    hGauss->Fill(rand3->Gaus(0.82,0.63));
  }

  double mean = -3.5;
  double sigma = 1.23;

  double meanStep = 0.01;
  double sigmaStep = 0.01;

  double changeRate = 0.000002;

  double chi2current, chi2meanUp, chi2meanDown, chi2sigmaUp, chi2sigmaDown;
  double devmeanUp, devmeanDown, devsigmaUp, devsigmaDown;

  double chi2manual, ndfcount;

  int times = 150;

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 800, 600);

  TH1D* hSave = new TH1D();

  TGraph* grPath = new TGraph();
  grPath->GetXaxis()->SetTitle("Mean");
  grPath->GetYaxis()->SetTitle("Sigma");

  for (int k = 1; k <= times; k++){
    grPath->AddPoint(mean, sigma);

    TH1D* hChi2 = new TH1D("hChi2", "hChi2", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2meanUp = new TH1D("hChi2meanUp", "hChi2meanUp", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2meanDown = new TH1D("hChi2meanDown", "hChi2meanDown", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2sigmaUp = new TH1D("hChi2sigmaUp", "hChi2sigmaUp", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());
    TH1D* hChi2sigmaDown = new TH1D("hChi2sigmaDown", "hChi2sigmaDown", hGauss->GetNbinsX(), hGauss->GetXaxis()->GetXmin(), hGauss->GetXaxis()->GetXmax());

    for (int i = 1; i <=events; i++){
      hChi2->Fill(rand3->Gaus(mean,sigma));
    }

    double meanUp = mean + meanStep;
    for (int i = 1; i <=events; i++){
      hChi2meanUp->Fill(rand3->Gaus(meanUp,sigma));
    }

    double meanDown = mean - meanStep;
    for (int i = 1; i <=events; i++){
      hChi2meanDown->Fill(rand3->Gaus(meanDown,sigma));
    }
    
    double sigmaUp = sigma + sigmaStep;
    for (int i = 1; i <=events; i++){
      hChi2sigmaUp->Fill(rand3->Gaus(mean,sigmaUp));
    }

    double sigmaDown = sigma - sigmaStep;
    for (int i = 1; i <=events; i++){
      hChi2sigmaDown->Fill(rand3->Gaus(mean,sigmaDown));
    }

    chi2manual = 0.;
    ndfcount = 0.;

    //! loop for chi-squared
    for (int i = 1; i <= hGauss->GetNbinsX(); i++){
      if(hGauss->GetBinContent(i) + hChi2->GetBinContent(i) != 0.){
        chi2manual += pow( hChi2->GetBinContent(i) - hGauss->GetBinContent(i), 2 ) / (hGauss->GetBinContent(i) + hChi2->GetBinContent(i));
        ndfcount ++;
      }
    }

    chi2manual = chi2manual/(ndfcount-1); //! formula for reduced chi-squared, NDF = n-m = bins-1

    chi2current = hGauss->Chi2Test(hChi2, "UU CHI2/NDF");
    chi2meanUp = hGauss->Chi2Test(hChi2meanUp, "UU CHI2");
    chi2meanDown = hGauss->Chi2Test(hChi2meanDown, "UU CHI2");
    chi2sigmaUp = hGauss->Chi2Test(hChi2sigmaUp, "UU CHI2");
    chi2sigmaDown = hGauss->Chi2Test(hChi2sigmaDown, "UU CHI2");

    std::cout << "step = " << k << "\n";
    std::cout << "current chi2 = " << chi2current << "\n";
    std::cout << "manual chi2 = " << chi2manual << "\n";

    // devmeanUp = (chi2meanUp-chi2current)/(meanRate*mean);
    // devmeanDown = (chi2meanDown-chi2current)/(meanRate*mean);
    // devsigmaUp = (chi2sigmaUp-chi2current)/(sigmaRate*sigma);
    // devsigmaDown =  (chi2sigmaDown-chi2current)/(sigmaRate*sigma);

    double devMean = (chi2meanUp - chi2meanDown)/(2*meanStep);
    double devSigma = (chi2sigmaUp - chi2sigmaDown)/(2*sigmaStep);

    std::cout << "dev mean = " << devMean << "; dev sigma = " << devSigma << "\n";

    mean = mean - changeRate*devMean;
    sigma = sigma - changeRate*devSigma;

    std::cout << "new mean = " << mean << "; new sigma = " << sigma << "\n";

    hSave = (TH1D*)hChi2->Clone();


    delete hChi2;
    delete hChi2meanUp;
    delete hChi2meanDown;
    delete hChi2sigmaUp;
    delete hChi2sigmaDown;

    canvas2->cd();
    hGauss->Draw();
    hSave->SetLineColor(kRed);
    hSave->Draw("same");
    canvas2->cd();
    canvas2->Modified();
    canvas2->Update();
    canvas2->Print("output.gif+5");
  }

  // std::cout << "minimum chi2 = " << chi2current << "\n";
  // std::cout << "final mean = " << mean << "; final sigma = " << sigma << "\n";

  // TLegend *leg1 = new TLegend(0.75, 0.6, 0.98, 0.75);
  // leg1->SetHeader("Gaussians", "C");
  // leg1->SetBorderSize(2);
  // leg1->AddEntry(hGauss, "Model", "l");
  // leg1->AddEntry(hSave, "Gradient descent", "l");
  // leg1->Draw();

  std::cout << "time = " << timer->RealTime() << " (s)\n";

}
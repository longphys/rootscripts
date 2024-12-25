#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TStopwatch.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TLine.h"

//! Range of experimental Cs137 and Na22 histograms (Bins)
double xminExp = 100;
double xmaxExp = 800;
int binExp = xmaxExp-xminExp;

//! Threshold parameters
// double para_k1 = 0.1;
double para_k1 = 0.1;
double para_c1 = 100.0;

double para_k2 = 0.1;
double para_c2 = 100.0;

//! Compton edge position range (Bins)
double minCompt511 = 100;
double maxCompt511 = 150;

double minCompt1274 = 320;
double maxCompt1274 = 370;

double minCompt662 = 130;
double maxCompt662 = 160;

//! Derivative = 0 position range (Bins)
double min0Der511 = 60;
double max0Der511 = 90;

double min0Der1274 = 150;
double max0Der1274 = 250;

double min0Der662 = 80;
double max0Der662 = 110;

//! Compton edge by energy
double physE_c[3] = {340.6, 1061.7, 477.3};

//! tanTheta on derivative fiting range (Bins)
double minThetafit511 = 170;
double maxThetafit511 = 200;

double minThetafit1274 = 550;
double maxThetafit1274 = 700;

double minThetafit662 = 245;
double maxThetafit662 = 290;

//! Initial relative resolution for fitting
double relRes662 = 0.1;
double relRes1274 = 0.06;

//! Energy spectrum fitting range (keV)
double minEFit1274 = 550.;
double maxEFit1274 = 1300.;

double minEFit662 = 270.;
double maxEFit662 = 700.;

TF1* fit = new TF1("fit", "[0]*x+[1]");

TH1D* hChannel1 = new TH1D();
TH1D* hChannel2 = new TH1D();

TH1D* hFilDiff5Channel1 = new TH1D("hFilDiff5Channel1", "Na22 5 points derivative", binExp, xminExp, xmaxExp);
TH1D* hFilDiff5Channel2 = new TH1D("hFilDiff5Channel2", "Cs137 5 points derivative", binExp, xminExp, xmaxExp);

TH1D* hFiltered1 = new TH1D("hFiltered1", "Exp data Diff filtered", binExp, xminExp, xmaxExp);
TH1D* hFiltered2 = new TH1D("hFiltered2", "Exp data Diff filtered", binExp, xminExp, xmaxExp);

int id[3] = {0, 0};
int id2[3] = {0, 0};

double func(double a, double b, int i){
  double xminCal = xminExp*a+b;
  double xmaxCal = xmaxExp*a+b;
  // std::cout << "xminCal = " << xminCal <<"; xmaxCal = " << xmaxCal << "\n";

  TH1D* hCal1 = new TH1D("hcal1", "Na22 Calibrated Spectrum", binExp, xminCal, xmaxCal);
  for (int i = 1; i <= binExp; i++){
    hCal1->SetBinContent(i, hChannel1->GetBinContent(i));
  }
  hCal1->SetXTitle("Energy(keV)");
  hCal1->SetYTitle("Count");
  hCal1->Draw();

  // TMarker* mark;
  // for (int i = 0; i <= 1; i++){
  //   mark = new TMarker((id[i]+100)*a_cal+b_cal, hCal1->GetBinContent(id[i]), 8);
  //   mark->SetMarkerColor(kRed);
  //   mark->Draw();
  // }

  // canvas5->cd(3);
  TH1D* hCal2 = new TH1D("hcal2", "Cs137 Calibrated Spectrum", binExp, xminCal, xmaxCal);
  for (int i = 1; i <= binExp; i++){
    hCal2->SetBinContent(i, hChannel2->GetBinContent(i));
  }
  hCal2->SetXTitle("Energy(keV)");
  hCal2->SetYTitle("Count");
  hCal2->Draw();

  // mark = new TMarker((id[2]+100)*a_cal+b_cal, hCal2->GetBinContent(id[2]), 8);
  // mark->SetMarkerColor(kRed);
  // mark->Draw();
  
  // canvas5->cd(4);
  TH1D* hCalDiff1 = new TH1D("hCalDiff1", "Diff1", binExp, xminCal, xmaxCal);
  for (int i = 1; i <= binExp; i++){
    hCalDiff1->SetBinContent(i, hFilDiff5Channel1->GetBinContent(i));
  }

  //! tanTheta
  double tanTheta[3];
  // TF1* pol = new TF1("pol", "");
  // hCalDiff5->Fit("pol", "", "", 170, 200);
  hCalDiff1->Fit("fit", "", "", minThetafit511, maxThetafit511);
  tanTheta[0] = fit->GetParameter(0);

  // canvas5->cd(5);
  TH1D* hCalDiff2 = (TH1D*)hCalDiff1->Clone();
  hCalDiff2->Fit("fit", "", "", minThetafit1274, maxThetafit1274);
  tanTheta[1] = fit->GetParameter(0);

  // canvas5->cd(6);
  TH1D* hCalDiff3 = new TH1D("hCalDiff3", "Diff3", binExp, xminCal, xmaxCal);
  for (int i = 1; i <= binExp; i++){
    hCalDiff3->SetBinContent(i, hFilDiff5Channel2->GetBinContent(i));
  }
  hCalDiff3->Fit("fit", "", "", minThetafit662, maxThetafit662);
  tanTheta[2] = fit->GetParameter(0);

  // std::cout << "\ntanTheta[0] = " << tanTheta[0] << "; tanTheta[1] = " << tanTheta[1] << "; tanTheta[2] = " << tanTheta[2] << "\n";

  double A[3], B[3], C[3], D[3], E_c[3];
  // //! 511
  A[0] = hCal1->GetBinContent(id2[1]);
  B[0] = hCal1->GetBinCenter(id2[0]);
  C[0] = hCal1->GetBinContent(id[1]);
  D[0] = tanTheta[0];
  E_c[0] = hCal1->GetBinCenter(id[0]);

  //! 1274
  A[1] = hCal1->GetBinContent(id2[1]);
  B[1] = hCal1->GetBinCenter(id2[1]);
  C[1] = hCal1->GetBinContent(id[1]);
  D[1] = tanTheta[1];
  E_c[1] = hCal1->GetBinCenter(id[1]);
  
  //! 662
  A[2] = hCal2->GetBinContent(id2[2]);
  B[2] = hCal2->GetBinCenter(id2[2]);
  C[2] = hCal2->GetBinContent(id[2]);
  D[2] = tanTheta[2];
  E_c[2] = hCal2->GetBinCenter(id[2]);

  // std::cout << "\nDiff=0 height A = " << A[1] << "\n";
  // std::cout << "Diff=0 position B = " << B[1] << " (keV)\n";
  // std::cout << "Compton edge height C = " << C[1] << "\n";
  // std::cout << "tanTheta D = " << D[1] << "\n";
  // std::cout << "Compton edge position E_c = " << E_c[1] << " (keV)\n";

  // TCanvas* canvas6 = new TCanvas("canvas6", "test", 1400, 700);
  // canvas6->Divide(2, 1);

  double a_res[3], b_res[3], c_res[3];
  for(int i = 0; i <= 2; i++){
    a_res[i] = D[i];
    b_res[i] = -2*B[i]*D[i];
    c_res[i] = A[i] + pow(B[i],2)*D[i];    
  }

  // std::cout << "a_res = " << a_res << "; b_res = " << b_res << "; c_res = " << c_res << "\n";
  
  TF1* newfit1 = new TF1("newfit1", "[5] + [6]*x + [7]*x*x + (1./2.)*([0]*(x*x + [4]*[4]) + [1]*x + [2])*erfc((x - [3])/(sqrt(2.)*[4])) + ((-[4]/sqrt(2.*TMath::Pi()))*[0]*(x + [3]) + [1])*exp(-pow((x - [3]),2)/(2.*[4]*[4]))", xminCal, xmaxCal);
  newfit1->SetParameter(0, a_res[1]);
  newfit1->SetParameter(1, b_res[1]);
  newfit1->SetParameter(2, c_res[1]);
  newfit1->SetParameter(3, E_c[1]);
  newfit1->SetParameter(4, E_c[1]*relRes1274);
  // newfit1->SetParameter(5, 200.);
  newfit1->SetNpx(10000);

  TF1* newfit2 = new TF1("newfit2", "[5] + [6]*x + [7]*x*x + (1./2.)*([0]*(x*x + [4]*[4]) + [1]*x + [2])*erfc((x - [3])/(sqrt(2.)*[4])) + ((-[4]/sqrt(2.*TMath::Pi()))*[0]*(x + [3]) + [1])*exp(-pow((x - [3]),2)/(2.*[4]*[4]))", xminCal, xmaxCal);
  newfit2->SetParameter(0, a_res[2]);
  newfit2->SetParameter(1, b_res[2]);
  newfit2->SetParameter(2, c_res[2]);
  newfit2->SetParameter(3, E_c[2]);
  newfit2->SetParameter(4, E_c[2]*relRes662);
  // newfit2->SetParameter(5, 200.);
  newfit2->SetNpx(10000);

  // canvas6->cd();
  TH1D* hFit1 = new TH1D("hFit1", "hFit1", binExp, xminCal, xmaxCal);
  for(int i = 1; i <= binExp; i++){
    hFit1->SetBinContent(i, hFiltered1->GetBinContent(i));
  }
  hFit1->SetXTitle("Energy(keV)");
  hFit1->SetYTitle("Count");

  TH1D* hFit2 = new TH1D("hFit2", "hFit2", binExp, xminCal, xmaxCal);
  for(int i = 1; i <= binExp; i++){
    hFit2->SetBinContent(i, hFiltered2->GetBinContent(i));
  }
  hFit2->SetXTitle("Energy(keV)");
  hFit2->SetYTitle("Count");

  // canvas6->cd(1);
  std::cout << "\nNa22 R FUNCTION FITTING\n";
  hFit1->Fit("newfit1", "", "", minEFit1274, maxEFit1274);
  // hCal1->Fit("newfit1", "", "", minEFit1274, maxEFit1274);

  // canvas6->cd(2);
  std::cout << "\nCs137 R FUNCTION FITTING\n";
  hFit2->Fit("newfit2", "", "", minEFit662, maxEFit662);
  // hCal2->Fit("newfit2", "", "", minEFit662, maxEFit662);
  // newfit->Draw();

  // std::cout << "\nInitial parameters\na = " << a_res[1] << "; " << a_res[2]
  // << "\nb = " << b_res[1] << "; " << b_res[2] 
  // << "\nc = " << c_res[1] << "; " << c_res[2]
  // << "\nE_c = " << E_c[1] << "; " << E_c[2];

  // TCanvas* canvas7 = new TCanvas("canvas7", "Fit", 1200, 600);
  // canvas7->Divide(2,2);

  TH1D* hFitSigma = new TH1D("hFitSigma", "Fit Sigma", binExp, xminCal/1000., xmaxCal/1000.);

  double newE_c[3], newSigma[3];
  newE_c[1] = newfit1->GetParameter(3);
  newE_c[2] = newfit2->GetParameter(3);
  newSigma[1] = newfit1->GetParameter(4);
  newSigma[2] = newfit2->GetParameter(4);

  // newE_c[0] = 340.;
  // newSigma[0] = 38.5;

  TF1* fRes = new TF1("fRes", "sqrt([0]*[0]/x + [1]*[1]/(x*x))");
  // TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");
  // hFitSigma->SetBinContent(hCal1->FindBin(newE_c[0]), pow(newSigma[0],2) );
  hFitSigma->SetBinContent(hCal1->FindBin(newE_c[1]), newSigma[1]/newE_c[1] );
  hFitSigma->SetBinContent(hCal2->FindBin(newE_c[2]), newSigma[2]/newE_c[2] );
  // canvas7->cd(1);
  hFitSigma->Fit("fRes", "", "", 0.4, 1.1);
  // hFitSigma->Draw();

  std::cout << "\nFIRST RESOLUTION PARAMETERS:\n"; 
  std::cout << "RESco_B = " << fRes->GetParameter(0) << "; RESco_C = " << fRes->GetParameter(1) << "\n\n";
  // std::cout << "RESco_A = " << fRes->GetParameter(0) << "; RESco_C = " << fRes->GetParameter(1) << "\n";

  TH1D* hFitnewE_c = new TH1D("hFitnewE_c", "Fit new E_c", binExp, xminCal, xmaxCal);
  hFitnewE_c->SetBinContent(hCal1->FindBin(newE_c[1]), physE_c[1]);
  hFitnewE_c->SetBinContent(hCal2->FindBin(newE_c[2]), physE_c[2]);
  // canvas7->cd(2);
  hFitnewE_c->Fit("fit", "", "", 400, 1100);
  hFitnewE_c->Draw();

  double newa_cal = a*fit->GetParameter(0), newb_cal = b*fit->GetParameter(0) + fit->GetParameter(1);

  // std::cout << "(MeV)new a_cal = " << newa_cal << "; new b_cal = " << newb_cal << "\n";
  std::cout << "\n(keV)new a_cal = " << newa_cal/1000. << "; new b_cal = " << newb_cal/1000. << "\n";

  switch(i){
    case 0:
      return newa_cal;
    case 1:
      return newb_cal;
    case 2:
      return fRes->GetParameter(0);
    case 3:
      return fRes->GetParameter(1);
    default:
      return 0;
  }
  
  delete hCal1; 
  delete hCal2;
  delete hCalDiff1;
  delete hCalDiff2;
  delete hCalDiff3;
  delete newfit1;
  delete newfit2;
  delete hFit1;
  delete hFit2;
  delete hFitSigma;
  delete fRes;
  delete hFitnewE_c;
}
void fft1()
{
  auto timer = new TStopwatch();
  timer->Start();

  //! Canvas 1
  TCanvas *canvas1 = new TCanvas("canvas1", "Exp data show", 1500, 1000);
  canvas1->cd();
  canvas1->Divide(4, 2);

  //! Measurement
  canvas1->cd(1);
  TFile* fileOpen1 = new TFile("../../expfiles/new/save_stilbene_dividers_Plastic_ch2_Na22_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
  // TFile* fileOpen1 = new TFile("./expfiles/usable/save_stilbene_divider_Na22_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9841_run_0_0001.root", "read");

  //!
  TFile* file2 = new TFile("../../expfiles/time/6/plastic_detectors_Na22_ch0HV1899_ch1HV1915_ch2HV1950_ch0_9668_ch1_9806_ch2_9852_run_0_0001.root", "read");

  TTree* tree2 =  (TTree*) file2->Get("AnalysisxTree");

  UShort_t value_exp[48];
  tree2->SetBranchAddress("NeEvent.neutAmp[48]", &value_exp);

  TH1D* hist2 = new TH1D("hist2", "Na22 Measurement", binExp, xminExp, xmaxExp);

  // for (int i = 0; i<tree2->GetEntriesFast(); i++){
  for (int i = 0; i<300000; i++){
    tree2->GetEntry(i);
    if(value_exp[0]){hist2->Fill(value_exp[0]);}
  }

  hChannel1 = (TH1D*)hist2->Clone();
  // hChannel1 = (TH1D*)fileOpen1->Get("hCh1");
  hChannel1->SetLineColor(kRed);
  hChannel1->SetXTitle("Channel");
  hChannel1->SetYTitle("Count");
  hChannel1->Draw();

  canvas1->cd(5);
  TFile* fileOpen2 = new TFile("../../expfiles/new/save_stilbene_dividers_Plastic_ch2_Cs137_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_0001.root", "read");
  
  hChannel2 = (TH1D*)fileOpen2->Get("hCh1");
  hChannel2->SetLineColor(kRed);
  hChannel2->SetXTitle("Channel");
  hChannel2->SetYTitle("Count");
  hChannel2->Draw();

  //! Pad3: Even function
  int binExpEven = 2*binExp ;

  double binContent, binContentPlus, binLength, newBinContent;
  binLength = (xmaxExp - xminExp)/ binExp;

  canvas1->cd(2);
  TH1D* hEvenChannel1 = new TH1D("hEvenChannel1", "Transformed to even function", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);

  for(int i = 1; i <= binExpEven/2; i++)
  {
    binContent = hChannel1->GetBinContent(i);
    // binContent = hDiffChannel->GetBinContent(i);
    hEvenChannel1->SetBinContent(i, binContent);
    hEvenChannel1->SetBinContent(binExpEven - i, binContent);
  }
  hEvenChannel1->SetBinContent(binExp, hEvenChannel1->GetBinContent(binExp + 1));
  hEvenChannel1->Draw();

  canvas1->cd(6);
  TH1D* hEvenChannel2 = new TH1D("hEvenChannel2", "Transformed to even function", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);

  for(int i = 1; i <= binExpEven/2; i++)
  {
    binContent = hChannel2->GetBinContent(i);
    // binContent = hDiffChannel->GetBinContent(i);
    hEvenChannel2->SetBinContent(i, binContent);
    hEvenChannel2->SetBinContent(binExpEven - i, binContent);
  }
  hEvenChannel2->SetBinContent(binExp, hEvenChannel2->GetBinContent(binExp + 1));
  hEvenChannel2->Draw();  

  TF1* fLogis1 = new TF1("fLogis1", "1./(1.+exp([0]*(x-[1])))", 0, binExpEven);
  fLogis1->SetParameter(0, para_k1);
  fLogis1->SetParameter(1, para_c1);
  fLogis1->SetNpx(10000);

  TF1* fLogis2 = new TF1("fLogis2", "1./(1.+exp([0]*(x-[1])))", 0, binExpEven);
  fLogis2->SetParameter(0, para_k2);
  fLogis2->SetParameter(1, para_c2);
  fLogis2->SetNpx(10000);
  // fLogis->Draw();

  //! Pad4: Magnitude
  canvas1->cd(3);
  TH1 *hm1 =nullptr;
  TVirtualFFT::SetTransform(nullptr);
  hm1 = hEvenChannel1->FFT(hm1, "MAG");
  hm1->SetTitle("Magnitude of the FFT transformation");
  hm1->Draw();
  
  TH1D* newhm1 = new TH1D("newhm1", "newhm1", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    newhm1->SetBinContent(i, hm1->GetBinContent(i)/binExpEven);
  }
  // newhm->Draw();

  TH1D* hLogis1 = new TH1D("hLogis1", "Logis", binExpEven, 0, binExpEven);
  for(int i = 1; i <= binExpEven/2; i++)
  {
    hLogis1->SetBinContent(i, fLogis1->Eval(i));
    hLogis1->SetBinContent(binExpEven-i, fLogis1->Eval(i));
  }

  hLogis1->SetLineColor(kRed);
  hLogis1->SetLineWidth(3);
  hLogis1->Draw("same");

  TCanvas* CanvasFilter = new TCanvas("CanvasFilter", "CanvasFilter", 800, 600);
  CanvasFilter->cd()->SetLogy();
  CanvasFilter->cd()->SetLogx();
  CanvasFilter->cd()->SetTicks();
  hm1->Draw();

  hm1->Rebin(5);
  hm1->SetLineWidth(3);
  hm1->SetStats(0);
  hm1->GetXaxis()->SetRangeUser(0.,700.);

  hm1->GetXaxis()->SetTitle("Frequency (a. units)");
  hm1->GetXaxis()->SetLabelFont(42);
  hm1->GetXaxis()->SetTitleFont(52);
  hm1->GetXaxis()->SetTitleSize(0.04);
  hm1->GetXaxis()->CenterTitle(true);

  hm1->GetYaxis()->SetTitle("Magnitude");
  hm1->GetYaxis()->SetLabelFont(42);
  hm1->GetYaxis()->SetTitleFont(52);
  hm1->GetYaxis()->SetTitleSize(0.04);
  hm1->GetYaxis()->CenterTitle(true);

  // hLogis1->Scale(300000, "noSW2");
  hLogis1->Draw("same");

  TLegend *legend = new TLegend(0.4, 0.55, 0.8, 0.85);
  legend->SetBorderSize(0);
  legend->SetLineWidth(2);
  legend->SetTextSize(0.06);
  legend->AddEntry(hm1, "Fourier image magnitude", "l");
  legend->AddEntry(hLogis1, "Logistic function", "l");

  legend->Draw();

  //! Pad8: Apply threshold
  canvas1->cd(4);

  TVirtualFFT *fft1 = TVirtualFFT::GetCurrentTransform();  
  //Use the following method to get the full output:
  Double_t *re_full1 = new Double_t[binExpEven];
  Double_t *im_full1 = new Double_t[binExpEven];
 
  fft1->GetPointsComplex(re_full1, im_full1);

  for(int i = 1; i <= binExpEven; i++)
  {
    re_full1[i] = re_full1[i]*hLogis1->GetBinContent(i);
    im_full1[i] = im_full1[i]*hLogis1->GetBinContent(i);
    // im_full1[i] = im_full1[i];
  }

  //Now let's make a backward transform:
  TVirtualFFT *fft_back1 = TVirtualFFT::FFT(1, &binExpEven, "C2R M K");
  fft_back1->SetPointsComplex(re_full1,im_full1);
  fft_back1->Transform();
  TH1 *hb1 = nullptr;
  //Let's look at the output
  hb1 = TH1::TransformHisto(fft_back1,hb1,"RE");
  // hb1 = TH1::TransformHisto(fft_back1,hb1,"MAG");

  TH1D* newhb1 = new TH1D("newhb1", "newhb1", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    // std::cout << "value = " << hb->GetBinContent(i)/binExpEven << "\n";
    newhb1->SetBinContent(i, hb1->GetBinContent(i)/binExpEven);
  }
  newhb1->SetTitle("The backward transform result");
  newhb1->Draw();

  //! Pad4: Magnitude
  canvas1->cd(7);
  TH1 *hm2 =nullptr;
  TVirtualFFT::SetTransform(nullptr);
  hm2 = hEvenChannel2->FFT(hm2, "MAG");
  hm2->SetTitle("Magnitude of the transform");
  hm2->Draw();
  
  TH1D* newhm2 = new TH1D("newhm2", "newhm2", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    newhm2->SetBinContent(i, hm2->GetBinContent(i)/binExpEven);
  }
  // newhm->Draw();

  TH1D* hLogis2 = new TH1D("hLogis", "Logis", binExpEven, 0, binExpEven);
  for(int i = 1; i <= binExpEven/2; i++)
  {
    hLogis2->SetBinContent(i, fLogis2->Eval(i));
    hLogis2->SetBinContent(binExpEven-i, fLogis2->Eval(i));
  }

  hLogis2->SetLineColor(kRed);
  hLogis2->Draw("same");

  //! Pad8: Apply threshold
  canvas1->cd(8);

  TVirtualFFT *fft2 = TVirtualFFT::GetCurrentTransform();  
  //Use the following method to get the full output:
  Double_t *re_full2 = new Double_t[binExpEven];
  Double_t *im_full2 = new Double_t[binExpEven];
 
  fft2->GetPointsComplex(re_full2, im_full2);

  for(int i = 1; i <= binExpEven; i++)
  {
    re_full2[i] = re_full2[i]*hLogis2->GetBinContent(i);
    im_full2[i] = im_full2[i]*hLogis2->GetBinContent(i);
  }

  //Now let's make a backward transform:
  TVirtualFFT *fft_back2 = TVirtualFFT::FFT(1, &binExpEven, "C2R M K");
  fft_back2->SetPointsComplex(re_full2,im_full2);
  fft_back2->Transform();
  TH1 *hb2 = nullptr;
  //Let's look at the output
  hb2 = TH1::TransformHisto(fft_back2,hb2,"RE");

  TH1D* newhb2 = new TH1D("newhb2", "newhb2", binExpEven, -(xmaxExp - xminExp), xmaxExp - xminExp);
  for(int i = 1; i <= binExpEven; i++)
  {
    // std::cout << "value = " << hb->GetBinContent(i)/binExpEven << "\n";
    newhb2->SetBinContent(i, hb2->GetBinContent(i)/binExpEven);
  }
  newhb2->SetTitle("The backward transform result");
  newhb2->Draw();

  //! Canvas 2
  TCanvas* canvas2 = new TCanvas("canvas2", "Compare", 800, 600);
  canvas2->Divide(2, 2);

  //! Pad 1: Measurement
  canvas2->cd(1);
  hChannel1->Draw();

  canvas2->cd(3);
  hChannel2->Draw();

  //! Pad 2: Filtered
  canvas2->cd(2);
  
  for(int i = 1; i <= binExp; i++)
  {
    hFiltered1->SetBinContent(i, newhb1->GetBinContent(i));
  }

  hFiltered1->SetLineColor(kRed);
  hFiltered1->SetXTitle("Channel");
  hFiltered1->SetYTitle("Count");
  hFiltered1->Draw();

  //! RESULT CANVAS
  // TCanvas* CanvasFilter_result = new TCanvas("CanvasFilter_result", "CanvasFilter_result", 600, 900);
  // CanvasFilter_result->Divide(1,2);

  TCanvas* CanvasFilter_result = new TCanvas("CanvasFilter_result", "CanvasFilter_result", 1500, 700);
  CanvasFilter_result->Divide(2,1);
  CanvasFilter_result->cd(1)->SetTicks();
  CanvasFilter_result->cd(1);
  CanvasFilter_result->cd(1)->SetLeftMargin(0.05);
  CanvasFilter_result->cd(1)->SetRightMargin(0.05);
  hChannel1->Draw();

  hChannel1->SetLineWidth(3);
  hChannel1->SetTitle("");
  hChannel1->SetLineColor(kBlue);
  hChannel1->SetStats(0);
  // hChannel1->GetXaxis()->SetRangeUser(100., 500.);
  hChannel1->GetYaxis()->SetRangeUser(0.,1500.);

  hChannel1->GetXaxis()->SetTitle("Channel");
  hChannel1->GetXaxis()->SetLabelFont(42);
  hChannel1->GetXaxis()->SetTitleFont(52);
  hChannel1->GetXaxis()->SetTitleSize(0.04);
  hChannel1->GetXaxis()->CenterTitle(true);

  hChannel1->GetYaxis()->SetTitle("Events");
  hChannel1->GetYaxis()->SetLabelFont(42);
  hChannel1->GetYaxis()->SetTitleFont(52);
  hChannel1->GetYaxis()->SetTitleSize(0.04);
  hChannel1->GetYaxis()->CenterTitle(true);

  CanvasFilter_result->cd(2);
  CanvasFilter_result->cd(2)->SetTicks();
  CanvasFilter_result->cd(2)->SetLeftMargin(0.05);
  CanvasFilter_result->cd(2)->SetRightMargin(0.15);
  hFiltered1->Draw();

  hFiltered1->SetLineWidth(3);
  hFiltered1->SetTitle("");
  hFiltered1->SetStats(0);
  // hFiltered1->GetXaxis()->SetRangeUser(100., 500.);
  hFiltered1->GetYaxis()->SetRangeUser(0.,1500.);

  hFiltered1->GetXaxis()->SetTitle("Channel");
  hFiltered1->GetXaxis()->SetLabelFont(42);
  hFiltered1->GetXaxis()->SetTitleFont(52);
  hFiltered1->GetXaxis()->SetTitleSize(0.04);
  hFiltered1->GetXaxis()->CenterTitle(true);

  hFiltered1->GetYaxis()->SetTitle("Events");
  hFiltered1->GetYaxis()->SetLabelFont(42);
  hFiltered1->GetYaxis()->SetTitleFont(52);
  hFiltered1->GetYaxis()->SetTitleSize(0.04);
  hFiltered1->GetYaxis()->CenterTitle(true);

  TLegend *legend2 = new TLegend(0.45, 0.45, 0.8, 0.8);
  legend2->SetBorderSize(0);
  legend2->SetLineWidth(2);
  legend2->SetTextSize(0.04);
  legend2->AddEntry(hChannel1, "Sample {}^{137}Cs Spectrum", "l");
  legend2->AddEntry(hFiltered1, "Filtered {}^{137}Cs Spectrum", "l");

  legend2->Draw();

  canvas2->cd(4);
  
  for(int i = 1; i <= binExp; i++)
  {
    hFiltered2->SetBinContent(i, newhb2->GetBinContent(i));
  }

  hFiltered2->SetLineColor(kRed);
  hFiltered2->SetXTitle("Channel");
  hFiltered2->SetYTitle("Count");
  hFiltered2->Draw();

  //! Canvas 3
  TCanvas* canvas3 = new TCanvas("canvas3", "Compare Diff", 1500, 600);
  canvas3->Divide(2, 1);

  canvas3->cd(1);
  canvas3->cd(1)->SetGridx();
  canvas3->cd(1)->SetGridy();
  for(int i = 3; i <= binExp - 2; i++)
  {
    newBinContent = (-hFiltered1->GetBinContent(i+2) 
    + 8*hFiltered1->GetBinContent(i+1) 
    - 8*hFiltered1->GetBinContent(i-1) 
    + hFiltered1->GetBinContent(i-2))/12*binLength;
    hFilDiff5Channel1->SetBinContent(i, newBinContent);
  } 

  hFilDiff5Channel1->SetLineColor(kBlue);
  hFilDiff5Channel1->SetXTitle("Channel");
  hFilDiff5Channel1->Draw();
  hFiltered1->Draw("same");

  //! RESULT CANVAS
  TCanvas* CanvasDiff = new TCanvas("CanvasDiff", "CanvasDiff", 1200, 1000);

  CanvasDiff->cd()->SetTicks();
  CanvasDiff->cd()->SetGrid();
  CanvasDiff->cd()->SetLeftMargin(0.15);
  hFilDiff5Channel1->Draw();
  // hFiltered1->Scale(0.1, "noSW2");

  hFiltered1->SetLineColor(kBlack);
  hFiltered1->SetLineWidth(5);
  hFiltered1->Draw("same");

  hFilDiff5Channel1->SetLineWidth(5);
  hFilDiff5Channel1->SetTitle("");
  hFilDiff5Channel1->SetLineColor(kRed);
  hFilDiff5Channel1->SetStats(0);

  hFilDiff5Channel1->GetXaxis()->SetTitle("Channel");
  hFilDiff5Channel1->GetXaxis()->SetLabelFont(42);
  hFilDiff5Channel1->GetXaxis()->SetTitleFont(52);
  hFilDiff5Channel1->GetXaxis()->SetTitleSize(0.04);
  hFilDiff5Channel1->GetXaxis()->CenterTitle(true);

  TLegend *legendDiff = new TLegend(0.55, 0.55, 0.85, 0.85);
  legendDiff->SetBorderSize(0);
  legendDiff->SetLineWidth(3);
  legendDiff->SetTextSize(0.03);
  legendDiff->AddEntry(hFiltered1, "Filtered {}^{137}Cs Spectrum", "l");
  legendDiff->AddEntry(hFilDiff5Channel1, "Corresponding first derivative", "l");

  legendDiff->Draw();

  TLine* dotted_line = new TLine(496, -5.5, 496, 15);
  dotted_line->SetLineColor(kRed);
  dotted_line->SetLineStyle(9);
  dotted_line->SetLineWidth(2);
  dotted_line->Draw();

  canvas3->cd(2);
  canvas3->cd(2)->SetGridx();
  canvas3->cd(2)->SetGridy();
  for(int i = 3; i <= binExp - 2; i++)
  {
    newBinContent = (-hFiltered2->GetBinContent(i+2) 
    + 8*hFiltered2->GetBinContent(i+1) 
    - 8*hFiltered2->GetBinContent(i-1) 
    + hFiltered2->GetBinContent(i-2))/12*binLength;
    hFilDiff5Channel2->SetBinContent(i, newBinContent);
  } 

  hFilDiff5Channel2->SetLineColor(kBlue);
  hFilDiff5Channel2->SetXTitle("Channel");
  hFilDiff5Channel2->Draw();
  // hFiltered2->Scale(0.04);
  hFiltered2->Draw("same");

  // double a = -0.5;
  // double b = 0.5;
  // double c = (a+b)/2;
  // double f_a = TMath::Erfc(a) - theta/tanTheta;
  // double f_b = TMath::Erfc(b) - theta/tanTheta;
  // double f_c = TMath::Erfc(c) - theta/tanTheta;

  // if (f_a*f_b < 0.){
  //   while (abs(f_c) > 0.0001){
  //     if (f_a*f_c < 0){
  //       b = c;
  //       c=(a+b)/2;
  //     }
  //     else if (f_b*f_c < 0){
  //       a = c;
  //       c=(a+b)/2;
  //     }
  //     else if(f_a==0){
  //       c = a;
  //     }
  //     else if(f_b==0){
  //       c = b;
  //     }
  //     f_a = TMath::Erfc(a) - theta/tanTheta;
  //     f_b = TMath::Erfc(b) - theta/tanTheta;
  //     f_c = TMath::Erfc(c) - theta/tanTheta;
  //   }
  // }
  // else if(f_a == 0.){
  //   c = a;
  // }
  // else if(f_b == 0.){
  //   c = b;
  // }
  // else{
  //   std::cout << "f(a) and f(b) don't have opposite signs. \n";
  // }

  // std::cout << "solution for x = " << c << "\n";
  // std::cout << "y(x) = " << TMath::Erfc(c) << "\n";
  
  // TCanvas* canvas4 = new TCanvas("canvas4", "erfc function", 800, 600);
  // TF1* ferfc = new TF1("ferfc", "erfc(x)", -TMath::Pi(), TMath::Pi());
  // TF1* f1 = new TF1("f1", "x/tan(x)", -TMath::Pi(), TMath::Pi());
  // canvas4->cd();
  // ferfc->SetNpx(10000);
  // ferfc->Draw();
  // f1->SetNpx(10000);
  // f1->Draw("same");

  // std::cout << 1000/(sqrt(2)*c) << "\n";

  TCanvas* canvas5 = new TCanvas("canvas5", "Calibrate", 1500, 700);
  canvas5->Divide(3,2);

  canvas5->cd(1);

  //! Compton edge position
  //! 511 keV 
  double min = 0.;
  for(int i = minCompt511; i <= maxCompt511; i++){
    if(hFilDiff5Channel1->GetBinContent(i) < min){
      min = hFilDiff5Channel1->GetBinContent(i);
      id[0] = i;
    }
  }
  std::cout << "\nCompt edge 340 keV\nDiff channel[" << id[0] << "] = " << min << "\n";
  std::cout << "Channel[" << id[0] << "] = " << hChannel1->GetBinContent(id[0]) << "\n";

  //! 1274 keV
  min = 0.;
  for(int i = minCompt1274; i <= maxCompt1274; i++){
    if(hFilDiff5Channel1->GetBinContent(i) < min){
      min = hFilDiff5Channel1->GetBinContent(i);
      id[1] = i;
    }
  }
  std::cout << "\nCompt edge 1061 keV\nDiff channel[" << id[1] << "] = " << min << "\n";
  std::cout << "Channel[" << id[1] << "] = " << hChannel1->GetBinContent(id[1]) << "\n";

  //! 662 keV 
  min = 0.;
  for(int i = minCompt662; i <= maxCompt662; i++){
    if(hFilDiff5Channel2->GetBinContent(i) < min){
      min = hFilDiff5Channel2->GetBinContent(i);
      id[2] = i;
    }
  }
  std::cout << "\nCompt edge 477 keV\nDiff channel[" << id[2] << "] = " << min << "\n";
  std::cout << "Channel[" << id[2] << "] = " << hChannel1->GetBinContent(id[2]) << "\n";

  //! Diff = 0 position
  //! 511 keV
  min = 100.;
  for(int i = min0Der511; i <= max0Der511; i++){
    if(abs(hFilDiff5Channel1->GetBinContent(i)) < min){
      min = abs(hFilDiff5Channel1->GetBinContent(i));
      id2[0] = i;
    }
  }
  std::cout << "\nCompt edge 340 keV\nDiff channel[" << id2[0] << "] = " << hFilDiff5Channel1->GetBinContent(id2[0]) << "\n";
  std::cout << "Channel [" << id2[0] << "] = " << hChannel1->GetBinContent(id2[0]) << "\n"; //!
  
  //! 1274 keV position
  min = 100.;
  for(int i = min0Der1274; i <= max0Der1274; i++){
    if(abs(hFilDiff5Channel1->GetBinContent(i)) < min){
      min = abs(hFilDiff5Channel1->GetBinContent(i));
      id2[1] = i;
    }
  }
  std::cout << "\nCompt edge 1061 keV\nDiff channel[" << id2[1] << "] = " << hFilDiff5Channel1->GetBinContent(id2[1]) << "\n";
  std::cout << "Channel [" << id2[1] << "] = " << hChannel1->GetBinContent(id2[1]) << "\n"; //!

  //! 662 keV position
  min = 100.;
  for(int i = min0Der662; i <= max0Der662; i++){
    if(abs(hFilDiff5Channel2->GetBinContent(i)) < min){
      min = abs(hFilDiff5Channel2->GetBinContent(i));
      id2[2] = i;
    }
  }
  std::cout << "\nCompt edge 477 keV\nDiff channel[" << id2[2] << "] = " << hFilDiff5Channel2->GetBinContent(id2[2]) << "\n";
  std::cout << "Channel [" << id2[2] << "] = " << hChannel2->GetBinContent(id2[2]) << "\n\n"; //!

  //! Energy calibration fit
  TH1D* cal = new TH1D("cal", "test", binExp, xminExp, xmaxExp);

  // cal->SetBinContent(id[0], physE_c[0]);
  cal->SetBinContent(id[1], physE_c[1]);
  cal->SetBinContent(id[2], physE_c[2]);
  cal->Draw();
  std::cout << "\nENERGY CALIBRATION\n";
  cal->Fit("fit", "", "", id[0]+100 -20 , id[1]+100 +20);

  double a_cal = fit->GetParameter(0);
  double b_cal = fit->GetParameter(1);

  std::cout << "\nFIRST CALIBRATION:\n";
  // std::cout << "Calibration(MeV): a = " << a_cal << "; b = " << b_cal << "\n";
  std::cout << "Calibration(keV): a = " << a_cal/1000. << "; b = " << b_cal/1000. << "\n\n";

  // canvas5->cd(2);
  // double xminCal = xminExp*a_cal+b_cal;
  // double xmaxCal = xmaxExp*a_cal+b_cal;
  // std::cout << "xminCal = " << xminCal <<"; xmaxCal = " << xmaxCal << "\n";

  // TH1D* hCal1 = new TH1D("hcal1", "Na22 Calibrated Spectrum", binExp, xminCal, xmaxCal);
  // for (int i = 1; i <= binExp; i++){
  //   hCal1->SetBinContent(i, hChannel1->GetBinContent(i));
  // }
  // hCal1->SetXTitle("Energy(keV)");
  // hCal1->SetYTitle("Count");
  // hCal1->Draw();

  // TMarker* mark;
  // for (int i = 0; i <= 1; i++){
  //   mark = new TMarker((id[i]+100)*a_cal+b_cal, hCal1->GetBinContent(id[i]), 8);
  //   mark->SetMarkerColor(kRed);
  //   mark->Draw();
  // }

  // canvas5->cd(3);
  // TH1D* hCal2 = new TH1D("hcal2", "Cs137 Calibrated Spectrum", binExp, xminCal, xmaxCal);
  // for (int i = 1; i <= binExp; i++){
  //   hCal2->SetBinContent(i, hChannel2->GetBinContent(i));
  // }
  // hCal2->SetXTitle("Energy(keV)");
  // hCal2->SetYTitle("Count");
  // hCal2->Draw();

  // mark = new TMarker((id[2]+100)*a_cal+b_cal, hCal2->GetBinContent(id[2]), 8);
  // mark->SetMarkerColor(kRed);
  // mark->Draw();
  
  // canvas5->cd(4);
  // TH1D* hCalDiff1 = new TH1D("hCalDiff1", "Diff1", binExp, xminCal, xmaxCal);
  // for (int i = 1; i <= binExp; i++){
  //   hCalDiff1->SetBinContent(i, hFilDiff5Channel1->GetBinContent(i));
  // }

  // //! tanTheta
  // double tanTheta[3];
  // // TF1* pol = new TF1("pol", "");
  // // hCalDiff5->Fit("pol", "", "", 170, 200);
  // hCalDiff1->Fit("fit", "", "", minThetafit511, maxThetafit511);
  // tanTheta[0] = fit->GetParameter(0);

  // canvas5->cd(5);
  // TH1D* hCalDiff2 = (TH1D*)hCalDiff1->Clone();
  // hCalDiff2->Fit("fit", "", "", minThetafit1274, maxThetafit1274);
  // tanTheta[1] = fit->GetParameter(0);

  // canvas5->cd(6);
  // TH1D* hCalDiff3 = new TH1D("hCalDiff3", "Diff3", binExp, xminCal, xmaxCal);
  // for (int i = 1; i <= binExp; i++){
  //   hCalDiff3->SetBinContent(i, hFilDiff5Channel2->GetBinContent(i));
  // }
  // hCalDiff3->Fit("fit", "", "", minThetafit662, maxThetafit662);
  // tanTheta[2] = fit->GetParameter(0);

  // // std::cout << "\ntanTheta[0] = " << tanTheta[0] << "; tanTheta[1] = " << tanTheta[1] << "; tanTheta[2] = " << tanTheta[2] << "\n";

  // double A[3], B[3], C[3], D[3], E_c[3];
  // // //! 511
  // A[0] = hCal1->GetBinContent(id2[1]);
  // B[0] = hCal1->GetBinCenter(id2[0]);
  // C[0] = hCal1->GetBinContent(id[1]);
  // D[0] = tanTheta[0];
  // E_c[0] = hCal1->GetBinCenter(id[0]);

  // //! 1274
  // A[1] = hCal1->GetBinContent(id2[1]);
  // B[1] = hCal1->GetBinCenter(id2[1]);
  // C[1] = hCal1->GetBinContent(id[1]);
  // D[1] = tanTheta[1];
  // E_c[1] = hCal1->GetBinCenter(id[1]);
  
  // //! 662
  // A[2] = hCal2->GetBinContent(id2[2]);
  // B[2] = hCal2->GetBinCenter(id2[2]);
  // C[2] = hCal2->GetBinContent(id[2]);
  // D[2] = tanTheta[2];
  // E_c[2] = hCal2->GetBinCenter(id[2]);

  // // std::cout << "\nDiff=0 height A = " << A[1] << "\n";
  // // std::cout << "Diff=0 position B = " << B[1] << " (keV)\n";
  // // std::cout << "Compton edge height C = " << C[1] << "\n";
  // // std::cout << "tanTheta D = " << D[1] << "\n";
  // // std::cout << "Compton edge position E_c = " << E_c[1] << " (keV)\n";

  // double L[3], M[3], N[3];
  // for(int i = 0; i <= 2; i++){
  //   L[i] = D[i]/2;
  //   M[i] = -sqrt(2/TMath::Pi())*D[i]*E_c[i];
  //   N[i] = (D[i]*pow(E_c[i],2) + A[i])/2 - 2*B[i]*D[i]*( E_c[i]/2 + B[i]*D[i]/4 + 1) - C[i];
  // }

  // // std::cout << "\nL = " << L[1] << "; M = " << M[1] << "; N = " << N[1] << "\n";
  // double sigma1[3], sigma2[3];
  // for(int i = 0; i <= 2; i++){
  //   sigma1[i] = (-M[i] + sqrt(M[i]*M[i] - 4.*L[i]*N[i]))/(2.*L[i]);
  //   sigma2[i] = (-M[i] - sqrt(M[i]*M[i] - 4.*L[i]*N[i]))/(2.*L[i]);
  // }
    
  // // std::cout << "\nsigma = " << sigma1[1] << " keV or sigma = " << sigma2[1] << " keV\n";
  // // std::cout << "relative sigma = " << 100*sigma1/E_c << "% or relative sigma = " << 100*sigma2/E_c << "%\n";

  // TCanvas* canvas6 = new TCanvas("canvas6", "test", 1400, 700);
  // canvas6->Divide(2, 1);

  // // canvas6->cd(1);
  // // TH1D* hDiff1 = new TH1D("hDiff1", "compare filter", binExp, xminExp, xmaxExp);
  // // for(int i = 3; i <= binExp - 2; i++)
  // // {
  // //   newBinContent = (-hChannel1->GetBinContent(i+2) 
  // //   + 8*hChannel1->GetBinContent(i+1) 
  // //   - 8*hChannel1->GetBinContent(i-1) 
  // //   + hChannel1->GetBinContent(i-2))/12*binLength;
  // //   hDiff1->SetBinContent(i, newBinContent);
  // // }
  // // hDiff1->Draw();
  // // hDiff1->SetLineColor(kRed + 2);
  // // hFilDiff5Channel1->Draw("same");
  
  // // canvas6->cd(2);
  // // TH1D* hDiff2 = new TH1D("hDiff2", "compare filter", binExp, xminExp, xmaxExp);
  // // for(int i = 3; i <= binExp - 2; i++)
  // // {
  // //   newBinContent = (-hChannel2->GetBinContent(i+2) 
  // //   + 8*hChannel2->GetBinContent(i+1) 
  // //   - 8*hChannel2->GetBinContent(i-1) 
  // //   + hChannel2->GetBinContent(i-2))/12*binLength;
  // //   hDiff2->SetBinContent(i, newBinContent);
  // // }
  // // hDiff2->Draw();
  // // hDiff2->SetLineColor(kRed + 2);
  // // hFilDiff5Channel2->Draw("same");

  // double a_res[3], b_res[3], c_res[3];
  // for(int i = 0; i <= 2; i++){
  //   a_res[i] = D[i];
  //   b_res[i] = -2*B[i]*D[i];
  //   c_res[i] = A[i] + pow(B[i],2)*D[i];    
  // }

  // // std::cout << "a_res = " << a_res << "; b_res = " << b_res << "; c_res = " << c_res << "\n";
  
  // TF1* newfit1 = new TF1("newfit1", "[5] + [6]*x + [7]*x*x + (1./2.)*([0]*(x*x + [4]*[4]) + [1]*x + [2])*erfc((x - [3])/(sqrt(2.)*[4])) + ((-[4]/sqrt(2.*TMath::Pi()))*[0]*(x + [3]) + [1])*exp(-pow((x - [3]),2)/(2.*[4]*[4]))", xminCal, xmaxCal);
  // newfit1->SetParameter(0, a_res[1]);
  // newfit1->SetParameter(1, b_res[1]);
  // newfit1->SetParameter(2, c_res[1]);
  // newfit1->SetParameter(3, E_c[1]);
  // newfit1->SetParameter(4, E_c[1]*relRes1274);
  // // newfit1->SetParameter(5, 200.);
  // newfit1->SetNpx(10000);

  // TF1* newfit2 = new TF1("newfit2", "[5] + [6]*x + [7]*x*x + (1./2.)*([0]*(x*x + [4]*[4]) + [1]*x + [2])*erfc((x - [3])/(sqrt(2.)*[4])) + ((-[4]/sqrt(2.*TMath::Pi()))*[0]*(x + [3]) + [1])*exp(-pow((x - [3]),2)/(2.*[4]*[4]))", xminCal, xmaxCal);
  // newfit2->SetParameter(0, a_res[2]);
  // newfit2->SetParameter(1, b_res[2]);
  // newfit2->SetParameter(2, c_res[2]);
  // newfit2->SetParameter(3, E_c[2]);
  // newfit2->SetParameter(4, E_c[2]*relRes662);
  // // newfit2->SetParameter(5, 200.);
  // newfit2->SetNpx(10000);

  // canvas6->cd();
  // TH1D* hFit1 = new TH1D("hFit1", "hFit1", binExp, xminCal, xmaxCal);
  // for(int i = 1; i <= binExp; i++){
  //   hFit1->SetBinContent(i, hFiltered1->GetBinContent(i));
  // }
  // hFit1->SetXTitle("Energy(keV)");
  // hFit1->SetYTitle("Count");

  // TH1D* hFit2 = new TH1D("hFit2", "hFit2", binExp, xminCal, xmaxCal);
  // for(int i = 1; i <= binExp; i++){
  //   hFit2->SetBinContent(i, hFiltered2->GetBinContent(i));
  // }
  // hFit2->SetXTitle("Energy(keV)");
  // hFit2->SetYTitle("Count");

  // canvas6->cd(1);
  // std::cout << "\nNa22 R FUNCTION FITTING\n";
  // hFit1->Fit("newfit1", "", "", minEFit1274, maxEFit1274);
  // // hCal1->Fit("newfit1", "", "", minEFit1274, maxEFit1274);

  // canvas6->cd(2);
  // std::cout << "\nCs137 R FUNCTION FITTING\n";
  // hFit2->Fit("newfit2", "", "", minEFit662, maxEFit662);
  // // hCal2->Fit("newfit2", "", "", minEFit662, maxEFit662);
  // // newfit->Draw();

  // std::cout << "\nInitial parameters\na = " << a_res[1] << "; " << a_res[2]
  // << "\nb = " << b_res[1] << "; " << b_res[2] 
  // << "\nc = " << c_res[1] << "; " << c_res[2]
  // << "\nE_c = " << E_c[1] << "; " << E_c[2];

  TCanvas* canvas7 = new TCanvas("canvas7", "Fit", 1200, 600);
  canvas7->Divide(2,2);

  // TH1D* hFitSigma = new TH1D("hFitSigma", "Fit Sigma", binExp, xminCal/1000., xmaxCal/1000.);

  // double newE_c[3], newSigma[3];
  // newE_c[1] = newfit1->GetParameter(3);
  // newE_c[2] = newfit2->GetParameter(3);
  // newSigma[1] = newfit1->GetParameter(4);
  // newSigma[2] = newfit2->GetParameter(4);

  // // newE_c[0] = 340.;
  // // newSigma[0] = 38.5;

  // TF1* fRes = new TF1("fRes", "sqrt([0]*[0]/x + [1]*[1]/(x*x))");
  // // TF1* fRes = new TF1("fRes", "sqrt([0]*[0] + [1]*[1]/(x*x))");
  // // hFitSigma->SetBinContent(hCal1->FindBin(newE_c[0]), pow(newSigma[0],2) );
  // hFitSigma->SetBinContent(hCal1->FindBin(newE_c[1]), newSigma[1]/newE_c[1] );
  // hFitSigma->SetBinContent(hCal2->FindBin(newE_c[2]), newSigma[2]/newE_c[2] );
  // canvas7->cd(1);
  // hFitSigma->Fit("fRes", "", "", 0.4, 1.1);
  // hFitSigma->Draw();

  // std::cout << "\nFIRST RESOLUTION PARAMETERS:\n"; 
  // std::cout << "RESco_B = " << fRes->GetParameter(0) << "; RESco_C = " << fRes->GetParameter(1) << "\n\n";
  // // std::cout << "RESco_A = " << fRes->GetParameter(0) << "; RESco_C = " << fRes->GetParameter(1) << "\n";

  // TH1D* hFitnewE_c = new TH1D("hFitnewE_c", "Fit new E_c", binExp, xminCal, xmaxCal);
  // hFitnewE_c->SetBinContent(hCal1->FindBin(newE_c[1]), physE_c[1]);
  // hFitnewE_c->SetBinContent(hCal2->FindBin(newE_c[2]), physE_c[2]);
  // canvas7->cd(2);
  // hFitnewE_c->Fit("fit", "", "", 400, 1100);
  // hFitnewE_c->Draw();

  // double newa_cal = a_cal*fit->GetParameter(0), newb_cal = b_cal*fit->GetParameter(0) + fit->GetParameter(1);

  // std::cout << "(MeV)new a_cal = " << newa_cal << "; new b_cal = " << newb_cal << "\n";
  // std::cout << "\n(keV)new a_cal = " << newa_cal/1000. << "; new b_cal = " << newb_cal/1000. << "\n";
  
  // std::cout << "\n(keV)new a_cal = " << func(a_cal, b_cal, 0)/1000. << "; new b_cal = " << func(a_cal, b_cal, 1)/1000. << "\n";
  func(a_cal, b_cal, 0);

  double a_change, b_change, coB_change, coC_change;

  for (int i = 1; i <= 1; i++){
    func(a_cal, b_cal, 4);

    a_change = func(a_cal, b_cal, 0);
    b_change = func(a_cal, b_cal, 1);
    coB_change = func(a_cal, b_cal, 2);
    coC_change = func(a_cal, b_cal, 3);

    a_cal = a_change;
    b_cal = b_change;
  }

  std::cout << "final a_cal = " << a_change/1000 << "; final b_cal = " << b_change/1000 << "\n";
  std::cout << "final coB = " << coB_change << "; final coC = " << coC_change << "\n";
  // TH1D* hhCal1 = new TH1D("hhcal1", "Na22 Calibrated Spectrum", binExp, newa_cal*xminExp + newb_cal, newa_cal*xmaxExp + newb_cal);
  // for(int i = 1; i<=binExp; i++){
  //   hhCal1->SetBinContent(i, hChannel1->GetBinContent(i));
  // }
  // canvas7->cd(3);
  // hhCal1->Draw();

  // TH1D* hhCal2 = new TH1D("hhcal2", "Cs137 Calibrated Spectrum", binExp, newa_cal*xminExp + newb_cal, newa_cal*xmaxExp + newb_cal);
  // for(int i = 1; i<=binExp; i++){
  //   hhCal2->SetBinContent(i, hChannel2->GetBinContent(i));
  // }
  // canvas7->cd(4);
  // hhCal2->Draw();

  std::cout << "\ntime: " << timer->RealTime() << " seconds\n";
}

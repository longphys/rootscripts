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
#include "TChain.h"

//! histogram options
double xminSim = 0.;
double xmaxSim = 1.3;
double xminExp = 100;
double xmaxExp = 1500;
int binSim = (xmaxSim-xminSim)*1000.;
int binExp = xmaxExp-xminExp;

double E1 = 0.477;
double E2 = 1.061;
double E3 = 1.117;

//! Zoomed in Histograms options
double minTestCs137 = 0.3; //MeV
double maxTestCs137 = 0.7;
double minTestNa22 = 1.0; //MeV
double maxTestNa22 = 1.5;
double minTestCo60 = 0.9; //MeV
double maxTestCo60 = 1.5;

//! ENERGY CALIBRATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// double Ch1 = 315;
// double Ch1 = 410;
// double Ch2 = 800;
//? Now determined by FFT
int minRangeNa22 = 300;
int minRangeCo60 = 300;

//! RESOLUTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double dSig1 = 8.;
double dSig2 = 5.5;
double dSig3 = 5.;
//? In percentage

//! SCALING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ScaleCs137 = 8.;
double ScaleNa22 = 2.2;
double ScaleCo60 = 3.;

//! CHANNELS
// int channel = 0;
// int channel = 1;
// int channel = 2;

//! FFT Threshold
double threshExp = 100.;
double threshSim = 100.;

//! Tuning Rate
double tuningRateCh1 = 2.;
double tuningRateCh2 = 2.;
double tuningRateCh3 = 2.;
double tuningRateSig1 = 0.05;
double tuningRateSig2 = 0.05;
double tuningRateSig3 = 0.05;

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

TH1D* Diff(TH1D* hChannel, std::string TH1name, std::string TH1title, int bin, double xmin, double xmax)
{
  TH1namechar = TH1name.c_str();
  TH1titlechar = TH1title.c_str();
  double newBinContent;
  double binLength = (xmax-xmin)/bin;

  TH1D* hDiff = new TH1D(TH1namechar, TH1titlechar, bin, xmin, xmax);
  for(int i = 3; i <= binExp - 2; i++)
  {
    newBinContent = (-hChannel->GetBinContent(i+2) 
    + 8*hChannel->GetBinContent(i+1) 
    - 8*hChannel->GetBinContent(i-1) 
    + hChannel->GetBinContent(i-2))/12*binLength;
    hDiff->SetBinContent(i, newBinContent);
  }
  return hDiff;
}

void gainCo60(){
  TCanvas* c1 = new TCanvas("c1", "c1", 1500, 1000);
  c1->Divide(3,1);

  TCanvas* c2 = new TCanvas("c2", "c2", 1500, 1000);
  c2->Divide(4,2);

  const char* number_char;
  const char* hist = "Histogram ";
  const char* graph = "Channel ";
  const char* name = "./expfiles/new/stilbene_dividers_Plastic_ch2_Co60_no_coinc_ch0HV2100_ch1HV2300_ch2HV1900_ch0_9751_ch1_9800_ch2_9421_run_0_000";
  const char* extension = ".root";

  int maxfilenum = 5;

  TGraph* gr[3];
  for (int channel = 0; channel <= 2; channel++){

    number_char = std::to_string(channel).c_str();
    char* graphname = new char[strlen(graph)+1];
    strcpy(graphname, graph);
    strcat(graphname, number_char);

    gr[channel] = new TGraph();
    gr[channel]->SetTitle(graphname);
    gr[channel]->GetXaxis()->SetTitle("File number");
    gr[channel]->GetYaxis()->SetTitle("Compton edge channel");

    gr[channel]->SetMarkerStyle(kFullCircle);
    gr[channel]->SetMarkerSize(1);
    for (int i = 1; i <= maxfilenum; i++){
      number_char = std::to_string(i).c_str();

      char* fullname = new char[strlen(name)+1+strlen(extension)];
      strcpy(fullname, name); /* copy name into the new var */
      strcat(fullname, number_char); /* add the extension */
      strcat(fullname, extension); /* add the extension */

      TFile* fexpCo60 = new TFile(fullname, "read");
      TTree* texpCo60 =  (TTree*) fexpCo60->Get("AnalysisxTree");

      UShort_t ch[48];
      texpCo60->SetBranchAddress("NeEvent.neutAmp[48]", ch);
      Long64_t entriesExpCo60 = texpCo60->GetEntries();

      // std::cout << "entries = " << entriesExpCo60 << "\n";

      TH1D* hExpCo60 = new TH1D("hExpCo60", "", binExp, xminExp, xmaxExp);
      char* histname = new char[strlen(hist)+1];
      strcpy(histname, hist);
      strcat(histname, number_char);

      hExpCo60->SetTitle(histname);

      for(int i = 0; i < entriesExpCo60; i++){
        texpCo60->GetEntry(i);
        hExpCo60->Fill(ch[channel]);
      }

      TH1D* hExpCo60Filtered = fft(hExpCo60, 0.1, 75, "hExpCo60Filtered", "Co60 Spectrum Filtered", binExp, xminExp, xmaxExp);
      TH1D* hExpCo60Diff = Diff(hExpCo60Filtered, "hExpCo60Diff", "Co60 Spectrum Derivative", binExp, xminExp, xmaxExp);

      hExpCo60Diff->GetXaxis()->SetRangeUser(minRangeCo60, xmaxExp);
      int edge = hExpCo60Diff->GetMinimumBin() + xminExp;
      hExpCo60Diff->GetXaxis()->UnZoom();

      gr[channel]->AddPoint(i, edge);
      c1->cd(channel+1);
      gr[channel]->Draw("ALP");
      c1->Modified();
      c1->Update();

      c2->cd(i);
      hExpCo60Diff->Draw();
      c2->Modified();
      c2->Update();
    }
  }
}
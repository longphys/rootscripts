#ifndef DATAANALYSER_HH
#define DATAANALYSER_HH

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TVector3.h"
#include "TVirtualFFT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TLegend.h"

#include "ArgumentParser.hh"

class DataAnalyser{
public:
  TH1D* fft(TH1D* h_channel, double para_k, double para_c, std::string name_h_channel, std::string title_h_channel);
  TH1D* Diff(TH1D* h_channel, std::string name_h_channel, std::string title_h_channel);
  void Analyze();
};
#endif
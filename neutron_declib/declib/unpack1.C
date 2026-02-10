#include "DecLibTest.hh"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <sstream>
#include <fstream>
#include <iostream>


int startingChannelNum = 0;
int endingChannelNum = 31;
int nChannels = endingChannelNum - startingChannelNum + 1;

int64_t EventNr = -1;
vector <double> Area;
vector <double> Time;
vector <double> Width;
vector <double> Height;
vector <unsigned int> Chan;

// void unpack1(string fileIn="../test06_2_Plast_BC404_14x_14y.dec")
// void unpack1(string fileIn="../test07_2_Plast_BC404_14x_14y.dec")
void unpack1(string fileIn="/home/long/data/neutron_eff_plastic_ing27/test08_2_Plast_BC404_14x_14y_Pos_x0_y0.dec")
// void unpack1(string fileIn="../test09_2_Plast_BC404_14x_14y_Pos_x20_y-20.dec")
// void unpack1(string fileIn="../test10_2_Plast_BC404_14x_14y_Pos_x0_y-92.dec")
// void unpack1(string fileIn="../test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1.dec")
// void unpack1(string fileIn="../test12_2_Plast_BC404_Cs137.dec")
// void unpack1(string fileIn="../test13_2_Plast_BC404_Na22.dec")
// void unpack1(string fileIn="../test14_2_Plast_BC404_Background.dec")
{
  auto timer = new TStopwatch();
  timer->Start();

	ROOT::EnableThreadSafety();
	// TH1::AddDirectory(false);
  
  DecManager &dm=DecManager::Instance();

  dm.SetNBuffers(8);
  dm.SetNThreads(8);
  dm.SetBufferSize(1e7);
  dm.AddFile(fileIn);
  Event ev;
  dm.Start();

  // TFile* fileOut = new TFile("./output/test06_2_Plast_BC404_14x_14y_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test07_2_Plast_BC404_14x_14y_out.root","recreate");
  TFile* fileOut = new TFile("/home/long/data/neutron_eff_plastic_ing27/output/test08_2_Plast_BC404_14x_14y_Pos_x0_y0_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test09_2_Plast_BC404_14x_14y_Pos_x20_y-20_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test10_2_Plast_BC404_14x_14y_Pos_x0_y-92_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test12_2_Plast_BC404_Cs137.root","recreate");
  // TFile* fileOut = new TFile("./output/test13_2_Plast_BC404_Na22.root","recreate");
  // TFile* fileOut = new TFile("./output/test14_2_Plast_BC404_Background.root","recreate");

  TTree* tree = new TTree("tree", "tree");
  tree->Branch("EventID", &EventNr, "EventID/I");
  tree->Branch("Area", &Area);
  tree->Branch("Time", &Time);
  tree->Branch("Width", &Width);
  tree->Branch("Height", &Height);
  tree->Branch("Channel", &Chan);

  while(dm.GetNextEvent(&ev,0))
  // for (int counter = 0; counter < (int)100000000; counter++)
  {
    Area.clear();
    Time.clear();
    Chan.clear();
    Width.clear();
    Height.clear();

    Area.resize(nChannels,0.);
    Time.resize(nChannels,0.);
    Chan.resize(nChannels,0.);
    Width.resize(nChannels,0.);
    Height.resize(nChannels,0.);

    for(int i = startingChannelNum; i <= endingChannelNum; i++){
      Pulse* pulse=ev.GetPulse(i);
      if(pulse){
        Area[i-startingChannelNum] = pulse->Area;
        Time[i-startingChannelNum] = pulse->Time;
        Chan[i-startingChannelNum] = pulse->Chan;
        Width[i-startingChannelNum] = pulse->Width;
        Height[i-startingChannelNum] = pulse->Height;

        // std::cout << "Channel=" << pulse->Chan << "\n";
      }
    }
    tree->Fill();
    EventNr++;
    if (EventNr%100000 == 0)
    {
      std::cout << "EventNr = " << EventNr << std::endl;
    }
  }
  tree->Write();
  fileOut->Close();

  std::cout << "\ntime: " << timer->RealTime() << " (s)\n";
}

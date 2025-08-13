#include "DecLibTest.hh"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <sstream>
#include <fstream>
#include <iostream>


int startingChannelNum = 0;
int endingChannelNum = 29;

int64_t EventNr = -1;
vector <double> dE;
vector <double> dTime;
vector <unsigned int> channel;

// void unpack1(string fileIn="../test06_2_Plast_BC404_14x_14y.dec")
// void unpack1(string fileIn="../test07_2_Plast_BC404_14x_14y.dec")
// void unpack1(string fileIn="../test08_2_Plast_BC404_14x_14y_Pos_x0_y0.dec")
// void unpack1(string fileIn="../test09_2_Plast_BC404_14x_14y_Pos_x20_y-20.dec")
// void unpack1(string fileIn="../test10_2_Plast_BC404_14x_14y_Pos_x0_y-92.dec")
// void unpack1(string fileIn="../test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1.dec")
// void unpack1(string fileIn="../test12_2_Plast_BC404_Cs137.dec")
// void unpack1(string fileIn="../test13_2_Plast_BC404_Na22.dec")
void unpack1(string fileIn="../test14_2_Plast_BC404_Background.dec")
{
  auto timer = new TStopwatch();
  timer->Start();


	ROOT::EnableThreadSafety();
	TH1::AddDirectory(false);

  // TFile* fileOut = new TFile("./output/test06_2_Plast_BC404_14x_14y_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test07_2_Plast_BC404_14x_14y_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test08_2_Plast_BC404_14x_14y_Pos_x0_y0_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test09_2_Plast_BC404_14x_14y_Pos_x20_y-20_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test10_2_Plast_BC404_14x_14y_Pos_x0_y-92_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test11_2_Plast_BC404_14x_14y_Pos_x0_y-92_wo_upper_detector-ch1_out.root","recreate");
  // TFile* fileOut = new TFile("./output/test12_2_Plast_BC404_Cs137.root","recreate");
  // TFile* fileOut = new TFile("./output/test13_2_Plast_BC404_Na22.root","recreate");
  TFile* fileOut = new TFile("./output/test14_2_Plast_BC404_Background.root","recreate");


  TTree* tree = new TTree("ETree", "ETree");
  tree->Branch("EventID", &EventNr, "EventID/I");
  tree->Branch("EDep", &dE);
  tree->Branch("Time", &dTime);
  tree->Branch("Channel", &channel);
  
  DecManager &dm=DecManager::Instance();

  dm.SetNBuffers(2);//количество буферов
  dm.SetNThreads(1);//количество потоков. В простейшем случае, 1
  dm.SetBufferSize(1e7);//размер буфера с событиями
  dm.AddFile(fileIn);
  Event ev;
  dm.Start();

  while(dm.GetNextEvent(&ev,0))
  // for (int counter = 0; counter < (int)100000000; counter++)
  {
    // dm.GetNextEvent(&ev,0);
    dE.clear();
    dTime.clear();
    channel.clear();

    dE.resize(endingChannelNum-startingChannelNum+1,0.);
    dTime.resize(endingChannelNum-startingChannelNum+1);
    channel.resize(endingChannelNum-startingChannelNum+1);

    for(int i = startingChannelNum; i <= endingChannelNum; i++){
      Pulse* pulse=ev.GetPulse(i);//импульс с hpge детектора
      if(pulse){
        dE[i] = pulse->Area;
        dTime[i] = pulse->Time;
        channel[i] = pulse->Chan;
        // std::cout << "Channel=" << pulse->Chan << "\n";
        // dE.push_back(pulse->Area);
        // dTime.push_back(pulse->Time);
        // channel.push_back(pulse->Chan);
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

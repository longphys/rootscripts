#include "DecLibTest.hh"
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <sstream>
#include <fstream>
#include <iostream>

int64_t EventNr = -1;
vector <double> dE;
vector <double> dTime;
vector <unsigned int> channel;

int startingChannelNum = 0;
int endingChannelNum = 1;

// void unpack0(string fileIn="../test01_2_Plast_BC404_Co60.dec")
// void unpack0(string fileIn="../test02_2_Plast_BC404_Mn54.dec")
void unpack0(string fileIn="../test04_2_Plast_BC404_Bskgr.dec")
{
	ROOT::EnableThreadSafety();
	TH1::AddDirectory(false);

  // TFile* fileOut = new TFile("out_test01_2_Plast_BC404_Co60.root","recreate");
  // TFile* fileOut = new TFile("out_test02_2_Plast_BC404_Mn54.root","recreate");
  TFile* fileOut = new TFile("out_test04_2_Plast_BC404_Bskgr.root","recreate");
  TTree* tree = new TTree("ETree", "ETree");
  tree->Branch("EventID", &EventNr, "EventID/I");
  tree->Branch("EDep", &dE);
  tree->Branch("Time", &dTime);
  tree->Branch("Channel", &channel);
  
  DecManager &dm=DecManager::Instance();


  dm.SetNBuffers(1);//количество буферов
  dm.SetNThreads(1);//количество потоков. В простейшем случае, 1
  dm.SetBufferSize(1e7);//размер буфера с событиями
  dm.AddFile(fileIn);
  Event ev;
  dm.Start();

  while(dm.GetNextEvent(&ev,0))
  {
    dE.clear();
    dTime.clear();
    channel.clear();

    for(int i = startingChannelNum; i <= endingChannelNum; i++){
      Pulse* pulse=ev.GetPulse(i);//импульс с hpge детектора
      if(pulse){
        dE.push_back(pulse->Area);
        dTime.push_back(pulse->Time);
        channel.push_back(pulse->Chan);
      }
    }
    EventNr++;
    tree->Fill();
  }
  tree->Write();
  fileOut->Close();
}

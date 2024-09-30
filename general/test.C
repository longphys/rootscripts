#include "TRandom3.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"

#include <iostream>

void test()
{
  std::vector <double> num;
  int size = 10;
  num.resize(size, 0.);

  std::cout << "First loop\n";
  for(int i = 0; i < num.size(); i++){
    std::cout << "num[" << i << "] = " << num[i] << "\n";
  }
  std::cout << "\n";

  num[num.size()-1] = 10.;

  std::cout << "Second loop\n";
  for(int i = 0; i < num.size(); i++){
    std::cout << "num[" << i << "] = " << num[i] << "\n";
  }
  std::cout << "\n";

  num.clear();
  // size += 2;
  //num.resize(0);
  //num.shrink_to_fit();
  //num.push_back(42);
  //num.push_back(-12);

  //num.erase(num.begin(), num.end());

  std::cout << "Final loop\n";
  for(int i = 0; i < num.size(); i++){
  //for(int i = 0; i < 10; i++){
    std::cout << "num[" << i << "] = " << num[i] << "\n";
  }
  std::cout << "\n";
}

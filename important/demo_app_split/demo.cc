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
#include "MacroReader.hh"

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <getopt.h>
#include <cstdlib>  // For atof, atoi
#include <cstring>  // For strcmp

bool CheckRootFile(const char* filename) {
    TFile file(filename);
    return !file.IsZombie();  // Check if the file exists
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "Error: Expected exactly one argument (macro file) or --help.\n";
        return 1;
    }

    // If the user asks for help, display usage info and exit
    if (std::string(argv[1]) == "--help") {
        std::cout << "Usage: ./demo macro.mac or ./demo --help\n"
                  << "  macro.mac: A script file containing arguments\n"
                  << "  --help   : Show this help message\n";
        return 0;
    }

    // Otherwise, assume the argument is a macro file and process it
    MacroReader aMacroReader;
    std::vector<std::string> args = aMacroReader.Read(argv[1]);

    // Convert vector of strings to array of C-style strings
    std::vector<char*> c_args;
    c_args.push_back(argv[0]); // Program name
    for (std::string& arg : args) {
        c_args.push_back(&arg[0]);
    }
    c_args.push_back(nullptr); // Null termination

    int new_argc = c_args.size() - 1;
    char** new_argv = c_args.data();

    ArgumentParser anArgumentParser;
    if (anArgumentParser.Parse(new_argc, new_argv) != 0) return 1;

    // Initialize application
    TApplication app("app", &argc, argv);

	std::cerr << "\nMain is running...\n";

	// Get Files and Trees
	GetTrees();
	
    // Create canvas and Handle window close event
    //TCanvas* c = new TCanvas("c", "Plot", 0, 0, 800, 600);

    //TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
    //rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	
    // Run the application
    std::cerr << "\nRun success.\n";
    app.Run();

    return 0;
}


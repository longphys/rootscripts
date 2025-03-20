#ifndef ARGUMENTPARSER_HH
#define ARGUMENTPARSER_HH

#include "Config.hh"

#include "TFile.h"

#include <iostream>
#include <getopt.h>
#include <cstdlib>  // For atof, atoi

class ArgumentParser{
public:

  bool CheckRootFile(const char* filename) {
    TFile file(filename);
    return !file.IsZombie();  // Check if the file exists
  }

  enum OptionIDs {
    opt_name_f_sim_1 = 1000,
    opt_name_f_sim_2,
    opt_ascii,
    opt_name_f_mea_1,
    opt_name_f_mea_2,
    opt_entries_sim_fft,
    opt_entries_mea_fft,
    opt_entries_sim_descent,
    opt_entries_mea_descent,
    opt_x_min_sim,
    opt_x_max_sim,
    opt_x_min_mea,
    opt_x_max_mea,
    opt_energy_1,
    opt_energy_2,
    opt_res_1,
    opt_res_2,
    opt_x_min_fft_1,
    opt_x_max_fft_1,
    opt_x_min_fft_2,
    opt_x_max_fft_2,
    opt_x_min_descent_1,
    opt_x_max_descent_1,
    opt_x_min_descent_2,
    opt_x_max_descent_2,
    opt_rate_fft_ini_1,
    opt_rate_fft_ini_2,
    opt_thresh_fft_ini_1,
    opt_thresh_fft_ini_2,
    opt_rate_fft_descent_1,
    opt_rate_fft_descent_2,
    opt_thresh_fft_descent_1,
    opt_thresh_fft_descent_2,
    opt_learning_rate_channel_1,
    opt_learning_rate_channel_2,
    opt_learning_rate_res_1,
    opt_learning_rate_res_2,
    opt_scale_sim_1,
    opt_scale_sim_2,
    opt_scale_mea_1,
    opt_scale_mea_2,
    opt_step_channel_1,
    opt_step_channel_2,
    opt_step_res_1,
    opt_step_res_2
  };

  template <typename T>
  void WriteOutput(const char* optarg, T& variable, T defaultValue, const std::string& name, const std::string& unit = ""){
    // if (optarg == nullptr || optarg[0] == '-') {
    if (optarg == nullptr) {
      std::cerr << "Default " << name << " used: " << defaultValue << " " << unit << "\n";
      variable = defaultValue;
    } else {
      variable = atof(optarg);
      std::cerr << name << " used: " << variable << " " << unit << "\n";
    }
  } 
  
  int Parse(int argc, char** argv, Config& config);
};
#endif
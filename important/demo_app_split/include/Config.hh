#ifndef CONFIG_HH
#define CONFIG_HH

#include <iostream>
#include <string>

struct Config{
  // Choice of using an ascii measurement file
  int ascii = 0;

  // Name of input file
  std::string name_f_sim_1 = "";
  std::string name_f_sim_2 = "";
  std::string name_f_mea_1 = "";
  std::string name_f_mea_2 = "";

  // Default number of entries for FFT
  int entries_sim_fft = 100000;
  int entries_mea_fft = 100000;

  // Default number of entries for gradient descent
  int entries_sim_descent = 500000;
  int entries_mea_descent = 500000;

  // Default options for all histograms
  double x_min_sim = 0.;
  double x_max_sim = 1.3; // in MeV
  int bin_sim = (x_max_sim - x_min_sim)*1000.;

  double x_min_mea = 0.;
  double x_max_mea = 1000.; // in channels
  int bin_mea = (x_max_mea - x_min_mea);

  // Default energies of Compton edges in MeV (Cs_137 and Na_22)
  double energy_1 = 0.477;
  double energy_2 = 1.061;

  // Default first value for energy resolution in percentage
  double res_1 = 8.;
  double res_2 = 5.5;

  // Default channel window to perform FFT
  double x_min_fft_1 = 200;
  double x_max_fft_1 = x_max_mea;
  double x_min_fft_2 = 400;
  double x_max_fft_2 = x_max_mea;

  // Default energy window to perform gradient descent
  double x_min_descent_1 = 0.35;
  double x_max_descent_1 = 0.7;
  double x_min_descent_2 = 0.8;
  double x_max_descent_2 = 1.3;

  // Default FFT cut-off growth rate for estimating Compton edge
  double rate_fft_ini_1 = 0.1;
  double rate_fft_ini_2 = 0.1;

  // Default FFT threshold for estimating Compton edge
  double thresh_fft_ini_1 = 50.;
  double thresh_fft_ini_2 = 50.;

  // Default FFT cut-off growth rate for performing gradient descent
  double rate_fft_descent_1 = 0.1;
  double rate_fft_descent_2 = 0.1;

  // Default FFT threshold for performing gradient descent
  double thresh_fft_descent_sim = 800.;
  double thresh_fft_descent_mea = 800.;

  // Default learning rate for gradient descent
  double learning_rate_channel_1 = 2.;
  double learning_rate_channel_2 = 2.; 
  double learning_rate_res_1 = 0.05;
  double learning_rate_res_2 = 0.05;

  // Default step for calculating gradient
  double step_channel_1 = 1.;
  double step_channel_2 = 1.;
  double step_res_1;
  double step_res_2;

  // Default scaling factor of histograms
  double scale_sim_1 = 1.;
  double scale_sim_2 = 1.;
  double scale_mea_1 = 1.;
  double scale_mea_2 = 1.;
};
#endif
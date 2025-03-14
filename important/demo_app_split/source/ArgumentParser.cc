#include "ArgumentParser.hh"
ArgumentParser::ArgumentParser() {
  // Default values
  name_f_sim_1 = "";
  name_f_sim_2 = "";
  name_f_mea_1 = "";
  name_f_mea_2 = "";

  channel = 0;
  entries_sim_fft = 100000;
  entries_mea_fft = 100000;
  entries_sim_descent = 500000;
  entries_mea_descent = 500000;

  x_min_sim = 0.;
  x_max_sim = 1.3;
  bin_sim = (x_max_sim - x_min_sim) * 1000.;

  x_min_mea = 0.;
  x_max_mea = 1000.;
  bin_mea = (x_max_mea - x_min_mea);

  energy_1 = 0.477;
  energy_2 = 1.061;

  res_1 = 8.;
  res_2 = 5.5;

  x_min_fft_1 = 200;
  x_max_fft_1 = x_max_mea;
  x_min_fft_2 = 400;
  x_max_fft_2 = x_max_mea;

  x_min_descent_1 = 0.35;
  x_max_descent_1 = 0.7;
  x_min_descent_2 = 0.8;
  x_max_descent_2 = 1.3;

  rate_fft_ini_1 = 0.1;
  rate_fft_ini_2 = 0.1;

  thresh_fft_ini_1 = 50.;
  thresh_fft_ini_2 = 50.;

  rate_fft_descent_1 = 0.1;
  rate_fft_descent_2 = 0.1;

  thresh_fft_descent_sim = 800.;
  thresh_fft_descent_mea = 800.;

  learning_rate_channel_1 = 2.;
  learning_rate_channel_2 = 2.;
  learning_rate_res_1 = 0.05;
  learning_rate_res_2 = 0.05;

  step_channel_1 = 1.;
  step_channel_2 = 1.;
  step_res_1 = 0.;
  step_res_2 = 0.;

  scale_sim_1 = 1.;
  scale_sim_2 = 1.;
  scale_mea_1 = 1.;
  scale_mea_2 = 1.;
}


  // Name of input file
  std::string name_f_sim_1 = "";
  std::string name_f_sim_2 = "";
  std::string name_f_mea_1 = "";
  std::string name_f_mea_2 = "";

  // Default channel of the measurement file
  int channel = 0;

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


int ArgumentParser::Parse(int argc, char** argv){
  static struct option longOptions[] = {
    {"name_f_sim_1", required_argument, 0, opt_name_f_sim_1},
    {"name_f_sim_2", required_argument, 0, opt_name_f_sim_2},
    {"name_f_mea_1", required_argument, 0, opt_name_f_mea_1},
    {"name_f_mea_2", required_argument, 0, opt_name_f_mea_2},
    {"channel", required_argument, 0, opt_channel},
    {"entries_sim_fft", required_argument, 0, opt_entries_sim_fft},
    {"entries_mea_fft", required_argument, 0, opt_entries_mea_fft},
    {"entries_sim_descent", required_argument, 0, opt_entries_sim_descent},
    {"entries_mea_descent", required_argument, 0, opt_entries_mea_descent},
    {"x_min_sim", required_argument, 0, opt_x_min_sim},
    {"x_max_sim", required_argument, 0, opt_x_max_sim},
    {"x_min_mea", required_argument, 0, opt_x_min_mea},
    {"x_max_mea", required_argument, 0, opt_x_max_mea},
    {"energy_1", required_argument, 0, opt_energy_1},
    {"energy_2", required_argument, 0, opt_energy_2},
    {"res_1", required_argument, 0, opt_res_1},
    {"res_2", required_argument, 0, opt_res_2},
    {"x_min_fft_1", required_argument, 0, opt_x_min_fft_1},
    {"x_max_fft_1", required_argument, 0, opt_x_max_fft_1},
    {"x_min_fft_2", required_argument, 0, opt_x_min_fft_2},
    {"x_max_fft_2", required_argument, 0, opt_x_max_fft_2},
    {"x_min_descent_1", required_argument, 0, opt_x_min_descent_1},
    {"x_max_descent_1", required_argument, 0, opt_x_max_descent_1},
    {"x_min_descent_2", required_argument, 0, opt_x_min_descent_2},
    {"x_max_descent_2", required_argument, 0, opt_x_max_descent_2},
    {"rate_fft_ini_1", required_argument, 0, opt_rate_fft_ini_1},
    {"rate_fft_ini_2", required_argument, 0, opt_rate_fft_ini_2},
    {"thresh_fft_ini_1", required_argument, 0, opt_thresh_fft_ini_1},
    {"thresh_fft_ini_2", required_argument, 0, opt_thresh_fft_ini_2},
    {"rate_fft_descent_1", required_argument, 0, opt_rate_fft_descent_1},
    {"rate_fft_descent_2", required_argument, 0, opt_rate_fft_descent_2},
    {"thresh_fft_descent_sim", required_argument, 0, opt_thresh_fft_descent_1},
    {"thresh_fft_descent_mea", required_argument, 0, opt_thresh_fft_descent_2},
    {"learning_rate_channel_1", required_argument, 0, opt_learning_rate_channel_1},
    {"learning_rate_channel_2", required_argument, 0, opt_learning_rate_channel_2},
    {"learning_rate_res_1", required_argument, 0, opt_learning_rate_res_1},
    {"learning_rate_res_2", required_argument, 0, opt_learning_rate_res_2},
    {"step_channel_1", required_argument, 0, opt_step_channel_1},
    {"step_channel_2", required_argument, 0, opt_step_channel_2},
    {"step_res_1", required_argument, 0, opt_step_res_1},
    {"step_res_2", required_argument, 0, opt_step_res_2},
    {"scale_sim_1", required_argument, 0, opt_scale_sim_1},
    {"scale_sim_2", required_argument, 0, opt_scale_sim_2},
    {"scale_mea_1", required_argument, 0, opt_scale_mea_1},
    {"scale_mea_2", required_argument, 0, opt_scale_mea_2},
    {nullptr, 0, nullptr, 0}
  };

  int opt;
  while ((opt = getopt_long(argc, argv, "", longOptions, 0)) != -1){
    switch(opt) {
    case opt_name_f_sim_1:
      if (optarg == nullptr || optarg[0] == '-'){
        std::cerr << "Error: --opt_name_f_sim_1 requires an argument but none was provided.\n";
        return -1;
      }
      name_f_sim_1 = optarg;
      if (!CheckRootFile(name_f_sim_1.c_str())){
        std::cerr << "Error: Invalid simulated file 1: " << name_f_sim_1 << std::endl;
        return -1;
      }
      break;

    case opt_name_f_sim_2:
      if (optarg == nullptr || optarg[0] == '-'){
        std::cerr << "Error: --opt_name_f_sim_2 requires an argument but none was provided.\n";
        return -1;
      }
      name_f_sim_2 = optarg;
      if (!CheckRootFile(name_f_sim_2.c_str())){
        std::cerr << "Error: Invalid simulated file 2: " << name_f_sim_2 << std::endl;
        return -1;
      }
      break;
        
    case opt_name_f_mea_1:
      if (optarg == nullptr || optarg[0] == '-'){
        std::cerr << "Error: --opt_name_f_mea_1 requires an argument but none was provided.\n";
        return -1;
      }
      name_f_mea_1 = optarg;
      if (!CheckRootFile(name_f_mea_1.c_str())){
        std::cerr << "Error: Invalid measurement file 1: " << name_f_mea_1 << std::endl;
        return -1;
      }
      break;

    case opt_name_f_mea_2:
      if (optarg == nullptr || optarg[0] == '-'){
        std::cerr << "Error: --opt_name_f_mea_2 requires an argument but none was provided.\n";
        return -1;
      }
      name_f_mea_2 = optarg;
      if (!CheckRootFile(name_f_mea_2.c_str())){
        std::cerr << "Error: Invalid measurement file 2: " << name_f_mea_2 << std::endl;
        return -1;
      }
      break;
            
    case opt_channel:
      WriteOutput(optarg, channel, channel, "channel");
      break;
      
    case opt_entries_sim_fft:
      WriteOutput(optarg, entries_sim_fft, entries_sim_fft, "number of simulation entries for FFT");
      break;

    case opt_entries_mea_fft:
      WriteOutput(optarg, entries_mea_fft, entries_mea_fft, "number of measurement entries for FFT");
      break;
      
    case opt_entries_sim_descent:
      WriteOutput(optarg, entries_sim_descent, entries_sim_descent, "number of simulation entries for gradient descent");
      break;

    case opt_entries_mea_descent:
      WriteOutput(optarg, entries_mea_descent, entries_mea_descent, "number of measurement entries for gradient descent");
      break;

    case opt_x_min_sim:
      WriteOutput(optarg, x_min_sim, x_min_sim, "x_min for simulation");
      bin_sim = (x_max_sim - x_min_sim)*1000.;
      break;

    case opt_x_max_sim:
      WriteOutput(optarg, x_max_sim, x_max_sim, "x_max for simulation");
      bin_sim = (x_max_sim - x_min_sim)*1000.;
      break;

    case opt_x_min_mea:
      WriteOutput(optarg, x_min_mea, x_min_mea, "x_min for measurement");
      bin_mea = (x_max_mea - x_min_mea);
      break;

    case opt_x_max_mea:
      WriteOutput(optarg, x_max_mea, x_max_mea, "x_max for measurement");
      bin_mea = (x_max_mea - x_min_mea);
      break;

    case opt_energy_1:
      WriteOutput(optarg, energy_1, energy_1, "Energy of the first Compton edge", "(MeV)");
      break;

    case opt_energy_2:
      WriteOutput(optarg, energy_2, energy_2, "Energy of the second Compton edge", "(MeV)");
      break;
      
    case opt_res_1:
      WriteOutput(optarg, res_1, res_1, "Energy resolution of the first Compton edge", "(%)");
      break;
      
    case opt_res_2:
      WriteOutput(optarg, res_2, res_2, "Energy resolution of the second Compton edge", "(%)");
      break;
              
    case opt_x_min_fft_1:
      WriteOutput(optarg, x_min_fft_1, x_min_fft_1, "Minimum channel for FFT 1");
      break;
      
    case opt_x_max_fft_1:
      WriteOutput(optarg, x_max_fft_1, x_max_fft_1, "Maximum channel for FFT 1");
      break;
      
    case opt_x_min_fft_2:
      WriteOutput(optarg, x_min_fft_2, x_min_fft_2, "Minimum channel for FFT 2");
      break;
      
    case opt_x_max_fft_2:
      WriteOutput(optarg, x_max_fft_2, x_max_fft_2, "Maximum channel for FFT 2");
      break;
      
    case opt_x_min_descent_1:
      WriteOutput(optarg, x_min_descent_1, x_min_descent_1, "Minimum channel for gradient descent 1");
      break;
      
    case opt_x_max_descent_1:
      WriteOutput(optarg, x_max_descent_1, x_max_descent_1, "Maximum channel for gradient descent 1");
      break;
      
    case opt_x_min_descent_2:
      WriteOutput(optarg, x_min_descent_2, x_min_descent_2, "Minimum channel for gradient descent 2");
      break;
      
    case opt_x_max_descent_2:
      WriteOutput(optarg, x_max_descent_2, x_max_descent_2, "Maximum channel for gradient descent 2");
      break;
      
    case opt_rate_fft_ini_1:
      WriteOutput(optarg, rate_fft_ini_1, rate_fft_ini_1, "FFT cut-off growth rate 1 (estimation)");
      break;
      
    case opt_rate_fft_ini_2:
      WriteOutput(optarg, rate_fft_ini_2, rate_fft_ini_2, "FFT cut-off growth rate 2 (estimation)");
      break;
      
    case opt_thresh_fft_ini_1:
      WriteOutput(optarg, thresh_fft_ini_1, thresh_fft_ini_1, "FFT cut-off threshold 1 (estimation)");
      break;
      
    case opt_thresh_fft_ini_2:
      WriteOutput(optarg, thresh_fft_ini_2, thresh_fft_ini_2, "FFT cut-off threshold 2 (estimation)");
      break;
      
    case opt_rate_fft_descent_1:
      WriteOutput(optarg, rate_fft_descent_1, rate_fft_descent_1, "FFT cut-off growth rate 1 (descent)");
      break;
      
    case opt_rate_fft_descent_2:
      WriteOutput(optarg, rate_fft_descent_2, rate_fft_descent_2, "FFT cut-off growth rate 2 (descent)");
      break;
      
    case opt_thresh_fft_descent_1:
      WriteOutput(optarg, thresh_fft_descent_sim, thresh_fft_descent_sim, "FFT cut-off threshold 1 (descent)");
      break;
      
    case opt_thresh_fft_descent_2:
      WriteOutput(optarg, thresh_fft_descent_mea, thresh_fft_descent_mea, "FFT cut-off threshold 2 (descent)");
      break;
      
    case opt_learning_rate_channel_1:
      WriteOutput(optarg, learning_rate_channel_1, learning_rate_channel_1, "Gradient descent learning rate for estimating channel 1");
      break;
      
    case opt_learning_rate_channel_2:
      WriteOutput(optarg, learning_rate_channel_2, learning_rate_channel_2, "Gradient descent learning rate for estimating channel 2");
      break;
      
    case opt_learning_rate_res_1:
      WriteOutput(optarg, learning_rate_res_1, learning_rate_res_1, "Gradient descent learning rate for estimating energy resolution 1");
      break;
      
    case opt_learning_rate_res_2:
      WriteOutput(optarg, learning_rate_res_2, learning_rate_res_2, "Gradient descent learning rate for estimating energy resolution 2");
      break;
      
    case opt_step_channel_1:
      WriteOutput(optarg, step_channel_1, step_channel_1, "Step of channel 1 for calculating gradient");
      break;
      
    case opt_step_channel_2:
      WriteOutput(optarg, step_channel_2, step_channel_2, "Step of channel 2 for calculating gradient");
      break;
      
    case opt_step_res_1:
      WriteOutput(optarg, step_res_1, step_res_1, "Step of energy resolution 1 for calculating gradient");
      break;
      
    case opt_step_res_2:
      WriteOutput(optarg, step_res_2, step_res_2, "Step of energy resolution 2 for calculating gradient");
      break;
      
    case opt_scale_sim_1:
      WriteOutput(optarg, scale_sim_1, scale_sim_1, "Scaling factor for simulation histogram 1");
      break;
      
    case opt_scale_sim_2:
      WriteOutput(optarg, scale_sim_2, scale_sim_2, "Scaling factor for simulation histogram 2");
      break;
      
    case opt_scale_mea_1:
      WriteOutput(optarg, scale_mea_1, scale_mea_1, "Scaling factor for measurement histogram 1");
      break;
      
    case opt_scale_mea_2:
      WriteOutput(optarg, scale_mea_2, scale_mea_2, "Scaling factor for measurement histogram 2");
      break;

    case '?': // Handle unknown options
      return -1;
    }
  }

  // Ensure both simulated files are provided
  if (name_f_sim_1.empty() || name_f_sim_2.empty() || name_f_mea_1.empty() || name_f_mea_2.empty()) {
      std::cerr << "Error: --opt_name_f_sim_1 --opt_name_f_sim_2 --opt_name_f_mea_1 --opt_name_f_mea_2 must be specified.\n";
      return -1;
  }
  
  return 0;
}
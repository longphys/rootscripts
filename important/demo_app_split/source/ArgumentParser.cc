#include "ArgumentParser.hh"

int ArgumentParser::Parse(int argc, char** argv, Config& config){
  static struct option longOptions[] = {
    {"name_f_sim_1", required_argument, 0, opt_name_f_sim_1},
    {"name_f_sim_2", required_argument, 0, opt_name_f_sim_2},
    {"ascii", required_argument, 0, opt_ascii},
    {"name_f_mea_1", required_argument, 0, opt_name_f_mea_1},
    {"name_f_mea_2", required_argument, 0, opt_name_f_mea_2},
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
      config.name_f_sim_1 = optarg;
      if (!CheckRootFile(config.name_f_sim_1.c_str())){
        std::cerr << "Error: Invalid simulated file 1: " << config.name_f_sim_1 << std::endl;
        return -1;
      }
      break;

    case opt_ascii:
    WriteOutput(optarg, config.ascii, config.ascii, "ASCII option");
    break;

    case opt_name_f_sim_2:
      if (optarg == nullptr || optarg[0] == '-'){
        std::cerr << "Error: --opt_name_f_sim_2 requires an argument but none was provided.\n";
        return -1;
      }
      config.name_f_sim_2 = optarg;
      if (!CheckRootFile(config.name_f_sim_2.c_str())){
        std::cerr << "Error: Invalid simulated file 2: " << config.name_f_sim_2 << std::endl;
        return -1;
      }
      break;
        
    case opt_name_f_mea_1:
      if (optarg == nullptr || optarg[0] == '-'){
        std::cerr << "Error: --opt_name_f_mea_1 requires an argument but none was provided.\n";
        return -1;
      }
      config.name_f_mea_1 = optarg;
      if (!CheckRootFile(config.name_f_mea_1.c_str())){
        std::cerr << "Error: Invalid measurement file 1: " << config.name_f_mea_1 << std::endl;
        return -1;
      }
      break;

    case opt_name_f_mea_2:
      if (optarg == nullptr || optarg[0] == '-'){
        std::cerr << "Error: --opt_name_f_mea_2 requires an argument but none was provided.\n";
        return -1;
      }
      config.name_f_mea_2 = optarg;
      if (!CheckRootFile(config.name_f_mea_2.c_str())){
        std::cerr << "Error: Invalid measurement file 2: " << config.name_f_mea_2 << std::endl;
        return -1;
      }
      break;
            
    case opt_entries_sim_fft:
      WriteOutput(optarg, config.entries_sim_fft, config.entries_sim_fft, "number of simulation entries for FFT");
      break;

    case opt_entries_mea_fft:
      WriteOutput(optarg, config.entries_mea_fft, config.entries_mea_fft, "number of measurement entries for FFT");
      break;
      
    case opt_entries_sim_descent:
      WriteOutput(optarg, config.entries_sim_descent, config.entries_sim_descent, "number of simulation entries for gradient descent");
      break;

    case opt_entries_mea_descent:
      WriteOutput(optarg, config.entries_mea_descent, config.entries_mea_descent, "number of measurement entries for gradient descent");
      break;

    case opt_x_min_sim:
      WriteOutput(optarg, config.x_min_sim, config.x_min_sim, "x_min for simulation");
      config.bin_sim = (config.x_max_sim - config.x_min_sim)*1000.;
      break;

    case opt_x_max_sim:
      WriteOutput(optarg, config.x_max_sim, config.x_max_sim, "x_max for simulation");
      config.bin_sim = (config.x_max_sim - config.x_min_sim)*1000.;
      break;

    case opt_x_min_mea:
      WriteOutput(optarg, config.x_min_mea, config.x_min_mea, "x_min for measurement");
      config.bin_mea = (config.x_max_mea - config.x_min_mea);
      break;

    case opt_x_max_mea:
      WriteOutput(optarg, config.x_max_mea, config.x_max_mea, "x_max for measurement");
      config.bin_mea = (config.x_max_mea - config.x_min_mea);
      break;

    case opt_energy_1:
      WriteOutput(optarg, config.energy_1, config.energy_1, "Energy of the first Compton edge", "(MeV)");
      break;

    case opt_energy_2:
      WriteOutput(optarg, config.energy_2, config.energy_2, "Energy of the second Compton edge", "(MeV)");
      break;
      
    case opt_res_1:
      WriteOutput(optarg, config.res_1, config.res_1, "Energy resolution of the first Compton edge", "(%)");
      break;
      
    case opt_res_2:
      WriteOutput(optarg, config.res_2, config.res_2, "Energy resolution of the second Compton edge", "(%)");
      break;
              
    case opt_x_min_fft_1:
      WriteOutput(optarg, config.x_min_fft_1, config.x_min_fft_1, "Minimum channel for FFT 1");
      break;
      
    case opt_x_max_fft_1:
      WriteOutput(optarg, config.x_max_fft_1, config.x_max_fft_1, "Maximum channel for FFT 1");
      break;
      
    case opt_x_min_fft_2:
      WriteOutput(optarg, config.x_min_fft_2, config.x_min_fft_2, "Minimum channel for FFT 2");
      break;
      
    case opt_x_max_fft_2:
      WriteOutput(optarg, config.x_max_fft_2, config.x_max_fft_2, "Maximum channel for FFT 2");
      break;
      
    case opt_x_min_descent_1:
      WriteOutput(optarg, config.x_min_descent_1, config.x_min_descent_1, "Minimum channel for gradient descent 1");
      break;
      
    case opt_x_max_descent_1:
      WriteOutput(optarg, config.x_max_descent_1, config.x_max_descent_1, "Maximum channel for gradient descent 1");
      break;
      
    case opt_x_min_descent_2:
      WriteOutput(optarg, config.x_min_descent_2, config.x_min_descent_2, "Minimum channel for gradient descent 2");
      break;
      
    case opt_x_max_descent_2:
      WriteOutput(optarg, config.x_max_descent_2, config.x_max_descent_2, "Maximum channel for gradient descent 2");
      break;
      
    case opt_rate_fft_ini_1:
      WriteOutput(optarg, config.rate_fft_ini_1, config.rate_fft_ini_1, "FFT cut-off growth rate 1 (estimation)");
      break;
      
    case opt_rate_fft_ini_2:
      WriteOutput(optarg, config.rate_fft_ini_2, config.rate_fft_ini_2, "FFT cut-off growth rate 2 (estimation)");
      break;
      
    case opt_thresh_fft_ini_1:
      WriteOutput(optarg, config.thresh_fft_ini_1, config.thresh_fft_ini_1, "FFT cut-off threshold 1 (estimation)");
      break;
      
    case opt_thresh_fft_ini_2:
      WriteOutput(optarg, config.thresh_fft_ini_2, config.thresh_fft_ini_2, "FFT cut-off threshold 2 (estimation)");
      break;
      
    case opt_rate_fft_descent_1:
      WriteOutput(optarg, config.rate_fft_descent_1, config.rate_fft_descent_1, "FFT cut-off growth rate 1 (descent)");
      break;
      
    case opt_rate_fft_descent_2:
      WriteOutput(optarg, config.rate_fft_descent_2, config.rate_fft_descent_2, "FFT cut-off growth rate 2 (descent)");
      break;
      
    case opt_thresh_fft_descent_1:
      WriteOutput(optarg, config.thresh_fft_descent_sim, config.thresh_fft_descent_sim, "FFT cut-off threshold 1 (descent)");
      break;
      
    case opt_thresh_fft_descent_2:
      WriteOutput(optarg, config.thresh_fft_descent_mea, config.thresh_fft_descent_mea, "FFT cut-off threshold 2 (descent)");
      break;
      
    case opt_learning_rate_channel_1:
      WriteOutput(optarg, config.learning_rate_channel_1, config.learning_rate_channel_1, "Gradient descent learning rate for estimating channel 1");
      break;
      
    case opt_learning_rate_channel_2:
      WriteOutput(optarg, config.learning_rate_channel_2, config.learning_rate_channel_2, "Gradient descent learning rate for estimating channel 2");
      break;
      
    case opt_learning_rate_res_1:
      WriteOutput(optarg, config.learning_rate_res_1, config.learning_rate_res_1, "Gradient descent learning rate for estimating energy resolution 1");
      break;
      
    case opt_learning_rate_res_2:
      WriteOutput(optarg, config.learning_rate_res_2, config.learning_rate_res_2, "Gradient descent learning rate for estimating energy resolution 2");
      break;
      
    case opt_step_channel_1:
      WriteOutput(optarg, config.step_channel_1, config.step_channel_1, "Step of channel 1 for calculating gradient");
      break;
      
    case opt_step_channel_2:
      WriteOutput(optarg, config.step_channel_2, config.step_channel_2, "Step of channel 2 for calculating gradient");
      break;
      
    case opt_step_res_1:
      WriteOutput(optarg, config.step_res_1, config.step_res_1, "Step of energy resolution 1 for calculating gradient");
      break;
      
    case opt_step_res_2:
      WriteOutput(optarg, config.step_res_2, config.step_res_2, "Step of energy resolution 2 for calculating gradient");
      break;
      
    case opt_scale_sim_1:
      WriteOutput(optarg, config.scale_sim_1, config.scale_sim_1, "Scaling factor for simulation histogram 1");
      break;
      
    case opt_scale_sim_2:
      WriteOutput(optarg, config.scale_sim_2, config.scale_sim_2, "Scaling factor for simulation histogram 2");
      break;
      
    case opt_scale_mea_1:
      WriteOutput(optarg, config.scale_mea_1, config.scale_mea_1, "Scaling factor for measurement histogram 1");
      break;
      
    case opt_scale_mea_2:
      WriteOutput(optarg, config.scale_mea_2, config.scale_mea_2, "Scaling factor for measurement histogram 2");
      break;

    case '?': // Handle unknown options
      return -1;
    }
  }

  // Ensure both simulated files are provided
  if (config.name_f_sim_1.empty() || config.name_f_sim_2.empty() || config.name_f_mea_1.empty() || config.name_f_mea_2.empty()) {
      std::cerr << "Error: --opt_name_f_sim_1 --opt_name_f_sim_2 --opt_name_f_mea_1 --opt_name_f_mea_2 must be specified.\n";
      return -1;
  }
  
  return 0;
}
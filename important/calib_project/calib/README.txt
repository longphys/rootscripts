# This application is used to perform energy calibration based on 2 measurements of calibration sources.
# Requirements:

  CMake ≥ 3.10
  C++ Compiler (GCC/Clang/MSVC) supporting C++17
  ROOT Framework version 6.30.02

# Build the Application:

  mkdir build && cd build
  cmake ..
  make

# Run the Program:

  ./demo macro.file

# The result of calibration coefficients will be showed on the terminal and saved in a "coefficients.txt" file as default.

# The application accepts input arguments through a macro file.
# All arguments must be filled for the application to run.
# Lines that start with the symbol "--" are treated as arguments.
# Lines that start with the symbol "#" are treated as comments.

# The example ASCII files represent the accquired energy spectra from a measurement.
# The first columnn represents the bin corresponds to amplitude or channel, and the second column represents the height of the bin.
# Examples of measurement files (ROOT and ASCII type), simulation files and a macro file are included in the root folder of the application.

# Below is a list of the arguments and their descriptions:
  --name_f_sim_1 <path> - Path to the first simulation ROOT file.
  --name_f_sim_2 <path> - Path to the second simulation ROOT file.
  --ascii <int> - The choice to use a ROOT type or ASCII type for the measurement file. (0 for ROOT, other for ASCII)
  --name_f_mea_1 <path> - Path to the first measurement ROOT file.
  --name_f_mea_2 <path> - Path to the second measurement ROOT file.
  --entries_sim_fft <int> - Number of entries from the simulation file for quick Fast Fourier Transform (FFT) to perform the first calibration.
  --entries_mea_fft <int> - Number of entries from the measurement file for quick FFT to perform the first calibration.
  --entries_sim_descent <int> - Number of entries from the simulation file for Gradient descent to perform up to the final calibration.
  --entries_mea_descent <int> - Number of entries from the measurement file for Gradient descent to perform up to the final calibration.
# when "entries_sim_descent" and "entries_sim_descent" arguments are used with variables below 0, all events available in the files will be used.

### THESE PARAMETERS ONLY NEED TO BE DEFINED ONCE (THE RANGE TO DRAW THE HISTOGRAMS) ###
  --x_min_sim <double> - Lower limit of the simulation file histogram. (in MeV)
  --x_max_sim <double> - Upper limit of the simulation file histogram. (in MeV)
  --x_min_mea <double> - Lower limit of the measurement file histogram. (in channels)
  --x_max_mea <double> - Upper limit of the measurement file histogram. (in channels)
########################################################################################

### THESE PARAMETERS ONLY NEED TO BE DEFINED ONCE (CHARACTERISTICS OF THE COMPTON EDGES) ###
  --energy_1 <double> - Energy value of the first Compton edge. (in MeV)
  --energy_2 <double> - Energy value of the second Compton edge. (in MeV)
  --res_1 8. <double> - Energy resolution at the first Compton edge. (in %)
  --res_2 5. <double> - Energy resolution at the second Compton edge. (in %)
############################################################################################

### THESE PARAMETERS ONLY NEED TO BE DEFINED ONCE (THE RANGE TO FIND THE FIRST POSITION OF THE COMPTON EDGES) ###
  --x_min_fft_1 200 <double> - Lower limit to perform quick FFT for first calibration (first file). (in channels)
  --x_max_fft_1 1000 <double> - Upper limit to perform quick FFT for first calibration (first file). (in channels)
  --x_min_fft_2 300 <double> - Lower limit to perform quick FFT for second calibration (first file). (in channels)
  --x_max_fft_2 1000 <double> - Upper limit to perform quick FFT for second calibration (first file). (in channels)
#################################################################################################################

### THESE PARAMETERS ONLY NEED TO BE DEFINED ONCE (NOISE FILTER TO FIND THE FIRST POSITION OF THE COMPTON EDGES) ###
  --rate_fft_ini_1 0.1 <double> - First parameter to perform quick FFT for first calibration (first file).
  --rate_fft_ini_2 0.1 <double> - First parameter to perform quick FFT for first calibration (second file).
  --thresh_fft_ini_1 50 <double> - Second parameter to perform quick FFT for first calibration (first file).
  --thresh_fft_ini_2 50 <double> - Second parameter to perform quick FFT for first calibration (second file).
####################################################################################################################
  
### THESE ARE IMPORTANT PARAMETERS TO MODIFY MANY TIMES (THE WINDOW OF THE GRADIENT DESCENT LOOP) ###
  --x_min_descent_1 0.35 <double> - Lower limit to perform gradient descent up to last calibration (first file). (in MeV)
  --x_max_descent_1 0.6 <double> - Upper limit to perform gradient descent up to last calibration (first file). (in MeV)
  --x_min_descent_2 0.8 <double> - Lower limit to perform gradient descent up to last calibration (second file). (in MeV)
  --x_max_descent_2 1.2 <double> - Upper limit to perform gradient descent up to last calibration (second file). (in MeV)
#####################################################################################################
  
### THESE ARE IMPORTANT PARAMETERS TO MODIFY MANY TIMES (NOISE FILTER OF THE DATA) ###
  --rate_fft_descent_1 0.5 <double> - First parameter to perform FFT up to last calibration (first file).
  --rate_fft_descent_2 0.5 <double> - First parameter to perform FFT up to last calibration (second file).
  --thresh_fft_descent_sim 1000 <double> - First parameter to perform FFT up to last calibration (first file).
  --thresh_fft_descent_mea 1000 <double> - Second parameter to perform FFT up to last calibration (second file).
######################################################################################
  
### THESE ARE IMPORTANT PARAMETERS TO MODIFY MANY TIMES (LEARNING RATE OF THE GRADIENT DESCENT) ###
  --learning_rate_channel_1 1. <double> - Learning rate of the Compton edge position of the first file.
  --learning_rate_channel_2 1. <double> - Learning rate of the Compton edge position of the second file.
  --learning_rate_res_1 0.01 <double> - Learning rate of the resolution of the Compton edge position of the first file.
  --learning_rate_res_2 0.01 <double> - Learning rate of the resolution of the Compton edge position of the second file.
  --step_channel_1 1 <double> - Step size of the Compton edge position of the first file.
  --step_channel_2 1 <double> - Step size of the Compton edge position of the second file.
  --step_res_1 0.05 <double> - Step size of the resolution of the Compton edge position of the first file.
  --step_res_2 0.05 <double> - Step size of the resolution of the Compton edge position of the second file.
###################################################################################################

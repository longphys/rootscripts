Requirements:

  CMake ≥ 3.10
  C++ Compiler (GCC/Clang/MSVC) supporting C++17
  ROOT Framework (Download ROOT)

Build the Application:

  mkdir build && cd build
  cmake ..
  make

Run the Program:

  ./demo macro.file

The application accepts input arguments through a macro file.
Below is a list of available arguments and their descriptions:

  --name_f_sim_1 <path> – Path to the first simulation ROOT file.
  --name_f_sim_2 <path> – Path to the second simulation ROOT file.
  --name_f_mea_1 <path> – Path to the first measurement ROOT file.
  --name_f_mea_2 <path> – Path to the second measurement ROOT file.
  --channel 2 <int> - Channel of interest of the measurement file.
  --entries_sim_fft <int> - Number of entries from the simulation file for quick Fast Fourier Transform (FFT) to perform the first calibration.
  --entries_mea_fft <int> - Number of entries from the measurement file for quick FFT to perform the first calibration.
  --entries_sim_descent <int> - Number of entries from the simulation file for Gradient descent to perform up to the final calibration.
  --entries_mea_descent <int> - Number of entries from the measurement file for Gradient descent to perform up to the final calibration.
  --x_min_sim <double> - Lower limit of the simulation file histogram. (in MeV)
  --x_max_sim <double> - Upper limit of the simulation file histogram. (in MeV)
  --x_min_mea <double> - Lower limit of the measurement file histogram. (in channels)
  --x_max_mea <double> - Upper limit of the measurement file histogram. (in channels)
  --energy_1 <double> - Energy value of the first Compton edge. (in MeV)
  --energy_2 <double> - Energy value of the second Compton edge. (in MeV)
  --res_1 8.
  --res_2 5.
  --x_min_fft_1 200
  --x_max_fft_1 1000
  --x_min_fft_2 300
  --x_max_fft_2 1000
  --x_min_descent_1 0.35
  --x_max_descent_1 0.6
  --x_min_descent_2 0.8
  --x_max_descent_2 1.2
  --rate_fft_ini_1 0.1
  --rate_fft_ini_2 0.1
  --thresh_fft_ini_1 50
  --thresh_fft_ini_2 50
  --rate_fft_descent_1 0.5
  --rate_fft_descent_2 0.5
  --thresh_fft_descent_sim 1000
  --thresh_fft_descent_mea 1000
  --learning_rate_channel_1 1.
  --learning_rate_channel_2 1.
  --learning_rate_res_1 0.01
  --learning_rate_res_2 0.01
  --step_channel_1 1
  --step_channel_2 1
  --step_res_1 0.05
  --step_res_2 0.05
  --scale_sim_1 8.3
  --scale_sim_2 4.2
  --scale_mea_1 1.
  --scale_mea_2 1.

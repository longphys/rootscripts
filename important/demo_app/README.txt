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
  --channel 2 
  --entries_sim_fft 100000
  --entries_mea_fft 100000
  --entries_sim_descent 500000
  --entries_mea_descent 500000
  --x_min_sim 0.0
  --x_max_sim 1.3
  --x_min_mea 100
  --x_max_mea 1000
  --energy_1 0.477
  --energy_2 1.061
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
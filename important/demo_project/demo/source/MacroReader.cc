#include "MacroReader.hh"
// Function to read the macro file and simulate command-line arguments
std::vector<std::string> MacroReader::Read(const char* filename) {
  std::vector<std::string> args;
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Error: Cannot open macro file: " << filename << std::endl;
    exit(1);
  }
  
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string word;
    while (iss >> word) {
        args.push_back(word);  // Store each word as an argument
    }
  }
  
  return args;
}
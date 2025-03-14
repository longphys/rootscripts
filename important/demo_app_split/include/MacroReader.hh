#ifndef MACROREADER_HH
#define MACROREADER_HH

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

class MacroReader{
public:
  std::vector<std::string> Read(const char* filename);
};
#endif
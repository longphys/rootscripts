#ifndef ARGUMENTPARSER_HH
#define ARGUMENTPARSER_HH

#include <string>

class ArgumentParser {
public:
    std::string file1, file2;
    int channel = 0;
    
    void Parse(int argc, char** argv);
};

#endif
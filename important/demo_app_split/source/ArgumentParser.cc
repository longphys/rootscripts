// ArgumentParser.cc
#include "ArgumentParser.hh"
#include <iostream>
#include <getopt.h>
#include <cstdlib>

void ArgumentParser::Parse(int argc, char** argv) {
    struct option long_options[] = {
        {"file1", required_argument, nullptr, 'a'},
        {"file2", required_argument, nullptr, 'b'},
        {"channel", required_argument, nullptr, 'c'},
        {nullptr, 0, nullptr, 0}
    };
    
    int opt;
    while ((opt = getopt_long(argc, argv, "a:b:c:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'a':
                file1 = optarg;
                break;
            case 'b':
                file2 = optarg;
                break;
            case 'c':
                channel = std::stoi(optarg);
                break;
            default:
                std::cerr << "Usage: --file1 <file> --file2 <file> --channel <num>" << std::endl;
                exit(EXIT_FAILURE);
        }
    }
}
#include "load_reads2.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <cctype>
#include "utility.h"
#include "common.h"

const int MAX_STR = 1024;



void load_reads(std::vector<std::string>& Read, std::string read_left_file, std::string read_right_file) {

    time_t beg = time(NULL);

    Read.reserve(1000000);
    std::fstream left_in, right_in;
    
    left_in.open(read_left_file.c_str(), std::fstream::in);
    right_in.open(read_right_file.c_str(), std::fstream::in);
    
    if (!left_in.is_open()||!right_in.is_open()) {
        std::cerr << "[error] File " << read_left_file << "or" << read_right_file << " can't be opened." << std::endl;
        exit(1);
    }

    std::cerr << "Begin loading reads ..." << std::endl;
    int line_count_1 = 0, line_count_2 = 0;

    char c_line[MAX_STR];
    while (!left_in.eof()) {
        left_in.getline(c_line, MAX_STR);
        line_count_1++;
        if (line_count_1 % 4 == 2) {
            std::string sequence(c_line);
            std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                [](unsigned char c) { return std::toupper(c); });
            /*
            if (sequence.find("N") != std::string::npos) {
                sequence = "";
                Read.push_back(sequence);
                continue;
            }
            */
            Read.push_back(sequence);
        }
    }

    while (!right_in.eof()) {
        right_in.getline(c_line, MAX_STR);
        line_count_2++;

        if (line_count_2 % 4 == 2 ) {
            std::string sequence(c_line);
            std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                [](unsigned char c) { return std::toupper(c); });
            /*
            if (sequence.find("N") != std::string::npos) {
                sequence = "";
                Read.push_back(sequence);
                continue;
            }
            */
            if(g_reads_direction == 1)
                sequence = revcomp(sequence);
            Read.push_back(sequence);
        }
    }

    time_t end = time(NULL);

    std::cerr << Read.size() << " reads have been loaded ! (elapsed time : " << (end - beg) << " s)" << std::endl;
    right_in.close();
    /*
    std::ofstream outfile("/home/bioinfo/limh/paper/data/SRR1510027/Read.txt"); // 创建输出文件流
    if (outfile.is_open()) {
        for (const std::string& line : Read) {
            outfile << line << std::endl;
        }
        outfile.close();
    }
   */ 
}
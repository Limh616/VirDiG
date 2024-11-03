#include "load_reads2.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <cctype>
#include <thread>
#include <mutex>
#include "utility.h"
#include "common.h"

const int MAX_STR = 1024;
std::mutex read_mutex;
std::vector<std::string> right_reads;


void load_reads_from_file(const std::string& filename, std::vector<std::string>& local_reads, bool is_left) {
    std::fstream in_file;
    in_file.open(filename.c_str(), std::fstream::in);
    if (!in_file.is_open()) {
        std::cerr << "[error] File " << filename << " can't be opened." << std::endl;
        exit(1);
    }

    char c_line[MAX_STR];
    int line_count = 0;
    while (!in_file.eof()) {
        in_file.getline(c_line, MAX_STR);
        line_count++;
        if (line_count % 4 == 2) {
            std::string sequence(c_line);
            std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                [](unsigned char c) { return std::toupper(c); });
            if (!is_left && g_reads_direction == 1) {
                sequence = revcomp(sequence);
            }
         
            std::lock_guard<std::mutex> guard(read_mutex);
            local_reads.push_back(sequence);
        }
    }
    in_file.close();
}


void load_reads(std::vector<std::string>& Read, std::string read_left_file, std::string read_right_file) {

    time_t beg = time(NULL);

    std::cerr << "Begin loading reads ..." << std::endl;

   
    std::thread left_thread(load_reads_from_file, read_left_file, std::ref(Read), true);
    std::thread right_thread(load_reads_from_file, read_right_file, std::ref(right_reads), false);

   
    left_thread.join();
    right_thread.join();

    
    {
        std::lock_guard<std::mutex> guard(read_mutex); 
        Read.insert(Read.end(), right_reads.begin(), right_reads.end());
    }

    right_reads.clear();
    time_t end = time(NULL);
    std::cerr << Read.size() << " reads have been loaded ! (elapsed time : " << (end - beg) << " s)" << std::endl;
}



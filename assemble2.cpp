// assemble2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "load_reads2.h"
//#include "common.h"
#include "utility.h"
#include "transcript.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <errno.h>
#include <ctime>
#include <getopt.h>

#include <boost/unordered_map.hpp>


std::string reads_left_file, reads_right_file;
///*
int g_kmer_length = 31;
#define map_weight		301
#define non_canonical  302
bool g_help;
std::string out_dir;
int g_reads_direction = 1;
float g_map_weight = 0.7;
int g_non_canonical = 0;
int g_threads = 6;
// */

static const char* short_options = "l:r:k:o:d:h:t:n";

static struct option long_options[] = {
    // general options

    {"reads_left",                  required_argument,      0,      'l'},
    {"reads_right",                 required_argument,      0,      'r'},
    {"kmer_length",                 required_argument,      0,      'k'},
    {"out_dir",                     required_argument,      0,      'o'},
    {"reads_direction",             required_argument,      0,      'd'},
    {"threads",                     required_argument,      0,      't'},
    {"non_canonical",               required_argument,      0,      non_canonical},
    {"map_weight",                  required_argument,      0,      map_weight},
    {"help",                        no_argument,            0,      'h'},

    //...
      {0,0,0,0} // terminator

};

int parse_options(int argc, char* argv[]);
std::string usage();

static size_t rg_index = 0;

std::string base_name() {

    std::stringstream idx;
    idx << "comp" << rg_index;
    rg_index++;

    return idx.str();
}

int parse_options(int argc, char* argv[]) {

    int option_index = 0;
    int next_option;

    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index); //解析命令行参数
        switch (next_option) {
        case -1:     /* Done with options. */
            break;
        case 'k':
            g_kmer_length = atoi(optarg);
            break;
        case 'r':
            reads_right_file = optarg;
            break;
        case 'l':
            reads_left_file = optarg;
            break;
        case 'h':
            g_help = true;
            break;
        case 'd':
            g_reads_direction = atoi(optarg);
            break;
        case 't':
            g_threads = atoi(optarg);
            break;
        case  non_canonical:
            g_non_canonical = atoi(optarg);
            break;
        case  map_weight:
            g_map_weight = atof(optarg);
            break;
        case 'o':
            out_dir = optarg;
            break;
        default:
            std::cout << usage();
            exit(1);
        }
    } while (next_option != -1);

    if (g_help) {
        std::cout << usage();
        exit(0);
    }

    if (reads_left_file == "" || reads_right_file == "") {
        std::cerr << "Error : --input option needs an argument!! " << std::endl;
        std::cout << usage();
        exit(1);
    }
    if (g_reads_direction != 1 && g_reads_direction != 2) {
        std::cout << "Error: --fr can only be 1 or 2" << std::endl;
        exit(1);
    }


    if (g_kmer_length > 32) {
        errAbort(const_cast<char*>("Length of kmer can not be excess 32!\n"));
    }

    return 0;
}


std::string usage() {

    std::stringstream usage_info;
    usage_info
        << std::endl
        << "===============================================================================" << std::endl
        << "Usage " << std::endl
        << "===============================================================================" << std::endl
        << " ** Options: **" << std::endl
        << "  -k <int>: length of kmer, default 31. " << std::endl
        << "  -o <string>: output directory. " << std::endl
        << "  -h : help information. " << std::endl
        << "  Pair end reads: " << std::endl
        << "  -l <string>: left reads file name (.fasta). " << std::endl
        << "  -r <string>: right reads file name (.fasta). " << std::endl
        << "  -d <int>: pair-end reads directions can be defined, 1: opposite directions  2: same direction. default: 1. " << std::endl
        << "  -t <int>: number of threads, default 6. " << std::endl
        << "  --non_canonical <int>: whether to generate non-standard transcripts, 1 : true, 0 : false, default 0. " << std::endl
        << "  --map_weight <float>: paired-end reads are assigned paired node weights, recommended to be in the range of 0 to 1, default 0.7. " << std::endl
        << "===============================================================================" << std::endl
        << std::endl;

    return (usage_info.str());
}




int main(int argc, char* argv[]) {


    time_t beg = time(NULL);
    // process command line arguments，处理命令行参数
    int parse_ret = parse_options(argc, argv);
    if (parse_ret)
        return parse_ret;

    time_t begin = time(NULL);

    // load data
    std::vector<std::string> Read;
    if (reads_left_file != "" && reads_right_file != "")
        load_reads(Read, reads_left_file, reads_right_file);

    Construct_splicing_graphs(Read);
    time_t end = time(NULL);
    std::cout << "ALL finished  (elapsed time: " << (end - beg) << " s)" << std::endl;
    return 0;
}


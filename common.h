#ifndef COMMON_H
#define COMMON_H

/*
 *  common.h
 */

#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>


extern int first_last_k;
 // general options
extern int g_kmer_length ;
extern bool g_help;

extern bool g_is_paired_end;
extern std::string out_dir;
extern int g_reads_direction  ;
extern int g_ref_genome_len ;
extern float g_map_weight;

#endif


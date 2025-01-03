# VirDiG

VirDiG is an efficient de novo transcriptome assembler for coronavirus. It can assemble transcripts from 
short reads (paired) without using a reference. The software expects as input RNA-Seq reads in fastq format
and outputs all assembled candidate transcripts in fasta format.


## Requirements

- **Operating System:** Linux
- **C++ Standard:** C++11 or higher
- **Boost**
  ```
   wget http://downloads.sourceforge.net/project/boost/boost/1.80.0/boost_1_80_0.tar.gz
   tar xfz boost_1_80_0.tar.gz
   rm boost_1_80_0.tar.gz
   cd boost_1_80_0
   ./bootstrap.sh --prefix=/usr/local --with-libraries=program_options,regex,filesystem,system
   export
   ./b2 install
   cd /home
   rm -rf boost_1_80_0
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
  ```


## Data are obtained and preprocessed

1. See [NCBI SRA Tools](https://github.com/ncbi/sra-tools.git) and are required to build this sra-tools.

2. Download and build fastp:
   ```
   $ cd ./bin
   $ ./prefetch SRR1942956                     # (Accession code)
   $ ./fastq-dump --split-files SRR1942956 
   ```
3. Download hisat2 and sequence alignment:
	
	a) download hisat2.
	   				
	   $ wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download						
	   $ unzip download
	   
	b) create index file:
	   	
	   $ hisat2-build ref_genome.fasta /index_dir/genome_index
	   
	c) sequence alignment:
	   	
	   $ hisat2 -x index_dir/genome_index -1 SRR1942956_1.fastq -2 SRR1942956_2.fastq -S  SRR1942956.sam --no-unal
	   
4. Download samtools and Format conversion
	
	See [Samtools](https://github.com/samtools/samtools/) and build it
	```
	$ samtools fastq -1 forward.fastq -2 ./reverse.fastq -0 /dev/null -s /dev/null -n SRR1942956.sam 
	```


## Installation

### 1. Install via Bioconda
```
conda install virdig
```


### 2. Compile from source

1. Building  VirDiG
	```
	$ cd VirDiG
	$ g++ -o virdig assemble2.cpp load_reads2.cpp transcript.cpp utility.cpp 
	```
2. Usage
	
** Required :

    Paired_end reads:
      -l <string>  : left reads.
      -r <string>  : right reads.

    -o <string>  : name of directory for output.

A typical command of VirDiG might be:

    $ ./virdig -l forward.fq -r reverse.fq -o output
  	
** Options :
    	
	-h                     : help information.
	-k <int>               : length of kmer, default: 31.
	-d <int>               : pair-end reads directions can be defined, 1 indicates that the directions of reads at both ends of paired-end are inversely complementary, and 2 indicates that the directions of reads at both ends are the same. default: 1.
	-t <int>               : number of threads, default 6.

    --non_canonical <int>  : whether to generate non-standard transcripts, 1 : true, 0 : false, default 0.
    --map_weight <float>   : paired-end reads are assigned paired node weights, recommended to be in the range of 0 to 1, default: 0.7.


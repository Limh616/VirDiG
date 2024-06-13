#include "transcript.h"
#include <vector>
#include "utility.h"
#include "common.h"
#include <fstream>
#include <boost/unordered_map.hpp>
#include <map>

// serialize/restore data as a text stream
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// headers privating serialize() function for many c++ classes
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

typedef std::pair<int, int> pair_t;
typedef std::pair<kmer_int_type_t, size_t>  kmer_occurence_pair_t;
typedef typename boost::unordered_map<kmer_int_type_t, size_t> Kmer_Hash;
typedef typename boost::unordered_map<std::string, size_t> Read_Hash;
typedef typename boost::unordered_map<std::string, std::vector<size_t>> Kmer_Reads_Hash;
typedef typename boost::unordered_map<kmer_int_type_t, size_t>::iterator kmer_hash_iterator_t;

typedef int node_idx_t;
class Node {
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
        ar& sequence;
        ar& parents;
        ar& children;
    }

public:
    Node() : sequence("") {};
    Node(const std::string& mysequence) {
        sequence = mysequence;
    }

    Node(const Node& node) {
        sequence = node.sequence;
        children = node.children;
        parents = node.parents;
    }

    void set_sequence(const std::string& myseq) {
        sequence = myseq;
    }

    std::string get_sequence() {
        return sequence;
    }

    bool add_child(node_idx_t child) {
        if (child < 0)
            return false;
        if (!children.empty()) {
            for (size_t i = 0; i < children.size(); ++i) {
                if (children[i] == child) // if exist already
                    return false;
            }
        }
        this->children.push_back(child);
        return true;
    }

    bool add_parent(node_idx_t parent) {
        if (parent < 0)
            return false;
        if (!parents.empty()) {
            for (size_t i = 0; i < parents.size(); ++i) {
                if (parents[i] == parent)
                    return false;
            }
        }
        this->parents.push_back(parent);
        return true;
    }

    bool is_child(node_idx_t child) {
        if (child < 0)
            return false;
        std::vector<node_idx_t>::iterator it = children.begin();
        for (; it != children.end(); ++it) {
            if (*it == child)
                return true;
        }
        return false;
    }

    bool is_parent(node_idx_t parent) {
        if (parent < 0)
            return false;
        std::vector<node_idx_t>::iterator it = parents.begin();
        for (; it != parents.end(); ++it) {
            if (*it == parent)
                return true;
        }
        return false;
    }

    bool delete_child(node_idx_t child) {
        if (child < 0)
            return false;
        std::vector<node_idx_t>::iterator it = children.begin();
        for (; it != children.end(); ++it) {
            if (*it == child)
                break;
        }
        if (it != children.end()) {
            children.erase(it);
            return true;
        }
        else {
            return false;
        }
    }

    bool delete_parent(node_idx_t parent) {
        if (parent < 0)
            return false;
        std::vector<node_idx_t>::iterator it = parents.begin();
        for (; it != parents.end(); ++it) {
            if (*it == parent)
                break;
        }
        if (it != parents.end()) {
            parents.erase(it);
            return true;
        }
        else {
            return false;
        }
    }

    void clear_children() {
        children.clear();
    }

    void clear_parents() {
        parents.clear();
    }

    void clear() {
        sequence.clear();
        children.clear();
        parents.clear();
    }

public:
    std::string  sequence;
    std::vector<node_idx_t> parents;
    std::vector<node_idx_t> children;

};

std::vector<size_t> used_reads;

size_t node_size = 0;  // the number of nodes in this graph
int L = 101;  // g_kmer_length = 31;
int min_error_count = 3;
float min_error_rate = 0.05;
int max_rep_trunk_threshold = g_kmer_length + L;
float min_ratio_non_error = 0.05f;
float min_seed_entropy = 1.5f;
int trunk_head_threshoid = 300;


bool sortBySecond(const kmer_occurence_pair_t& pair1, const kmer_occurence_pair_t& pair2) {
    return pair1.second > pair2.second;
}



bool read_is_used(size_t i) {
    return (std::find(used_reads.begin(), used_reads.end(), i) != used_reads.end());
}

bool add_used_read(size_t i) {

    used_reads.push_back(i);
    return true;
}


size_t get_kmer_count(Kmer_Hash& kmer_hash, kmer_int_type_t kmer_val) {

    kmer_hash_iterator_t it = kmer_hash.find(kmer_val);

    if (it != kmer_hash.end()) {

        return(it->second);
    }
    else
        return 0;
}

bool exists(Kmer_Hash& kmer_hash, const kmer_int_type_t kmer_val) {

    kmer_hash_iterator_t it = kmer_hash.find(kmer_val);

    if (it != kmer_hash.end())
        return true;
    else
        return false;
}

bool exists(Kmer_Hash& kmer_hash, const std::string& kmer) {

    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    return (exists(kmer_hash, kmer_val));
}


void get_read_hash(Read_Hash& read_hash, std::vector<std::string>& Read, std::vector<std::string>& read_list) {

    std::cerr << "Beginning read hash ..." << std::endl;
    time_t beg = time(NULL);
    if (Read.empty())
        return;
   
    std::vector<std::string> f_reads_list;
    std::vector<std::string> r_reads_list;
    int reads_size = Read.size() / 2;

    size_t data_size = Read.size();
    for (size_t i = 0; i < data_size; ++i) {
        const std::string& sequence = Read[i];
        if (sequence == "")
            continue;
        if (read_hash.count(sequence) == 0) {
            if (i < reads_size) {
                r_reads_list.push_back(sequence);
                f_reads_list.push_back(Read[reads_size + i]);
            }
            else {
                f_reads_list.push_back(sequence);
                r_reads_list.push_back(Read[i - reads_size]);
            }
        }
        read_hash[sequence]++;
    }
    
    for (int i = 0; i < r_reads_list.size(); i++) {
        read_list.push_back(r_reads_list[i]);
    }
    for (int j = 0; j < f_reads_list.size(); j++) {
        read_list.push_back(f_reads_list[j]);
    }

    time_t end = time(NULL);
   

    std::cerr << "read hash finished, get " << read_hash.size() << " reads!"  << "read list finished, get " << read_list.size()
        << "paired reads!  (elapsed time: " << (end - beg) << " s)" << std::endl;
}

void get_kmer_hash(Kmer_Hash& kmer_hash, Read_Hash& read_hash) {

    std::cerr << "Beginning kmer hash ..." << std::endl;
    time_t beg = time(NULL);

    boost::unordered_map<std::string, size_t>::iterator it;
    for (it = read_hash.begin(); it != read_hash.end(); ++it) {
        const std::string& sequence = it->first;
        
       // if (sequence == "")
       //     continue;
        for (size_t j = 0; j <= sequence.length() - g_kmer_length; ++j) {
            const std::string& kmer = sequence.substr(j, g_kmer_length);
            if (kmer.find("N") != std::string::npos)
                continue;
            kmer_int_type_t kmer_val = kmer_to_intval(kmer, g_kmer_length);
            kmer_hash[kmer_val] += it->second;
        }
    }
    time_t end = time(NULL);
    if (kmer_hash.empty()) {
        std::cout << "kmer_hash is empty." << std::endl;
    }
    /*
    std::string fileName = out_dir + "kmer_hash.txt";
    std::ofstream outfile(fileName);
    // std::ofstream outfile("/home/bioinfo/limh/paper/data/SRR1942956/kmer_hash.txt");
    for (const auto& pair : kmer_hash) {
        outfile << pair.first << " " << pair.second << std::endl;
    }
    outfile.close();
   //   */
    std::cerr << "Kmer hash finished, get " << kmer_hash.size()
        << " kmer! (elapsed time: " << (end - beg) << " s)" << std::endl;
}

void get_forward_candidates(Kmer_Hash& kmer_hash, kmer_int_type_t seed_kmer, std::vector<kmer_occurence_pair_t>& candidates) {

    candidates.clear(); // clear the vector
    kmer_int_type_t forward_prefix
        = (seed_kmer << (33 - g_kmer_length) * 2) >> (32 - g_kmer_length) * 2;

    for (kmer_int_type_t i = 0; i < 4; ++i) {

        kmer_occurence_pair_t candidate;
        candidate.first = forward_prefix | i;
        candidate.second = get_kmer_count(kmer_hash, candidate.first);

        if (candidate.second) {
            candidates.push_back(candidate);
        }
    }
}

void get_reverse_candidates(Kmer_Hash& kmer_hash, kmer_int_type_t seed_kmer, std::vector<kmer_occurence_pair_t>& candidates) {

    candidates.clear();
    kmer_int_type_t reverse_suffix = seed_kmer >> 2;

    for (kmer_int_type_t i = 0; i < 4; ++i) {

        kmer_occurence_pair_t candidate;
        candidate.first = (i << (g_kmer_length * 2 - 2)) | reverse_suffix;
        candidate.second = get_kmer_count(kmer_hash, candidate.first);

        if (candidate.second) {
            candidates.push_back(candidate);
        }
    }
}

void remove_erroneous(Kmer_Hash& kmer_hash, float min_ratio_non_error) {
    std::cout << "Beginning remove erroneous ..." << std::endl;
    time_t beg = time(NULL);
    kmer_hash_iterator_t it;
    Kmer_Hash deletion_hash;
   
    for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {

        kmer_int_type_t kmer_val = it->first;
        std::vector<kmer_occurence_pair_t> f_candidates;
        get_forward_candidates(kmer_hash, kmer_val, f_candidates);
        std::vector<kmer_occurence_pair_t> r_candidates;
        get_reverse_candidates(kmer_hash, kmer_val, r_candidates);
        sort(f_candidates.begin(), f_candidates.end(), sortBySecond);
        sort(r_candidates.begin(), r_candidates.end(), sortBySecond);
        int f_dominant_count = 0;
        for (unsigned int i = 0; i < f_candidates.size(); ++i) {
            if (f_candidates[i].second == 1) {
               // deletion_hash[f_candidates[i].first]++;
                kmer_hash.erase(f_candidates[i].first);

            }

            if (f_candidates[i].second > 1) {
                int candidate_count = f_candidates[i].second;
                if (f_dominant_count == 0) {
                    f_dominant_count = candidate_count;   //同组最大的kmer的丰度
                }
                else if ((float)candidate_count / f_dominant_count < min_ratio_non_error) {
                    deletion_hash[f_candidates[i].first]++;
                    kmer_hash.erase(f_candidates[i].first);
                }
            }
        }
        int r_dominant_count = 0;
        for (unsigned int i = 0; i < r_candidates.size(); ++i) {
            if (r_candidates[i].second == 1) {
              //  deletion_hash[r_candidates[i].first]++;
                kmer_hash.erase(r_candidates[i].first);
              
            }

            if (r_candidates[i].second > 1) {
                int candidate_count = r_candidates[i].second;
                if (r_dominant_count == 0) {
                    r_dominant_count = candidate_count;   //同组最大的kmer的丰度
                }
                else if ((float)candidate_count / r_dominant_count < min_ratio_non_error) {
                    deletion_hash[r_candidates[i].first]++;
                    kmer_hash.erase(r_candidates[i].first);
                }
            }
        }

    }
    
    
    time_t end = time(NULL);
    std::cout << "kmer count after errors deletion : " << kmer_hash.size()  
       << " (elapsed time: " << (end - beg) << " s)"  <<std::endl;

}

float compute_entropy(kmer_int_type_t kmer) {
    char counts[] = { 0, 0, 0, 0 };
    for (unsigned int i = 0; i < g_kmer_length; ++i) {
        int c = kmer & 3;
        kmer = kmer >> 2;
        counts[c]++;
    }
    float entropy = 0;
    for (unsigned int i = 0; i < 4; i++) {
        float prob = (float)counts[i] / g_kmer_length;
        if (prob > 0) {
            float val = prob * log(1.f / prob) / log(2.f);
            entropy += val;
        }
    }
    return(entropy);
}

void get_seed_kmer(Kmer_Hash& kmer_hash, std::vector<kmer_occurence_pair_t>& seed_hash) {
    std::cout << "Beginning get_seed_kmer ..." << std::endl;
    time_t beg = time(NULL);
    kmer_hash_iterator_t it;
    for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {
        if (compute_entropy(it->first) < min_seed_entropy)
            continue;
        seed_hash.push_back(std::make_pair(it->first,it->second));
    }
    
   // std::sort(seed_hash.begin(), seed_hash.end(), sortBySecond) ;
    time_t end = time(NULL);
    std::cerr << "get seed kmer finished, get " << seed_hash.size()
        << " seed kmers! (elapsed time: " << (end - beg) << " s)" << std::endl;
}



std::string forward_extend_contig(Kmer_Hash& kmer_hash, kmer_int_type_t kmer_val, Kmer_Hash& used_kmers_) {

    kmer_int_type_t intval = kmer_val;
    std::string str = intval_to_kmer(intval, g_kmer_length);
    std::vector<kmer_occurence_pair_t> candidates;
  
    while (1) {
        get_forward_candidates(kmer_hash, intval, candidates);
        if (candidates.empty()) break;  

        std::sort(candidates.begin(), candidates.end(), sortBySecond); //降序

        if (used_kmers_.find(candidates[0].first) != used_kmers_.end()) {  
            break;
        }
        kmer_int_type_t candidate = candidates[0].first;
        size_t count = candidates[0].second;
        used_kmers_[candidate] = count;
        
        int base_num = candidate & 3ll;
        char base = int_to_base(base_num);
        str += base;
        intval = candidate;
    }
 
    return str;
}

std::string forward_extend_contig(Kmer_Hash& kmer_hash, std::string contig) {

    std::string rkmer = contig.substr(contig.length() - g_kmer_length);
    kmer_int_type_t intval = kmer_to_intval(rkmer);
 
    std::string str = contig;
    std::vector<kmer_occurence_pair_t> candidates;
  
    while (1) {
        get_forward_candidates(kmer_hash, intval, candidates);
        if (candidates.empty()) break;

        std::sort(candidates.begin(), candidates.end(), sortBySecond); //降序
        kmer_int_type_t candidate = candidates[0].first;
        std::string kmer = intval_to_kmer(candidate, g_kmer_length);
        if (str.find(kmer) != std::string::npos) {
            break;
        }

        size_t count = candidates[0].second;
        
        int base_num = candidate & 3ll;
        char base = int_to_base(base_num);
        str += base;
        intval = candidate;
    }
    return str;
}


std::string reverse_extend_contig(Kmer_Hash& kmer_hash, kmer_int_type_t kmer_val, Kmer_Hash& used_kmers_) {

    kmer_int_type_t intval = kmer_val;
    std::string str = intval_to_kmer(intval, g_kmer_length);
    std::vector<kmer_occurence_pair_t> candidates;
 
    while (true) {
        get_reverse_candidates(kmer_hash, intval, candidates);
        if (candidates.empty()) 
            break;

        std::sort(candidates.begin(), candidates.end(), sortBySecond); //降序
        if (used_kmers_.find(candidates[0].first) != used_kmers_.end()) { 
            break;
        }
        kmer_int_type_t candidate = candidates[0].first;
        size_t count = candidates[0].second;
        used_kmers_[candidate] = count; 
        int base_num = (candidate >> (g_kmer_length * 2 - 2)) & 3ll;
        char base = int_to_base(base_num);
        str = base + str;
        intval = candidate;
    }
    return str;
}
std::string reverse_extend_contig(Kmer_Hash& kmer_hash, std::string contig) {
    std::string lkmer = contig.substr(0, g_kmer_length);
    kmer_int_type_t intval = kmer_to_intval(lkmer);
    std::string str = contig;

    std::vector<kmer_occurence_pair_t> candidates;
    while (true) {
        get_reverse_candidates(kmer_hash, intval, candidates);
        if (candidates.empty())
            break;

        std::sort(candidates.begin(), candidates.end(), sortBySecond); //降序
        kmer_int_type_t candidate = candidates[0].first;
        std::string kmer = intval_to_kmer(candidate,g_kmer_length);
        if (str.find(kmer)!=std::string::npos) { 
            break;
        }
        size_t count = candidates[0].second;
        int base_num = (candidate >> (g_kmer_length * 2 - 2)) & 3ll;
        char base = int_to_base(base_num);
        str = base + str;
        intval = candidate;

    }
    return str;
}



std::string refine_contig(std::string contig, Kmer_Hash& kmer_hash, Kmer_Hash& used_kmers_) {
    
    int length = contig.length();
    if (length < 3 * g_kmer_length)
        return contig;
    // reverse
    std::vector<kmer_occurence_pair_t> candidates;
    for (int i = 0; i < 2 * g_kmer_length; ++i) {
        const std::string& kmer = contig.substr(i, g_kmer_length);
        kmer_int_type_t intval = kmer_to_intval(kmer);
        get_reverse_candidates(kmer_hash,intval, candidates);
        if (candidates.size() <= 1)
            continue;
        bool hit_point = false;
        for (size_t j = 0; j < candidates.size(); ++j) {
            if (used_kmers_.find(intval) == used_kmers_.end()) {
                hit_point = true;
                break;
            }
        }
        if (hit_point) {
            const std::string& sequence = reverse_extend_contig(kmer_hash, intval,used_kmers_);
            if ((int)sequence.length() > 3 * g_kmer_length) {
                contig = sequence + contig.substr(i + g_kmer_length);
                break;
            }
        }
    }

    // forward
    length = contig.length();
    for (int i = length - 3 * g_kmer_length; i < length - g_kmer_length; ++i) {
        const std::string& kmer = contig.substr(i, g_kmer_length);
        kmer_int_type_t intval = kmer_to_intval(kmer);
        get_forward_candidates(kmer_hash,intval, candidates);
        if (candidates.size() <= 1)
            continue;
        bool hit_point = false;
        for (size_t j = 0; j < candidates.size(); ++j) {
            if (used_kmers_.find(intval) == used_kmers_.end()) {
                hit_point = true;
                break;
            }
        }
        if (hit_point) {
            const std::string& sequence = forward_extend_contig(kmer_hash, intval,used_kmers_);
            if ((int)sequence.length() > 3 * g_kmer_length) {
                contig = contig.substr(0, i) + sequence;
                break;
            }
        }
    }

    return contig;
}


void extend_contigs(Kmer_Hash kmer_hash,  std::vector<std::string>& contigs_list) {

    time_t beg = time(NULL);
    std::cerr << "extend_contigs ..." << std::endl;
    
    Kmer_Hash used_kmers_;
    std::vector<kmer_occurence_pair_t>  seed_kmers;
    get_seed_kmer(kmer_hash, seed_kmers);
    for (int i = 0; i < seed_kmers.size(); i++) {
        kmer_int_type_t kmer_val = seed_kmers[i].first;

        if (used_kmers_.find(kmer_val) == used_kmers_.end() ) { 
            used_kmers_.emplace(kmer_val, kmer_hash[kmer_val]);
            std::string left = reverse_extend_contig(kmer_hash, kmer_val, used_kmers_);
            std::string right = forward_extend_contig(kmer_hash, kmer_val, used_kmers_);
            std::string contig = left + right.substr(g_kmer_length);
           // contig = refine_contig(contig, kmer_hash,used_kmers_);
            if (contig.length() > 3 * g_kmer_length) {
                ///*
                left = reverse_extend_contig(kmer_hash, contig);
                right = forward_extend_contig(kmer_hash,contig);
                if (right.length() == contig.length())
                    contig = left;
                else
                    contig = left + right.substr(contig.length());
               // */
                contigs_list.push_back(contig);

            }
            /*
            else {
                for (int i = 0; i < contig.length(); i++) {
                    std::string kmer = contig.substr(i, g_kmer_length);
                    kmer_int_type_t k_value = kmer_to_intval(kmer);
                    used_kmers_.erase(k_value);
                }
            }
            */
        }

    }
    

    time_t end = time(NULL);
  ///*
    int i = 0;
    std::string fileName = out_dir + "contigs.fasta";
    std::ofstream outputFile(fileName);
    
    if (outputFile.is_open()) {
        for (const auto& pair : contigs_list) {
            outputFile << "> " << i++ << std::endl;
            outputFile << pair << "\n";
        }
    }
    outputFile.close();
//  */
    std::cerr << "contigs  finished, get " << contigs_list.size()
        << " contigs! (elapsed time: " << (end - beg) << " s)" << std::endl;
}

bool is_similar( const std::string& str1, const std::string& str2) {
    if (str1.length() == 0 || str2.length() == 0) {
        return false;
    }
    int mismatch = 0;
    int kmer_length = g_kmer_length;
    std::string first_kmer = str2.substr(0, kmer_length);
    size_t start = str1.find(first_kmer);
    std::string last_kmer = str2.substr(str2.length() - kmer_length);
    size_t end = str1.find(last_kmer);
   //   std::cout << "569" << std::endl;

    if (start == std::string::npos && end == std::string::npos) {
        kmer_length = g_kmer_length - 4;
        first_kmer = str2.substr(0, kmer_length);
        start = str1.find(first_kmer);
        last_kmer = str2.substr(str2.length() - kmer_length);
        end = str1.find(last_kmer);
        if (start == std::string::npos && end == std::string::npos)
            return false;
    }
       
    if (start + kmer_length == str1.length() || end == 0)
        return true;
    //   std::cout << "738" << std::endl;
    if (start != std::string::npos) {
        const std::string& str3 = str1.substr(start + kmer_length);
        const std::string& str4 = str2.substr(kmer_length);
        int length = (str3.length() < str4.length()) ? str3.length() : str4.length();

        for (int i = 0; i < length; ++i) {
            if (str3[i] != str4[i])
                mismatch++;
        }
        if (static_cast<float>(mismatch) / length <= 0.05 || mismatch <= 2) {
            return true;
        }
        else {
            return false;
        }
    }
    //  std::cout << "754" << std::endl;
    if (end != std::string::npos) {
        const std::string& str3 = str1.substr(0, end);
        const std::string& str4 = str2.substr(0, str2.length() - kmer_length);
        int length_1 = str3.length() - 1;
        int length_2 = str4.length() - 1;
        int length = (length_1 < length_2) ? length_1 : length_2;
        while (length_1 >= 0 && length_2 >= 0) {
            if (str1[length_1] != str2[length_2])
                mismatch++;
            length_1--;
            length_2--;
        }
        if (static_cast<float>(mismatch) / (length + 1) <= 0.05 || mismatch <= 2) {
            return true;
        }
        else {
            return false;
        }
    }

}

float get_str_coverage(std::string str, Kmer_Hash& kmer_hash,bool flag) {
    int cov_size = str.length() - g_kmer_length + 1;
    std::vector<int> cov_v(cov_size);
    if (!flag)
        str = revcomp(str);
    for (int i = 0; i < cov_size; i++) {
        std::string kmer = str.substr(i, g_kmer_length);
        kmer_int_type_t intval = kmer_to_intval(kmer);
        cov_v[i] = kmer_hash[intval];
    }
    int first = 0;
    int last = cov_size - 1;
   // /*
    sort(cov_v.begin(), cov_v.end());
    int quantile = static_cast<int>(0.05 * cov_size + 0.5);
    
    if (quantile > 0) {
        first = quantile;
        last = cov_size - quantile;
    }
   // */
    int sum = 0;
    // if (first > last)
    //     break;
    for (; first <= last; ++first) {
        sum += cov_v[first];
    }
    float cov = sum * 1.0 / (cov_size - 2 * quantile);
  //  std::cout << cov << std::endl;
  //  float cov = static_cast<float>(sum * 1.0 / cov_size );
    return cov;
}





node_idx_t add_node(Node& node, std::vector<Node>& node_set) {

    node_set.push_back(node);
    return (node_size++);
}



void reverse_similar(const std::string& str1, const std::string& str2, int& i, int& j) {
    i = str1.length() - 1;
    j = str2.length() - 1;
    while (i >= 0 && j >= 0 && str1[i] == str2[j]) {
        i--;
        j--;
    }
    if (str1.length() - i < L) {
        i = -1;
    }
}

int editDistance(const std::string& str1, const std::string& str2) {
    int m = str1.length();
    int n = str2.length();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            if (i == 0) {
                dp[i][j] = j;
            }
            else if (j == 0) {
                dp[i][j] = i;
            }
            else if (str1[i - 1] == str2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            }
            else {
                dp[i][j] = 1 + std::min({ dp[i][j - 1], dp[i - 1][j], dp[i - 1][j - 1] });
            }
        }
    }
    return dp[m][n];
}

bool forward_similar(const std::string& str1, const std::string& str2 ,int n ) {
    int i = 0;
    int j = 0;
    int len = (str1.length() < str2.length()) ? str1.length() : str2.length();
    int mismatch = n;
    std::string str3 = str1.substr(0, len);
    std::string str4 = str2.substr(0, len);
    int distance = editDistance(str3, str4);
    return distance <= mismatch;

}

bool reverse_different2(const std::string& str1, const std::string& str2) {
    int i = 0;
    int j = 0;
    int len = (str1.length() < str2.length()) ? str1.length() : str2.length();
    int mismatch = (int(len * 0.5 + 0.5) > 20) ? int(len * 0.5 + 0.5) : 20;
    std::string str3 = str1.substr(str1.length() - len);
    std::string str4 = str2.substr(str2.length() - len);
    int distance = editDistance(str3, str4);
    return distance > mismatch;

}


bool reverse_different(const std::string& str1, const std::string& str2) {
    int length_1 = str1.length() - 1;
    int length_2 = str2.length() - 1;
    int match = 0;
    int length = (length_1 < length_2) ? length_1 : length_2;
    while (length_1 >= 0 && length_2 >= 0) {
        if (str1[length_1] == str2[length_2])
            match++;
        length_1--;
        length_2--;
    }
    if (static_cast<float>(match) / (length + 1) > 0.7 || match > 40 ) {
        return true;
    }
    else {
        return false;
    }

}

std::vector<size_t> mismatch_pos(const std::string& str1, const std::string& str2) {
    std::vector<size_t> mismatch_set;
    for (size_t i = 0; i < str1.length(); i++) {
        if (str1[i] != str2[i]) {
            mismatch_set.push_back(i);
            if (mismatch_set.size() > static_cast<int>(str1.length() * 0.01 + 0.5))
                break;
        }
    }
    return mismatch_set;
}


int read_num(std::string read, std::vector<std::string> Read) {
    int i;
    for (i = 0; i < Read.size(); ++i) {
        if (Read[i] == read)
            break;
    }
    if (i < Read.size())
        return i;
    else
        return -1;
}

bool list_contain_str(const std::string& str, const std::vector<std::string>& list){
    int i;
    for ( i = 0; i < list.size(); i++) {
        if (list[i] == str)
            break;
    }
    if (i < list.size())
        return true;
    else
        return false;
}


size_t find_longest_read(std::vector<std::string>& list) {
    if (list.empty()) {
        return -1;
    }
    int max_langth = 0;
    int num = -1;
    size_t date_size = list.size();
    for (int i = 0; i < date_size; i++) {
        if (read_is_used(i))
            continue;
        std::string str = list[i];
        if (str.length() > max_langth) {
            num = i;
            max_langth = str.length();
        }
    }
    return num;
}

void find_read_contain_kmer(std::vector<std::string>& read_list, std::string candidate, std::vector<size_t>& candidate_reads) {
    for (size_t i = 0; i < read_list.size(); ++i) {
        if (read_list[i].find(candidate) != std::string::npos)
            candidate_reads.push_back(i);
    }
}

void set_child(std::vector<Node>& node_set) {

    for (size_t i = 0; i < node_set.size(); ++i) {
        if (!node_set[i].children.empty())
            node_set[i].children.clear();
    }

    // reset parents of each node 
    for (size_t i = 0; i < node_set.size(); ++i) {
        std::vector<node_idx_t>::const_iterator it;
        for (it = node_set[i].parents.begin();
            it != node_set[i].parents.end(); ++it) {
            node_set[*it].add_child(i);
        }
    }
}

void set_parents(std::vector<Node>& node_set) {

    for (size_t i = 0; i < node_set.size(); ++i) {
        if (!node_set[i].parents.empty())
            node_set[i].parents.clear();
    }

    // reset parents of each node 
    for (size_t i = 0; i < node_set.size(); ++i) {
        std::vector<node_idx_t>::const_iterator it;
        for (it = node_set[i].children.begin();
            it != node_set[i].children.end(); ++it) {
            node_set[*it].add_parent(i);
        }
    }
}

//是否存在多条reads支持分支
bool support_edge_decision(std::string str, std::vector<size_t> reads, std::vector<std::string> read_list) {

    std::vector<int> similarity;
    std::string read;
    for (int i = 0; i < reads.size(); ++i) {
        read = read_list[reads[i]];

        int kmer_len = g_kmer_length;
        std::string first_kmer = read.substr(0, kmer_len);
        std::string last_kmer = read.substr(read.length() - kmer_len);
        int start = str.find(first_kmer);
        int end = str.find(last_kmer);
        if (start == std::string::npos && end == std::string::npos) {
            kmer_len = g_kmer_length - 6;
            first_kmer = read.substr(0, kmer_len);
            last_kmer = read.substr(read.length() - kmer_len);
            start = str.find(first_kmer);
            end = str.find(last_kmer);
            if (start == std::string::npos && end == std::string::npos)
                continue;
        }
        int match = 0;
        if (start != std::string::npos) {
            std::string new_str = str.substr(start);
            int j = 0;
            int len = (read.length() < new_str.length()) ? read.length() : new_str.length();
            while (j < len) {
                if (read[j] == new_str[j])
                    match++;
                j++;
            }
        }
        else if (end != std::string::npos) {
            std::string new_str = str.substr(0, end + kmer_len);
            int l = read.length() - 1, k = new_str.length() - 1;
            while (l >= 0 && k >= 0) {
                if (read[l] == new_str[k])
                    match++;
                l--;
                k--;
            }
        }
        similarity.push_back(match);

    }
    std::sort(similarity.begin(), similarity.end(), [](int a, int b) { return a > b; });
    /*
    for (size_t i = 0; i < 4; ++i) {
        std::cout << similarity[i] << " ";
    }
    std::cout << std::endl;
   // */
    if (similarity.size() >= 4 && similarity[3] > read.length() * 0.95)
        return true;
    else
        return false;

}


void  add_reverse_branch(int p, std::vector<std::string>& read_list, Kmer_Hash kmer_hash, std::vector<Node>& node_set,bool& flag) {

    std::cout << "Beginning add branch ..." << std::endl;
    time_t beg = time(NULL);
    std::string fileName = out_dir + "branch.fasta";
    std::ofstream ofile(fileName);
    std::string trunk = node_set[p].sequence;
    
    for (int i = 0 ; i <= trunk.length() - g_kmer_length; i++) {
        std::string kmer = trunk.substr(i, g_kmer_length);
        if (kmer.substr(0, 3) != "ATG") //
            continue;

        std::vector<size_t> candidate_reads;
        find_read_contain_kmer(read_list, kmer, candidate_reads);
        for (int j = 0; j < candidate_reads.size(); j++) {
            std::string candidate_read = read_list[candidate_reads[j]];
            int pos = candidate_read.find(kmer);
            if (pos < g_kmer_length || pos > candidate_read.length() - g_kmer_length)
                continue;
            std::string first_kmer = candidate_read.substr(0, g_kmer_length);
            std::string last_kmer = candidate_read.substr(candidate_read.length() - g_kmer_length);
            int start = trunk.find(first_kmer);
            int end = trunk.find(last_kmer);
            if ( end == std::string::npos  || end - start < candidate_read.length() + 3)  //
                continue;
            int pos2 = node_set[p].sequence.find(kmer);
            if (pos2 == std::string::npos || pos2 <= 3)
                continue;
            std::string r_str = candidate_read.substr(0, pos);
            std::string f_str = candidate_read.substr(pos);
            int len2 = (i < r_str.length()) ? i : r_str.length();
            bool r_flag = reverse_different2(r_str, trunk.substr(i - len2, len2));
            if (!r_flag)
                continue;
            int len1 = (trunk.length() - i < f_str.length()) ? trunk.length() - i : f_str.length();
            bool f_flag = forward_similar(f_str, trunk.substr(i , len1), 2);
            if (!f_flag)
                continue;
            int junction = pos - 0.5 * (g_kmer_length + 1);
            std::string jun_kmer = candidate_read.substr(junction, g_kmer_length);
            if (!flag)
                jun_kmer = revcomp(jun_kmer);
            kmer_int_type_t jun_val = kmer_to_intval(jun_kmer);
           
            if ( kmer_hash[jun_val] < 3)
                continue;

            if (f_flag && r_flag) {
                std::vector<size_t> jun_reads;
                find_read_contain_kmer(read_list, jun_kmer, jun_reads);
                std::string ref_str =  candidate_read + trunk.substr(i + f_str.length(), 20);
                if (!support_edge_decision(ref_str, jun_reads, read_list))
                    continue;

                // std::cout  << i << " kmer: " << kmer << " jun_kmer_num: " << kmer_hash[jun_val] << std::endl;
                ofile << ">" << "pos:  " << i << std::endl;
                ofile << candidate_read << std::endl;
                // std::cout << candidate_read << std::endl;

                Node node1, node2;
                node1 = node_set[p].sequence.substr(0, pos2);
                node_set[p].sequence = node_set[p].sequence.substr(pos2);
                node_idx_t q1 = add_node(node1, node_set);
                node2 = candidate_read.substr(0, pos);
                node_idx_t q2 = add_node(node2, node_set);
                node_set[q1].parents = node_set[p].parents;
                node_set[p].parents.clear();
                node_set[p].add_parent(q1);
                node_set[p].add_parent(q2);

                break;
            }
        }
    }
    Node node;
    node = node_set[p].sequence;
    node_idx_t q = add_node(node, node_set);
    node_set[q].parents = node_set[p].parents;
    node_set[p].clear();

    set_child(node_set);
    time_t end = time(NULL);
    std::cout << "Finished add branch. get " << node_set.size() << " nodes. (elapsed time : " << (end - beg) << " s)" << std::endl;
}


int different_point(std::string str1, std::string str2) {
    int i = 0;
    while (i < str1.length()-1) {
        if (str1[i] != str2[i] && str1[i+1] != str2[i+1])
            break;
        i++;
    }
    return i;
}

void get_read_list(Read_Hash& read_hash, std::vector<std::string>& read_list,bool flag) {
    boost::unordered_map<std::string, size_t>::iterator it;
    if(flag)
        for (it = read_hash.begin(); it != read_hash.end(); ++it) {
            read_list.push_back(it->first);
        }
    else if (!flag)
        for (it = read_hash.begin(); it != read_hash.end(); ++it) {
            read_list.push_back( revcomp(it->first));
        }
}


void check_reverse_branch(int p, std::vector<std::string>& read_list, Kmer_Hash kmer_hash, std::vector<Node>& node_set,bool flag) {

    std::cout << "Beginning add branch ..." << std::endl;
    time_t beg = time(NULL);

    boost::unordered_map<size_t, size_t> used_read_hash;
    std::string trunk = node_set[p].sequence;
   
    std::string fileName = out_dir + "branch.fasta";
    std::ofstream ofile(fileName);

    for (int i = 0.6 * trunk.length(); i <= trunk.length() - g_kmer_length; i++) {
        std::string kmer = trunk.substr(i, g_kmer_length);
        if (kmer.substr(0, 3) != "ATG") //
            continue;

        std::vector<size_t> candidate_reads;
        find_read_contain_kmer(read_list, kmer, candidate_reads);
        for (int j = 0; j < candidate_reads.size(); j++) {

            std::string candidate_read = read_list[candidate_reads[j]];
            int pos1 = candidate_read.find(kmer);
            int len = 0.3 * candidate_read.length();
            if (pos1 < len || pos1 > candidate_read.length() - len)
                continue;
            int kmer_len = g_kmer_length;
            std::string first_kmer = candidate_read.substr(0, kmer_len);
            std::string last_kmer = candidate_read.substr(candidate_read.length() - kmer_len);
            int start = trunk.find(first_kmer);
            int end = trunk.find(last_kmer);
            if (start == std::string::npos || end == std::string::npos || start > 300 || end - start < candidate_read.length() + 3)  //
                continue;
            // /*
            int pos3 = different_point(trunk.substr(start, candidate_read.length()), candidate_read);
            if (pos3 > pos1)
                pos3 = pos1;
            // */
            int pos2 = node_set[p].sequence.find(kmer);
            int pos4 = node_set[p].sequence.find(candidate_read.substr(pos3, g_kmer_length));

            if (pos4 == std::string::npos || pos4 <= 2 || pos2 <= 2)
                continue;

            std::string f_str = candidate_read.substr(pos3);
            int len1 = (trunk.length() - i < f_str.length()) ? trunk.length() - i : f_str.length();
            bool f_flag = forward_similar(f_str, trunk.substr(i - (pos1 - pos3), len1), 2);  // 956 ：0

            int junction = pos3 - 0.5 * (g_kmer_length + 1);
            std::string jun_kmer = candidate_read.substr(junction, g_kmer_length);
            if (!flag)
                jun_kmer = revcomp(jun_kmer);
            kmer_int_type_t jun_val = kmer_to_intval(jun_kmer);

            
            // std::cout << f_flag << std::endl; 
            
        //    std::cout << i << " kmer: " << kmer << " jun_kmer_num: " << kmer_hash[jun_val] << " " << f_flag << std::endl;
            if (f_flag && kmer_hash[jun_val] > 3) {  //&& l_flag1 && !l_flag2  

               // std::cout << i << " ";
                std::vector<size_t> jun_reads;
                if (!flag)
                    jun_kmer = revcomp(jun_kmer);
                find_read_contain_kmer(read_list, jun_kmer, jun_reads);
                std::string ref_str = trunk.substr(0, start) + candidate_read + trunk.substr(end + kmer_len, 40);
                if (!support_edge_decision(ref_str, jun_reads, read_list))
                    continue;

              //   std::cout  << i << " kmer: " << kmer << " jun_kmer_num: " << kmer_hash[jun_val] << std::endl;

                // std::cout << candidate_read << std::endl;
                ofile << ">" << "pos:  " << i << std::endl;
                ofile << candidate_read << std::endl;

                Node node1, node2;
                node1 = node_set[p].sequence.substr(0, pos2);
                node_set[p].sequence = node_set[p].sequence.substr(pos2);
                node_idx_t q1 = add_node(node1, node_set);
                node2 = candidate_read.substr(0, pos1);
                node_idx_t q2 = add_node(node2, node_set);
                node_set[q1].parents = node_set[p].parents;
                node_set[p].parents.clear();
                node_set[p].add_parent(q1);
                node_set[p].add_parent(q2);

                break;
            }

        }

    }
    Node node;
    node = node_set[p].sequence;
    node_idx_t q = add_node(node, node_set);
    node_set[q].parents = node_set[p].parents;
    node_set[p].clear();

    set_child(node_set);
    time_t end = time(NULL);
    std::cout << "Finished add branch  (elapsed time: " << (end - beg) << " s)" << std::endl;


}

void printNodes(std::vector<Node> node_set,int p, int depth = 0) {
    
    for (int i = 0; i < depth; ++i) {
        std::cout << "  ";
    }
    std::cout << node_set[p].sequence << std::endl;
    
    for (int j = 0; j < node_set[p].children.size(); j++) {
        printNodes(node_set,j, depth + 1);
    }
    
}

int find_node_contain_kmer(std::vector<Node>& node_set, std::string kmer) {
    int i;
    for ( i = 0; i < node_set.size(); ++i) {
        std::string sequence = node_set[i].sequence;
        if (sequence.find(kmer) != std::string::npos)
            break;
    }
    if (i < node_set.size())
        return i;
    else
        return -1;
}

int find_jun_contain_kmer(std::vector<Node>& node_set, std::string kmer) {
    int i;
    for (i = 0; i < node_set.size(); ++i) {
        int len = (node_set[i].sequence.length() < g_kmer_length - 1) ? 0 : g_kmer_length - 1;
        std::string sequence = node_set[i].sequence.substr(node_set[i].sequence.length()-len );
        while (sequence.length() >= 2 * g_kmer_length - 2) {
            if (node_set[i].children.size() == 0)
                break;
            int j = node_set[i].children[0];
            sequence += node_set[j].sequence;
            i = j;
        }
        sequence = sequence.substr(0, 2 * g_kmer_length - 2);
        if (sequence.find(kmer) != std::string::npos)
            break;
    }
    if (i < node_set.size())
        return i;
    else
        return -1;
}

int find_node_map_read(std::vector<Node>& node_set, std::string read) {
    std::string kmer = read.substr(0, g_kmer_length);
    int node_num ;
    for (int i = 0; i < node_set.size(); ++i) {
        std::string sequence = node_set[i].sequence;
        if (sequence.find(kmer) != std::string::npos) {
            node_num = i;
            std::string sequence = node_set[node_num].sequence;
            int pos = sequence.find(kmer);
            sequence = sequence.substr(pos);
           
            if (sequence.length() < read.length()) {
                while (sequence.length() < read.length()) {
                    if (node_set[node_num].children.size() == 0)
                        return -1;
                    int child = node_set[node_num].children[0];
                    sequence = sequence + node_set[child].sequence;
                    node_num = child;
                }
                sequence = sequence.substr(0, read.length());
                if (forward_similar(sequence, read,1))
                    return node_num;
            }
        }
    }
    for (int i = 0; i < node_set.size(); ++i) {
        int len = (node_set[i].sequence.length() < g_kmer_length - 1) ? 0 : g_kmer_length - 1;
        std::string sequence = node_set[i].sequence.substr(node_set[i].sequence.length() - len);
        while (sequence.length() >= 2 * g_kmer_length - 2) {
            if (node_set[i].children.size() == 0)
                break;
            int j = node_set[i].children[0];
            sequence += node_set[j].sequence;
            i = j;
        }
        node_num = i;
        sequence = sequence.substr(0, 2 * g_kmer_length - 2);
        if (sequence.find(kmer) != std::string::npos) {
            int pos = sequence.find(kmer);
            sequence = sequence.substr(pos);
            if (sequence.length() < read.length()) {
                while (sequence.length() < read.length()) {
                    if (node_set[node_num].children.size() == 0)
                        return -1;
                    int child = node_set[node_num].children[0];
                    sequence = sequence + node_set[child].sequence;
                    node_num = child;
                }
                sequence = sequence.substr(0, read.length());
                if (forward_similar(sequence, read,1))
                    return node_num;
            }
        }
               
    }
    return -1;
}

std::vector<int> findKeysWithValuesContaining(const boost::unordered_map<int, std::vector<size_t>>& node_read_hash, int value) {
    std::vector<int> result;

    for (const auto& pair : node_read_hash) {
        const std::vector<size_t>& vec = pair.second;

        for (const auto& element : vec) {
            if (element == value) {
                result.push_back(pair.first);
                break;
            }
        }
    }

    return result;
}

bool pairIsContained(const std::pair<int, int>& pair1, const std::vector<std::pair<int, int>>& transcript_path) {
    for (const auto& pair : transcript_path)
        if (pair1 != pair && pair1.first >= pair.first && pair1.second <= pair.second)
            return true;
    return false;
}


void get_revcomp_paired_reads_list(std::vector<std::string>& paired_read_list) {
    
    int size = paired_read_list.size() / 2;
    for (int i = 0; i < size; i++) {
        std::string read = paired_read_list[i];
        paired_read_list[i] = revcomp(paired_read_list[i + size]);
        paired_read_list[i + size] = revcomp(read);
    }
   
}


void find_path(std::vector<Node>& node_set, std::vector<std::string> paired_read_list, Kmer_Hash& kmer_hash,bool& flag, std::vector <std::string>& transcript_list) {
    std::cout << "Beginning find_path ..." << std::endl;
    time_t beg = time(NULL);
  
    std::string fileName1 = out_dir + "mate_node_list.txt";
    std::ofstream ofile1(fileName1);
 

    boost::unordered_map<std::string,std::vector<size_t>>  kmer_read_hash;
    boost::unordered_map<int, std::vector<size_t>>  node_read_hash;
    boost::unordered_map<int, int> read_to_node_hash;
    
    int read_len = paired_read_list[0].length();
    int L = read_len - g_kmer_length;
    std::cout << "paired_reads_map..." << std::endl;
    // first kmer -> read num
    for (int i = 0; i < paired_read_list.size(); i++) {
        std::string read = paired_read_list[i];
        std::string f_kmer = read.substr(0, g_kmer_length);
        kmer_read_hash[f_kmer].push_back(i);
    }
   // std::cout << kmer_read_hash.size() << std::endl;
    std::vector<size_t> vector;

   
    for(int j = 1; j < node_set.size(); ++j) {  
       // std::cout << j << std::endl;
        std::string sequence = node_set[j].sequence;
        if (sequence.length() < read_len) {
            int l = j;
            while (sequence.length() < read_len) {
                if (node_set[l].children.size() == 0)
                    break;
                int child = node_set[l].children[0];
                sequence = sequence + node_set[child].sequence;
                l = child;
            }
            if (sequence.length() < read_len)
                break;
            sequence = sequence.substr(0, read_len);
        }

        boost::unordered_map<size_t, int>  read_set;
        for (int k = 0; k <= sequence.length() - g_kmer_length; ++k) {;
            std::string kmer = sequence.substr(k,g_kmer_length);
            vector = kmer_read_hash[kmer];
            for (int ii = 0; ii < vector.size(); ii++) {
                std::string read = paired_read_list[vector[ii]];
                if (is_similar(sequence, read)) {
                    node_read_hash[j].push_back(vector[ii]);
                    read_to_node_hash[vector[ii]] = j;
                }
            }
            vector.clear();
        }  
        /*
        for (auto it = read_set.begin(); it != read_set.end(); ) {
            if (it->second < read_len-g_kmer_length-2) {
                it = read_set.erase(it);
            }
            else {
                ++it;
            }
        }
        
        
        for (auto it = read_set.begin(); it != read_set.end(); it++) {
            std::string read = paired_read_list[it->first];
            if (is_similar( sequence, read)) {
                node_read_hash[j].push_back(it->first);
                read_to_node_hash[it->first] = j;
            }
        } 
        */
    }

    /*
    std::ofstream file("/home/bioinfo/limh/paper/data/SRR1942956/node_read_list.txt");

    for (const auto& entry : node_read_hash) {
        int key = entry.first;
        const std::vector<size_t>& value = entry.second;
        file << "Key: " << key << std::endl;
        file << "Values:";
        for (const auto& elem : value) {
            file << " " << elem;
        }
        file << std::endl << std::endl;
    }

    file.close();
    */

    std::cout << "Begins calculating the abundance of each node... " << std::endl;
    std::vector<std::pair<int, int>> transcript_path;
    boost::unordered_map<int, int > mate_node_hash;
    boost::unordered_map<int, int > mate_node_cov_hash;
    boost::unordered_map<int, int > node_cov_hash;
    for (int j = 1; j < node_set.size(); ++j) {    //计算各节点的丰度，如果是1或者偶数节点，计算连接处的丰度，非1的结束节点计算节点平均丰度
        if (node_set.size() == 2)
            continue;
        if (j == 1 || j % 2 == 0) {
            int chil = node_set[j].children[0];
            int l = (node_set[chil].sequence.length() < g_kmer_length - 1) ? node_set[chil].sequence.length() : g_kmer_length - 1;
            std::string con_seq = node_set[j].sequence.substr(node_set[j].sequence.length() - g_kmer_length + 1) + node_set[chil].sequence.substr(0, l);
            int cov_con_seq = get_str_coverage(con_seq, kmer_hash,flag);
            node_cov_hash[j] = cov_con_seq;
        }
        else {
            int l = j;
            std::string str = node_set[l].sequence;
            int len_threshold = str.length();
            len_threshold = str.length() + g_kmer_length - 1;
            while (str.length() < len_threshold) {
                if (node_set[l].children.size() == 0)
                    break;
                int child = node_set[l].children[0];
                str += node_set[child].sequence;
                l = child;
            }
            str = str.substr(0, len_threshold);
            int cov = get_str_coverage(str, kmer_hash,flag);
            node_cov_hash[j] = cov;
        }
    }
    /*
    for (auto it = node_cov_hash.begin(); it != node_cov_hash.end(); ++it) {
        std::cout << "Key: " << it->first << ", Value: " << it->second << std::endl;
    }
    */
    //重建转录本路径

    std::cout << "Start building transcript paths..." << std::endl;
    for (int j = 1; j < node_set.size(); ++j) {   //遍历所有的没有父节点的节点
        if (!node_set[j].parents.empty() || node_cov_hash[j] < 0.01 * node_cov_hash[j + 1]) {
            continue;
        }
           
        std::vector<size_t> read_set = node_read_hash[j];
        boost::unordered_map<int, int> node_hash;;
        for (int i = 0; i < read_set.size(); i++) {
         
            int mate_read_num;
            if (read_set[i] >= paired_read_list.size() / 2)
                mate_read_num = read_set[i] - paired_read_list.size() / 2;
            else
                mate_read_num = read_set[i] + paired_read_list.size() / 2;
            if (read_to_node_hash.count(mate_read_num) == 0 || read_to_node_hash[mate_read_num] <= j || 
                read_to_node_hash[mate_read_num] % 2 == 0 || read_to_node_hash[mate_read_num] - j >= 10 )
                continue;
            
          //  ofile << j <<  "  " << read_to_node_hash[mate_read_num]  << " " << read_set[i]  << " " << mate_read_num << std::endl;
          //  ofile << paired_read_list[read_set[i]] << "    " << paired_read_list[mate_read_num] << std::endl;

            node_hash[read_to_node_hash[mate_read_num]]++;  //配对节点的个数
        }

        ofile1 << j  << " " << node_cov_hash[j] << std::endl;  // j节点连接处的丰度
        mate_node_hash[j] = 0;
 
        for (auto it = node_hash.begin(); it != node_hash.end(); ++it) {
            if (it->first < j || it->second <= 1)
                continue;
            if (it->first > mate_node_hash[j])
                mate_node_hash[j] = it->first;
     
            ofile1 << it->first << " " << it->second  << " str.cov: " << node_cov_hash[it->first] << std::endl;  // 
        }

  //      /*
        int l = j , sum = 0;
        float cov = 0, cov_p = 0, cov_v = 0 , mw = 0;
        while (!node_set[l].children.empty() && mate_node_hash[j] != 0) { //&& node_set[l].children[0] <= mate_node_hash[j]

            if (cov == 0) {
                cov = node_cov_hash[node_set[l].children[0]];
                cov_p = cov;
                sum++;
            }

            else {
                cov_v = node_cov_hash[node_set[l].children[0]];
                sum++;
                if (l <= mate_node_hash[j])
                    mw = g_map_weight;
                else
                    mw = 0;
                float cov_w = (-log2(fabs(1 - cov_v / cov_p)));
                float score = mw + cov_w;
                ofile1 << node_set[l].children[0] << " " << score << " " << mw << " " << cov_w << " " << cov_v / cov_p << std::endl;
                if (score <= 1)
                    break;
               
                cov += node_cov_hash[node_set[l].children[0]];
                cov_p = cov / sum;
            }
            l = node_set[l].children[0];

        }
        /*
        if(!node_set[l].children.empty())
            mate_node_hash[j] = node_set[l].children[0];
        else
        */
            mate_node_hash[j] = l;


     //   */
        node_hash.clear();
        if (mate_node_hash[j] == j) {
            mate_node_hash.erase(j);
            continue;
        }
        transcript_path.emplace_back(j, mate_node_hash[j]);
        ofile1 << " **************************************************************" << std::endl;

    }

   /*
    for (const auto& pair : transcript_path) {
        std::cout << pair.first << " " << pair.second << std::endl;
    }
  //  
    /*
    // 去除了被包含的路径
    std::vector<std::pair<int, int>> filtered_transcript_path;
    for (const auto& pair1 : transcript_path) {
        if (pairIsContained(pair1, transcript_path)) {
            filtered_transcript_path.push_back(pair1);
        }
    }

    for (const auto& pair : filtered_transcript_path) {
        mate_node_hash.erase(pair.first);
    }
   // /*
    for (auto it = mate_node_hash.begin(); it != mate_node_hash.end(); it++) {
        std::cout << it->first << " " << it->second << std::endl;
    }
   // */



    //根据转录本路径重建转录本
    std::cout << "Begin rebuilding the transcript..." << std::endl;
    //int k = 1;
    
    for (int j = 1; j < node_set.size(); ++j) {   //遍历所有的没有父节点的节点
        if (!node_set[j].parents.empty() || mate_node_hash[j] == 0)
            continue;
        std::string transcript;
        if(j == 1)
            transcript  = node_set[j].sequence;

        int l = j;
        int cov = 0, i = 0;
        while (!node_set[l].children.empty() && node_set[l].children[0] <=mate_node_hash[j]) { //
   
            transcript += node_set[node_set[l].children[0]].sequence;
            l = node_set[l].children[0];
        }

        if (j == 1) {
            int start = 0, end1 = 0, end2 = 0;
            for (int i = 0; i < transcript.length() - 3; i++) {
                std::string s = transcript.substr(i, 3);
                if (s == "ATG" && i > 200) {
                    start = i;
                    break;
                }
            }
            for (int i = start; i < transcript.length() - 3; i+=3) {
                std::string s = transcript.substr(i, 3);

                if (s == "TAA" ) {
                    end1 = i+3;
                    break;
                }
            }
            for (int i = start; i < transcript.length() - 3; i++) {
                std::string s = transcript.substr(i, 3);
                if (s == "TAA" && (i - start) % 3 == 2 && (i - start) * 1.0 / (end1 - start) > 1.5) {
                    end2 = i + 3;
                    break;
                }
            }
            if (end1 != 0) {
                transcript_list.push_back(transcript.substr(start,end1-start));
             //   ofile2 << ">" << "transcript_" << k++ << "__" << j << std::endl;
             //   ofile2 << transcript.substr(start, end1 - start) << std::endl;
            }
            if (end2 != 0) {
                transcript_list.push_back(transcript.substr(start, end2 - start));
               // ofile2 << ">" << "transcript_" << k++ << "__" << j << std::endl;
               // ofile2 << transcript.substr(start, end2 - start) << std::endl;
            }
            if (end2 == 0) {
                transcript_list.push_back(transcript);

             //   ofile2 << ">" << "transcript_" << k++ << "__" << j << std::endl;
             //   ofile2 << transcript << std::endl;
            }
        }
       
       // for (int i = transcript.length()-3 ; i >=0; i--) {
        else {
            for (int i = 0; i <= transcript.length() - 3; i += 3) {
                std::string stop_codon = transcript.substr(i, 3);


                if (stop_codon == "TAA" || stop_codon == "TAG" || stop_codon == "TGA") {
                    transcript = transcript.substr(0, i + 3);
                    if (transcript.length() <= 100)
                        continue;
                    transcript_list.push_back(transcript);
                  //  ofile2 << ">" << "transcript_" << k++ << "__" << j << std::endl;
                  //  ofile2 << transcript << std::endl;
                    break;
                }

            }
        }
    }

    time_t end = time(NULL);
    std::cout << "Finished find_path "  <<  transcript_list.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;
}

void get_branch_reads(std::vector<std::string>& contigs_list, Read_Hash& read_hash, std::vector<std::string>& branch_read_list,bool& flag) {
    boost::unordered_map<std::string, size_t>::iterator it;
    it = read_hash.begin();
    int len = it->first.length() + 10;
    for (; it != read_hash.end(); ++it) {
        std::string read = it->first;
        std::string fkmer = read.substr(0, g_kmer_length);
        std::string lkmer = read.substr(read.length() - g_kmer_length);
        bool map = false;
        for (int i = 0; i < contigs_list.size(); i++) {
            std::string contig = contigs_list[i];
            int first = contig.find(fkmer);
            int last = contig.find(lkmer);
            if (first != std::string::npos && last != std::string::npos && last - first < len) {
                map = true;
                break;
            }
        }
        if (!map) {
            if(flag)
                branch_read_list.push_back(read);
            else
                branch_read_list.push_back(revcomp(read));
        }
    }
}


void construct_graph(std::vector<std::string>& contigs_list, Kmer_Hash& kmer_hash, Read_Hash& read_hash, std::vector<std::string>& paired_read_list) {
    
    std::cout << "Beginning get_transcript ..." << std::endl;
    time_t beg = time(NULL);
    bool flag = true;
    std::sort(contigs_list.begin(), contigs_list.end(), [](std::string& str1, std::string& str2) { return str1.length() > str2.length(); });
    
   
    std::string trunk = contigs_list[0];
    std::vector<std::string> read_list;

    if (trunk.length() >= 0.9 * g_ref_genome_len) {
        std::string fileName = out_dir + "node_set.fasta";
        std::ofstream outfile(fileName);
        std::string fileName2 = out_dir + "transcript.fasta";
        std::ofstream ofile2(fileName2);
        std::vector <std::string> transcript_list;

        std::vector<Node> node_set;
        for (int i = 0; i < 10; i++) {
            if (trunk.substr(i, 5) == "TTTTT") {
                flag = false;
                trunk = revcomp(trunk);
                break;
            }
        }
        transcript_list.push_back(trunk);

        get_read_list(read_hash, read_list, flag);
        Node node(trunk);
        int p = add_node(node, node_set);
        check_reverse_branch(p, read_list, kmer_hash, node_set, flag);
        //add_reverse_branch(p, read_list, kmer_hash, node_set);
        if (!flag) 
            get_revcomp_paired_reads_list(paired_read_list);
      
        find_path(node_set, paired_read_list, kmer_hash, flag, transcript_list);

        for (int i = 0; i < transcript_list.size(); i++) {
            ofile2 << ">" << "transcript_" << i << std::endl;
            ofile2 << transcript_list[i] << std::endl;
        }
        ofile2.close();

        for (int i = 1; i < node_set.size(); i++) {
            outfile << "> " << i << std::endl;
            outfile << node_set[i].sequence << std::endl;
        }
        outfile.close();

    }
    else {
        int k = 0;
        std::vector<std::string> transcript_list;
        std::string fileName2 = out_dir + "transcript.fasta";
        std::ofstream ofile2(fileName2);
        
        for (int i = 0; i < contigs_list.size(); i++) {
            std::string trunk = contigs_list[i];
            if (trunk.length() > 0.65 * g_ref_genome_len) {
                int start = 0, end1 = 0;
                std::string fkmer = trunk.substr(100, g_kmer_length);
                kmer_int_type_t fkmer_val = kmer_to_intval(fkmer);
                int cov1 = kmer_hash[fkmer_val];
                std::string lkmer = trunk.substr(trunk.length() - 100, g_kmer_length);
                kmer_int_type_t lkmer_val = kmer_to_intval(lkmer);
                int cov2 = kmer_hash[lkmer_val];
                if (cov1 < cov2  ) 
                    trunk = revcomp(trunk);

                for (int i = 0.5 * trunk.length();  i < 0.7 * trunk.length() - 3; i ++) {
                    std::string kmer = trunk.substr(i, g_kmer_length);

                    if (kmer.substr(g_kmer_length-3) == "TAA") {
                        end1 = i + g_kmer_length;
                        kmer_int_type_t kmer_value = kmer_to_intval(kmer);
                        int cov = kmer_hash[kmer_value];
                        std::string kmer2 = trunk.substr(i + g_kmer_length, g_kmer_length);
                        kmer_int_type_t kmer_value2 = kmer_to_intval(kmer2);
                        int cov2 = kmer_hash[kmer_value2];
                     //   std::cout << end1 << " " << cov << " " << cov2 << " " << cov2 * 1.0 / cov << std::endl;
                        if (fabs(1 - cov2 * 1.0 / cov) < 0.3 || fabs(cov-cov2) < 10)
                            continue;
                        else
                            break;
                    }
                }

                transcript_list.push_back(trunk.substr(0, end1 ));
                ofile2 << ">" << "trunk_" << k++  << std::endl;
                ofile2 << trunk.substr(0, end1) << std::endl;

                transcript_list.push_back(trunk);
                ofile2 << ">" << "trunk_" << k++ << std::endl;
                ofile2 << trunk << std::endl;
               
            }
            else {
                transcript_list.push_back(trunk);
                ofile2 << ">" << "trunk_" << k++ << std::endl;
                ofile2 << trunk << std::endl;
            }     
        }
        std::cout << "Finished get " << transcript_list.size() << " transcripts " << std::endl;
       
    }

    time_t end = time(NULL);
    std::cout << "Finished ,(elapsed time: " << (end - beg) << " s)" << std::endl;
}


std::vector<std::string> remove_isContained(  std::vector<std::string>& strings) {
    
    std::vector<std::string> filteredStrings;
    std::sort(strings.begin(), strings.end(), [](const std::string& str1, const std::string& str2) {return str1.length() > str2.length(); });
    for (const std::string& str : strings) {
        bool isContained = false;
        for (const std::string& otherStr : filteredStrings) {
            if (str != otherStr && str <= otherStr && is_similar(otherStr,str)|| is_similar(otherStr,revcomp(str))) {
                isContained = true;
                break;
            }
        }

        if (!isContained) {
            filteredStrings.push_back(str);
        }
       
    }
    return filteredStrings;

}
/*
void merge_extension(std::vector<std::string>& contigs_list) {

    size_t i = find_longest_read(contigs_list);
    std::string trunk = contigs_list[i];
    int kmer_len = g_kmer_length - 1;
    std::string first_kmer = trunk.substr(0, kmer_len);
    std::string last_kmer = trunk.substr(trunk.length() - kmer_len);
    for (int j = 0; j < contigs_list.size(); j++) {
        if (j == i)
            continue;
        std::string extend_contig = contigs_list[j];
        int start = extend_contig.find(first_kmer);
        int end = extend_contig.find(last_kmer);
        if (start != std::string::npos && end != std::string::npos) {
            trunk = extend_contig.substr(0, start) + trunk + extend_contig.substr(end + kmer_len);
            contigs_list[i] = trunk;
            contigs_list.erase(contigs_list.begin() + j);
            merge_extension(contigs_list);
        }
    }
}

*/

void merge_extension(std::vector<std::string>& contigs_list) {
    size_t i = find_longest_read(contigs_list);
    std::string trunk = contigs_list[i];
    int kmer_len = g_kmer_length - 1;
    std::string first_kmer = trunk.substr(0, kmer_len);
    std::string last_kmer = trunk.substr(trunk.length() - kmer_len);
    bool merged = false; 
    for (int j = 0; j < contigs_list.size(); j++) {
        if (j == i)
            continue;

        std::string extend_contig = contigs_list[j];
        int start = extend_contig.find(first_kmer);
        int end = extend_contig.find(last_kmer);
        if (start != std::string::npos && end != std::string::npos) {
            std::string new_trunk = extend_contig.substr(0, start) + trunk + extend_contig.substr(end + kmer_len);
            contigs_list[i] = new_trunk;
            contigs_list.erase(contigs_list.begin() + j);
            merged = true; 
            merge_extension(contigs_list);
            break;
        }
    }
    if (!merged) {
        return;
    }
}

int haveCommonSequence(const std::string& str1, const std::string& str2, int commonLength) {
    if (str1.length() < commonLength || str2.length() < commonLength) {
        return false;
    }
    std::string kmer = str1.substr( str1.length() - commonLength);
    int pos = str2.find(kmer);
    if (pos == std::string::npos || pos > g_kmer_length)
        return -1;
    if (str2.substr(0,pos+commonLength) == str1.substr(str1.length()-pos-commonLength))
        return pos+commonLength;
    else
        return -1;
}

void head_tail_connection(std::vector<std::string>& strings) {
    int commonLength = 7;

    for (size_t i = 0; i < strings.size(); ++i) {
        bool merged = false;
        for (size_t j = 0; j < strings.size(); ++j) {
            if (i == j) continue;
            int L = haveCommonSequence(strings[i], strings[j], commonLength);
            if (L != -1) {
                std::string mergedString = strings[i] + strings[j].substr(L);
                
                strings.erase(strings.begin() + i);
                strings.erase(strings.begin() + j);
                strings.push_back(mergedString);
                head_tail_connection(strings);
            }
        }
        
    }
}

void remove_revcomp(std::vector<std::string>& contigs_list) {

    std::cerr << "Beginning remove_revcomp  ..." << std::endl;
    time_t beg = time(NULL);

   
    contigs_list = remove_isContained(contigs_list);
   // std::cout << contigs_list.size() << std::endl;
   // merge_extension(contigs_list);
   
   // std::cout << contigs_list.size() << std::endl;
   
    head_tail_connection(contigs_list);

    time_t end = time(NULL);
    int j = 0;
    std::string fileName = out_dir + "contigs2.fasta";
    std::ofstream outfile(fileName);
   // std::ofstream outputFile("/home/bioinfo/limh/paper/data/SRR1942956/contigs2.fasta");
   // std::ofstream outputFile("C:/Users/lmh/Desktop/contigs2.fasta");
    if (outfile.is_open()) {
        for (const auto& pair : contigs_list) {
            outfile << "> " << j++ << std::endl;
            outfile << pair << "\n";
        }
    }
    outfile.close();
    
    std::cout << "contigs_list count after remove_revcomp : " << contigs_list.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;
}





void remove_revread(std::vector<std::string>& contigs_list, Read_Hash& read_hash, std::vector<std::string>& read_list) {

    std::cerr << "Beginning remove_revread  ..." << std::endl;
    time_t beg = time(NULL);
    std::string trunk = contigs_list[0];
    for (auto it = read_hash.begin(); it != read_hash.end(); ++it) {
        std::string read = it->first;
        if (!is_similar(trunk, revcomp(read))) {
            read_list.push_back(read);
        }
    }
   
    time_t end = time(NULL);
    std::cout << "Finished remove_revread : " << read_list.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;
}

void get_contigs_cov(std::vector<std::string>& contigs_list, Kmer_Hash& kmer_hash) {
   
    for (int i = 0; i < contigs_list.size(); i++) {
       // size_t n = find_longest_read(contigs_list);
        std::string trunk = contigs_list[i];

        ///*
        std::vector<size_t> trunk_cov;
        for (int i = 0; i <= trunk.length() - g_kmer_length; i++) {
            std::string kmer = trunk.substr(i, g_kmer_length);
            kmer_int_type_t val = kmer_to_intval(kmer);
            size_t cov = kmer_hash[val];
            trunk_cov.push_back(cov);
        }
        std::string fileName = out_dir + "trunk_cov_" + std::to_string(i) + ".txt";
        //std::string fileName = "/home/bioinfo/limh/paper/data/SRR1942956/transcript_cov_"  + std::to_string(t) + ".txt";
        std::ofstream of(fileName);

        for (int i = 0; i < trunk_cov.size(); i++) {
            of << i << " : " << trunk[i] << " : " << trunk_cov[i] << std::endl;
        }

        of.close();
    }
}


void Construct_splicing_graphs(std::vector<std::string>& Read) {

    Read_Hash read_hash;
    std::vector<std::string> paired_read_list;
    get_read_hash(read_hash, Read, paired_read_list);

    Kmer_Hash kmer_hash;
    get_kmer_hash(kmer_hash, read_hash);
    remove_erroneous(kmer_hash, min_ratio_non_error);



    std::vector<std::string> contigs_list;
    extend_contigs(kmer_hash, contigs_list);

    remove_revcomp(contigs_list);

 
    
    
     construct_graph(contigs_list, kmer_hash, read_hash,paired_read_list);
    

}

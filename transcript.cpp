#include "transcript.h"
#include <vector>
#include "utility.h"
#include "common.h"
#include <fstream>
#include <boost/unordered_map.hpp>
#include <map>
#include <thread>
#include <mutex>

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
typedef typename boost::unordered_map<std::string, std::string> Paired_Read_Hash;
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
Read_Hash read_hash;
Kmer_Hash kmer_hash;
std::vector<std::string> read_list;

size_t node_size = 0;  // the number of nodes in this graph
int L = 101;  // g_kmer_length = 31;
int min_error_count = 3;
float min_error_rate = 0.05;
int max_rep_trunk_threshold = g_kmer_length + L;
float min_ratio_non_error = 0.05f;
float min_seed_entropy = 1.5f;
int trunk_head_threshoid = 300;
std::mutex hash_mutex;


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


void reads_hash_process_chunk(const std::vector<std::string>& reads, size_t start, size_t end) {
    boost::unordered_map<std::string, size_t> local_hash;

    for (size_t i = start; i < end; ++i) {
        const std::string& sequence = reads[i];
        if (sequence.empty()) continue;
        local_hash[sequence]++;
    }

    std::lock_guard<std::mutex> lock(hash_mutex);
    for (const auto& entry : local_hash) {
        read_hash[entry.first] += entry.second;
    }
}

void get_read_hash(Read_Hash& read_hash, std::vector<std::string>& Read) {

    std::cerr << "Beginning read hash ..." << std::endl;
    time_t beg = time(NULL);
    if (Read.empty()) return;
    size_t data_size = Read.size();
    size_t num_threads = g_threads;
    size_t chunk_size = data_size / num_threads;
    std::vector<std::thread> threads;

    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? data_size : start + chunk_size;
        threads.emplace_back(reads_hash_process_chunk, std::cref(Read), start, end);
    }

    for (auto& th : threads) {
        th.join();
    }
    time_t end = time(NULL);

    std::cerr << "read hash finished, get " << read_hash.size() << " reads! (elapsed time: " << (end - beg) << " s)" << std::endl;
}

void process_kmer(const boost::unordered_map<std::string, size_t>& read_hash_chunk) {
    boost::unordered_map<kmer_int_type_t, size_t> local_kmer_hash;
    for (const auto& entry : read_hash_chunk) {
        const std::string& sequence = entry.first;
        for (size_t j = 0; j <= sequence.length() - g_kmer_length; ++j) {
            const std::string& kmer = sequence.substr(j, g_kmer_length);
            if (kmer.find("N") != std::string::npos)
                continue;
            kmer_int_type_t kmer_val = kmer_to_intval(kmer, g_kmer_length);
            local_kmer_hash[kmer_val] += entry.second;
        }
    }

    std::lock_guard<std::mutex> lock(hash_mutex);
    for (const auto& entry : local_kmer_hash) {
        kmer_hash[entry.first] += entry.second;
    }
}

void get_kmer_hash(Kmer_Hash& kmer_hash, Read_Hash& read_hash) {

    std::cerr << "Beginning kmer hash ..." << std::endl;
    time_t beg = time(NULL);
    size_t num_threads = g_threads;
    std::vector<std::thread> threads;
    std::vector<boost::unordered_map<std::string, size_t>> read_hash_chunks(num_threads);

    size_t index = 0;
    for (const auto& entry : read_hash) {
        read_hash_chunks[index % num_threads][entry.first] = entry.second;
        ++index;
    }

    for (size_t i = 0; i < num_threads; ++i) {
        threads.emplace_back(process_kmer, std::ref(read_hash_chunks[i]));
    }

    for (auto& th : threads) {
        th.join();
    }
    read_hash.clear();

    time_t end = time(NULL);
    if (kmer_hash.empty()) {
        std::cout << "kmer_hash is empty." << std::endl;
    }
    std::cerr << "Kmer hash finished, got " << kmer_hash.size()
        << " kmers! (elapsed time: " << (end - beg) << " s)" << std::endl;
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

                kmer_hash.erase(f_candidates[i].first);

            }

            if (f_candidates[i].second > 1) {
                int candidate_count = f_candidates[i].second;
                if (f_dominant_count == 0) {
                    f_dominant_count = candidate_count;   //同组最大的kmer的丰度
                }
                else if ((float)candidate_count / f_dominant_count < min_ratio_non_error) {
                    kmer_hash.erase(f_candidates[i].first);
                }
            }
        }
        int r_dominant_count = 0;
        for (unsigned int i = 0; i < r_candidates.size(); ++i) {
            if (r_candidates[i].second == 1) {

                kmer_hash.erase(r_candidates[i].first);

            }

            if (r_candidates[i].second > 1) {
                int candidate_count = r_candidates[i].second;
                if (r_dominant_count == 0) {
                    r_dominant_count = candidate_count;   //同组最大的kmer的丰度
                }
                else if ((float)candidate_count / r_dominant_count < min_ratio_non_error) {
                    kmer_hash.erase(r_candidates[i].first);
                }
            }
        }

    }

    time_t end = time(NULL);
    std::cout << "kmer count after errors deletion : " << kmer_hash.size()
        << " (elapsed time: " << (end - beg) << " s)" << std::endl;

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
        seed_hash.push_back(std::make_pair(it->first, it->second));
    }

    time_t end = time(NULL);
    std::cerr << "get seed kmer finished, get " << seed_hash.size()
        << " seed kmers! (elapsed time: " << (end - beg) << " s)" << std::endl;
}



std::string forward_extend_kunitig(Kmer_Hash& kmer_hash, kmer_int_type_t kmer_val, Kmer_Hash& used_kmers_) {

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

std::string forward_extend_kunitig(Kmer_Hash& kmer_hash, std::string kunitig) {

    std::string rkmer = kunitig.substr(kunitig.length() - g_kmer_length);
    kmer_int_type_t intval = kmer_to_intval(rkmer);

    std::string str = kunitig;
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


std::string reverse_extend_kunitig(Kmer_Hash& kmer_hash, kmer_int_type_t kmer_val, Kmer_Hash& used_kmers_) {

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


std::string reverse_extend_kunitig(Kmer_Hash& kmer_hash, std::string kunitig) {
    std::string lkmer = kunitig.substr(0, g_kmer_length);
    kmer_int_type_t intval = kmer_to_intval(lkmer);
    std::string str = kunitig;

    std::vector<kmer_occurence_pair_t> candidates;
    while (true) {
        get_reverse_candidates(kmer_hash, intval, candidates);
        if (candidates.empty())
            break;

        std::sort(candidates.begin(), candidates.end(), sortBySecond); //降序
        kmer_int_type_t candidate = candidates[0].first;
        std::string kmer = intval_to_kmer(candidate, g_kmer_length);
        if (str.find(kmer) != std::string::npos) {
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




void extend_kunitigs(Kmer_Hash kmer_hash, std::vector<std::string>& kunitigs_list) {

    time_t beg = time(NULL);
    std::cerr << "extend_kunitigs ..." << std::endl;

    Kmer_Hash used_kmers_;
    std::vector<kmer_occurence_pair_t>  seed_kmers;
    get_seed_kmer(kmer_hash, seed_kmers);
    for (int i = 0; i < seed_kmers.size(); i++) {
        kmer_int_type_t kmer_val = seed_kmers[i].first;

        if (used_kmers_.find(kmer_val) == used_kmers_.end()) {
            used_kmers_.emplace(kmer_val, kmer_hash[kmer_val]);
            std::string left = reverse_extend_kunitig(kmer_hash, kmer_val, used_kmers_);
            std::string right = forward_extend_kunitig(kmer_hash, kmer_val, used_kmers_);
            std::string kunitig = left + right.substr(g_kmer_length);

            if (kunitig.length() > 3 * g_kmer_length) {
                ///*
                left = reverse_extend_kunitig(kmer_hash, kunitig);
                right = forward_extend_kunitig(kmer_hash, kunitig);
                if (right.length() == kunitig.length())
                    kunitig = left;
                else
                    kunitig = left + right.substr(kunitig.length());
                // */
                kunitigs_list.push_back(kunitig);

            }

        }

    }


    time_t end = time(NULL);
    /*
    int i = 0;
    std::string fileName = out_dir + "kunitigs.fasta";
    std::ofstream outputFile(fileName);

    if (outputFile.is_open()) {
        for (const auto& pair : kunitigs_list) {
            outputFile << "> " << i++ << std::endl;
            outputFile << pair << "\n";
        }
    }
    outputFile.close();
      */
    std::cerr << "kunitigs  finished, get " << kunitigs_list.size()
        << " kunitigs! (elapsed time: " << (end - beg) << " s)" << std::endl;
}

bool is_similar(std::string str1, std::string str2) {
    if (str1.length() == 0 || str2.length() == 0) {
        return false;
    }

    if (str1.length() < str2.length()) {
        std::string t = str2;
        str2 = str1;
        str1 = t;
    }

    int mismatch = 0;
    int kmer_length = g_kmer_length;

    std::string first_kmer = str2.substr(0, kmer_length);
    size_t start = str1.find(first_kmer);
    std::string last_kmer = str2.substr(str2.length() - kmer_length);
    size_t end = str1.find(last_kmer);


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

float get_str_coverage(std::string str, Kmer_Hash& kmer_hash, bool flag) {
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

    for (; first <= last; ++first) {
        sum += cov_v[first];
    }
    float cov = sum * 1.0 / (cov_size - 2 * quantile);

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

bool forward_similar(const std::string& str1, const std::string& str2, int n) {
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
    if (static_cast<float>(match) / (length + 1) > 0.7 || match > 40) {
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

    if (similarity.size() >= 4 && similarity[3] > read.length() * 0.95)
        return true;
    else
        return false;

}


int different_point(std::string str1, std::string str2) {
    int i = 0;
    while (i < str1.length() - 1) {
        if (str1[i] != str2[i] && str1[i + 1] != str2[i + 1])
            break;
        i++;
    }
    return i;
}

void get_paired_read_list(std::vector<std::string>& paired_reads, std::vector<std::string>& Read) {

    std::cout << "Beginning to get paired reads..." << std::endl;
    time_t beg = time(NULL);

    boost::unordered_map<std::string, std::string> pair_reads_hash;
    int size = Read.size();

    std::vector<std::thread> threads;
    for (int i = 0; i < size; i++) {
        if (i < size / 2)
            pair_reads_hash[Read[i]] = Read[i + size / 2];
        else
            pair_reads_hash[Read[i]] = Read[i - size / 2];
    }

    Read.clear();

    size = pair_reads_hash.size();
    std::vector<std::string> mate_reads;
    for (auto it = pair_reads_hash.begin(); it != pair_reads_hash.end(); it++) {
        paired_reads.push_back(it->first);
        mate_reads.push_back(it->second);
    }
    pair_reads_hash.clear();

    paired_reads.insert(paired_reads.end(), mate_reads.begin(), mate_reads.end());
    mate_reads.clear();
    time_t end = time(NULL);
    std::cout << "Finished, get paired_read " << paired_reads.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;

}

bool is_ATG_pos(int pos, std::string trunk, int& ATG_pos) {

    for (int i = 0; i < 15; i++) {
        if (trunk.substr(pos + i, 3) == "ATG") {
            ATG_pos = pos + i;
            return true;
        }

        if (trunk.substr(pos - i, 3) == "ATG") {
            ATG_pos = pos - i;
            return true;
        }
    }
    return false;

}


void check_branch(int p, std::vector<std::string>& read_list, Kmer_Hash kmer_hash, std::vector<Node>& node_set, bool flag, std::vector <int>& jun_poss) {

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



            if (f_flag && kmer_hash[jun_val] > 3) {


                std::vector<size_t> jun_reads;
                if (!flag)
                    jun_kmer = revcomp(jun_kmer);
                find_read_contain_kmer(read_list, jun_kmer, jun_reads);
                std::string ref_str = trunk.substr(0, start) + candidate_read + trunk.substr(end + kmer_len, 40);
                if (!support_edge_decision(ref_str, jun_reads, read_list))
                    continue;
                if (g_non_canonical == 1) {

                    jun_poss.push_back(i);
                }

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

bool is_jun_pos(int jum_pos, std::vector <int>& jun_poss) {
    for (int i = 0; i < jun_poss.size(); i++) {
        if (jun_poss[i] - 20 < jum_pos && jum_pos < jun_poss[i] + 20) {
            return true;
        }
    }
    return false;
}

void check_non_canonical_branch(int p, std::vector<std::string>& read_list, Kmer_Hash kmer_hash, std::vector<Node>& node_set, bool flag, std::vector <int>& jun_poss) {

    std::cout << "Beginning add branch ..." << std::endl;
    time_t beg = time(NULL);

    std::string trunk = node_set[p].sequence;

    std::string fileName = out_dir + "branch2.fasta";
    std::ofstream ofile(fileName);

    boost::unordered_map<int, int> jump_pos_hash;
    std::map<int, std::string> jump_reads_hash;

    int read_len = read_list[0].length();
    for (int i = 0; i < read_list.size(); i++) {
        std::string read = read_list[i];
        std::string first_kmer = read.substr(0, g_kmer_length);
        std::string last_kmer = read.substr(read_len - g_kmer_length);
        int start = trunk.find(first_kmer);
        int end = trunk.find(last_kmer);
        int pos = different_point(trunk.substr(start, read_len), read);
        int jum_pos = end + g_kmer_length + pos - read_len;

        bool f_flag = forward_similar(read.substr(pos), trunk.substr(jum_pos, read_len - pos), 2);
        if (!f_flag) { continue; }

        jump_pos_hash[jum_pos]++;
        if (jump_pos_hash[jum_pos] > 5 && !is_jun_pos(jum_pos, jun_poss)) {
            jump_reads_hash[jum_pos] = read.substr(0, pos);
            ofile << ">" << "pos:  " << jum_pos << std::endl;
            ofile << read << std::endl;
        }
    }

    int pre_pos = 0;
    for (auto it = jump_reads_hash.begin(); it != jump_reads_hash.end(); it++) {

        Node node1, node2;
        node1 = node_set[p].sequence.substr(0, it->first - pre_pos);
        node_set[p].sequence = node_set[p].sequence.substr(it->first - pre_pos);
        node_idx_t q1 = add_node(node1, node_set);
        node2 = it->second;
        node_idx_t q2 = add_node(node2, node_set);
        node_set[q1].parents = node_set[p].parents;
        node_set[p].parents.clear();
        node_set[p].add_parent(q1);
        node_set[p].add_parent(q2);
        pre_pos = it->first;
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



int find_node_contain_kmer(std::vector<Node>& node_set, std::string kmer) {
    int i;
    for (i = 0; i < node_set.size(); ++i) {
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
        std::string sequence = node_set[i].sequence.substr(node_set[i].sequence.length() - len);
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
    int node_num;
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
                if (forward_similar(sequence, read, 1))
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
                if (forward_similar(sequence, read, 1))
                    return node_num;
            }
        }

    }
    return -1;
}



void find_path(std::vector<Node>& node_set, std::vector<std::string> paired_read_list, Kmer_Hash& kmer_hash, bool& flag, std::vector <std::string>& transcript_list, bool create_non_canonical) {
    std::cout << "Beginning find_path ..." << std::endl;
    time_t beg = time(NULL);

    std::string fileName1 = out_dir + "mate_node_list.txt";
    std::ofstream ofile1(fileName1);
    std::string fileName2 = out_dir + "non_tran_cov_change.txt";
    std::ofstream ofile2(fileName2);


    boost::unordered_map<std::string, std::vector<size_t>>  kmer_read_hash;
    boost::unordered_map<int, std::vector<size_t>>  node_read_hash;
    boost::unordered_map<int, int> read_to_node_hash;

    int read_len = paired_read_list[0].length();
    int L = read_len - g_kmer_length;
    std::cout << "paired_reads_map..." << std::endl;
    // first kmer -> read num
    if (flag) {
        for (int i = 0; i < paired_read_list.size(); i++) {
            std::string read = paired_read_list[i];
            std::string f_kmer = read.substr(0, g_kmer_length);
            kmer_read_hash[f_kmer].push_back(i);
        }
    }
    else {
        for (int i = 0; i < paired_read_list.size(); i++) {
            std::string read = paired_read_list[i];
            std::string l_kmer = revcomp(read.substr(read_len - g_kmer_length));
            kmer_read_hash[l_kmer].push_back(i);
        }
    }
    // std::cout << kmer_read_hash.size() << std::endl;

    std::vector<size_t> vector;
    for (int j = 1; j < node_set.size(); ++j) {

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
        for (int k = 0; k <= sequence.length() - g_kmer_length; ++k) {

            std::string kmer = sequence.substr(k, g_kmer_length);
            vector = kmer_read_hash[kmer];
            for (int ii = 0; ii < vector.size(); ii++) {
                std::string read = paired_read_list[vector[ii]];
                if (!flag)
                    read = revcomp(read);
                if (is_similar(sequence, read)) {
                    node_read_hash[j].push_back(vector[ii]);
                    read_to_node_hash[vector[ii]] = j;
                }
            }
            vector.clear();
        }

    }



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
            int start = (0 < int(node_set[j].sequence.length() - g_kmer_length + 1)) ? node_set[j].sequence.length() - g_kmer_length + 1 : 0;
            std::string con_seq = node_set[j].sequence.substr(start) + node_set[chil].sequence.substr(0, l);

            int cov_con_seq = get_str_coverage(con_seq, kmer_hash, flag);
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
            int cov = get_str_coverage(str, kmer_hash, flag);
            node_cov_hash[j] = cov;
        }
    }

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
                read_to_node_hash[mate_read_num] % 2 == 0 || read_to_node_hash[mate_read_num] - j >= 10)
                continue;



            node_hash[read_to_node_hash[mate_read_num]]++;  //配对节点的个数
        }

        ofile1 << j << " " << node_cov_hash[j] << std::endl;  // j节点连接处的丰度
        mate_node_hash[j] = 0;

        for (auto it = node_hash.begin(); it != node_hash.end(); ++it) {
            if (it->first < j || it->second <= 1)
                continue;
            if (it->first > mate_node_hash[j])
                mate_node_hash[j] = it->first;

            ofile1 << it->first << " " << it->second << " reads_support: " << node_cov_hash[it->first] << std::endl;  // 
        }

        //      /*
        int l = j, sum = 0;
        float cov = 0, cov_p = 0, cov_v = 0, mw = 0, C = 1;
        if (create_non_canonical)
            C += 0.5;
        while (!node_set[l].children.empty() && mate_node_hash[j] != 0) {

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
                if (score <= C)
                    break;

                cov += node_cov_hash[node_set[l].children[0]];
                cov_p = cov / sum;
            }
            l = node_set[l].children[0];

        }

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

    //根据转录本路径重建转录本
    std::cout << "Begin rebuilding the transcript..." << std::endl;
    //int k = 1;

    for (int j = 1; j < node_set.size(); ++j) {   //遍历所有的没有父节点的节点
        if (!node_set[j].parents.empty() || mate_node_hash[j] == 0)
            continue;
        std::string transcript;
        if (j == 1 && !create_non_canonical)
            transcript = node_set[j].sequence;

        int l = j;
        int cov = 0, i = 0;
        while (!node_set[l].children.empty() && node_set[l].children[0] <= mate_node_hash[j]) { //

            transcript += node_set[node_set[l].children[0]].sequence;
            l = node_set[l].children[0];
        }

        ofile2 << "nt_" << j << std::endl;
        if (create_non_canonical && j != 1) {
            int current_kmer_cov, sum_kmer_cov = 0;
            float cov_p = kmer_hash[kmer_to_intval(transcript.substr(0, g_kmer_length))];
            int i = 0;
            for (; i < transcript.length() - g_kmer_length; i++) {
                kmer_int_type_t kmerint = kmer_to_intval(transcript.substr(i, g_kmer_length));
                current_kmer_cov = kmer_hash[kmerint];
                ofile2 << i + 1 << " current_kmer_cov " << current_kmer_cov << " cov_p " << cov_p << " " << fabs(1 - current_kmer_cov / cov_p) << std::endl;
                if (fabs(1 - current_kmer_cov / cov_p) > 0.3) {

                    break;
                }
                sum_kmer_cov += current_kmer_cov;
                cov_p = sum_kmer_cov / (i + 1);
            }
            transcript_list.push_back(transcript.substr(0, i + 1));
        }
        else {
            if (j == 1) {
                int start = 0, end1 = 0, end2 = 0;
                for (int i = 0; i < transcript.length() - 3; i++) {
                    std::string s = transcript.substr(i, 3);
                    if (s == "ATG" && i > 200) {
                        start = i;
                        break;
                    }
                }
                for (int i = start; i < transcript.length() - 3; i += 3) {
                    std::string s = transcript.substr(i, 3);

                    if (s == "TAA") {
                        end1 = i + 3;
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
                    transcript_list.push_back(transcript.substr(start, end1 - start));

                }
                if (end2 != 0) {
                    transcript_list.push_back(transcript.substr(start, end2 - start));

                }
                if (end2 == 0) {
                    transcript_list.push_back(transcript);


                }
            }

            else {
                for (int i = 0; i <= transcript.length() - 3; i += 3) {
                    std::string stop_codon = transcript.substr(i, 3);


                    if (stop_codon == "TAA" || stop_codon == "TAG" || stop_codon == "TGA") {
                        transcript = transcript.substr(0, i + 3);
                        if (transcript.length() <= 100)
                            continue;
                        transcript_list.push_back(transcript);

                        break;
                    }

                }
            }
        }
    }

    time_t end = time(NULL);
    std::cout << "Finished find_transcripts " << transcript_list.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;
}


void find_jump_reads(std::vector<std::string>& read_list, std::string trunk, std::vector<std::string>& jump_reads, bool flag) {

    std::cout << "Beginning find jump read..." << std::endl;
    time_t beg = time(NULL);
    int read_len = read_list[0].length();
    boost::unordered_map<std::string, int> kmer_pos_hash;
    for (int i = 0; i <= trunk.length() - g_kmer_length; i++) {
        std::string kmer = trunk.substr(i, g_kmer_length);
        kmer_pos_hash[kmer] = i;
    }

    for (int i = 0; i < read_list.size(); i++) {
        std::string read = read_list[i];
        if (!flag)
            read = revcomp(read);
        std::string first_kmer = read.substr(0, g_kmer_length);
        std::string last_kmer = read.substr(read_len - g_kmer_length);
        auto it_first = kmer_pos_hash.find(first_kmer);
        auto it_last = kmer_pos_hash.find(last_kmer);
        if (it_first == kmer_pos_hash.end() || it_last == kmer_pos_hash.end())
            continue;
        int first_pos = it_first->second;
        int last_pos = it_last->second;
        if (last_pos < first_pos || last_pos - first_pos < 3 + read_len)
            continue;
        jump_reads.push_back(read);
    }
    kmer_pos_hash.clear();
    time_t end = time(NULL);
    std::cout << "Finished find jump reads " << jump_reads.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;
}

void process_chunk(const Read_Hash& read_hash, int start, int end) {
    std::vector<std::string> local_list;
    auto it = std::next(read_hash.begin(), start);
    for (int i = start; i < end && it != read_hash.end(); ++i, ++it) {
        local_list.push_back(it->first);
    }

    {
        std::lock_guard<std::mutex> lock(hash_mutex);
        read_list.insert(read_list.end(), local_list.begin(), local_list.end());
    }
}

void get_reads(Read_Hash& read_hash) {

    std::cout << "Begining get reads..." << std::endl;
    time_t beg = time(NULL);
    int size = read_hash.size();
    int num_threads = std::thread::hardware_concurrency();
    int chunk_size = size / num_threads;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; i++) {
        int start = i * chunk_size;
        int end = (i == num_threads - 1) ? size : start + chunk_size;

        threads.emplace_back(process_chunk, std::ref(read_hash), start, end);
    }

    for (auto& thread : threads) {
        thread.join();
    }
    time_t end = time(NULL);
    std::cout << "Finished, get " << read_list.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;
}


void construct_graph(std::vector<std::string>& kunitigs_list, Kmer_Hash& kmer_hash, std::vector<std::string>& read_list, std::vector<std::string>& paired_read_list) {

    std::cout << "Beginning get_transcript ..." << std::endl;
    time_t beg = time(NULL);
    bool flag = true;
    bool create_non_canonical = false;
    std::sort(kunitigs_list.begin(), kunitigs_list.end(), [](std::string& str1, std::string& str2) { return str1.length() > str2.length(); });

    std::string trunk = kunitigs_list[0];
    int second_length = 0, other_length_sum = 0;
    int used_kmer_num = trunk.length() - g_kmer_length + 1;
    Kmer_Hash used_kmer_hash;

    if (kunitigs_list.size() > 1) {
        second_length = kunitigs_list[1].length();

        for (int i = 0; i < kunitigs_list.size(); ++i) {

            std::string sequence = kunitigs_list[i];

            for (size_t j = 0; j <= sequence.length() - g_kmer_length; ++j) {
                const std::string& kmer = sequence.substr(j, g_kmer_length);
                kmer_int_type_t kmer_val = kmer_to_intval(kmer, g_kmer_length);
                used_kmer_hash[kmer_val] += 1;
            }
        }

        for (int j = 1; j < kunitigs_list.size(); ++j) {
            other_length_sum += kunitigs_list[j].length();
        }
    }


    std::vector <std::string> transcript_list;
    if (trunk.length() - g_kmer_length + 1 > 0.9 * used_kmer_hash.size() || trunk.length() > 20 * second_length || trunk.length() > 5 * other_length_sum) {
        std::string fileName = out_dir + "node_set.fasta";
        std::ofstream outfile(fileName);
        std::string fileName2 = out_dir + "transcript.fasta";
        std::ofstream ofile2(fileName2);

        std::vector<Node> node_set;
        for (int i = 0; i < 10; i++) {
            if (trunk.substr(trunk.length() - i - 7, 7) == "AAAAAAA")
                break;
            if (trunk.substr(i, 5) == "TTTTT") {
                flag = false;
                trunk = revcomp(trunk);
                break;
            }
        }

        transcript_list.push_back(trunk);
        Node node(trunk);
        int p = add_node(node, node_set);

        std::vector<std::string> jump_reads;
        find_jump_reads(read_list, trunk, jump_reads, flag);
        read_list.clear();

        std::vector <int> jun_poss;
        check_branch(p, jump_reads, kmer_hash, node_set, flag, jun_poss);

        // /*
        for (int i = 1; i < node_set.size(); i++) {
            outfile << "> " << i << std::endl;
            outfile << node_set[i].sequence << std::endl;
        }
        outfile.close();
        // */

        find_path(node_set, paired_read_list, kmer_hash, flag, transcript_list, create_non_canonical);

        for (int i = 0; i < transcript_list.size(); i++) {
            ofile2 << ">" << "transcript_" << i << std::endl;
            ofile2 << transcript_list[i] << std::endl;
        }

        node_set.clear();
        transcript_list.clear();

        node_size = 0;
        if (g_non_canonical == 1) {
            create_non_canonical = true;
            Node node(trunk);
            p = add_node(node, node_set);

            check_non_canonical_branch(p, jump_reads, kmer_hash, node_set, flag, jun_poss);
            find_path(node_set, paired_read_list, kmer_hash, flag, transcript_list, create_non_canonical);

            for (int i = 0; i < transcript_list.size(); i++) {
                ofile2 << ">" << "non_canonical_transcript_" << i << std::endl;
                ofile2 << transcript_list[i] << std::endl;
            }

        }
        ofile2.close();


    }

    else {
        int k = 0;
        std::string fileName2 = out_dir + "transcript.fasta";
        std::ofstream ofile2(fileName2);

        int genome_len = 0.9 * used_kmer_hash.size() + g_kmer_length;

        for (int i = 0; i < kunitigs_list.size(); i++) {
            std::string trunk = kunitigs_list[i];
            if (trunk.length() > 0.65 * genome_len) {
                int start = 0, end1 = 0;
                std::string fkmer = trunk.substr(100, g_kmer_length);
                kmer_int_type_t fkmer_val = kmer_to_intval(fkmer);
                int cov1 = kmer_hash[fkmer_val];
                std::string lkmer = trunk.substr(trunk.length() - 100, g_kmer_length);
                kmer_int_type_t lkmer_val = kmer_to_intval(lkmer);
                int cov2 = kmer_hash[lkmer_val];
                if (cov1 < cov2)
                    trunk = revcomp(trunk);

                for (int i = 0.5 * trunk.length(); i < 0.7 * trunk.length() - 3; i++) {
                    std::string kmer = trunk.substr(i, g_kmer_length);

                    if (kmer.substr(g_kmer_length - 3) == "TAA") {
                        end1 = i + g_kmer_length;
                        kmer_int_type_t kmer_value = kmer_to_intval(kmer);
                        int cov = kmer_hash[kmer_value];
                        std::string kmer2 = trunk.substr(i + g_kmer_length, g_kmer_length);
                        kmer_int_type_t kmer_value2 = kmer_to_intval(kmer2);
                        int cov2 = kmer_hash[kmer_value2];
                        //   std::cout << end1 << " " << cov << " " << cov2 << " " << cov2 * 1.0 / cov << std::endl;
                        if (fabs(1 - cov2 * 1.0 / cov) < 0.3 || fabs(cov - cov2) < 10)
                            continue;
                        else
                            break;
                    }
                }

                transcript_list.push_back(trunk.substr(0, end1));
                ofile2 << ">" << "transcript_" << k++ << std::endl;
                ofile2 << trunk.substr(0, end1) << std::endl;

                transcript_list.push_back(trunk);
                ofile2 << ">" << "transcript_" << k++ << std::endl;
                ofile2 << trunk << std::endl;

            }
            else {
                transcript_list.push_back(trunk);
                ofile2 << ">" << "transcript_" << k++ << std::endl;
                ofile2 << trunk << std::endl;
            }
        }
        std::cout << "Finished get " << transcript_list.size() << " transcripts " << std::endl;
    }

    time_t end = time(NULL);
    std::cout << "Finished ,(elapsed time: " << (end - beg) << " s)" << std::endl;
}



bool forward_similar(const std::string& str1, const std::string& str2) {
    int i = 0;
    int j = 0;
    int len = (str1.length() < str2.length()) ? str1.length() : str2.length();
    int mismatch = 0.02 * len;

    int distance = editDistance(str1, str2);
    return distance <= mismatch;

}


std::vector<std::string> remove_isContained(std::vector<std::string>& strings) {

    std::vector<std::string> filteredStrings;
    std::sort(strings.begin(), strings.end(), [](const std::string& str1, const std::string& str2) {return str1.length() > str2.length(); });
    for (int i = 0; i < strings.size(); i++) {
        std::string str = strings[i];
        bool isContained = false;
        for (const std::string& otherStr : filteredStrings) {
            if (str == otherStr || is_similar(otherStr, str) || is_similar(otherStr, revcomp(str))) {
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

std::vector<std::string> remove_similar(std::vector<std::string>& strings) {

    std::vector<std::string> filteredStrings;
    std::sort(strings.begin(), strings.end(), [](const std::string& str1, const std::string& str2) {return str1.length() > str2.length(); });
    for (const std::string& str : strings) {
        bool isContained = false;
        for (const std::string& otherStr : filteredStrings) {
            if (forward_similar(otherStr, str) || forward_similar(otherStr, revcomp(str))) {    //  
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



int haveCommonSequence(const std::string& str1, const std::string& str2, int commonLength) {
    if (str1.length() < commonLength || str2.length() < commonLength) {
        return false;
    }
    std::string kmer = str1.substr(str1.length() - commonLength);
    int pos = str2.find(kmer);
    if (pos == std::string::npos || pos > g_kmer_length)
        return -1;
    if (str2.substr(0, pos + commonLength) == str1.substr(str1.length() - pos - commonLength))
        return pos + commonLength;
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

void remove_revcomp(std::vector<std::string>& kunitigs_list) {

    std::cerr << "Beginning remove_revcomp  ..." << std::endl;
    time_t beg = time(NULL);


    kunitigs_list = remove_isContained(kunitigs_list);

    std::cout << "kunitigs_list count after remove_isContained : " << kunitigs_list.size() << std::endl;

    if (kunitigs_list.size() > 1)
        kunitigs_list = remove_similar(kunitigs_list);

    head_tail_connection(kunitigs_list);

    time_t end = time(NULL);
    /*
    int j = 0;
    std::string fileName = out_dir + "kunitigs2.fasta";
    std::ofstream outfile(fileName);

    if (outfile.is_open()) {
        for (const auto& pair : kunitigs_list) {
            outfile << "> " << j++ << std::endl;
            outfile << pair << "\n";
        }
    }
    outfile.close();
    */
    std::cout << "kunitigs_list count after remove_revcomp : " << kunitigs_list.size() << " (elapsed time: " << (end - beg) << " s)" << std::endl;
}


void Construct_splicing_graphs(std::vector<std::string>& Read) {

    get_read_hash(read_hash, Read);
    std::vector<std::string> paired_read_list;
    get_paired_read_list(paired_read_list, Read);

    get_reads(read_hash);
    get_kmer_hash(kmer_hash, read_hash);

    remove_erroneous(kmer_hash, min_ratio_non_error);

    std::vector<std::string> kunitigs_list;
    extend_kunitigs(kmer_hash, kunitigs_list);
    remove_revcomp(kunitigs_list);

    construct_graph(kunitigs_list, kmer_hash, read_list, paired_read_list);

}
// kmc_histogram_fixed.cpp
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include "kmc_file.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <kmc_database> [output.csv]" << std::endl;
        return 1;
    }

    std::string db_path = argv[1];
    std::string output_file = argc > 2 ? argv[2] : "histogram.csv";

    std::cout << "Opening database: " << db_path << std::endl;

    CKMCFile kmc_file;
    
    if (!kmc_file.OpenForListing(db_path)) {
        std::cerr << "Error: Cannot open KMC database: " << db_path << std::endl;
        return 1;
    }

    std::cout << "Database opened successfully" << std::endl;

    uint32 kmer_length = kmc_file.KmerLength();
    std::cout << "K-mer length: " << kmer_length << std::endl;

    uint64 mode, counter_size, lut_prefix_length, signature_len, min_count, max_count, total_kmers;

    std::unordered_map<uint32, uint64_t> hist;
    CKmerAPI kmer(kmer_length);
    uint32 count;
    uint64_t total = 0;


    kmc_file.RestartListing();

    while (kmc_file.ReadNextKmer(kmer, count)) {
        hist[count]++;
        total++;
    }

    kmc_file.Close();

    std::cout << "\nFinished! Total k-mers read: " << total << std::endl;

    if (total == 0) {
        std::cerr << "Error: No k-mers were read from the database!" << std::endl;
        return 1;
    }

    std::ofstream out(output_file);
    if (!out) {
        std::cerr << "Cannot create output file: " << output_file << std::endl;
        return 1;
    }

    out << "count,num_kmers\n";
    
    
    for (const auto& p : hist) {
        out << p.first << "," << p.second << "\n";
    }
    out.close();

    std::cout << "Success!" << std::endl;

    return 0;
}
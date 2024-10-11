#include <algorithm>
#include <cerrno>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>

namespace bio = boost::iostreams;

std::map<char, char> nucleotide_complements = {
    {'A', 'T'},
    {'T', 'A'},
    {'C', 'G'},
    {'G', 'C'},
    {'N', 'N'},
    {'a', 't'},
    {'t', 'a'},
    {'c', 'g'},
    {'g', 'c'},
    {'n', 'n'},
};

std::string complement(const std::string& seq) {
    std::string complemented;
    for (const auto it : seq) {
        complemented += nucleotide_complements.at(it);
    }
    return complemented;
}

std::string reverse_complement(const std::string& seq) {
    std::string rc;
    for (int i = seq.length() - 1; i >= 0; i--) {
        rc += nucleotide_complements.at(seq[i]);
    }
    return rc;
}

class FileException: public std::runtime_error {
public:
    FileException() : std::runtime_error("") { }
    explicit FileException(std::string msg) : std::runtime_error(msg) { }
};

bool is_gzipped_file(const std::string& filename) {
    bool gzipped = false;
    FILE* f = fopen(filename.c_str(), "rb");

    if (f == NULL) {
        throw FileException("Could not open file \"" + filename + "\": " + strerror(errno));
    } else {
        if (fgetc(f) == 0x1f && fgetc(f) == 0x8b) {
            gzipped = true;
        }
    }
    fclose(f);
    return gzipped;
}

bool is_gzipped_filename(const std::string& filename) {
    std::string ext = ".gz";
    return std::equal(ext.rbegin(), ext.rend(), filename.rbegin());
}

void log_message(const std::string& message) {
    std::time_t now = std::time(0);
    std::tm* local = std::localtime(&now);
    char buffer[80];
    std::strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", local);
    std::cerr << buffer << " " << message << std::endl;
}

class FASTQRecord {
public:
    std::string name;
    std::string comment;
    std::string sequence;
    std::string quality;

    void transform(int offset, int length, bool rc) {
        sequence = sequence.substr(offset, length);
        quality = quality.substr(offset, length);
        if (rc) {
            sequence = reverse_complement(sequence);
            std::reverse(quality.begin(), quality.end());
        }
    }

};

std::ostream& operator<<(std::ostream& os, const FASTQRecord& record) {
    os << record.name << '\n'
       << record.sequence << '\n'
       << record.comment << '\n'
       << record.quality << '\n';
    return os;
}

std::istream& operator>>(std::istream& is, FASTQRecord& record) {
    getline(is, record.name);
    getline(is, record.sequence);
    getline(is, record.comment);
    getline(is, record.quality);
    return is;
}

std::vector<FASTQRecord> read_fastq(const std::string& filename, long long int first_n = -1) {

    bio::file_source input_file(filename);
    bio::stream<bio::file_source> input_stream(input_file);

    bio::filtering_stream<bio::input> in;
    if (is_gzipped_file(filename)) {
        in.push(bio::gzip_decompressor());
    }
    in.push(input_stream);

    FASTQRecord record;
    std::vector<FASTQRecord> records;
    long long int record_count = 0;

    while (in >> record) {
        record_count++;
        if (first_n != -1 && record_count > first_n) {
            break;
        }
        records.push_back(record);
    }

    bio::close(in);

    return records;
}


std::unordered_set<std::string> read_whitelist(const std::string& whitelist_filename) {
    std::ifstream whitelist_file(whitelist_filename);
    std::unordered_set<std::string> whitelist_barcodes;
    std::string line;

    while (std::getline(whitelist_file, line)) {
        whitelist_barcodes.insert(line);
    }

    whitelist_file.close();

    return whitelist_barcodes;
}


std::vector<int> determine_transform(const std::string& whitelist_filename, const std::string& barcode_fastq_file, const int check_first, bool verbose = false) {

    if (verbose) {
        log_message("Determining transform");
        log_message("Reading whitelist from " + whitelist_filename + "...");
    }
    std::unordered_set<std::string> whitelist_barcodes = read_whitelist(whitelist_filename);
    std::unordered_set<std::string> whitelist_barcodes_rc;

    for (const auto& barcode : whitelist_barcodes) {
        whitelist_barcodes_rc.insert(reverse_complement(barcode));
    }
    
    if (verbose) {
        log_message("Whitelist size: " + std::to_string(whitelist_barcodes.size()));
    }

    // determine barcode lengths (from whitelist_barcodes)
    int barcode_length;
    std::unordered_set<int> barcode_lengths;

    for (const auto& barcode : whitelist_barcodes) {
        barcode_lengths.insert(barcode.size());
    }

    if (barcode_lengths.size() > 1) {
        throw std::runtime_error("Barcodes in the whitelist are not all the same length.");
    } else {
        barcode_length = *barcode_lengths.begin();
    }

    if (verbose) {
        log_message("Inferred barcode length: " + std::to_string(barcode_length));
    }
    
    // Determine the transformation
    if (verbose) {
        log_message("Reading the first " + std::to_string(check_first) + " records from " + barcode_fastq_file + "...");
    }
    std::vector<FASTQRecord> records = read_fastq(barcode_fastq_file, check_first);

    std::map<int, int> match_counts; // offset --> match_count
    std::map<int, int> rc_match_counts; // offset --> match_count
    int best_offset = 0;
    bool best_rc = false;
    int best_match_count = 0;

    for (const auto& record : records) {
        for (size_t offset = 0; offset < (record.sequence.size() - barcode_length + 1); offset++) {
            std::string subsequence = record.sequence.substr(offset, barcode_length);
            if (whitelist_barcodes.find(subsequence) != whitelist_barcodes.end()) {
                match_counts[offset]++;
            }
            if (whitelist_barcodes_rc.find(subsequence) != whitelist_barcodes_rc.end()) {
                rc_match_counts[offset]++;
            }
        }
    }

    for (const auto& it : match_counts) {
        if (it.second > best_match_count) {
            best_match_count = it.second;
            best_offset = it.first;
            best_rc = false;
        }
    }

    for (const auto& it : rc_match_counts) {
        if (it.second > best_match_count) {
            best_match_count = it.second;
            best_offset = it.first;
            best_rc = true;
        }
    }

    size_t record_count = records.size();

    if (verbose) {
        log_message("Best offset: " + std::to_string(best_offset) + ", best rc: " + std::to_string(best_rc) + ", best match count: " + std::to_string(best_match_count) + " out of " + std::to_string(record_count) + " records" + " (" + std::to_string(100*(double)best_match_count / record_count) + "%)");
    }

    int best_rc_bool = best_rc ? 1 : 0;
    
    return std::vector<int> {best_offset, barcode_length, best_rc_bool};

}

void print_usage() {
    std::cout << "Transform and count cell barcodes from a 10X snATAC-seq library.\n\n"

              << "Usage:\n\n" << "transform-and-count-barcodes [options] input_file barcode_whitelist output_fastq output_counts\n\n"
              << "where:\n"
              << "    input_file is the fastq file of barcode reads\n"
              << "    barcode_whitelist is the barcode whitelist\n"
              << "    output_fastq will be the fastq file of barcodes\n"
              << "    output_counts will be the file of barcode counts\n\n"

              << "Options may include:\n\n"

              << "-h|--help: show this usage message.\n\n"
              << "-v|--verbose: show more details and progress updates." << std::endl << std::endl;

}

int main(int argc, char **argv)
{
    int c, option_index = 0;
    bool verbose = 0;

    std::string input_filename;
    std::string whitelist_filename;
    std::string output_fastq;
    std::string output_counts;

    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", no_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // parse the command line arguments
    while ((c = getopt_long(argc, argv, "vh?", long_options, &option_index)) != -1) {
        switch (c) {
        case 'h':
        case '?':
            print_usage();
            exit(1);
        case 'v':
            verbose = true;
            break;
        case 0:
            if (long_options[option_index].flag != 0){
                break;
            }

        default:
            print_usage();
            exit(1);
        }
    }

    if (argc != (optind + 4)) {
        print_usage();
        exit(1);
    }

    input_filename = argv[optind];
    whitelist_filename = argv[optind + 1];
    output_fastq = argv[optind + 2];
    output_counts = argv[optind + 3];

    std::vector<int> best_transform = determine_transform(whitelist_filename, input_filename, 10000);

    int best_offset = best_transform[0];
    int barcode_length = best_transform[1];
    bool best_rc = best_transform[2];

    std::unordered_set<std::string> whitelist = read_whitelist(whitelist_filename);

    if (verbose) {
        log_message("Transforming records...");
    }

    // input stream
    bio::file_source input_file(input_filename);
    bio::stream<bio::file_source> input_stream(input_file);

    bio::filtering_stream<bio::input> in;
    if (is_gzipped_file(input_filename)) {
        in.push(bio::gzip_decompressor());
    }
    in.push(input_stream);


    //output fastq
    std::ofstream output_file(output_fastq, std::ios_base::out | std::ios_base::binary);

    bio::filtering_stream<bio::output> out;
    if (is_gzipped_filename(output_fastq)) {
        out.push(bio::gzip_compressor(bio::gzip_params(bio::gzip::best_speed)));
    }
    out.push(output_file);


    // transform
    FASTQRecord record;
    std::unordered_map<std::string, long long int> counts;
    long long int record_count = 0;

    while (in >> record) {
        record_count++;
        if (record_count % 1000000 == 0 && verbose) {
            log_message("Processed " + std::to_string(record_count) + " records so far...");
        }
        record.transform(best_offset, barcode_length, best_rc);
        counts[record.sequence]++;
        out << record;
    }

    bio::close(in);
    bio::close(out);

    // output counts
    std::ofstream counts_file(output_counts);
    for (const auto& it : counts) {
        counts_file << it.first << "\t" << it.second << std::endl;
    }

    counts_file.close();

    if (verbose) {
        log_message("Transformed " + std::to_string(record_count) + " records.");
        log_message("Done.");
    }

}

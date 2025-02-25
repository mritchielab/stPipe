//' @name RunDemultiplex
 //' @title Demultiplex Stereo-seq data
 //' @description This function processes sequencing data for spatial transcriptomics.
 //' @param read_1_fq Path to the first FASTQ file.
 //' @param read_2_fq Path to the second FASTQ file.
 //' @param h5_mapping Path to the HDF5 barcode mapping file.
 //' @param output_fq Path to the output FASTQ file.
 //' @param n_reads Number of reads to process.
 //' @param bc_start Start position of barcode.
 //' @param bc_len Length of barcode.
 //' @param umi_start Start position of UMI.
 //' @param umi_len Length of UMI.
 //' @param bin_size Binning size of n * n.
 //' @return None. Writes the demultiplexed FASTQ file to the specified path.
 //' @export

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <iterator>
#include "H5Cpp.h"
#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include "progressbar.hpp"
#include <zlib.h>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include "Rcpp.h"
using namespace Rcpp;

using std::cout;
using std::endl;
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

// uncomment next line to enable debugging:
// #define _DEBUGGING true

#ifndef _DEBUGGING
#define _DEBUGGING false
#endif

// uncomment next line to enable progress bar (ONLY USE IF SINGLE THREADED):
//#define _PROGRESS true

#ifndef _PROGRESS
#define _PROGRESS false
#endif

unsigned int CoordPairToInt(unsigned int a, unsigned int b) {
    return (a << 16) + b;
}

KSEQ_INIT(gzFile, gzread);

unsigned long int seq_to_int(std::string seq) {
    std::map<char, unsigned long int> base_to_int_map{{'A', 0UL}, {'C', 1UL}, {'T', 2UL}, {'G', 3UL}};
    unsigned long int acc = 0;
    for (int i = seq.length() - 1; i >= 0; i--) {
        acc += base_to_int_map[seq[i]] << (i * 2);
    }
    return acc;
}

class ReadStats {
    public:
        int total, clean_hit, miss, err_hit;
        ReadStats();
        double getHitAcc();
        void reportStats();
        void reportStatsR();
};

ReadStats::ReadStats() {
    total = 0;
    clean_hit = 0;
    miss = 0;
    err_hit = 0;
}

double ReadStats::getHitAcc() {
    return (1.0 - miss/total);
}

void ReadStats::reportStats() {
    cout << "Total number of reads: " << total << endl;
    cout << "Clean hits: " << clean_hit << endl;
    cout << "Single Error hits: " << err_hit << endl;
    cout << "Misses: " << miss << endl;
    cout << "Total hits: " << (clean_hit+err_hit) << endl;
}

void ReadStats::reportStatsR() {
    Rcpp::Rcout << "Total number of reads: " << total << endl;
    Rcpp::Rcout << "Clean hits: " << clean_hit << endl;
    Rcpp::Rcout << "Single Error hits: " << err_hit << endl;
    Rcpp::Rcout << "Misses: " << miss << endl;
    Rcpp::Rcout << "Total hits: " << (clean_hit+err_hit) << endl;
}

/*
Rcpp exposed function to run entire method
IO:
    In:
    - fastq read 1 path (string)
    - fastq read 2 path (string)
    - h5 barcode mapping file path (string)
    - fastq output path (string)
    - number of reads in fastq (int)
    - start position of spatial barcode in each sequence (int)
    - length of spatial barcode (int)
    - start position of UMI in each sequence (int)
    - length of UMI (int)
    - bin size (int, default = 1)

    Out:
    - writes fastq file to output path (*.fq.gz)
        - equivalent to read 2 but with demultiplexed information in header:
            - CID (region of barcode defined by CID start pos and length)
            - UMI (region of barcode defined by UMI start pos and length)
            - x coordinate (from demultiplexing based on mapping file)
            - y coordinate (from demultiplexing based on mapping file)
*/
//' @export
// [[Rcpp::export]]
void RunDemultiplex(const char* read_1_fq_path, const char* read_2_fq_path, const char* h5_mapping_path, const char* output_fq_path, int n_reads, int coord_bc_start, int coord_bc_len, int umi_start, int umi_len, int bin_size = 1) {
    gzFile _fp = gzopen(read_1_fq_path, "r");
    gzFile _fp_r2 = gzopen(read_2_fq_path, "r");
    gzFile fp_write = gzopen(output_fq_path, "w");
    kseq_t* _seq = kseq_init(_fp);
    kseq_t* _seq_r2 = kseq_init(_fp_r2);

    H5File *file = new H5File(h5_mapping_path, H5F_ACC_RDONLY);

    if (_DEBUGGING) {
        Rcpp::Rcout << "file pointers:" << endl;
        Rcpp::Rcout << &_fp << endl;
        Rcpp::Rcout << &_fp_r2 << endl;
        Rcpp::Rcout << &fp_write << endl;
        Rcpp::Rcout << &file << endl;
    }

    DataSet dataset = file->openDataSet("bpMatrix_1");
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims_out[3];
    int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

    if (_DEBUGGING) {
        Rcpp::Rcout << "rank " << rank << ", dimensions " <<
                (unsigned long)(dims_out[0]) << " x " <<
                (unsigned long)(dims_out[1]) << endl;
    }

    unsigned long dim_x = (unsigned long)(dims_out[0]);
    unsigned long dim_y = (unsigned long)(dims_out[1]);
    unsigned long dim_z = (unsigned long)(dims_out[2]);
    DataSpace memspace(3, dims_out);
    Rcpp::Rcout << "x:" << dim_x << " y:" << dim_y << " z:" << dim_z << endl;

    std::vector<unsigned long> buf(dim_x*dim_y*dim_z);
    unsigned long *buf_ptr = buf.data();
    dataset.read(buf_ptr, PredType::NATIVE_ULONG, memspace, dataspace);
    std::unordered_map<unsigned long, unsigned int> barcode_map;
    std::set<unsigned long> duplicate_barcodes;

    int c = 0;
    int collisions = 0;
    progressbar bar;
    Rcpp::Rcout << "building map..." << endl;
    if (_PROGRESS) {
        progressbar bar((int) dim_x*dim_y*dim_z/10000);
    }
    for (unsigned long b : buf) {
        if (b == 0) {
            if (_PROGRESS) {
                if (c % 10000 == 0) bar.update();
            }
            c++;
            continue;
        }
        if (_PROGRESS) {
            if (c % 10000 == 0) bar.update();
        }
        if (barcode_map.find(b) != barcode_map.end()) {
            barcode_map[b] = CoordPairToInt(0, 0);
            duplicate_barcodes.insert(b);
            collisions++;
        } else {
            barcode_map[b] = CoordPairToInt((unsigned int) (c % dim_x), (unsigned int) (c / dim_x));
        }
        c++;
    }
    Rcpp::Rcout << "done!" << endl;

    if (_DEBUGGING) {
        Rcpp::Rcout << endl;
        Rcpp::Rcout << "max size: " << barcode_map.max_size() << endl;
        Rcpp::Rcout << buf[0] << " == " << ((barcode_map.find(buf[0])->second & 0xFFFF0000) >> 16) << ", " << (barcode_map.find(buf[0])->second & 0x0000FFFF) << endl;
        Rcpp::Rcout << buf[10] << " == " << ((barcode_map.find(buf[10])->second & 0xFFFF0000) >> 16) << ", " << (barcode_map.find(buf[10])->second & 0x0000FFFF) << endl;
        Rcpp::Rcout << buf[100000] << " == " << ((barcode_map.find(buf[100000])->second & 0xFFFF0000) >> 16) << ", " << (barcode_map.find(buf[100000])->second & 0x0000FFFF) << endl;
        Rcpp::Rcout << "num collisions: " << collisions << endl;
        Rcpp::Rcout << "duplicate barcode number: " << duplicate_barcodes.size() << endl;
        if (!duplicate_barcodes.empty()) {
            for (auto b : duplicate_barcodes) {
                Rcpp::Rcout << b << endl;
            }
        }
    }

    ReadStats stats;
    int l = 0;
    int l_r2 = 0;
    progressbar bar_2;
    Rcpp::Rcout << "beginning deconvolution..." << endl;
    if (_PROGRESS) {
        progressbar bar_2((int) n_reads/10000);
        // check for interrupt
        Rcpp::checkUserInterrupt();
    }

    for (int i = 0; (i < n_reads && l >= 0); i++) {
        l = kseq_read(_seq);
        l_r2 = kseq_read(_seq_r2);
        if (_PROGRESS) {
            if (i % 10000 == 0) bar_2.update();
        }
        if (_DEBUGGING) {
            Rcpp::Rcout << _seq->seq.s << endl;
        }

        std::string seq_str = _seq->seq.s;
        std::string trimmed = seq_str.substr(coord_bc_start, coord_bc_len);
        std::string umi = seq_str.substr(umi_start, umi_len);
        if (_DEBUGGING) {
            Rcpp::Rcout << umi << endl;
            Rcpp::Rcout << trimmed << endl;
        }
        int bit_mask_index = 0;
        unsigned long trimmed_int = seq_to_int(trimmed);
        unsigned long trimmed_masked_int = trimmed_int;
        // number of bases * 2 to get binary seq len
        int seq_len = trimmed.size() * 2;

        while ((barcode_map.find(trimmed_masked_int) == barcode_map.end()) & (bit_mask_index < seq_len)){
            trimmed_masked_int = trimmed_int ^ (1 << bit_mask_index);
            bit_mask_index++;
        }

        if (bit_mask_index < seq_len) {
            if (bit_mask_index == 0) stats.clean_hit = stats.clean_hit + 1;
            if (bit_mask_index > 0) stats.err_hit = stats.err_hit + 1;
            std::string x = std::to_string(((barcode_map.find(trimmed_masked_int)->second & 0xFFFF0000) >> 16) / bin_size);
            std::string y = std::to_string((barcode_map.find(trimmed_masked_int)->second & 0x0000FFFF) / bin_size);
            // read_id structure:
            // @[barcode sequence]_[umi sequence]#[read name]:x[x coordinate]:y[y coordinate]
            std::string read_id = "@" + trimmed + "_" + umi + "#" + (_seq_r2->name.s) + ":x" + x + ":y" + y;
            std::string out_buf = read_id + "\n" + (_seq_r2->seq.s) + "\n+\n" + (_seq_r2->qual.s) + "\n";

            gzwrite(fp_write, out_buf.c_str(), sizeof(char)*out_buf.size());
        } else {
            stats.miss = stats.miss + 1;
        }
        stats.total = stats.total + 1;
    }
    stats.reportStats();

    Rcpp::Rcout << "done!\n closing files." << endl;
    kseq_destroy(_seq);
    kseq_destroy(_seq_r2);
    gzclose(_fp);
    gzclose(_fp_r2);
    gzclose(fp_write);

    Rcpp::Rcout << "Demultiplexing completed." << std::endl;

}

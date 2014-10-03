// This file is the main file for indexing all RNA signatures.

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "proto_data.h"
#include "fa_reader.h"
#include "proto/rnasigs.pb.h"
#include "rs_bloom.h"
#include "rs_common.h"
#include "rs_thread.h"

using std::string;
using std::pair;
using std::vector;
using std::fstream;
using std::ios;

DEFINE_string(transcript_fasta, "",
              "The fasta file for transcript sequences of every gene. This is suggested to use the output of rs_cluster. We use a specialized FASTA format, in which each line represents a gene and its transcripts.  A transcript sequence per line is not recommended. ");
DEFINE_string(index_file, "",
              "The path to the index file (output).");
DEFINE_int32(num_threads, -1,
           "The num of threads used in the program. "
           "[default: -1]: the num of CPUs in the machine.");
DEFINE_int32(rs_length, 40,
           "The length of the RS signature.");


namespace rs {

// This thread gets all duplicated substrings in a bloom.
class IndexThread : public ThreadInterface {
public:
  IndexThread(SingleFastaReader* reader, RSBloom* all_bloom, RSBloom* dup_bloom)
    : reader_(reader), all_bloom_(all_bloom), dup_bloom_(dup_bloom) {}

  void add(const string& seq, RSBloom* line_bloom) {
    for (int start = 0; start < (int)seq.size() - FLAGS_rs_length;
         start ++) {
      StringPiece temp(seq.data() + start, FLAGS_rs_length);
      // This is the first rs signature in this gene.

      if (line_bloom->add(temp)) {
        // This is NOT the first  rs signature in this scan.
        if (!all_bloom_->add(temp)) {
          // add to the duplicated set
          dup_bloom_->add(temp);
        }
      }
    }
  }

  void run() {
    vector<string> ids, seqs;

    while (reader_->read(&ids, &seqs) != 0) {
      for (size_t i = 0; i < ids.size(); i++) {
        string whole_seq = seqs[i];
        // TODO(zzj): memory leak, need to figure out why it crashes
        // when the whole_seq size is huge (e.g. 17M)
        RSBloom *rsb = new RSBloom(whole_seq.size() * 4, 0.001);
        vector<string> seqs = split_seq(whole_seq, '|');
        for (auto& seq : seqs) {

          // add four different possible sequences into bloom
          // normal
          add(seq, rsb);

          // complimentary
          compliment(&seq);
          add(seq, rsb);

          // complimentary and reversed
          reverse(seq.begin(), seq.end());
          add(seq, rsb);

          // reversed and normal
          compliment(&seq);
          add(seq, rsb);
        }
      }
    }
  }

private:
  SingleFastaReader* reader_;
  RSBloom* all_bloom_;
  RSBloom* dup_bloom_;
};

// This class is used for dumping the data into a file
// This is thread safe
// This should be used as a singleton.
class IndexDumper {
public:
  IndexDumper(const string& filename)
    : stream_(filename, ios::out | ios::binary | ios::trunc) {
  }
  void dump(const vector<GeneSignatures>& data) {
    std::lock_guard<std::mutex> lock(m_);
    for (auto& s : data) {
      write_protobuf_data(&stream_, &s);
    }
  }
private:
  fstream stream_;
  std::mutex m_;
};

// This thread runs after the IndexThread and will find out all unique
// strings of the transcripts.
class LookupUniqueStringThread : public ThreadInterface {
public:
  LookupUniqueStringThread(SingleFastaReader* reader, RSBloom* dup_bloom,
                           IndexDumper* dumper)
    : reader_(reader), dup_bloom_(dup_bloom), dumper_(dumper) {}

  void run() {
    vector<string> ids, seqs;
    while(reader_->read(&ids, &seqs) != 0) {
      for (size_t i = 0; i < ids.size(); i++) {
        string whole_seq = seqs[i];
        string header = ids[i];

        // first one is the gene id, the remaining ones are transcript ids.
        vector<string> ids = split_seq(header, '|');
        vector<string> seqs = split_seq(whole_seq, '|');
        if (ids.size() != seqs.size() + 1) {
          LOG(ERROR) << "Wrong format for header: " << header;
          LOG(ERROR) << ids.size() << " " << seqs.size();
        }

        GeneSignatures result;
        result.set_id(ids[0]);
        for (size_t i = 0; i < seqs.size(); i++) {
          auto& seq = seqs[i];
          auto transcript = result.add_transcripts();
          transcript->set_id(ids[i + 1]);
          transcript->set_length(seq.size());
          bool is_consecutive = false;
          int last_start = -1;
          for (int start = 0; start < (int)seq.size() - FLAGS_rs_length;
               start ++) {
            StringPiece temp (seq.data() + start, FLAGS_rs_length);
            if (!dup_bloom_->contain(temp)) {
              if (!is_consecutive) {
                // record the first position of the consecutive region
                last_start = start;
              }
              is_consecutive = true;
            } else {
              if (is_consecutive) {
                // save the current consecutive region
                auto meta = transcript->add_signatures();
                meta->set_seq(seq.substr(last_start,
                                         start - last_start + FLAGS_rs_length));
                meta->set_position(last_start);
              }
              is_consecutive = false;
            }
          }
          if (is_consecutive) {
            auto meta = transcript->add_signatures();
            meta->set_seq(seq.substr(last_start));
            meta->set_position(last_start);
          }
        }
        gene_signatures_.push_back(result);
        if (gene_signatures_.size() > 10) {
          dumper_->dump(gene_signatures_);
          gene_signatures_.clear();
        }
      }
    }
    dumper_->dump(gene_signatures_);
  }

private:
  SingleFastaReader* reader_;
  RSBloom* dup_bloom_;
  IndexDumper* dumper_;
  std::vector<GeneSignatures> gene_signatures_;
};

class IndexMain {
public:
  IndexMain(const string& transcript_fasta_filename,
            const string& output_file,
            int num_threads)
    : file_(transcript_fasta_filename),
      output_file_(output_file), num_threads_(num_threads) {
    const uint64_t total_length = get_file_size(transcript_fasta_filename);
    LOG(INFO) << "The estimated total number of k-mers is " << total_length;
    all_bloom_ = new RSBloom(total_length * 10, 0.001);
    dup_bloom_ = new RSBloom(total_length * 10, 0.001);
  }

  void run() {
    LOG(INFO) << "Loading k-mers and transcripts";
    SingleFastaReader* reader_ = new SingleFastaReader(file_, 1);
    IndexThread index_thread(reader_, all_bloom_, dup_bloom_);
    std::vector<std::thread> threads(num_threads_);
    for (int i = 0; i < num_threads_; i ++) {
      threads[i] = std::thread{RSThread(&index_thread)};
    }
    barrier(threads);
    LOG(INFO) << "all k-mers are loaded into memory";
    reader_->reset();
    IndexDumper dumper(output_file_);
    LookupUniqueStringThread lookup_thread(reader_, dup_bloom_,
                                           &dumper);
    vector<LookupUniqueStringThread> lookup_threads(num_threads_,
                                                    lookup_thread);
    for (int i = 0; i < num_threads_; i ++) {
      threads[i] = std::thread{RSThread(&lookup_threads[i])};
    }
    barrier(threads);
  }

private:
  std::string file_;

  // every rs signature will be added to this bloom
  RSBloom* all_bloom_;

  // only duplicated rs signature will be added to this bloom
  RSBloom* dup_bloom_;

  string output_file_;

  int num_threads_;
};

}  // namespace rs

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_num_threads == -1) {
    FLAGS_num_threads = std::thread::hardware_concurrency();
  }
  rs::IndexMain rsim(FLAGS_gene_fasta, FLAGS_index_file,
                     FLAGS_num_threads);
  rsim.run();
}

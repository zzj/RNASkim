// This file count the occurrence of the RNA signatures

#include <algorithm>
#include <chrono>
#include <condition_variable>
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
#include "rs_common.h"
#include "rs_thread.h"
#include "rs_estimate_lib.h"
#include "rolling_hash_counter.h"

using std::fstream;
using std::ios;
using std::map;
using std::string;
using std::vector;

DEFINE_string(selected_keys_file, "",
              "The path to the selected keys file (input).");
DEFINE_string(count_file, "",
              "The path to the file that contains all selected keys with their counts (output).");
DEFINE_int32(num_threads, -1,
           "The num of threads used in the program. "
           "[default: -1]: the num of CPUs in the machine.");
DEFINE_string(read_files1, "",
              "The fasta read files, splitted by ','");
DEFINE_string(read_files2, "",
              "The fasta read files, splitted by ','");
DEFINE_int32(rs_length, 40,
           "The length of the RS signature.");
DEFINE_bool(run_em, false,
           "Whether to run EM when counting.");
DEFINE_bool(fastq, false,
           "Whether the data is fastq format");

namespace rs {

void update_SelectedKey(const RollingHashCounter& counter, SelectedKey* sk) {
  for (int i = 0; i < sk->keys_size(); i++) {
    int count = 0;
    string key = sk->keys(i).key();

    count += counter.find(key);
    compliment(&key);
    count += counter.find(key);
    reverse(key.begin(), key.end());
    count += counter.find(key);
    compliment(&key);
    count += counter.find(key);
    sk->mutable_keys(i)->set_count(count);
  }
}

class CountThread : public ThreadInterface {
public:
  CountThread(RSPairReader* reader, RollingHashCounter* counter)
    : reader_(reader), counter_(counter) {}

  void run() {
    vector<string> reads1, reads2;
    int total = 0;
    while (reader_->read(&reads1, &reads2) != 0) {
      for (uint32_t i = 0; i < reads1.size(); i++) {
        // since the counter contains all four different keys,
        // here, we only need to process the sequence once.
        counter_->process(reads1[i]);
        counter_->process(reads2[i]);
      }
      total += reads1.size();
      //LOG(INFO) << "Processed " << total;
    }
  }

private:
  RSPairReader* reader_;
  RollingHashCounter* counter_;
};

// This should be only one thread
class EMThread : public ThreadInterface {
public:
  EMThread(const RollingHashCounter* counter,
           vector<SelectedKey>* selected_keys,
           std::atomic<bool>* is_running)
    : selected_keys_(selected_keys), counter_(counter),
      is_running_(is_running) {
    pi_.resize(selected_keys->size());
  }
  void run_em() {
    for (size_t i = 0; i < selected_keys_->size(); i++) {
      SignatureInfoDB db;
      SelectedKey& sk = selected_keys_->at(i);
      update_SelectedKey(*counter_, &sk);
      if (prepare_SignatureInfoDB(sk, &db)) {
        if (pi_[i].size() == 0) {
          pi_[i].resize(sk.tids_size(), 1.0 / sk.tids_size());
        }
        vector<double> count_per_tid;
        EM(sk.tids_size(), db, &count_per_tid, &pi_[i]);
        for (int j = 0; j < sk.tids_size(); j++) {
          profile_[sk.tids(j)] = count_per_tid[j] * sk.lengths(j);
        }
      }
      if (!(*is_running_)){
        break;
      }
    }
  }
  void run() {
    while (*is_running_) {
      int tick = 0;
      // wake up every second to see whether the counting threads are running.
      while (*is_running_) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        tick ++;
        if (tick > 10) break;
      }
      if (*is_running_)
        run_em();
    }
    run_em();
  }
  void dump_result() {
    map<string, int> tid2length;
    for (size_t i = 0; i < selected_keys_->size(); i++) {
      SelectedKey& sk = selected_keys_->at(i);
      tid2length[sk.tids(i)] = sk.lengths(i);
    }
    for (auto iter : profile_) {
      double rpkm = 0;
      if (tid2length.find(iter.first) != tid2length.end()) {
        rpkm = iter.second / tid2length[iter.first];
      }
      std::cout << iter.first << '\t' << iter.second
                << '\t' << rpkm << std::endl;
    }
  }
private:
  vector<SelectedKey>* selected_keys_;
  const RollingHashCounter* counter_;
  vector<vector<double> > pi_;
  map<string, double> profile_;
  std::atomic<bool>* is_running_;
};

class CountMain {
public:
  CountMain(const string& index_file, const string& read_files1,
            const string& read_files2, int num_threads)
    : num_threads_(num_threads), index_file_(index_file),
      read_files1_(read_files1), read_files2_(read_files2),
      stream_(index_file, ios::in | ios::binary) {}

  void run() {
    int buffer_size = 200000000;
    ::google::protobuf::uint8 * buffer =
        new ::google::protobuf::uint8[buffer_size];
    SelectedKey sk;
    vector<SelectedKey> selected_keys;
    vector<SelectedKey> selected_keys_for_em;
    vector<string> keys;
    LOG(INFO) << "Loading selected keys .. ";
    while(load_protobuf_data(&stream_, &sk, buffer, buffer_size)) {
      for (int i = 0; i < sk.keys_size(); i++) {
        string key = sk.keys(i).key();
        keys.push_back(key);
        // complimentary
        compliment(&key);
        keys.push_back(key);

        // complimentary and reversed
        reverse(key.begin(), key.end());
        keys.push_back(key);

        // reversed
        compliment(&key);
        keys.push_back(key);
      }
      selected_keys.push_back(sk);
      selected_keys_for_em.push_back(sk);
    }
    LOG(INFO) << "Building the index ...";
    LOG(INFO) << "There are totally " << keys.size() << " keys";
    RollingHashCounter counter(keys, 10);
    LOG(INFO) << "Counting the occurrences of the keys in the reads .. ";
    vector<string> fa_files1 = split_seq(read_files1_, ',');
    vector<string> fa_files2 = split_seq(read_files2_, ',');
    RSPairReader* reader_ = nullptr;
    if (FLAGS_fastq)
        reader_ = new RSFastqPairReader(fa_files1, fa_files2);
    else
        reader_ = new RSPairReader(fa_files1, fa_files2);
    std::vector<std::thread> threads(num_threads_);
    CountThread count_thread(reader_, &counter);
    for (int i = 0; i < num_threads_; i ++) {
      threads[i] = std::thread{RSThread(&count_thread)};
    }
    std::atomic<bool> is_running (true);
    EMThread em_thread(&counter, &selected_keys_for_em, &is_running);
    std::thread em;
    if (FLAGS_run_em) {
      em = std::thread{RSThread(&em_thread)};
    }
    barrier(threads);
    if (FLAGS_run_em) {
      is_running = false;
      LOG(INFO) << "Notify other thread the counting is done.";
    }
    LOG(INFO) << "Dumping the results ...";
    fstream ostream(FLAGS_count_file, ios::out | ios::binary | ios::trunc);
    for (auto& sk : selected_keys) {
      update_SelectedKey(counter, &sk);
      write_protobuf_data(&ostream, &sk);
    }
    if (FLAGS_run_em) {
      em.join();
      em_thread.dump_result();
    }
    counter.dump_info();
    delete[] buffer;
  }
private:
  int num_threads_;
  string index_file_;
  string read_files1_;
  string read_files2_;
  fstream stream_;
};

}  // namespace rs


int main(int argc, char *argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_num_threads == -1) {
    FLAGS_num_threads = std::thread::hardware_concurrency();
  }
  LOG(INFO) << "Using " << FLAGS_num_threads << " threads...";
  if (FLAGS_run_em) {
    LOG(INFO) << "Run EM algorithm at another thread";
  }
  rs::CountMain cm(FLAGS_selected_keys_file, FLAGS_read_files1,
                   FLAGS_read_files2, FLAGS_num_threads);
  cm.run();
}

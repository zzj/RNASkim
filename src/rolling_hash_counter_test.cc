#include <cstdio>
#include <thread>

#include "gtest/gtest.h"
#include "glog/logging.h"
#include "rs_thread.h"

#include "rolling_hash_counter.h"

namespace rs {
namespace {

void setup_rolling_hash_array(RollingHashArray *rha, uint32_t size) {
  vector<string> keys;
  // need to create keys first because every key needs to be out live
  // when query it.
  char temp[10];
  for (uint32_t i = 0; i < size; i++) {
    std::sprintf(temp, "%ud" , i);
    rha->insert(temp, i, 0);
  }
  ASSERT_EQ(size, rha->size());
}

TEST(RollingHashArray, test) {
  int test_size = 1000000;
  RollingHashArray rha(test_size);
  setup_rolling_hash_array(&rha, test_size);
  char temp[10];
  for (uint32_t i = 0; i < rha.size(); i++) {
    std::sprintf(temp, "%ud", i);
    ASSERT_TRUE(rha.increase(temp, i, 1));
  }
}

class RollingHashArrayTestThread : public ThreadInterface {
public:
  RollingHashArrayTestThread(RollingHashArray *hash_array, int rep)
    : hash_array_(hash_array), rep_(rep) {}

  void run() {
    char temp[10];
    for (uint32_t i = 0; i < hash_array_->size(); i++) {
      std::sprintf(temp, "%ud", i);
      for (int rerun = 0; rerun < rep_; rerun ++) {
        ASSERT_TRUE(hash_array_->increase(temp, i, 1));
      }
    }
  }
private:
  RollingHashArray* hash_array_;
  int rep_;
};


TEST(RollingHashArray, test_multi_thread) {
  int test_size = 100000;
  int num_threads = 10;
  int num_reps = 10;

  RollingHashArray rha(test_size + 100);
  setup_rolling_hash_array(&rha, test_size);
  vector<std::thread> threads;
  RollingHashArrayTestThread test_thread(&rha, num_reps);
  for (int i = 0; i < num_threads; i ++) {
    threads.emplace_back(std::thread{RSThread(&test_thread)});
  }
  barrier(threads);
  char temp[10];
  for (uint32_t i = 0; i < rha.size(); i++) {
    std::sprintf(temp, "%ud", i);
    RollingHashArray::iterator item = rha.find(temp, i);
    ASSERT_TRUE(nullptr != item);
    ASSERT_EQ(num_reps * num_threads, item->value);
  }
}

TEST(RollingHashCounter, test) {
  vector<string> keys = {"ATCG", "CGAT", "AAAA", "TTTT"};
  RollingHashCounter counter(keys, 10);
  string s = "ATCGATCGATCGATCGATCGATCG";
  counter.process(s);
  ASSERT_EQ(6, counter.find("ATCG"));
  ASSERT_EQ(5, counter.find("CGAT"));
  ASSERT_EQ(0, counter.find("AAAA"));
  ASSERT_EQ(0, counter.find("TTTT"));
}

TEST(RollingHashCounter, another_test) {
  vector<string> keys = {"TTGAAAGACTAAAAGCATTGATAAATCCAGCCAATGTAAC",
                         "TTCCCCGGGACATGGTGCTCGGGGTCTGGACAGAACGGAG"};
  RollingHashCounter counter(keys, 10);
  string s;
  s = "GCCATGGAGATTGTGACCCTTTAGTTCCCCTAATGTTTGGTTCTCTTACTGTTGAAAGACTAAAAGCATTGATAAATCCAGCCAATGTAACCTTCAAAAT";
  counter.process(s);
  s = "AAGGTGAGAAGGCCCATTTTCAAAACTGGCTACTCCTTCCCCGGGACATGGTGCTCGGGGTCTGGACAGAACGGAGAACGGCTCTGAGCAGTGGCACCCT";
  counter.process(s);
  ASSERT_EQ(1, counter.find("TTGAAAGACTAAAAGCATTGATAAATCCAGCCAATGTAAC"));
  ASSERT_EQ(1, counter.find("TTCCCCGGGACATGGTGCTCGGGGTCTGGACAGAACGGAG"));
}

class CounterThread : public ThreadInterface {
public:
  CounterThread(RollingHashCounter *counter, const string& seq)
    : counter_(counter), seq_(seq) {}

  void run() {
    for (int i = 0; i < 100000; i++) {
      counter_->process(seq_);
    }
  }
private:
  RollingHashCounter* counter_;
  string seq_;
};

TEST(RollingHashCounter, multi_thread_test) {
  vector<string> keys = {"TTGAAAGACTAAAAGCATTGATAAATCCAGCCAATGTAAC",
                         "TTCCCCGGGACATGGTGCTCGGGGTCTGGACAGAACGGAG"};
  RollingHashCounter counter(keys, 10);
  string seq1 = "GCCATGGAGATTGTGACCCTTTAGTTCCCCTAATGTTTGGTTCTCTTACTGTTGAAAGACTAAAAGCATTGATAAATCCAGCCAATGTAACCTTCAAAAT";
  string seq2 = "AAGGTGAGAAGGCCCATTTTCAAAACTGGCTACTCCTTCCCCGGGACATGGTGCTCGGGGTCTGGACAGAACGGAGAACGGCTCTGAGCAGTGGCACCCT";
  CounterThread ct1(&counter, seq1), ct2(&counter, seq2),
    ct3(&counter, seq1), ct4(&counter, seq2);

  std::thread t1{RSThread(&ct1)};
  std::thread t2{RSThread(&ct2)};
  std::thread t3{RSThread(&ct3)};
  std::thread t4{RSThread(&ct4)};
  t1.join();
  t2.join();
  t3.join();
  t4.join();

  ASSERT_EQ(200000, counter.find("TTGAAAGACTAAAAGCATTGATAAATCCAGCCAATGTAAC"));
  ASSERT_EQ(200000, counter.find("TTCCCCGGGACATGGTGCTCGGGGTCTGGACAGAACGGAG"));

}

}  // namespace
}  // namespace rs

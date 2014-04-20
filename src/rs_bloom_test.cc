#include <thread>
#include "gtest/gtest.h"

#include "fa_reader.h"
#include "rs_bloom.h"
#include "rs_thread.h"

namespace rs {
namespace {

class reader_thread : public ThreadInterface {
public:
  reader_thread(SingleFastaReader* reader, RSBloom* bloom)
    : reader_(reader), bloom_(bloom) {}

  void run() {
    vector<string> ids, seqs;
    while(reader_->read(&ids, &seqs) != 0) {
      for (string seq : seqs) {
        bloom_->add(seq);
        ASSERT_TRUE(bloom_->contain(seq));
      }
    }
  }
private:
  SingleFastaReader* reader_;
  RSBloom* bloom_;
};

TEST(RSBloom, single_thread_test) {
  RSBloom bloom(1000, 0.0001);
  std::string key = "abc";
  ASSERT_TRUE(bloom.add(key));
  ASSERT_FALSE(bloom.add(key));
  ASSERT_TRUE(bloom.contain(key));
}

TEST(RSBloom, multi_thread_test) {
  std::vector<std::string> test_filenames = {
    "test_data/fa_reader_test.fasta.1",
    "test_data/fa_reader_test.fasta.2",
  };
  SingleFastaReader fa_reader(test_filenames[0]);
  RSBloom bloom(10000, 0.001);

  reader_thread rt1(&fa_reader, &bloom), rt2(&fa_reader, &bloom);
  std::thread t1{RSThread(&rt1)};
  std::thread t2{RSThread(&rt2)};
  t1.join();
  t2.join();
  fa_reader.reset();
  vector<string> ids, seqs;
  while(fa_reader.read(&ids, &seqs) != 0) {
    for (string seq : seqs) {
      bloom.add(seq);
      ASSERT_TRUE(bloom.contain(seq));
    }
  }
  std::string nokey = "ABC";
  ASSERT_FALSE(bloom.contain(nokey));
}

}  // namespace
}  // namespace rs

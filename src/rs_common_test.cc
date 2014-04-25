#include <cstdio>
#include <string>

#include "gtest/gtest.h"

#include "rs_common.h"

using std::string;
using std::vector;
namespace rs {
namespace {

TEST(EvenlySample, simple_test) {
  string seq = "ABAAABAAABA";
  vector<string> samples = evenly_sample(seq, 1, 3);
  ASSERT_EQ(samples.size(), 3);
  ASSERT_EQ(samples[0], "B");
  ASSERT_EQ(samples[1], "B");
  ASSERT_EQ(samples[2], "B");

  samples = evenly_sample(seq, 2, 3);
  ASSERT_EQ(samples.size(), 3);
  for (int i = 0; i < samples.size(); i++) {
    ASSERT_EQ(2, samples[i].size());
  }
  samples = evenly_sample(seq, 7, 3);
  ASSERT_EQ(samples.size(), 3);
  for (int i = 0; i < samples.size(); i++) {
    ASSERT_EQ(7, samples[i].size());
  }
}

TEST(GetFileSize, simple_test) {
  ASSERT_EQ(160173, get_file_size("test_data/fa_reader_test.fasta.1"));
}

}  // namespace
}  // namespace rs

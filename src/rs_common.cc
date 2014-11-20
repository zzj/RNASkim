#include <algorithm>
#include <iostream>
#include <fstream>
#include "glog/logging.h"

#include "rs_common.h"
namespace rs {

void compliment(string* seq) {
  for (size_t i = 0; i < seq->size(); i++) {
    switch (seq->at(i)) {
    case 'A':
      seq->at(i) = 'T'; break;
    case 'T':
      seq->at(i) = 'A'; break;
    case 'G':
      seq->at(i) = 'C'; break;
    case 'C':
      seq->at(i) = 'G'; break;

    default:
      LOG(ERROR) << "Unexpected character: " << seq->at(i) << std::endl
                 << "seq: " << *seq;
    }
  }
}

vector<string> all_keys(const string& key) {
  vector<string> keys;
  string temp = key;
  keys.push_back(temp);
  reverse(temp.begin(), temp.end());
  keys.push_back(temp);
  compliment(&temp);
  keys.push_back(temp);
  reverse(temp.begin(), temp.end());
  keys.push_back(temp);
  return keys;
}

vector<string> split_seq(const string& whole_seq, char delimiter) {
  vector<string> seqs;
  for (size_t i = 0; i < whole_seq.size();) {
    size_t end = i;
    while (end < whole_seq.size() && whole_seq[end] != delimiter &&
           whole_seq[end] != '\n')
      end++;
    seqs.emplace_back(whole_seq.substr(i, end - i));
    i = end + 1;
  }
  return seqs;
}

vector<string> evenly_sample(const string&seq, int sample_len, int sample_size) {
  DLOG_IF(FATAL, sample_size <= 1)
    << "The sample_size cannot be less than 2, sample_size=" << sample_size;
  int step = (seq.size() - sample_len * 2 ) / (sample_size - 1);
  if (step <= 0) step = 0;
  // It is better if the first sample is not start at the leftmost position.
  size_t start = sample_len / 2 + 1;
  vector<string> results;
  for (int i = 0; i < sample_size; i++) {
    // TODO(zzj): temporarily fix the bug, but should work on a fix later.
    if (start > seq.size()) continue;
    results.push_back(seq.substr(start, sample_len));
    start += step;
  }
  return results;
}

string first_seq_in_order(string seq) {
  vector<string> seqs;
  seqs.push_back(seq);

  compliment(&seq);
  seqs.push_back(seq);

  reverse(seq.begin(), seq.end());
  seqs.push_back(seq);

  compliment(&seq);
  seqs.push_back(seq);

  sort(seqs.begin(), seqs.end());
  return seqs[0];
}

uint64_t get_file_size(const string& filename) {
  std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
  LOG_IF(FATAL, !in.is_open()) << "Failed to open " << filename;
  in.seekg(0, std::ifstream::end);
  return in.tellg();
}

}  // namespace rs

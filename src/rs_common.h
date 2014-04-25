#include <vector>
#include <string>
#include <iostream>
#include <cstdint>

using std::vector;
using std::string;
using std::cout;
using std::endl;

namespace rs {

template <class Vector>
void debug_vector(const Vector& vec) {
  for (int i = 0; i < vec.size(); i++) {
    std::cout << vec[i] << ' ';
  }
  std::cout << std::endl;
}

void compliment(string* seq);
vector<string> all_keys(const string& key);
vector<string> split_seq(const string& whole_seq, char delimiter);
vector<string> evenly_sample(const string&seq, int sample_len, int sample_size);
string first_seq_in_order(string seq);
uint64_t get_file_size(const string& fliename);

}

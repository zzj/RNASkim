// This is used for implementing the rolling hash pattern when
// searching the string.

#ifndef RS_HASHARRAY_H
#define RS_HASHARRAY_H

#include <atomic>
#include <string>
#include <vector>

#include "stringpiece.h"
#include "karp_robin_hash.h"

using std::vector;
using std::string;

namespace rs {

struct RollingHashItem {
  string key;
  std::atomic<int> value;
};

// The construction is not thread safe at all, only one thread (main) is
// allowed to construct the object and call the insert operation.
// if you only have find operations after, you could pass it to other threads.
// The class is a open addressed hash array. The user needs to provide
// the hash value for the function.
class RollingHashArray {
public:
  typedef RollingHashItem* iterator;
  RollingHashArray(uint32_t capacity);

  // This function is not thread safe, and should be called only in the
  // main thread
  bool insert(const StringPiece& key, uint32_t hashvalue, const int value);

  // This is thread safe
  // increase the counter by one of the key by one
  bool increase(const StringPiece& key, uint32_t hashvalue, int delta = 1);

  iterator find(const StringPiece& key, uint32_t hashvalue) const;
  iterator end();
  uint32_t size();
  uint32_t capacity();

  std::atomic<long long> hits;
  std::atomic<long long> empty_hits;
  std::atomic<long long> misses;

private:
  uint32_t index(uint32_t hashvalue);
  bool find_next(uint32_t start, const StringPiece& key, uint32_t* next) const;

  RollingHashItem* arena_;
  uint32_t capacity_;
  uint32_t size_;
};

// The construction is not thread safe at all, only one thread (main) is
// allowed to construct the object.
// This is thread safe after the construction
class RollingHashCounter {
public:
  RollingHashCounter(const vector<string>& keys, double factor);
  void process(const string& seq);
  uint32_t find(const string& key) const;
  void dump_info();
private:
  KarpRobinHash *hash_func_;
  RollingHashArray *hash_array_;
  uint32_t key_length_;
  uint32_t capacity_;
};  // namespace rs

}

#endif  // RS_HASHARRAY_H

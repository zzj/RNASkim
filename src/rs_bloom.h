// This is a wrapper of libbloomd in C++.
// libbloomd is a lock free implementation of bloom filter.

#ifndef RS_BLOOM_H
#define RS_BLOOM_H

#include <inttypes.h>
#include <string>

#include "libbloomd/bitmap.h"
#include "libbloomd/bloom.h"

namespace rs {

class StringPiece;

class RSBloom {
 public:
  RSBloom(uint64_t capacity, double fp_probability);
  ~RSBloom();
  bool add(const StringPiece &key);
  bool contain(const StringPiece &key);
 private:
  bloom_bitmap map_;
  bloom_bloomfilter filter_;
};

}  // namespace rs

#endif // RS_BLOOM_H

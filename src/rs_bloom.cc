#include "rs_bloom.h"
#include "stringpiece.h"

namespace rs {

RSBloom::RSBloom(uint64_t capacity, double fp_probability) {
  bloom_filter_params params = {0, 0, capacity, fp_probability};
  bf_params_for_capacity(&params);
  bitmap_from_file(-1, params.bytes, ANONYMOUS, &map_);
  bf_from_bitmap(&map_, params.k_num, 1, &filter_);
}

RSBloom::~RSBloom() {
  bf_close(&filter_);
  bitmap_close(&map_);
}

bool RSBloom::add(const StringPiece &key) {
  int ret = bf_add(&filter_, key.data(), key.size());
  return ret != 0;
}

bool RSBloom::contain(const StringPiece &key) {
  return bf_contains(&filter_, key.data(), key.size());
}

}  // namespace rs

#include "glog/logging.h"

#include "karp_robin_hash.h"

namespace rs {

KarpRobinHash::KarpRobinHash(int size)
  : B_(33), BtoN_(1), size_(size) {
  for (uint32_t i = 0; i < size_ ; ++i) {
    BtoN_ *= B_;
  }
}

KarpRobinHash::HashValueType KarpRobinHash::hash(const StringPiece& key) {
  reset();
  const char* p_limit = key.data() + key.size();
  const char* p = key.data();
  for (; p < p_limit; p++) {
    eat(*p);
  }
  return hashvalue_;
}

void KarpRobinHash::eat(char inchar) {
  hashvalue_ = B_ * hashvalue_ + char_hash(inchar);
}
KarpRobinHash::HashValueType KarpRobinHash::update(char inchar, char outchar) {
  hashvalue_ = B_ * hashvalue_ + char_hash(inchar) - BtoN_ * char_hash(outchar);
  return hashvalue_;
}

KarpRobinHash::HashValueType KarpRobinHash::hash_value() const {
  return hashvalue_;
}

void KarpRobinHash::reset() {
  hashvalue_ = 0;
}

KarpRobinHash::HashValueType KarpRobinHash::char_hash(char c) {
  // just use the character's value is very good. The ratio of
  // false_positive/total is about the given factor for the counter.
  // if we assign the character differnet value, the program runs longer.
  return c;
}

}  // namespace rs

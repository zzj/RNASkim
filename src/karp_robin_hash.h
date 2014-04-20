#include "stringpiece.h"

namespace rs {

class KarpRobinHash {
 public:
  typedef uint32_t HashValueType;
  explicit KarpRobinHash(int size);
  HashValueType hash(const StringPiece& key);
  void eat(char inchar);
  HashValueType update(char inchar, char outchar);
  HashValueType hash_value() const;
  void reset();
 private:
  // return the corresponding hash value for the character
  HashValueType char_hash(char c);
  const HashValueType B_;
  HashValueType BtoN_;
  uint32_t size_;
  HashValueType hashvalue_;
};

}  // namespace rs

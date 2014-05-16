#include <iostream>

#include "rolling_hash_counter.h"
#include "glog/logging.h"
#include "likely.h"

namespace rs {

RollingHashArray::RollingHashArray(uint32_t capacity)
  : capacity_(capacity), size_(0), hits(0), misses(0), empty_hits(0) {
  arena_ = new RollingHashItem[capacity];
  for (uint32_t i = 0; i < capacity_; i++) {
    arena_[i].value.store(0);
  }
}

RollingHashArray::iterator RollingHashArray::end() {
  return nullptr;
}

uint32_t RollingHashArray::size() {
  return size_;
}

uint32_t RollingHashArray::capacity() {
  return capacity_;
}

uint32_t RollingHashArray::index(uint32_t hashvalue) {
  return hashvalue % capacity_;
}

// return false if the key is found.
// return true if the key is empty.
// terminated if the space is not available.
bool RollingHashArray::find_next(uint32_t start, const StringPiece& key,
                                 uint32_t* next) const {
  // linear prob to find next avaiable one.
  uint32_t last = start;
  while (true) {
    if (LIKELY(arena_[start].key.size() == 0)) {
      empty_hits.fetch_add(1, std::memory_order_relaxed);
      *next = start;
      return true;
    }
    // The running time does not improve when we only compare the prefix
    //if (compare_prefix(arena_[start].key, key, 10)) {
    if (UNLIKELY(arena_[start].key == key)) {
      hits.fetch_add(1, std::memory_order_relaxed);
      *next = start;
      return false;
    }
    misses.fetch_add(1, std::memory_order_relaxed);
    start = (start + 1) % capacity_;
    if (UNLIKELY(last == start)) {
      LOG(FATAL) << "The RollingHashArray is full" ;
    }
  }
}

bool RollingHashArray::insert(const StringPiece& key, uint32_t hashvalue,
                              const int value) {
  uint32_t available_index;
  bool should_insert = find_next(index(hashvalue), key, &available_index);
  if (should_insert) {
    arena_[available_index].key = key.ToString();
    arena_[available_index].value.store(value);
    size_ ++;
  }
  return should_insert;
}

bool RollingHashArray::increase(const StringPiece& key, uint32_t hashvalue,
                                int delta) {
  uint32_t key_index;
  bool is_empty = find_next(index(hashvalue), key, &key_index);
  if (LIKELY(is_empty)) {
    return false;
  }
  arena_[key_index].value.fetch_add(delta, std::memory_order_relaxed);
  return true;
}

RollingHashArray::iterator RollingHashArray::find(const StringPiece& key,
                                                  uint32_t hashvalue) const {
  uint32_t key_index;
  bool is_not_found = find_next(index(hashvalue), key, &key_index);
  if (is_not_found) {
    return NULL;
  }
  return &arena_[key_index];
}

RollingHashCounter::RollingHashCounter(const vector<string>& keys, double factor)
  : capacity_(keys.size() * factor) {
  // Make sure there is no thread level variable here
  LOG_IF(FATAL, keys.size() == 0) << "The keys size is 0.";
  key_length_ = keys[0].size();
  hash_func_ = new KarpRobinHash(key_length_);
  hash_array_ = new RollingHashArray(capacity_);
  for (auto& key : keys) {
    auto hashvalue = hash_func_->hash(key);
    hash_array_->insert(key, hashvalue, 0);
  }
}

// This function should be thread safe.
void RollingHashCounter::process(const string& seq) {
  LOG_IF(ERROR, seq.size() < key_length_)
    << "The seq size is smaller than the key_length: seq:" << seq
    << " key_length_:" << key_length_;
  // we only want to look into the kmer without 'N'.
  // if we found a 'N' in seq, we will skip all kmers that cover the
  // 'N'
  const char* p_start = seq.data();
  const char* p_limit = seq.data() + seq.length();
  // We need to copy the hash func here, since this function may be
  // called by multiple threads
  KarpRobinHash hash_func = *hash_func_;

  while (p_start + key_length_ <= p_limit) {
    const char* p_end = p_start + key_length_ - 1;
    bool restart = false;
    for (const char* j = p_end; j >= p_start; j--) {
      if (UNLIKELY(*j == 'N')) {
        // skipped
        p_start = j + 1;
        restart = true;
        break;
      }
    }
    if (UNLIKELY(restart)) continue;
    StringPiece key(p_start, key_length_);
    auto hashvalue = hash_func.hash(key);
    hash_array_->increase(key, hashvalue, 1);
    p_end ++;
    for (; p_end < p_limit; p_start ++, p_end ++) {
      if (UNLIKELY(*p_end == 'N')) {
        p_start = p_end;
        break;
      }
      auto hashvalue = hash_func.update(*p_end,  // inchar
                                        *p_start); // outchar
      StringPiece key(p_start + 1, key_length_);
      hash_array_->increase(key, hashvalue, 1);
    }
    p_start ++;
  }
}

uint32_t RollingHashCounter::find(const string& key) const {
  auto hashvalue = hash_func_->hash(key);
  auto iter = hash_array_->find(key, hashvalue);
  if (UNLIKELY(iter == nullptr)) {
    LOG(ERROR) << "Misuse of the RollingHashCounter. "
               << "The key does not exist in the hash counter";
    return 0;
  }
  return iter->value.load(std::memory_order_relaxed);
}

void RollingHashCounter::dump_info() {
  LOG(INFO) << "Hits: " << hash_array_->hits;
  LOG(INFO) << "Misses: " << hash_array_->misses;
  LOG(INFO) << "Empty hits (last hit is empty item): " << hash_array_->empty_hits;
}

}  // namespace rs

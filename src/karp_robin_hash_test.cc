#include <cstdio>
#include <string>

#include "gtest/gtest.h"

#include "karp_robin_hash.h"

using std::string;
namespace rs {
namespace {

void karp_robin_hash_test_helper(string s1, char inchar) {
  KarpRobinHash hasher(s1.size());

  hasher.reset();
  hasher.hash(s1);
  hasher.update(inchar, s1[0]);
  string s2 = s1.substr(1) + inchar;
  auto hash2 = hasher.hash_value();
  hasher.reset();
  auto hash3 = hasher.hash(s2);
  ASSERT_EQ(hash2, hash3);
  
}

TEST(KarpRobinHash, test) {
  string s1 = "ATCGATCGATCGATCGATCGATCG";
  // add 'C' at the end, remove 'A' at the head
  string s2 = "TCGATCGATCGATCGATCGATCGC";
  karp_robin_hash_test_helper(s1, 'A');
  karp_robin_hash_test_helper(s1, 'T');
  karp_robin_hash_test_helper(s1, 'C');
  karp_robin_hash_test_helper(s1, 'G');
  s1[0] = 'T';
  karp_robin_hash_test_helper(s1, 'A');
  karp_robin_hash_test_helper(s1, 'T');
  karp_robin_hash_test_helper(s1, 'C');
  karp_robin_hash_test_helper(s1, 'G');
}

}  // namespace 
}  // namespace rs

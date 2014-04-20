#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "rs_common.h"
#include "proto/rnasigs.pb.h"
#include "proto_data.h"
#include "rs_estimate_lib.h"

using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::map;
using std::unordered_map;
using std::set;
using std::string;
using std::vector;

#define DEBUG 0

DECLARE_int32(rs_length);

namespace rs{
  double average(const vector<double>& values) {
    if (values.size() == 0) return 0;
    double sum = 0;
    int p = 0;
    for (size_t i = 0; i < values.size(); i++) {
      if (values[i] == 0) return 0;
      sum += values[i];
      p ++;
    }
    if (p == 0) return 0;
    return sum / p;
  }

  int left_possible_locations(int transcript_length,
                              int fragment_length,
                              int read_length,
                              int position) {
    if (position > transcript_length - fragment_length + read_length - FLAGS_rs_length)
      return 0;
    // the maximal number of  possible ways to cover the kmer
    int ret = read_length - FLAGS_rs_length + 1;
    ret = std::min<int>(position + 1, ret);
    ret = std::min<int>(transcript_length - fragment_length + read_length - position - FLAGS_rs_length + 1, ret);
    return ret;
  }

  double weight_of_kmer(int length, vector<int> positions) {
    double weight = 0;
    // why 150? this is a rough estimation, and I am not sure what
    // is the best value. Only a half part of the read will cover
    // the beginning or the end of the transcript.
    int fragment_length = 180;
    int read_length = 95;
    // total number of ways to cover a kmer (left side and right side
    // fragment)
    int total_num_ways = (read_length - FLAGS_rs_length + 1) * 2;
    for (int pos : positions) {
      // the number of ways that the kmer is covered by the left read
      // in the fragment
      int w = left_possible_locations(length, fragment_length,
                                      read_length, pos);
      // the number of ways that the kmer is covered by the right read
      // in the fragment. (just need to reverse the position)
      w += left_possible_locations(length, fragment_length,
                                   read_length, length - pos - 1);
      if (w < 0) {
        LOG(ERROR) << length << ' ' << pos;
        LOG(ERROR) << w;
      }
      if ( w < 1) {
        w = 1;
      }
      weight += w;
    }
    weight /= positions.size();
    return positions.size();
    return 1;
    return weight / total_num_ways;
  }

  vector<vector<int> > find_covered_transcript(const SignatureInfoDB& db,
                                               int num_tids) {
    vector<vector<int> > ret;
    ret.resize(num_tids);
    for (int i= 0; i < num_tids; i++) {
      for (int j = 0; j < num_tids; j ++) {
        if (i == j) continue;
        bool is_covered = true;
        for (auto& item : db) {
          vector<double> weights = item.first;
          for (size_t k = 0; k < weights.size(); k++) {
            if (weights[i] == 0 && weights[j] != 0) {
              is_covered = false;
            }
          }
          if (!is_covered) break;
        }
        if (is_covered) {
          // the ith transcript cover the jth transcript
          ret[i].push_back(j);
        }
      }
    }
    return ret;
  }

  vector<double> normalize(const vector<double> pi) {
    double sum = 0;
    for (size_t i = 0; i < pi.size(); i++) {
      sum += pi[i];
    }
    if (sum == 0) {
      return vector<double>(pi.size(), 0);
    }
    vector<double> new_pi = pi;
    for (size_t i = 0; i < pi.size(); i++) {
      new_pi[i] /= sum;
    }
    return new_pi;
  }

  double target_value(const SignatureInfoDB& db,
                      const vector<double>& pi,
                      const vector<double>& num_kmers) {
    double target = 0;
    for (const auto& item : db) {
      const vector<double>& weights = item.first;
      int count = item.second.total_counts;
      double sum = 0;
      for (size_t i = 0; i < weights.size(); i ++) {
        if (num_kmers[i] != 0) {
          sum += weights[i] * pi[i] * item.second.occurences / num_kmers[i];
        }
      }
      if (sum > 0) {
        target += count * log(sum);
      }
    }
    return target;
  }

  void EM(int num_tids, const SignatureInfoDB& db,
          vector<double>* count_per_tid, vector<double>* pi) {
    vector<double> new_pi;
    int niter = 0;
    count_per_tid->resize(num_tids, 0);
    vector<double> num_kmers(num_tids, 0);
    // if has_weight[key_id][0] is tid, it means the tid contains
    // the key_id. The spars
    vector<vector<int> > key2tids(db.size());

    int item_id = 0;
    for (const auto& item : db) {
      const vector<double>& weights = item.first;
      for (int i = 0; i < num_tids; i++) {
        if (weights[i] > 0) {
          key2tids[item_id].push_back(i);
          num_kmers[i] += weights[i] * item.second.occurences;
        }
      }
      item_id ++;
    }
    double target = target_value(db, *pi, num_kmers);

    while (true) {
      // resize only set default value for the new variables
      // These two lines must be present at the same time
      // need to write a test to avoid further errors. TODO(zzj)
      count_per_tid->resize(0);
      count_per_tid->resize(num_tids, 0);

      int item_id = 0;
      for (const auto& item : db) {
        const vector<double>& weights = item.first;
        vector<double> new_pi_tid(key2tids[item_id].size(), 0);

        for (size_t i = 0; i < new_pi_tid.size(); i++) {
          int tid = key2tids[item_id][i];
          if (num_kmers[tid] != 0) {
            new_pi_tid[i] = weights[tid] * pi->at(tid) * item.second.occurences / num_kmers[tid];
          }
        }
        new_pi_tid = normalize(new_pi_tid);
        int count = item.second.total_counts;
        for (size_t i = 0; i < new_pi_tid.size(); i++) {
          int tid = key2tids[item_id][i];
          count_per_tid->at(tid) += count * new_pi_tid[i];
        }
        if (DEBUG && niter == 1000) {
          cout << "count: " << count << '\n';
          for (size_t i = 0; i < new_pi_tid.size(); i++) {
            int tid = key2tids[item_id][i];
            cout << tid << ':' << weights[tid] << ' ' ;
          }
          cout << '\n';
          for (size_t i = 0; i < new_pi_tid.size(); i++) {
            int tid = key2tids[item_id][i];
            cout << tid << ':' << count * new_pi_tid[tid] << ' ' ;
          }
          cout << '\n';
        }
        item_id ++;
      }
      for (int i = 0; i < num_tids; i++) {
        if (num_kmers[i] == 0) {
          count_per_tid->at(i) = 0;
          // LOG(ERROR) << "error";
          continue;
        }
        count_per_tid->at(i) = count_per_tid->at(i);
        if (count_per_tid->at(i) < 0.0000001) {
          count_per_tid->at(i) = 0;
        }
      }

      new_pi = normalize(*count_per_tid);

      double new_target = 0;

      // new_target = target_value(db, new_pi, num_kmers);
      // *pi = new_pi;
      // if (new_target - target < 0) {
      //   LOG(ERROR) << "EM algorithm does not maximize the target function, iter=" << niter;
      //   LOG(ERROR) << new_target - target;
      //   LOG(ERROR) << num_tids;
      // }
      // // target is negative
      // // check relative change is smaller enough or not
      // if ((new_target - target) / (-target) < 1e20) {
      //   break;
      // }
      double diff = 0;
      for (size_t i = 0; i < pi->size(); i++) {
        diff = std::max<double>(diff, fabs(new_pi[i] - pi->at(i)));
      }
      *pi = new_pi;
      if (diff < 0.0000000001) {
        break;
      }
      target = new_target;
      niter ++;
      if (niter > 500000) break;
      // LOG_IF(ERROR, num_tids > 100) << niter;
    }
    for (int tid = 0; tid < num_tids; tid++) {
      if (num_kmers[tid] != 0)
        count_per_tid->at(tid) /= num_kmers[tid];
    }
    *pi = normalize(*count_per_tid);
    if (DEBUG) {
      for (const auto& kv : db) {
        int count = kv.second.total_counts;
        // const vector<double>& weights = kv.first;
        // debug_vector(weights);
        cout << count << endl;
      }
      debug_vector(new_pi);
    }
    return ;
  }

  bool prepare_SignatureInfoDB(const SelectedKey &sk, SignatureInfoDB *db) {
    bool run_em = false;
    for (int i = 0; i < sk.keys_size(); i++) {
      auto& key = sk.keys(i);
      set<string> tids;
      vector<double> weights(sk.tids_size(), 0);

      for (int j = 0; j < key.transcript_infos_size(); j ++) {
        auto & info = key.transcript_infos(j);
        tids.insert(sk.tids(info.tidx()));
        vector<int> pos;
        for (auto p : info.positions()) {
          pos.push_back(p);
        }
        weights[info.tidx()] +=
          weight_of_kmer(sk.lengths(info.tidx()), pos);

        // cout << sk.tids(info.tidx()) << endl;
        // cout << key.key() << endl;
        // cout << key.count() << endl;
      }
      auto iter = db->find(weights);
      if (iter != db->end()) {
        iter->second.total_counts += key.count();
        iter->second.occurences += 1;
      } else {
        SignatureInfo si;
        si.total_counts = key.count();
        si.occurences = 1;
        (*db)[weights] = si;
      }
      if (key.count() != 0) {
        run_em = true;
      }
      // if (DEBUG && weights[inspected_tid] != 0) {
      //   cout << weights[inspected_tid] << ' ' << key.count() << '\n';
      //   cout << key.key() << '\n';
      // }
    }
    return run_em;
  }
} // namespace rs

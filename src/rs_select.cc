// TODO(zzj): make this program multithread

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "proto/rnasigs.pb.h"
#include "proto_data.h"
#include "rs_common.h"

using std::fstream;
using std::ios;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

DEFINE_string(index_file, "",
              "The path to the index file (input).");
DEFINE_string(selected_keys_file, "",
              "The path to the selected keys file (output).");
DEFINE_int32(rs_length, 40,
             "The length of sig-mer.");
DEFINE_int32(num_kmer_per_region, 10,
             "The number of sig-mers selected from each sig-mer region.");
using namespace rs;

void find_all_positions(const string& seq, const string& key,
                        vector<int> *positions, int shift = 0) {
  int position = 0;
  while (true) {
    position = seq.find(key, position);
    if (position == string::npos) {
      break;
    }
    positions->push_back(shift + position);
    position ++;
    // LOG(ERROR) << i << ' ' << j << ' ' << key <<  ' ' <<  position;
  }
}

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);
  int buffer_size = 200000000;
  ::google::protobuf::uint8 * buffer =
      new ::google::protobuf::uint8[buffer_size];
  GeneSignatures gene;
  fstream istream(FLAGS_index_file, ios::in | ios::binary);
  fstream ostream(FLAGS_selected_keys_file, ios::out | ios::binary | ios::trunc);
  long long total_selected_keys = 0;
  long long total_keys = 0;
  while(load_protobuf_data(&istream, &gene, buffer, buffer_size)) {
    set<string> gene_keys;

    // select keys from every transcripts
    for (int i = 0; i < gene.transcripts_size(); i ++) {
      const auto& transcript = gene.transcripts(i);
      int step = std::min(50 + (10 - FLAGS_num_kmer_per_region) * 10,
                          (transcript.length() - 20) / FLAGS_num_kmer_per_region);
      if (step < 10) {
        step = 10;
      }
      set<string> signatures;
      for (int j = 0; j < transcript.signatures_size(); j++) {
        const auto& sigs = transcript.signatures(j);
        const string& seq = sigs.seq();
        int length = seq.length();
        total_keys += length - FLAGS_rs_length + 1;
        int start = std::min(20, length - FLAGS_rs_length);
        for (int p = start; p < length - FLAGS_rs_length; p += step) {
          string sig = seq.substr(p, FLAGS_rs_length);
          signatures.insert(first_seq_in_order(sig));
        }
      }
      for (const string& sig : signatures) {
        gene_keys.insert(sig);
      }
    }
    total_selected_keys += gene_keys.size();
    // LOG(ERROR) << gene_keys.size();
    // check whether the selected key occurs more than once in the
    // current gene's transcripts' sequences.
    SelectedKey sk;
    sk.set_gid(gene.id());
    for (int i = 0; i < gene.transcripts_size(); i ++) {
      const auto& transcript = gene.transcripts(i);
      sk.add_tids(transcript.id());
      sk.add_lengths(transcript.length());
    }

    map<string, vector<pair<int, int> > > key_index;
    for (int i = 0; i < gene.transcripts_size(); i ++) {
      const auto& transcript = gene.transcripts(i);
      vector<int> positions;
      for (int j = 0; j < transcript.signatures_size(); j++) {
        const auto& sigs = transcript.signatures(j);
        string seq = sigs.seq();
        for (int k = 0; k < seq.size() - FLAGS_rs_length; k++) {
          key_index[seq.substr(k, FLAGS_rs_length)]
            .push_back(std::make_pair(i, sigs.position() + k));
        }
      }
    }

    set<int> skipped_tids;
    // tid -> sig positions
    map<int, set<string> > tid2sigs;

    for (auto key : gene_keys) {
      vector<string> all_keys_result = all_keys(key);
      for (auto iter_key: all_keys_result) {
        if (key_index.find(iter_key) != key_index.end()) {
          for (auto& item : key_index[iter_key]) {
            tid2sigs[item.first].insert(key);
          }
        }
      }
    }

    for (auto& item : tid2sigs) {
      if (item.second.size() <= 2) {
        LOG(ERROR) << gene.transcripts(item.first).id() << " has "
                   << item.second.size() << " rna_signatures. Skipped.";
        skipped_tids.insert(item.first);
      }
    }

    for (auto key : gene_keys) {
      SelectedKey::Key* key_info = sk.add_keys();
      key_info->set_key(key);
      vector<string> all_keys_result = all_keys(key);
      map<int, vector<int> > tid2positions;
      for (auto key: all_keys_result) {
        if (key_index.find(key) != key_index.end()) {
          for (auto& item : key_index[key]) {
            if (skipped_tids.find(item.first) != skipped_tids.end()) continue;
            tid2positions[item.first].push_back(item.second);
          }
        }
      }

      for (auto &item : tid2positions) {
        SelectedKey::Key::TranscriptInfo* ti = key_info->add_transcript_infos();
        ti->set_tidx(item.first);
        for (int position : item.second) {
          LOG_IF(ERROR, position >= sk.lengths(item.first))
            << "The position is out of the boundary";
          ti->add_positions(position);
        }
        // LOG(ERROR) << item.first;
        // debug_vector(item.second);
      }

      // Do not delete this. This is for the future reference.
      // for (int i = 0; i < gene.transcripts_size(); i ++) {
      //   const auto& transcript = gene.transcripts(i);
      //   vector<int> positions;
      //   for (int j = 0; j < transcript.signatures_size(); j++) {
      //     const auto& sigs = transcript.signatures(j);
      //     const string& seq = sigs.seq();
      //     // performance alert!!!
      //     // should convert the transcript first.
      //     // TODO(zzj): fix this later
      //     for (auto key : all_keys_result) {
      //       find_all_positions(seq, key, &positions,
      //                          sigs.position());  // shifted position
      //     }
      //   }
      //   // LOG(ERROR) << i;
      //   // debug_vector(positions);
      //   if (positions.size() > 0) {
      //     SelectedKey::Key::TranscriptInfo* ti = key_info->add_transcript_infos();
      //     ti->set_tidx(i);
      //     for (auto position : positions) {
      //       if (position >= sk.lengths(i)) {
      //         LOG(ERROR) << "The position is out of the boundary";
      //       }
      //       ti->add_positions(position);
      //     }
      //   }
      // }
    }
    if (sk.keys_size() == 0) {
      for (int i = 0; i < sk.tids_size(); i++) {
        LOG(ERROR) << sk.tids(i) << " does not have any rna_signatures";
      }
    }
    write_protobuf_data(&ostream, &sk);
  }
  LOG(ERROR) << total_selected_keys << " sig-mers are selected";
  LOG(ERROR) << total_keys << " sig-mers are scanned";
  delete buffer;

}

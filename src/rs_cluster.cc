#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <iostream>
#include <set>
#include <map>
#include <cstdio>
#include <omp.h>

#include "gflags/gflags.h"
#include "glog/logging.h"

#include "fa_reader.h"
#include "rs_bloom.h"
#include "rs_common.h"

using namespace std;

DEFINE_string(transcript_fasta, "",
              "The fasta file for every transcript sequences.");
DEFINE_string(output, "clustered_gene.fa",
              "The output of the result");
DEFINE_int32(num_threads, -1,
             "The num of threads used in the program. "
             "[default: -1]: the num of CPUs in the machine.");
DEFINE_int32(rs_length, 40,
             "The length of the RS signature.");
DEFINE_double(threshold, 0.1,
              "The threshold of minimum similarity for adding edges to the graph");
DEFINE_string(map_file, "overlap_map",
              "The file for storing the graph inforamtion");

namespace rs {

struct FastaRecord {
public:
  // gene id
  string gid;
  // sequence information for every transcript
  vector<string> seqs;
  // transcript id
  vector<string> tids;
  // bloom filter for all subsequences in this one
  RSBloom* bloom_filter;
  // index in the vector of FastaRecords
  set<int> similar_genes;
};

void add_into_bloom_filter(string& seq, RSBloom* bloom) {
  if ((int) seq.size() < FLAGS_rs_length) return;
  for (int start = 0; start < (int)seq.size() - FLAGS_rs_length;
       start ++) {
    bloom->add(StringPiece(seq.data() + start, FLAGS_rs_length));
  }
}

void load_records(const string& filename, vector<FastaRecord> *records) {
  vector<string> ids, seqs;
  SingleFastaReader* reader = new SingleFastaReader(filename);

  int total = 0;
  while (reader->read(&ids, &seqs) != 0) {
    for (size_t i = 0; i < ids.size(); i++) {
      string whole_seq = seqs[i];
      string header = ids[i];

      // first one is the gene id, the remaining ones are transcript ids.
      vector<string> ids = split_seq(header, '|');
      vector<string> seqs = split_seq(whole_seq, '|');

      FastaRecord fr;
      fr.gid = ids[0];
      ids.erase(ids.begin());
      fr.seqs = seqs;
      fr.tids = ids;
      fr.bloom_filter = new RSBloom(whole_seq.size() * 4, 0.001);
      for (const string& seq : fr.seqs) {
        // add four different possible sequences into bloom
        // normal
        string temp = seq;
        add_into_bloom_filter(temp, fr.bloom_filter);

        // complimentary
        compliment(&temp);
        add_into_bloom_filter(temp, fr.bloom_filter);

        // complimentary and reversed
        reverse(temp.begin(), temp.end());
        add_into_bloom_filter(temp, fr.bloom_filter);

        // reversed and normal
        compliment(&temp);
        add_into_bloom_filter(temp, fr.bloom_filter);
      }
      records->push_back(fr);
    }
    total += ids.size();
    LOG(INFO) << "Loaded " << total << " genes/transcripts into the memory";
  }
  LOG(INFO) << "Data is loaded";
}

void calculate_similarity(vector<FastaRecord> *records) {
  FILE* fd = fopen(FLAGS_map_file.c_str(), "w+");
  #pragma omp parallel for
  for (size_t i = 0; i < records->size(); i++) {
    const FastaRecord& ri = records->at(i);
    vector<string> samples;
    for (const string& seq : ri.seqs) {
      vector<string> temp = evenly_sample(seq, FLAGS_rs_length, 20);
      samples.insert(samples.end(), temp.begin(), temp.end());
    }
    for (size_t j = 0; j < records->size(); j++) {
      if (i == j) continue;
      int total = 0;
      int dup = 0;
      const FastaRecord& rj = records->at(j);
      for (const string& sample : samples) {
        if (rj.bloom_filter->contain(sample)) {
          dup ++;
        }
        total ++;
      }
      if (dup >= 1) {
        // http://www.gnu.org/software/libc/manual/html_node/Streams-and-Threads.html#Streams-and-Threads
        // The POSIX standard requires that by default the stream
        // operations are atomic. I.e., issuing two stream operations
        // for the same stream in two threads at the same time will
        // cause the operations to be executed as if they were issued
        // sequentially.
        fprintf(fd, "%s\t%s\t%d\t%d\n", ri.gid.c_str(), rj.gid.c_str(), dup, total);
      }
    }
  }
  fclose(fd);

  fd = fopen(FLAGS_map_file.c_str(), "r");
  int dup, total;
  char g1[200], g2[200];
  map<string, int> gid2idx;
  for (size_t i = 0; i < records->size(); i++) {
    const FastaRecord& ri = records->at(i);
    gid2idx[ri.gid] = i;
  }
  int num_edges = 0;
  int num_loaded_edges = 0;
  while(fscanf(fd, "%s%s%d%d", g1, g2, &dup, &total) == 4) {
    if (dup * 1.0 / total >= FLAGS_threshold) {
      records->at(gid2idx[g1]).similar_genes.insert(gid2idx[g2]);
      records->at(gid2idx[g2]).similar_genes.insert(gid2idx[g1]);
      num_edges ++;
    }
    num_loaded_edges ++;
  }
  LOG(INFO) << "Map is loaded. "
             << num_loaded_edges << " loaded into the memory."
             << num_edges << " added to the map.";
  fclose(fd);
}

void dfs(const vector<FastaRecord>& records,
         int start,
         set<int>* connected_nodes,
         set<int>* visited_nodes) {

  if (visited_nodes->find(start) != visited_nodes->end()) return;
  connected_nodes->insert(start);
  visited_nodes->insert(start);
  for (int next_gidx : records[start].similar_genes) {
    dfs(records, next_gidx, connected_nodes, visited_nodes);
  }
}

void cluster(const vector<FastaRecord>& records,
             vector<vector<int> > *cluster_results) {
  cluster_results->resize(records.size());
  set<int> visited_nodes;
  for (size_t i = 0; i < records.size(); i ++) {
    if (visited_nodes.find(i) != visited_nodes.end()) {
      continue;
    }
    if (records[i].similar_genes.size() == 0) {
      cluster_results->at(i).push_back(i);
      continue;
    }
    set<int> connected_nodes;
    dfs(records, i, &connected_nodes, &visited_nodes);
    for (int gidx : connected_nodes) {
      cluster_results->at(i).push_back(gidx);
    }
  }
}

void dump_result(const vector<FastaRecord>& records,
                 const vector<vector<int> >& cluster_results) {
  FILE *fd = fopen(FLAGS_output.c_str(), "w+");
  for (size_t i = 0; i < cluster_results.size(); i++) {
    if (cluster_results[i].size() == 0) continue;
    string gid = records[cluster_results[i][0]].gid;
    vector<string> tids, seqs;
    for (int j : cluster_results[i]) {
      tids.insert(tids.end(), records[j].tids.begin(), records[j].tids.end());
      seqs.insert(seqs.end(), records[j].seqs.begin(), records[j].seqs.end());
    }
    fprintf(fd, ">%s", gid.c_str());
    for (size_t j = 0; j < tids.size(); j++) {
      fprintf(fd, "|%s", tids[j].c_str());
    }
    fprintf(fd, "\n");
    for (size_t j = 0; j < seqs.size(); j++) {
      fprintf(fd, "%s", seqs[j].c_str());
      if (j != seqs.size() - 1) {
        fprintf(fd, "|");
      }
    }
    fprintf(fd, "\n");
  }
  fclose(fd);
}

} // namespace rs

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_num_threads == -1) {
    FLAGS_num_threads = std::thread::hardware_concurrency();
  }
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(FLAGS_num_threads);

  vector<rs::FastaRecord> records;
  rs::load_records(FLAGS_transcript_fasta, &records);

  rs::calculate_similarity(&records);

  vector<vector<int> > cluster_results;
  rs::cluster(records, &cluster_results);
  LOG(INFO) << "There are " << cluster_results.size() << " clusters.";
  rs::dump_result(records, cluster_results);
}

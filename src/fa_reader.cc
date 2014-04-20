#include <time.h>
#include <sys/time.h>

#include "fa_reader.h"
#include "glog/logging.h"

namespace rs {

const int BufferSize = 10000;

SingleFastaReader::SingleFastaReader(const string& file,
                                     int buffer_size)
  : file_(file), fd_(nullptr), buffer_size_(buffer_size) {
  fd_.open(file_.c_str(), ios::in);
  LOG_IF(FATAL, !fd_.good()) << "Failed to open file " << file_;
  fd_.rdbuf()->pubsetbuf(buffer, 1024 * 1024 * 5);
}

int SingleFastaReader::read(vector<string>* ids, vector<string>* seqs) {
  ids->clear(), seqs->clear();
  string id;
  string read;
  int total_reads = 0;
  std::lock_guard<std::mutex> lock(m_);
  while(!fd_.eof()) {
    fd_ >> id >> read;
    // only add if the line is not empty
    if (read.size() > 2) {
      // remove the first letter in fasta file. (should be either '>'
      // or '@'.
      ids->push_back(id.substr(1));
      seqs->push_back(read);
      total_reads ++;
    } else {
      break;
    }
    if (total_reads >= buffer_size_) break;
  }
  return total_reads;
}

void SingleFastaReader::reset() {
  fd_.clear(); fd_.seekg(0);
}


RSPairReader::RSPairReader(const std::vector<std::string>& files1,
                           const std::vector<std::string>& files2,
                           int buffer_size)
  : files1_(files1), files2_(files2),
    fd1_(nullptr), fd2_(nullptr), buffer_size_(buffer_size) {

  LOG_IF(FATAL, files1.size() != files2.size())
    << "Different size of paired files";
  LOG_IF(FATAL, files1.size() != 1)
    << "Cannot support more than one fasta file yet";
  LOG_IF(FATAL, files1.size() == 0)
    << "No input files";
  current_file_idx_ = 0;
  fd1_.open(files1_[current_file_idx_].c_str(), ios::in);
  LOG_IF(FATAL, !fd1_.good()) << "Failed to open file "
                             << files1_[current_file_idx_];
  fd2_.open(files2_[current_file_idx_].c_str(), ios::in);
  LOG_IF(FATAL, !fd2_.good()) << "Failed to open file "
                             << files2_[current_file_idx_];
  fd1_.rdbuf()->pubsetbuf(buffer1, 1024 * 1024 * 20);
  fd2_.rdbuf()->pubsetbuf(buffer2, 1024 * 1024 * 20);
  total_time = 0;
}

// Return the number of reads
int RSPairReader::read(vector<string>* reads1, vector<string>* reads2) {
  struct timeval start_time, end_time;
  gettimeofday(&start_time, NULL);
  reads1->resize(buffer_size_); reads2->resize(buffer_size_);
  string id1, id2;
  string read1, read2;
  int total_reads = 0;
  std::lock_guard<std::mutex> lock(m_);
  while(!fd1_.eof() && !fd2_.eof()) {
    // ignore the id line
    fd1_.ignore(256 * 256,'\n');
    fd2_.ignore(256 * 256,'\n');
    fd1_ >> reads1->at(total_reads);
    fd2_ >> reads2->at(total_reads);
    // ignore the remaining new line character
    fd1_.ignore(256 * 256,'\n');
    fd2_.ignore(256 * 256,'\n');
    read_quality_score();
    // only add if the line is not empty
    if (reads1->at(total_reads).size() > 2) {
      total_reads ++;
    } else {
      break;
    }
    if (total_reads >= buffer_size_) break;
  }
  reads1->resize(total_reads), reads2->resize(total_reads);
  gettimeofday(&end_time, NULL);
  total_time += end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
  return total_reads;
}

void RSPairReader::read_quality_score() {
}

RSFastqPairReader::RSFastqPairReader(const std::vector<std::string>& files1,
                                     const std::vector<std::string>& files2,
                                     int buffer_size)
  : RSPairReader(files1, files2, buffer_size) {}

void RSFastqPairReader::read_quality_score() {
  fd1_.ignore(256 * 256,'\n');
  fd2_.ignore(256 * 256,'\n');
  fd1_.ignore(256 * 256,'\n');
  fd2_.ignore(256 * 256,'\n');
}

}  // namespace rs

// This file is a wrapper of read_parser in jellyfish
// This abstraction should not leak any stuff about jellyfish outside
// of this class

#ifndef FA_READER_H
#define FA_READER_H

#include <fstream>
#include <mutex>
#include <string>
#include <vector>

#include "stringpiece.h"

using std::fstream;
using std::ios;
using std::string;
using std::vector;

// This class is a thread safe class.
namespace rs {

  class SingleFastaReader {
  public:
    SingleFastaReader(const std::string& filename,
                      int buffer_size = 50000);
    int read(std::vector<std::string>* ids, std::vector<std::string>* seq);
    void reset();
  private:
    std::string file_;
    fstream fd_;
    char buffer [1024 * 1024 * 5];
    mutable std::mutex m_;
    int buffer_size_;
  };

  // TODO(zzj): support multiple files
  class RSPairReader {
  public:
    RSPairReader(const std::vector<std::string>& files1,
                 const std::vector<std::string>& files2,
                 int buffer_size = 50000);

    // This is thread safe
    int read(vector<string>* reads1, vector<string>* reads2);
  private:
    std::vector<std::string> files1_;
    std::vector<std::string> files2_;
  protected:
    virtual void read_quality_score();
    fstream fd1_;
    fstream fd2_;
  private:
    char buffer1 [1024 * 1024 * 5];
    char buffer2 [1024 * 1024 * 5];
    int current_file_idx_;
    mutable std::mutex m_;
    int buffer_size_;
    double total_time;
  };

  class RSFastqPairReader : public RSPairReader {
  public:
    RSFastqPairReader(const std::vector<std::string>& files1,
                 const std::vector<std::string>& files2,
                 int buffer_size = 50000);
  protected:
    virtual void read_quality_score();
  };
}  // namespace rs

#endif // FA_READER_H

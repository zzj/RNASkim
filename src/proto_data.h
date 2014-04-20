#ifndef _PROTO_DATA_H
#define _PROTO_DATA_H

#include <cstdio>
#include <fstream>
#include <google/protobuf/io/coded_stream.h>
#include "glog/logging.h"

#include <memory>
#include <string>

using namespace google::protobuf::io;
using std::string;
using std::fstream;
using std::fprintf;
using std::unique_ptr;

template<class T>
bool load_protobuf_data(fstream *file, T *data,
                        ::google::protobuf::uint8 *buffer,
                        uint32_t buffer_size) {
  int m;
  if (file->eof()) return false;
  file->read(reinterpret_cast<char *> (&m), sizeof(m));
  if (file->eof()) return false;
  if (buffer_size < m) {
    LOG(ERROR) << "The buffer size (" << buffer_size
               << ") is smaller than the size(" << m
               << "). Skipped the current one";
    file->seekg(m, std::ios_base::cur);
    return false;
  }
  file->read((char *) buffer, m);
  // must use assign because '\0' might be part of the message.
  // TODO: add unit test.
  unique_ptr<CodedInputStream> input(new CodedInputStream(buffer, m));
  input->SetTotalBytesLimit(1024 * 1024 * 1024, 10 * 1024 * 1024);
  if (!data->ParseFromCodedStream(input.get())) {
    fprintf(stderr, "Failed to parse the buffer. [fstream version]");
    return false;
  }
  return true;
}

// TODO(zzj): This function should not be copied. Refactor later.
template<class T>
bool load_protobuf_data(fstream *file, T *data) {
  int m;
  file->read(reinterpret_cast<char *> (&m), sizeof(m));
  if (file->eof()) return false;
  ::google::protobuf::uint8 * buffer = new ::google::protobuf::uint8[m + 1];
  file->read((char *) buffer, m);
  // must use assign because '\0' might be part of the message.
  // TODO: add unit test.
  unique_ptr<CodedInputStream> input(new CodedInputStream(buffer, m));
  input->SetTotalBytesLimit(1024 * 1024 * 1024, 10 * 1024 * 1024);
  if (!data->ParseFromCodedStream(input.get())) {
    fprintf(stderr, "Failed to parse the buffer. [fstream version]");
    delete buffer;
    return false;
  }
  delete buffer;
  return true;
}

template<class T>
bool write_protobuf_data(fstream *file, T *data) {
  string seq;
  data->SerializeToString(&seq);
  int m = seq.size();
  file->write(reinterpret_cast<char *> (&m), sizeof(m));
  file->write(seq.c_str(), seq.size());
  return true;
}


#endif  // _PROTO_DATA_H

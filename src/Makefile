# GTEST Starts
# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = ../lib/gtest-1.7.0/

# Where to find user code.
USER_DIR = .

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -static -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -static -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR)  $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^
# Remeber to add -I$(GTEST_DIR)/include to CPPFLAGS and add gtest_main.a
# when link to the target.
# GTEST ends

LIB =
# This is for release software, but this will cause test programs fail
# to compile, commet this line when you do not want to release the
# software
#STATIC = -static
LDFLAGS = -O3  ../lib/glog-0.3.3/.libs/libglog.a \
	../lib/gflags-2.0/.libs/libgflags.a -lpthread \
	 ../lib/protobuf-2.5.0/src/.libs/libprotobuf.a -fopenmp
LDFLAGS_WITH_STATIC = $(LDFLAGS) $(STATIC)

INC = -I. -I../lib/gflags-2.0/src/ \
	-I../lib/glog-0.3.3/src/ -I../lib/protobuf-2.5.0/src/

RS_COMMON_TEST_SRCS = rs_common.cc rs_common_test.cc
RS_COMMON_TEST_OBJECTS = $(RS_COMMON_TEST_SRCS:.cc=.o)
RS_COMMON_TEST_EXECUTABLE = rs_common_test

FA_READER_SRCS = fa_reader.cc
FA_READER_OBJECTS = $(FA_READER_SRCS:.cc=.o)
FA_READER_TEST_SRCS = fa_reader.cc fa_reader_test.cc
FA_READER_TEST_OBJECTS = $(FA_READER_TEST_SRCS:.cc=.o)
FA_READER_TEST_EXECUTABLE = fa_reader_test

RS_BLOOM_SRCS = fa_reader.cc  libbloomd/murmurhash/MurmurHash3.cc \
	libbloomd/spookyhash/spooky.cc \
	libbloomd/bitmap.cc  libbloomd/bloom.cc rs_bloom.cc
RS_BLOOM_OBJECTS = $(RS_BLOOM_SRCS:.cc=.o)
RS_BLOOM_TEST_SRCS = $(RS_BLOOM_SRCS) rs_bloom_test.cc
RS_BLOOM_TEST_OBJECTS = $(RS_BLOOM_TEST_SRCS:.cc=.o)
RS_BLOOM_TEST_EXECUTABLE = rs_bloom_test

# rolling hash counter
ROLLING_HASH_COUNTER_SRCS = rolling_hash_counter.cc karp_robin_hash.cc \
	stringpiece.cc
ROLLING_HASH_COUNTER_OBJECTS = $(ROLLING_HASH_COUNTER_SRCS:.cc=.o)
ROLLING_HASH_COUNTER_TEST_SRCS = $(ROLLING_HASH_COUNTER_SRCS) \
	rolling_hash_counter_test.cc
ROLLING_HASH_COUNTER_TEST_OBJECTS = $(ROLLING_HASH_COUNTER_TEST_SRCS:.cc=.o)
ROLLING_HASH_COUNTER_TEST_EXECUTABLE = rolling_hash_counter_test

RS_CLUSTER_SRCS = rs_cluster.cc proto/rnasigs.pb.cc rs_common.cc
RS_CLUSTER_OBJECTS = $(RS_CLUSTER_SRCS:.cc=.o)
RS_CLUSTER_EXECUTABLE = rs_cluster

RS_INDEX_SRCS = rs_index.cc proto/rnasigs.pb.cc rs_common.cc
RS_INDEX_OBJECTS = $(RS_INDEX_SRCS:.cc=.o)
RS_INDEX_EXECUTABLE = rs_index

RS_SELECT_SRCS = rs_select.cc proto/rnasigs.pb.cc rs_common.cc
RS_SELECT_OBJECTS = $(RS_SELECT_SRCS:.cc=.o)
RS_SELECT_EXECUTABLE = rs_select

RS_ESTIMATE_SRCS = rs_estimate_lib.cc rs_estimate.cc proto/rnasigs.pb.cc
RS_ESTIMATE_OBJECTS = $(RS_ESTIMATE_SRCS:.cc=.o)
RS_ESTIMATE_EXECUTABLE = rs_estimate

RS_COUNT_SRCS = rs_count.cc proto/rnasigs.pb.cc rs_common.cc \
	rs_estimate_lib.cc $(ROLLING_HASH_COUNTER_SRCS)
RS_COUNT_OBJECTS = $(RS_COUNT_SRCS:.cc=.o)
RS_COUNT_EXECUTABLE = rs_count

EXECUTABLES = $(RS_INDEX_EXECUTABLE) $(RS_CLUSTER_EXECUTABLE) \
	$(RS_SELECT_EXECUTABLE) $(RS_COUNT_EXECUTABLE) \
	$(RS_ESTIMATE_EXECUTABLE)
SOURSES = $(SOURCES)
OBJECTS = $(RS_INDEX_OBJECTS) $(FA_READER_TEST_OBJECTS) \
	$(RS_BLOOM_TEST_OBJECTS) $(ROLLING_HASH_COUNTER_TEST_OBJECTS)
TESTS = gtest.a  gtest_main.a $(FA_READER_TEST_EXECUTABLE) \
	$(RS_BLOOM_TEST_EXECUTABLE) $(ROLLING_HASH_COUNTER_TEST_EXECUTABLE) \
	$(RS_COMMON_TEST_EXECUTABLE) karp_robin_hash_test

all: proto/rnasigs.pb.h rs/rnasigs_pb2.py gtest_main.a $(SOURCES) \
	$(EXECUTABLES) $(TESTS)

# Flags passed to the preprocessor.
CPPFLAGS = $(INC) -I$(GTEST_DIR)/include -O3 -Wall -Wextra  -std=c++0x \
	-fpermissive $(STATIC) -D_GLIBCXX_USE_NANOSLEEP -fopenmp

$(RS_CLUSTER_EXECUTABLE): $(RS_CLUSTER_OBJECTS) $(RS_BLOOM_OBJECTS) $(FA_READER_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS_WITH_STATIC) -o $@

$(RS_INDEX_EXECUTABLE): $(RS_INDEX_OBJECTS) $(RS_BLOOM_OBJECTS) $(FA_READER_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS_WITH_STATIC) -o $@

$(RS_SELECT_EXECUTABLE): $(RS_SELECT_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS_WITH_STATIC) -o $@

$(RS_COUNT_EXECUTABLE): $(RS_COUNT_OBJECTS) $(FA_READER_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS_WITH_STATIC) -o $@

$(RS_ESTIMATE_EXECUTABLE): $(RS_ESTIMATE_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS_WITH_STATIC) -o $@

$(RS_COMMON_TEST_EXECUTABLE): $(RS_COMMON_TEST_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS) -o $@ gtest_main.a

$(FA_READER_TEST_EXECUTABLE): $(FA_READER_TEST_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS) -o $@ gtest_main.a

$(RS_BLOOM_TEST_EXECUTABLE): $(RS_BLOOM_TEST_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS) -o $@ gtest_main.a

$(ROLLING_HASH_COUNTER_TEST_EXECUTABLE): $(ROLLING_HASH_COUNTER_TEST_OBJECTS)
	$(CXX) $(LIB) $^ $(LDFLAGS) -o $@ gtest_main.a

karp_robin_hash_test: karp_robin_hash_test.cc karp_robin_hash.cc karp_robin_hash.h
	$(CXX) $(CPPFLAGS) $(LIB) $^ $(LDFLAGS) -o $@ gtest_main.a

.cc.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MD -c -o $@ $<
	cp $*.d $*.P; \
		sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
			-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
			rm -f $*.d

-include $(OBJECTS:.o=.P)

rs/rnasigs_pb2.py: proto/rnasigs.proto
	../third_party/bin/protoc --python_out=./ proto/rnasigs.proto
	mkdir -p python
	cp proto/rnasigs_pb2.py python/

proto/rnasigs.pb.h: proto/rnasigs.proto
	../third_party/bin/protoc --cpp_out=./ proto/rnasigs.proto


clean:
	rm -f $(sort $(OBJECTS)) $(EXECUTABLES) $(TESTS) gtest.a \
	gtest_main.a  *.P libbloomd/*.P libbloomd/murmurhash/*.P \
	libbloomd/spookyhash/*.P *.o

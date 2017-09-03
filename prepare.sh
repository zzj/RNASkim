# TODO: rewite everything in Python
# It is easier to maintain and incorporate with my other scripts.

export RS_HOME=`pwd`
export LDFLAGS="-L$RS_HOME/lib/gflags-2.0/.libs/ -L$RS_HOME/lib/glog-0.3.3/.libs/"
export CPPFLAGS="-I$RS_HOME/lib/gflags-2.0/src/  -I$RS_HOME/lib/glog-0.3.3/src/"

mkdir -p lib
cd lib
#wget --no-check-certificate  https://protobuf.googlecode.com/files/protobuf-2.5.0.tar.gz
wget --no-check-certificate https://github.com/google/protobuf/releases/download/v2.5.0/protobuf-2.5.0.tar.gz
tar zxf protobuf-2.5.0.tar.gz
cd protobuf-2.5.0
./configure --prefix=`pwd`/../../third_party/
make -j8
make install
cd python
python setup.py build
python setup.py install
cd ..
cd ..

#pip install leveldb protobuf pyfasta
pip install leveldb protobuf pyfasta --user

# install gtest
#wget http://googletest.googlecode.com/files/gtest-1.7.0.zip
#unzip gtest-1.7.0.zip
wget https://github.com/google/googletest/archive/release-1.7.0.zip
mv release-1.7.0.zip gtest-1.7.0.zip
unzip gtest-1.7.0.zip
mv googletest-release-1.7.0 gtest-1.7.0

#wget  --no-check-certificate https://gflags.googlecode.com/files/gflags-2.0.zip
wget  --no-check-certificate https://github.com/gflags/gflags/archive/v2.0.zip
mv v2.0.zip gflags-2.0.zip
unzip gflags-2.0.zip
cd gflags-2.0
./configure
make
cd ..

#wget  --no-check-certificate https://google-glog.googlecode.com/files/glog-0.3.3.tar.gz
#tar zxf glog-0.3.3.tar.gz
#cd glog-0.3.3
wget  --no-check-certificate https://github.com/google/glog/archive/v0.3.4.tar.gz
mv v0.3.4.tar.gz glog-0.3.4.tar.gz
tar zxf glog-0.3.4.tar.gz
mv glog-0.3.4 glog-0.3.3
cd glog-0.3.3
./configure
make

cd ..
cd ..

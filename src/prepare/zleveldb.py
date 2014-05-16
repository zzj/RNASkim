'''zleveldb library provides an interface to write protobuf messages
into a binary file.

'''

import struct
import leveldb

class ZLevelDB:
    def __init__(self, filename):
        self.db = leveldb.LevelDB(filename)

    def put(self, key, data):
        buf = data.SerializeToString()
        self.db.Put(key, buf, sync=True)

    def get(self, key, data):
        buf = None
        try:
            buf = self.db.Get(key)
        except KeyError, e:
            data = None
            return None
        data.ParseFromString(buf)
        return None

    def iter(self, cls):
        for k, v in self.db.RangeIter(include_value=True):
            data = cls()
            data.ParseFromString(v)
            yield k, data


def load_protobuf_data(fd, data):
    k = fd.read(4)
    if (len(k) != 4):
        return False
    m = struct.unpack("i", k)[0]
    buf = fd.read(m)
    if (len(buf) != m):
        return False
    data.ParseFromString(buf)
    return True

def load_all_protobuf_data(fd, cls):
    data = []
    d = cls()
    while load_protobuf_data(fd, d):
        data.append(d)
        d = cls()
    return data

def write_protobuf_data(fd, data):
    buf = data.SerializeToString()
    fd.write(struct.pack("i",len(buf)))
    fd.write(buf)

def load_single_protobuf_data(filename, data):
    fd = open(filename, "rb")
    data.ParseFromString(fd.read())
    fd.close()

def write_single_protobuf_data(filename, data):
    fd = open(filename, "wb")
    fd.write(data.SerializeToString())
    fd.close()

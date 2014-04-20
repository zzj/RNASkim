// The thread library in C++11 is passed by value.
// Therefore, it is tricky to pass in a thread class with members to
// the thread class, because the thread class holds the copy of the
// thread class and the members of the original object are not
// changed.

// This RSThread is designed to solve this problem by pass in a
// pointer of the thread class. In this way, the thread calls the
// member function of the original object, and no copy is done.

#ifndef RS_THREAD_H
#define RS_THREAD_H

#include <thread>
#include <vector>

namespace rs {

class ThreadInterface {
 public:
  virtual void run() = 0;
};

class RSThread {
 public :
  RSThread(ThreadInterface* thread) : thread_(thread) {}
  void operator() () {
    thread_->run();
  }
 private:
  ThreadInterface* thread_;
};

inline void barrier(std::vector<std::thread>& threads) {
  for (size_t i = 0; i < threads.size(); i ++) {
    threads[i].join();
  }
}


}  // namespace rs

#endif  // RS_THREAD_H

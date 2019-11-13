#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <queue>
#include <boost/bind.hpp>
#include <boost/thread.hpp>

namespace hlsutil {

class ThreadPool {
private:
	std::queue< boost::function< void() > > tasks_;
	boost::thread_group threads_;
	std::size_t available_;
	boost::mutex mutex_;
	boost::condition_variable condition_;
	bool running_;
public:
	ThreadPool(std::size_t pool_size);
	virtual ~ThreadPool();
	template < typename Task > void runTask( Task task );
	void runPool();
};

} /* namespace hlsutil */
#endif /* THREADPOOL_H_ */

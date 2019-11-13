// http://stackoverflow.com/questions/12215395/thread-pool-using-boost-asio/12267138#12267138

#include "ThreadPool.h"

namespace hlsutil {

/// @brief Constructor.
ThreadPool::ThreadPool( std::size_t pool_size ) : available_( pool_size ), running_( true ) {
	for ( std::size_t i = 0; i < pool_size; ++i ) {
		threads_.create_thread( boost::bind( &ThreadPool::runPool, this ) ) ;
	}
}

/// @brief Destructor.
ThreadPool::~ThreadPool() {
	// Set running flag to false then notify all threads.
	{
		boost::unique_lock< boost::mutex > lock( mutex_ );
		running_ = false;
		condition_.notify_all();
	} try {
		threads_.join_all();
	}
	// Suppress all exceptions.
	catch ( ... ) {}
}

/// @brief Add task to the thread pool if a thread is currently available.
template < typename Task >
void ThreadPool::runTask( Task task ) {
	boost::unique_lock< boost::mutex > lock( mutex_ );

	// If no threads are available, then return.
	if ( 0 == available_ ) {
		return;
	}

	// Decrement count, indicating thread is no longer available.
	--available_;

	// Set task and signal condition variable so that a worker thread will
	// wake up and use the task.
	tasks_.push( boost::function< void() >( task ) );
	condition_.notify_one();
}

/// @brief Entry point for pool threads.
void ThreadPool::runPool() {
	while( running_ ) {
		// Wait on condition variable while the task is empty and the pool is
		// still running.
		boost::unique_lock< boost::mutex > lock( mutex_ );
		while ( tasks_.empty() && running_ ) {
			condition_.wait( lock );
		}
		// If pool is no longer running, break out.
		if ( !running_ ) break;

		// Copy task locally and remove from the queue.  This is done within
		// its own scope so that the task object is destructed immediately
		// after running the task.  This is useful in the event that the
		// function contains shared_ptr arguments bound via bind.
		{
			boost::function< void() > task = tasks_.front();
			tasks_.pop();

			lock.unlock();

			// Run the task.
			try {
				task();
			}
			// Suppress all exceptions.
			catch ( ... ) {}
		}

		// Task has finished, so increment count of available threads.
		lock.lock();
		++available_;
	} // while running_
}

}

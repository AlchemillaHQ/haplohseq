/*
 * ThreadPoolTest.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#define BOOST_TEST_MODULE ThreadPoolTest
#include <boost/test/included/unit_test.hpp>

#include "ThreadPool.h"

BOOST_AUTO_TEST_SUITE(ThreadPoolTest)

//void work() {
//	std::cout << "work\n";
//};
//
//struct worker
//{
//  void operator()() {
//	  std::cout << "worker\n";
//  };
//};
//
//void moreWork( int ) {
//	std::cout << "more work\n";
//};

BOOST_AUTO_TEST_CASE(firstTest) {

//	hlsutil::ThreadPool pool(2);
//	pool.runTask(work);
//	pool.runTask(worker());
//	pool.runTask(boost::bind(moreWork, 5));
}

BOOST_AUTO_TEST_SUITE_END()

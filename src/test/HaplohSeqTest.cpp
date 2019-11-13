/*
 * HaplohSeqTest.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#define BOOST_TEST_MODULE HaplohSeqTest
#include <boost/test/included/unit_test.hpp>
#include "InputProcessor.h"
#include "StringUtil.h"

using namespace boost::assign; //bring += into scope

BOOST_AUTO_TEST_SUITE(HaplohSeqTest)

BOOST_AUTO_TEST_CASE(readHaplohAllelesTest) {
	hlsutil::StringUtil str;
	hlsutil::InputProcessor inputProc;
	std::vector< std::vector<char> > alleles;
	char delim = ' ';
	char colDelim = '\0';
	int allelesColNum = 0;
	inputProc.readAlleles("./resources/smalltest/inputs/short.hapguess", alleles, delim, colDelim, allelesColNum);
	BOOST_REQUIRE_EQUAL(alleles.size(), 2);
	BOOST_REQUIRE_EQUAL(alleles[0][0], '?');
	BOOST_REQUIRE_EQUAL(alleles[1][0], '?');
	BOOST_REQUIRE_EQUAL(alleles[0][139], 'A');
	BOOST_REQUIRE_EQUAL(alleles[1][139], 'B');
	BOOST_REQUIRE_EQUAL(alleles[0].size(), 140);
	BOOST_REQUIRE_EQUAL(alleles[1].size(), 140);
}

BOOST_AUTO_TEST_CASE(readHaplohseqAllelesTest) {
	hlsutil::StringUtil str;
	hlsutil::InputProcessor inputProc;
	std::vector< std::vector<char> > alleles;
	char delim = '\0';
	char colDelim = ' ';
	int allelesColNum = 2;

	// test use of indices
	std::vector<unsigned long> indices;
	indices += 1, 2, 3;
	inputProc.readAlleles("./resources/smalltest/inputs/mach.hapguess", alleles, indices, delim, colDelim, allelesColNum);
	BOOST_REQUIRE_EQUAL(alleles.size(), 2);
	BOOST_REQUIRE_EQUAL(alleles[0][0], 'G');
	BOOST_REQUIRE_EQUAL(alleles[0][1], 'G');
	BOOST_REQUIRE_EQUAL(alleles[0][2], 'C');
	BOOST_REQUIRE_EQUAL(alleles[1][0], 'G');
	BOOST_REQUIRE_EQUAL(alleles[1][1], 'T');
	BOOST_REQUIRE_EQUAL(alleles[1][2], 'C');

	// test ignoring indices
	std::vector< std::vector<char> > alleles2;
	inputProc.readAlleles("./resources/smalltest/inputs/mach.hapguess", alleles2, delim, colDelim, allelesColNum);
	BOOST_REQUIRE_EQUAL(alleles2.size(), 2);
	BOOST_REQUIRE_EQUAL(alleles2[0][0], 'G');
	BOOST_REQUIRE_EQUAL(alleles2[0][1], 'G');
	BOOST_REQUIRE_EQUAL(alleles2[0][2], 'G');
	BOOST_REQUIRE_EQUAL(alleles2[1][0], 'G');
	BOOST_REQUIRE_EQUAL(alleles2[1][1], 'G');
	BOOST_REQUIRE_EQUAL(alleles2[1][2], 'T');
	BOOST_REQUIRE_EQUAL(alleles2[0][187], 'C');
	BOOST_REQUIRE_EQUAL(alleles2[1][187], 'C');
	BOOST_REQUIRE_EQUAL(alleles2[0][189], 'G');
	BOOST_REQUIRE_EQUAL(alleles2[1][189], 'G');
	BOOST_REQUIRE_EQUAL(alleles2[0].size(), 190);
	BOOST_REQUIRE_EQUAL(alleles2[1].size(), 190);
}

BOOST_AUTO_TEST_SUITE_END()


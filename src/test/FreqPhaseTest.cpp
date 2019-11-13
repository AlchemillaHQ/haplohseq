/*
 * FreqPhaseTest.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#define BOOST_TEST_MODULE FreqPhaseTest
#include <algorithm>
#include <boost/assign/std/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include <sstream>
#include <vector>
#include "FreqPhase.h"

using namespace boost::assign; //bring += into scope

BOOST_AUTO_TEST_SUITE(FreqPhaseTest)

BOOST_AUTO_TEST_CASE(enumerateHetSwitchesHaplohTest) {
	std::vector<char> alleles0;
	std::vector<char> alleles1;
	alleles0 += 'B', 'A', 'A', 'B', 'B', 'B', 'A';
	alleles1 += 'B', 'A', 'B', 'A', 'A', 'A', 'B';

	std::vector< std::vector<char> > phasedAlleles;
	phasedAlleles += alleles0, alleles1;

	std::vector<double> freqs;
	freqs += 1.0, 0.0, 0.5, 0.4, 0.6, 0.7, 0.2;

	hlsutil::InputProcessor inputProc;
	haplohseq::FreqPhase fp(inputProc);
	unsigned int randomSeed = 0;
	std::vector< std::vector<char> > freqPhasedAlleles = fp.phaseAB(freqs, 0.5, phasedAlleles, randomSeed);

	std::vector<int> switches = fp.enumerateHetSwitchesHaploh(phasedAlleles, freqPhasedAlleles);
	BOOST_REQUIRE_EQUAL(switches.size(), 4);
	BOOST_REQUIRE_EQUAL(switches[0], 1);
	BOOST_REQUIRE_EQUAL(switches[1], 0);
	BOOST_REQUIRE_EQUAL(switches[2], 1);
	BOOST_REQUIRE_EQUAL(switches[3], 1);
}

BOOST_AUTO_TEST_CASE(phaseTest) {
	std::vector<char> alleles0;
	std::vector<char> alleles1;
	alleles0 += 'B', 'A', 'A', 'B', 'B';
	alleles1 += 'B', 'A', 'B', 'A', 'A';

	std::vector< std::vector<char> > alleles;
	alleles += alleles0, alleles1;

	std::vector<double> freqs;
	freqs += 1.0, 0.0, 0.5, 0.4, 0.6;

	hlsutil::InputProcessor inputProc;
	haplohseq::FreqPhase fp(inputProc);
	unsigned int randomSeed = 0;
	std::vector< std::vector<char> > phasedAlleles = fp.phaseAB(freqs, 0.5, alleles, randomSeed);

	BOOST_REQUIRE_EQUAL(phasedAlleles.size(), 2);
	BOOST_REQUIRE_EQUAL(phasedAlleles[0].size(), 5);
	BOOST_REQUIRE_EQUAL(phasedAlleles[1].size(), 5);

	std::string phAlleles0(phasedAlleles[0].begin(), phasedAlleles[0].end());
	BOOST_REQUIRE_EQUAL(phAlleles0, "BABAB");

	std::string phAlleles1(phasedAlleles[1].begin(), phasedAlleles[1].end());
	BOOST_REQUIRE_EQUAL(phAlleles1, "BAABA");
}

BOOST_AUTO_TEST_CASE(meanFreqTest) {
	std::vector<double> freqs;
	freqs += 1.0, 0.0, 0.5, 0.4, 0.6;

	hlsutil::InputProcessor inputProc;
	haplohseq::FreqPhase fp(inputProc);
	double mean = fp.meanValue(freqs);
	BOOST_REQUIRE_EQUAL(mean, 0.5);
	BOOST_REQUIRE_EQUAL(freqs[0], 1.0);
	BOOST_REQUIRE_EQUAL(freqs[4], 0.6);
}

BOOST_AUTO_TEST_CASE(medianFreqTest) {
	std::vector<double> freqs;
	freqs += 1.0, 0.0, 0.5, 0.4, 0.3;

	hlsutil::InputProcessor inputProc;
	haplohseq::FreqPhase fp(inputProc);
	double median = fp.medianValue(freqs);
	BOOST_REQUIRE_EQUAL(median, 0.4);
	BOOST_REQUIRE_EQUAL(freqs[0], 1.0);
	BOOST_REQUIRE_EQUAL(freqs[4], 0.3);
	freqs.clear();

	freqs += 1.0, 0.0, 0.5, 0.4, 0.3, 0.8;
	median = fp.medianValue(freqs);
	BOOST_REQUIRE_EQUAL(median, 0.45);
}

BOOST_AUTO_TEST_SUITE_END()


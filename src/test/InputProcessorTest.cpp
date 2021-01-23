/*
 * InputProcessorTest.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#define BOOST_TEST_MODULE InputProcessorTest
#include <boost/assign/std/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include "InputProcessor.h"

using namespace boost::assign; //bring += into scope

BOOST_AUTO_TEST_SUITE(InputProcessorTest)

BOOST_AUTO_TEST_CASE(extractTest) {
	std::vector<char> alleles0;
	std::vector<char> alleles1;
	alleles0 += 'A', 'C', 'G', 'T';
	alleles1 += 'A', 'T', 'T', '?';

	std::vector< std::vector<char> > alleles;
	alleles += alleles0, alleles1;

	std::vector<unsigned long> hetIndices;
	hetIndices += 1, 2;

	hlsutil::InputProcessor inputProc;

	std::vector<double> freqs;
	freqs += 1.0, 0.0, 0.5, 0.4;
	std::vector<double> hetFreqs = inputProc.extract(hetIndices, freqs);
	BOOST_REQUIRE_EQUAL(hetIndices.size(), hetFreqs.size());
	BOOST_REQUIRE_EQUAL(hetFreqs[0], 0.0);
	BOOST_REQUIRE_EQUAL(hetFreqs[1], 0.5);

	std::vector< std::vector<char> > hetAlleles = inputProc.extract(hetIndices, alleles);
	BOOST_REQUIRE_EQUAL(hetIndices.size(), hetAlleles[0].size());
	BOOST_REQUIRE_EQUAL(hetIndices.size(), hetAlleles[1].size());
	BOOST_REQUIRE_EQUAL(hetAlleles[0][0], 'C');
	BOOST_REQUIRE_EQUAL(hetAlleles[0][1], 'G');
	BOOST_REQUIRE_EQUAL(hetAlleles[1][0], 'T');
	BOOST_REQUIRE_EQUAL(hetAlleles[1][1], 'T');

}

BOOST_AUTO_TEST_CASE(hetIndicesTest) {
	std::vector<char> alleles0;
	std::vector<char> alleles1;
	alleles0 += 'A', 'C', 'G', 'T';
	alleles1 += 'A', 'T', 'T', '?';

	std::vector< std::vector<char> > alleles;
	alleles += alleles0, alleles1;

	hlsutil::InputProcessor inputProc;
	std::vector<unsigned long> hetIndices = inputProc.hetIndices(alleles);
	BOOST_REQUIRE_EQUAL(hetIndices.size(), 2);
	BOOST_REQUIRE_EQUAL(hetIndices[0], 1);
	BOOST_REQUIRE_EQUAL(hetIndices[1], 2);
}

BOOST_AUTO_TEST_CASE(isNoCallTest) {
	hlsutil::InputProcessor inputProc;
	BOOST_REQUIRE_EQUAL(inputProc.isNoCall('A'), false);
	BOOST_REQUIRE_EQUAL(inputProc.isNoCall('?'), true);

	std::vector<char> noCallDef;
	noCallDef += '.', ' ';
	hlsutil::InputProcessor inputProc2(noCallDef);
	BOOST_REQUIRE_EQUAL(inputProc2.isNoCall('?'), false);
}

BOOST_AUTO_TEST_CASE(isHetTest) {
	std::vector<char> noCallDef;
	noCallDef += '.', ' ', '?';
	hlsutil::InputProcessor inputProc(noCallDef);

	BOOST_REQUIRE_EQUAL(inputProc.isHet('A','T'), true);
	BOOST_REQUIRE_EQUAL(inputProc.isHet('T','T'), false);
	BOOST_REQUIRE_EQUAL(inputProc.isHet('A','?'), false);
	BOOST_REQUIRE_EQUAL(inputProc.isHet('A','.'), false);
	BOOST_REQUIRE_EQUAL(inputProc.isHet('A',' '), false);
	BOOST_REQUIRE_EQUAL(inputProc.isHet('C','\0'), true);
}

BOOST_AUTO_TEST_CASE(readAllelesTest) {
	hlsutil::StringUtil str;
	hlsutil::InputProcessor inputProc;
	std::vector< std::vector<char> > alleles;
	inputProc.readAlleles("./resources/smalltest/inputs/short.hapguess", alleles, ' ', '\t', 0);
	std::string actualAlleles0 =
			"?BABBABABBAABBBBABABBBBAAABBAABABBBBABBBBABABBBAAB"
			"BBBAABAABAAABABBABBBAABBBBBBAAABAABABABABABAABAABA"
			"BABABAABABAAABBBABAABABBBBAAABAAABBBBBBA";
	std::string actualAlleles1 =
			"?AAAABAAABABBBBBBABAAABBABBAAABBAABABBAAABABABABAA"
			"AABBBBBBBABBBABBBBBBAABBBABBBBABBBABBBABBBBAAAAABA"
			"BABABAABABAAABABABAABABBBBAAABAAAAABBAAB";
	std::string alleles0(alleles[0].begin(), alleles[0].end());
	std::string alleles1(alleles[1].begin(), alleles[1].end());

	BOOST_REQUIRE_EQUAL(alleles.size(), 2);
	BOOST_REQUIRE_EQUAL(alleles[0].size(), 140);
	BOOST_REQUIRE_EQUAL(alleles[1].size(), 140);
	BOOST_REQUIRE_EQUAL(alleles0, actualAlleles0);
	BOOST_REQUIRE_EQUAL(alleles1, actualAlleles1);
}

BOOST_AUTO_TEST_SUITE_END()

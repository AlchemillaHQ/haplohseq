/*
 * FreqPhase.h
 *
 * For B-allele frequencies:
 * (1)
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#ifndef FREQPHASE_H_
#define FREQPHASE_H_

#include <algorithm>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <numeric>
#include <vector>

#include "InputProcessor.h"

namespace haplohseq {

class FreqPhase {
public:
	FreqPhase(hlsutil::InputProcessor &inputProc);
	virtual ~FreqPhase();
	hlsutil::InputProcessor inputProc;

	std::vector<int> enumerateHetSwitchesHaploh(
			const std::vector< std::vector<char> > &phasedAlleles,
			const std::vector< std::vector<char> > &freqPhasedAlleles);

	std::vector< std::vector<char> > phase(
				const std::vector<double> &refFreqs,
				const double &refThreshold,
				const std::vector<char> &refAlleles,
				const std::vector<char> &altAlleles,
				const std::vector< std::vector<char> > &haplotypes,
				const unsigned int &randomSeed);

	std::vector< std::vector<char> > phaseAB(
				const std::vector<double> &bFreqs,
				const double &bThreshold,
				const std::vector< std::vector<char> > &haplotypes,
				const unsigned int &randomSeed);

	double meanValue(const std::vector<double>& values);
	double medianValue(const std::vector<double>& values);

	double meanValue(const std::vector<int>& values);
	double medianValue(const std::vector<int>& values);
};

} /* namespace haplohseq */
#endif /* FREQPHASE_H_ */

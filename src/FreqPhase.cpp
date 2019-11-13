/*
 * FreqPhase.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#include "FreqPhase.h"

namespace haplohseq {

FreqPhase::FreqPhase(hlsutil::InputProcessor &inputProc) {
	this->inputProc = inputProc;
}

FreqPhase::~FreqPhase() {}

/**
 * Vector of switches.  If freqPhasedAlleles is concordant with phasedAlleles
 * a '1' is put in switches.  If it is non-concordant, a '0' is put in switches.
 *
 * Determine switches for het sites only.
 * (1) Assign first het as concordant, keeping track
 *     of freq-phased haplotype1 (fph1).
 * (2) Keeping track of which phased haplotype fph1
 *     belongs to, if at allele i+1, fph1 is matching the same
 *     phased haplotype that allele i matched to,
 *     a '1' is added to the switches.  If a
 *     haplotype switch was made, a '0' is added to
 *     the switches.
 *
 * (Note) If fph1's allele doesn't match either phased
 * haplotype, throw an exception.
 */
std::vector<int> FreqPhase::enumerateHetSwitchesHaploh(
		const std::vector< std::vector<char> > &phasedAlleles,
		const std::vector< std::vector<char> > &freqPhasedAlleles) {

	std::vector<int> switches;
	std::vector<char> haplotype0 = freqPhasedAlleles[0];

	int concordant = 1;
	int discordant = 0;

	int previousState = -1;
	int currentState = -1;

	for (unsigned int i=0; i<haplotype0.size(); i++) {

		// only phase het sites
		bool het = this->inputProc.isHet(phasedAlleles[0][i], phasedAlleles[1][i]);
		if (!het) {
			continue;
		}

		if (haplotype0[i] == phasedAlleles[0][i]) {
			currentState = 0;
		} else {
			if (haplotype0[i] == phasedAlleles[1][i]) {
				currentState = 1;
			} else {
				currentState = -1;
			}
		}

		if (previousState != -1) {
			if (previousState == currentState) {
				switches.push_back(concordant);
			} else {
				switches.push_back(discordant);
			}
		}

		previousState = currentState;
	}

	return switches;
}

/**
 * Takes frequencies and alleles and performs phasing
 */
std::vector< std::vector<char> > FreqPhase::phase(
			const std::vector<double> &refFreqs,
			const double &refThreshold,
			const std::vector<char> &refAlleles,
			const std::vector<char> &altAlleles,
			const std::vector< std::vector<char> > &haplotypes,
			const unsigned int &randomSeed) {

//	std::cout << "Running frequency based phasing.\n";
//	std::cout << refFreqs.size() << "\n";
//	std::cout << refThreshold << "\n";
//	std::cout << refAlleles.size() << "\n";
//	std::cout << altAlleles.size() << "\n";
//	std::cout << haplotypes[0].size() << "\n";
//	std::cout << haplotypes[1].size() << "\n";
	std::vector< std::vector<char> > phasedAlleles;
	std::vector<char> alleles0, alleles1;

	// initialize random number generation from uniform (0,1) distribution
	boost::random::mt19937 gen((static_cast<unsigned int>(randomSeed)));
	boost::random::uniform_real_distribution<> unif(0,1);

	for (unsigned int i=0; i < refFreqs.size(); i++) {
		double refFreq = refFreqs[i];
		char alt = altAlleles[i];
		char ref = refAlleles[i];

		char allele0 = haplotypes[0][i];
		char allele1 = haplotypes[1][i];

		// if the site is not a het, just use the current
		// alleles and move to the next site.
		if (!this->inputProc.isHet(allele0, allele1)) {
			alleles0.push_back(allele0);
			alleles1.push_back(allele1);
			continue;
		}

		// if this is not a simple het (one allele
		// is ref and one allele is alt), then just randomly
		// assign the alleles.
		// as a side effect, anything that is not
		// A/B (for microarray data) or
		// A/C/G/T (for seq data), will be randomly assigned
		// and not phased.
		if ((std::toupper(allele0) != ref && std::toupper(allele1) != ref) ||
			(std::toupper(allele0) != alt && std::toupper(allele1) != alt)){
			alleles0.push_back(allele0);
			alleles1.push_back(allele1);
			continue;
		}

		// sites remaining can be phased.

		// In the infinitesimally small event that
		// freq == threshold, reset freq using a
		// uniform random number generator.
		while (refFreq == refThreshold) {
			refFreq = unif(gen);
		}
		if (refFreq > refThreshold) {
			alleles0.push_back(ref);
			alleles1.push_back(alt);
		} else {
			if (refFreq < refThreshold) {
				alleles0.push_back(alt);
				alleles1.push_back(ref);
			}
		}
	}

	phasedAlleles.push_back(alleles0);
	phasedAlleles.push_back(alleles1);
	std::cout << "Successfully completed frequency based phasing.\n";
	return phasedAlleles;
}

std::vector< std::vector<char> > FreqPhase::phaseAB(
			const std::vector<double> &bFreqs,
			const double &bThreshold,
			const std::vector< std::vector<char> > &haplotypes,
			const unsigned int &randomSeed) {

	std::vector<char> refAlleles;
	std::vector<char> altAlleles;

	for (unsigned int i = 0; i < bFreqs.size(); i++) {
		refAlleles.push_back('B');
		altAlleles.push_back('A');
	}
	return this->phase(bFreqs, bThreshold, refAlleles, altAlleles, haplotypes, randomSeed);
}

double FreqPhase::meanValue(const std::vector<double>& values) {
	return accumulate(values.begin(), values.end(), 0.0 )/values.size();
}

double FreqPhase::medianValue(const std::vector<double>& values) {
	double median;
	size_t size = values.size();
	std::vector<double> tempFreqs(values);
	sort(tempFreqs.begin(), tempFreqs.end());

	if (size  % 2 == 0) {
		median = (tempFreqs[size / 2 - 1] + tempFreqs[size / 2]) / 2;
	}
	else {
		median = tempFreqs[size / 2];
	}
	return median;
}

double FreqPhase::meanValue(const std::vector<int>& values) {
	return accumulate(values.begin(), values.end(), 0.0 )/values.size();
}

double FreqPhase::medianValue(const std::vector<int>& values) {
	double median;
	size_t size = values.size();
	std::vector<int> tempFreqs(values);
	sort(tempFreqs.begin(), tempFreqs.end());

	if (size  % 2 == 0) {
		median = (tempFreqs[size / 2 - 1] + tempFreqs[size / 2]) / 2;
	}
	else {
		median = tempFreqs[size / 2];
	}
	return median;
}

} /* namespace haplohseq */

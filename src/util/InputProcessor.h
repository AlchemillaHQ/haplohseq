/*
 * InputProcessor.h
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#ifndef INPUTPROCESSOR_H_
#define INPUTPROCESSOR_H_

#include <boost/assign/std/vector.hpp>
#include <boost/foreach.hpp>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>

#include "StringUtil.h"

namespace hlsutil {

class InputProcessor {
public:
	InputProcessor();
	InputProcessor(std::vector<char> noCallDef);
	virtual ~InputProcessor();

	static const std::string VCF_CHR, VCF_POS, VCF_REF, VCF_ALT, VCF_REF_DP, VCF_DP, VCF_ALT_DP, VCF_GT, VCF_REF_FREQ;
	static const std::string BAF_CHR, BAF_POS, BAF_REF, BAF_ALT, BAF_LOGRR, BAF_GT, BAF_RAF, BAF_SWITCH;

	std::vector<char> noCallDef;

	// These two methods do the same thing.
	// Can I use a generic vector for values?
	std::vector<double> extract(
			const std::vector<unsigned long> &indices,
			const std::vector<double> &values);

	std::vector<char> extract(
			const std::vector<unsigned long> &indices,
			const std::vector<char> &values);

	std::vector< std::vector<char> > extract(
			const std::vector<unsigned long> &indices,
			const std::vector< std::vector<char> > &values);

	std::vector<unsigned long> hetIndices(
			const std::vector< std::vector<char> > &alleles);

	bool isHet(const char &allele0, const char&allele1);

	bool isNoCall(const char &nucleotide);
	bool isNoCall(const std::string &value);

	/* Currently used for BAFs */
	unsigned long readAlleles(
			const std::string &filename,
			std::vector< std::vector<char> > &alleles,
			const char &delim='\0',
			const char &colDelim=' ',
			unsigned int allelesColNum=2);

	// allelescolnum of 2 is the default for mach
	// '\0' means the allele bases are not separated
	unsigned long readAlleles(
			const std::string &filename,
			std::vector< std::vector<char> > &alleles,
			std::vector<unsigned long> &indices,
			const char &delim='\0',
			const char &colDelim=' ',
			unsigned int allelesColNum=2);

	void readFrequencies(
			const std::string &filename,
			std::vector<double> &frequencies,
			const char &delim=' ');

	unsigned long readVcfInformativeGenotypes(
									const std::string &filename,
									unsigned int vcfMinDepth,
									std::vector<unsigned long> &indices,
									std::vector<double> &frequencies,
									std::vector<char> &informativeRef,
									std::vector<char> &informativeAlt,
									std::vector<unsigned int> &informativeDepth,
									std::vector<std::string> &informativePos);

	unsigned long readBafInformativeGenotypesChr(
									std::vector< std::vector<char> > &phasedAlleles,
									const std::string &filename,
									std::map<std::string, std::vector<std::string> >& bafStrFields,
									std::map<std::string, std::vector<unsigned int> >& bafIntFields,
									std::map<std::string, std::vector<double> >& bafDoubleFields,
									std::vector<std::string> &orderedChrs,
									std::vector<unsigned long> &bafInformativeIndices,
									std::map<std::string, std::vector<unsigned long> > &chrInformativeIndices,
									std::map<std::string, std::vector<double> > &chrInformativeBafFreqs,
									std::map<std::string, std::vector<char> > &chrInformativeB,
									std::map<std::string, std::vector<char> > &chrInformativeA,
									std::map<std::string, std::vector<double> > &chrInformativeLogRR,
									std::map<std::string, std::vector<std::string> > &chrInformativePos);


	unsigned long readVcfInformativeGenotypesChr(
									std::vector< std::vector<char> > &phasedAlleles,
									const std::string &filename,
									const std::string &vcfSampleName,
									std::map<std::string, std::vector<std::string> >& vcfStrFields,
									std::map<std::string, std::vector<unsigned int> >& vcfIntFields,
									std::map<std::string, std::vector<double> >& vcfDoubleFields,
									std::vector<std::string> &orderedChrs,
									unsigned int &vcfMinDepth,
									std::vector<unsigned long> &vcfIndices,
									std::map<std::string, std::vector<unsigned long> > &chrIndices,
									std::map<std::string, std::vector<double> > &chrFrequencies,
									std::map<std::string, std::vector<char> > &chrInformativeRef,
									std::map<std::string, std::vector<char> > &chrInformativeAlt,
									std::map<std::string, std::vector<unsigned int> > &chrInformativeDepth,
									std::map<std::string, std::vector<std::string> > &chrInformativePos);
};

} /* namespace hlsutil */
#endif /* INPUTPROCESSOR_H_ */

/*
 * InputProcessor.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#include "InputProcessor.h"

namespace hlsutil {

const std::string InputProcessor::VCF_CHR = "chr";
const std::string InputProcessor::VCF_POS = "pos";
const std::string InputProcessor::VCF_REF = "ref";
const std::string InputProcessor::VCF_ALT = "alt";
const std::string InputProcessor::VCF_DP = "depth";
const std::string InputProcessor::VCF_REF_DP = "refDepth";
const std::string InputProcessor::VCF_ALT_DP = "altDepth";
const std::string InputProcessor::VCF_GT = "genotype";
const std::string InputProcessor::VCF_REF_FREQ = "refFreq";

const std::string InputProcessor::BAF_CHR = "chr";
const std::string InputProcessor::BAF_POS = "pos";
const std::string InputProcessor::BAF_REF = "ref";
const std::string InputProcessor::BAF_ALT = "alt";
const std::string InputProcessor::BAF_LOGRR = "logRR";
const std::string InputProcessor::BAF_GT = "genotype";
const std::string InputProcessor::BAF_RAF = "raf";
const std::string InputProcessor::BAF_SWITCH = "switch";

// Set default no-call characters
InputProcessor::InputProcessor() {
	using namespace boost::assign; //bring += into scope
	this->noCallDef += '?', '.', ' ';
}

InputProcessor::InputProcessor(std::vector<char> noCallDef) {
	this->noCallDef = noCallDef;
}


InputProcessor::~InputProcessor() {

}

// These two methods do the same thing.
// Can I use a generic vector for values?
std::vector<double> InputProcessor::extract(
		const std::vector<unsigned long> &indices,
		const std::vector<double> &values) {
	std::vector<double> indexedValues;
	for (unsigned int i=0; i < indices.size(); i++) {
		indexedValues.push_back(values[indices[i]]);
	}
	return indexedValues;
}

std::vector<char> InputProcessor::extract(
		const std::vector<unsigned long> &indices,
		const std::vector<char> &values) {
	std::vector<char> indexedValues;
	for (unsigned int i=0; i < indices.size(); i++) {
		indexedValues.push_back(values[indices[i]]);
	}
	return indexedValues;
}

std::vector< std::vector<char> > InputProcessor::extract(
		const std::vector<unsigned long> &indices,
		const std::vector< std::vector<char> > &values) {
	std::vector< std::vector<char> > indexedValues;
	for (unsigned int i=0; i < values.size(); i++) {
		indexedValues.push_back(this->extract(indices, values[i]));
	}
	return indexedValues;
}

/**
 * Gets indices for all het sites.
 */
std::vector<unsigned long> InputProcessor::hetIndices(
						const std::vector< std::vector<char> > &alleles) {
	// TODO: consider throwing exceptions if alleles.size() != 2
	// or if alleles[0].size() != alleles[1].size()

	// loop through alleles and record indices for het sites
	std::vector<unsigned long> hetIndices;
	for (unsigned int i = 0; i < alleles[0].size(); i++) {
		char allele1 = alleles[0][i];
		char allele2 = alleles[1][i];

		if (this->isNoCall(allele1) or this->isNoCall(allele2)) {
			continue;
		}

		if (allele1 != allele2) {
			hetIndices.push_back(i);
		}
	}
	return hetIndices;
}

bool InputProcessor::isHet(const char &allele0, const char&allele1) {

	if (this->isNoCall(allele0) or this->isNoCall(allele1)) {
		return false;
	}

	if (allele0 != allele1) {
		return true;
	}
	return false;
}

bool InputProcessor::isNoCall(const char &nucleotide) {
	for (unsigned int i = 0; i < this->noCallDef.size(); i++) {
		if (nucleotide == this->noCallDef[i]) {
			return true;
		}
	}
	return false;
}

bool InputProcessor::isNoCall(const std::string &value) {
	for (unsigned int i = 0; i < this->noCallDef.size(); i++) {
		std::string noCall = "";
		noCall += noCallDef[i];
		if (value == noCall) {
			return true;
		}
	}
	return false;
}

/*
 * This method will read alleles from a file.  It assumes that each token in the file line
 * is of length 1.
 */
unsigned long InputProcessor::readAlleles(		const std::string &filename,
										std::vector< std::vector<char> > &alleles,
										const char &delim,
										const char &colDelim,
										unsigned int allelesColNum) {
	std::vector<unsigned long> emptyVector;	// an alternative is to use pointers, and set to null, but that logic would propagate.
	return readAlleles(filename, alleles, emptyVector, delim, colDelim, allelesColNum);
}

// if you don't want to use a delimiter,
// set the delimiter to '\0' = char(0)
unsigned long InputProcessor::readAlleles(	const std::string &filename,
									std::vector< std::vector<char> > &alleles,
									std::vector< unsigned long> &informativeIndices,
									const char &delim,
									const char &colDelim,
									unsigned int allelesColNum) {
	StringUtil str;
	std::ifstream infile(filename.c_str(), std::ifstream::in);
	std::string line;

	unsigned long numPhasedGenotypes = 0;

	while(std::getline(infile, line)) {
	    if (str.trim(line) == "") {
	    	continue;
	    }

	    std::vector<std::string> lineTokens;
	    std::vector<std::string> alleleTokens;
	    std::vector<char> lineAlleles;
	    str.split(line, colDelim, lineTokens);
	    if (delim != '\0') {
		    // get alleles column
		    if (lineTokens.size() < allelesColNum + 1) {
		    	str.split(lineTokens[0], delim, alleleTokens);
		    } else {
		    	str.split(lineTokens[allelesColNum], delim, alleleTokens);
		    }

		    if (!informativeIndices.empty()) {
				for (unsigned long i = 0; i < informativeIndices.size(); i++) {
					lineAlleles.push_back(alleleTokens[informativeIndices[i]][0]);	//TODO NOTE: this assumes each allele length is 1!!!!
				}
		    }
		    else {
				for (unsigned long i = 0; i < alleleTokens.size(); i++) {
					lineAlleles.push_back(alleleTokens[i][0]);	//TODO NOTE: this assumes each allele length is 1!!!!
				}
		    }

    		// extra logic to ensure the lengths of the haplotypes are equivalent - also allows the return of the haplotype
    		// length for comparison with the vcf file
    	    if (numPhasedGenotypes == 0) {
    	    	numPhasedGenotypes = alleleTokens.size();
    	    } else {
    	    	if (alleleTokens.size() != numPhasedGenotypes) {
    	    		std::cerr << "ERROR: length of haplotypes are not equivalent in " << filename << std::endl;
    	    	}
    	    }

	    } else {
	    	std::string haplotype;
    		if (lineTokens.size() < allelesColNum + 1) {
    			haplotype = lineTokens[0];
    		} else {
    			haplotype = lineTokens[allelesColNum];
    		}
	    	if (!informativeIndices.empty()) {
				for (unsigned long i = 0; i < informativeIndices.size(); i++) {
					lineAlleles.push_back(haplotype[informativeIndices[i]]);
				}
	    	} else {
	    		for (unsigned long i = 0; i < haplotype.size(); i++) {
	    			lineAlleles.push_back(haplotype[i]);
	    		}
	    	}

    		// extra logic to ensure the lengths of the haplotypes are equivalent - also allows the return of the haplotype
    		// length for comparison with the vcf file
    	    if (numPhasedGenotypes == 0) {
    	    	numPhasedGenotypes = haplotype.length();
    	    } else {
    	    	if (haplotype.length() != numPhasedGenotypes) {
    	    		std::cerr << "ERROR: length of haplotypes are not equivalent in " << filename << std::endl;
    	    	}
    	    }
	    }

    	if (lineAlleles.size() > 0) {
    		alleles.push_back(lineAlleles);
    	}
	}
	infile.close();

	return numPhasedGenotypes;
}

void InputProcessor::readFrequencies(const std::string &filename, std::vector<double> &frequencies, const char &delim) {
	StringUtil str;
	std::ifstream infile(filename.c_str(), std::ifstream::in);

	std::string line;

	std::getline(infile, line);
	std::vector<std::string> lineTokens;
	str.split(line, delim, lineTokens);
	for (unsigned int i = 0; i < lineTokens.size(); i++) {
		std::string value = lineTokens[i];
		if (this->isNoCall(value)) {
			frequencies.push_back(atof(value.c_str()));
		} else {
			frequencies.push_back(-1.0);  // for missing markers putting in an uninformative value.
		}
	}
	infile.close();
}

//unsigned long InputProcessor::readVcfInformativeGenotypes(
//								const std::string &filename,
//								unsigned int vcfMinDepth,
//								std::vector<unsigned long> &indices,
//								std::vector<double> &refFrequencies,
//								std::vector<char> &informativeRef,
//								std::vector<char> &informativeAlt,
//								std::vector<unsigned int> &informativeDepth,
//								std::vector<std::string> &informativePos) {
//
//	// VCF columns (standard)
//	// assumes VCF is only 1 chromosome!!!!!!!!
//	unsigned int CHR			= 0;
//	unsigned int POS			= 1;
//	unsigned int REF			= 3;
//	unsigned int ALT			= 4;
//	unsigned int FORMAT_FLAGS	= 8;
//	unsigned int FORMAT_VALUES	= 9;	// assuming 1 sample
//
//	std::string ALLELE_DEPTH	= "AD";
//	std::string GENOTYPE		= "GT";
//	std::string DEPTH			= "DP";
//
//	char delim = '\t';
//	char formatDelim = ':';
//	StringUtil str;
//	std::ifstream infile(filename.c_str(), std::ifstream::in);
//
//	std::string line;
//	unsigned long posIndex = 0;
//
//	std::stringstream buffer;
//	buffer << infile.rdbuf();
//	while(std::getline(buffer, line)) {
////	while(std::getline(infile, line)) {
//		if (str.trim(line)[0] == '#' || str.trim(line) == "") {
//			continue;	// comment line
//		}
//		std::vector<std::string> lineTokens;
//		str.split(line, delim, lineTokens);
//		std::vector<std::string> formatFlags;
//		std::vector<std::string> formatValues;
//
//		// split format flags and put in formatFlags vector
//		str.split(lineTokens[FORMAT_FLAGS], formatDelim, formatFlags);
//		// formatValues - order corresponds with formatFlags
//		str.split(lineTokens[FORMAT_VALUES], formatDelim, formatValues);
//		// put in a map for easy access
//		std::map<std::string, std::string> flagValues;
////		for (unsigned int i = 0; i < formatFlags.size(); i++) {
//		int i = 0;
//		BOOST_FOREACH(std::string flag, formatFlags) {
//			flagValues[flag] = formatValues[i++];
//		}
//
//		unsigned int depth = 0;
//		std::map<std::string, std::string>::iterator it = flagValues.find(GENOTYPE);
//		std::string gt = it->second;
//		if (gt != "0/1" && gt != "1/0") {
//			posIndex += 1;
//			continue;			// non-informative site
//		}
//
//		it = flagValues.find(DEPTH);
//		if (it != flagValues.end()) {
//			depth = atoi((it->second).c_str());	// element found
//			if (depth < vcfMinDepth) {
//				posIndex += 1;
//				continue;		// unreliable site because lack of coverage
//			}
//
//			std::string chr = lineTokens[CHR];
//			std::string pos = lineTokens[POS];
//			std::string ref = lineTokens[REF];
//			std::string alt = lineTokens[ALT];
//
//			if (ref.length() > 1 or alt.length() > 1) {
//					posIndex += 1;
//					continue;	// unreliable site because ref or alt is not a simple SNP (may not be dealt with correctly downstream)
//			}
//
//			std::string ad = flagValues.find(ALLELE_DEPTH)->second;
//			std::vector<std::string> alleleDepths;
//			str.split(ad, ',', alleleDepths);
//			unsigned int refDepth = atoi(alleleDepths[0].c_str());
//			unsigned int altDepth = atoi(alleleDepths[1].c_str());
//
//			double refFrequency = ((double) refDepth)/((double)(refDepth + altDepth));
//
//			// only 1 character refs and alts make it this far, so [0] access is ok
//			informativePos.push_back(chr + '\t' + pos + '\t' + str.intToStr(posIndex));
//			informativeRef.push_back(ref[0]);
//			informativeAlt.push_back(alt[0]);
//			refFrequencies.push_back(refFrequency);
//			informativeDepth.push_back(depth);
//			indices.push_back(posIndex);
//			posIndex += 1;
//
//		} else {
//			posIndex += 1;
//			continue;			// non-informative site
//		}
//	}
//	std::cout << "done loading " + filename + "\n";
//	infile.close();
//
//	return posIndex;
//}

unsigned long InputProcessor::readBafInformativeGenotypesChr(
								std::vector< std::vector<char> > &phasedAlleles,
								const std::string &filename,
								std::map<std::string, std::vector<std::string> >& bafStrFields,
								std::map<std::string, std::vector<unsigned int> >& bafIntFields,
								std::map<std::string, std::vector<double> >& bafDoubleFields,
								std::vector<std::string> &orderedChrs,
								std::vector<unsigned long> &bafInformativeIndices,
								std::map<std::string, std::vector<unsigned long> > &chrInformativeIndices,
								std::map<std::string, std::vector<double> > &chrInformativeRefFreqs,
								std::map<std::string, std::vector<char> > &chrInformativeRef,
								std::map<std::string, std::vector<char> > &chrInformativeAlt,
								std::map<std::string, std::vector<double> > &chrInformativeLogRR,
								std::map<std::string, std::vector<std::string> > &chrInformativePos) {

	// BAF columns (standard)
//	unsigned int CHR			= 2;
//	unsigned int POS			= 3;
//	unsigned int A_ALLELE		= 4;
//	unsigned int B_ALLELE		= 5;
//	unsigned int GENOTYPE		= 6;
//	unsigned int BAF			= 8;	// TODO needs to be switched with LOG_RR after I fix bug in Jerry's gsmdf report
//	unsigned int LOG_RR			= 7;

	unsigned int CHR			= 0;
	unsigned int POS			= 1;
	unsigned int REF			= 2;
	unsigned int ALT			= 3;
	unsigned int GENOTYPE		= 12;
	unsigned int RAF			= 10;
	unsigned int LOG_RR			= 11;

	char delim = '\t';
//	char formatDelim = ':';
	StringUtil str;
	std::ifstream infile(filename.c_str(), std::ifstream::in);

	std::string line;
	unsigned long posIndex = 0;
	std::string previousChr = "NA";

//	std::stringstream buffer;
//	buffer << infile.rdbuf();
//	while(std::getline(buffer, line)) {
	getline(infile,line);	// read past header (TODO: next time keep track of column headers)
	while(std::getline(infile, line)) {
//		std::map<std::string, std::string> vcfEntry;
		std::string trimmedLine = str.trim(line);
		if (trimmedLine[0] == '#' || trimmedLine == "") {
			continue;	// comment line
		}
		std::vector<std::string> lineTokens;
		str.split(line, delim, lineTokens);

		// [1] COLLECT VCF INFO
		std::string chr = lineTokens[CHR];
		std::string pos = lineTokens[POS];
		std::string ref = lineTokens[REF];
		std::string alt = lineTokens[ALT];
		std::string gt = lineTokens[GENOTYPE];

		double refFreq = atof(lineTokens[RAF].c_str());
		double logRR = atof(lineTokens[LOG_RR].c_str());

		// populate bafFields for later use
		bafStrFields[BAF_CHR].push_back(chr);
		bafIntFields[BAF_POS].push_back(atoi(pos.c_str()));
		bafStrFields[BAF_REF].push_back(ref);
		bafStrFields[BAF_ALT].push_back(alt);
		bafDoubleFields[BAF_LOGRR].push_back(logRR);
		bafDoubleFields[BAF_RAF].push_back(refFreq);
		bafStrFields[BAF_GT].push_back(gt);

		// [2] FILTER UNINFORMATIVE SITES

//		if (depth < vcfMinDepth) {  // this will capture depths of 0 (or missing depths)
//			posIndex += 1;
//			continue;		// unreliable site because lack of coverage
//		}

		// NOTE: a het call here should also be a het in the estimated phased genotypes
		// but if we eventually use paired samples, we should filter hets based on
		// the estimated phased genotypes that might be based on the "normal sample"
		// where it is easier to identify het sites.  If you can imagine, in a "tumor sample"
		// if there is a complete LOH, the site will be ignored by haplohseq which does
		// not currently use a paired sample.
		// TODO: GET GT FROM PHASED ALLELES (WILL BE FROM NORMAL IN PAIRED MODE OR THE TUMOR IN SINGLE-SAMPLE MODE)
		char allele1 = phasedAlleles[0][posIndex];
		char allele2 = phasedAlleles[1][posIndex];
		if (allele1 == allele2 || gt == "./." || islower(allele1) || islower(allele2) ) {
		//if (std::toupper(allele1) == std::toupper(allele2) || gt == "./.") { // || islower(allele1) || islower(allele2) ) {
			posIndex += 1;
			continue;			// non-informative site
		}

		if (ref.length() > 1 || alt.length() > 1) {
				posIndex += 1;
				continue;	// unreliable site because ref or alt is not a simple SNP (may not be dealt with correctly downstream)
		}

		// [3] CREATE CHR SPECIFIC VECTORS

		// usually previousChr will equal chr and this will get bypassed and the search through the orderedChrs set will get bypassed
		if (previousChr.compare(chr) != 0 && (std::find(orderedChrs.begin(), orderedChrs.end(), chr) == orderedChrs.end())) {

			previousChr = chr;
			std::cout << "adding " << chr << "\n";
			orderedChrs.push_back(chr);

			std::vector<unsigned long> indices;
			std::vector<double> informativeRefFreqs;
			std::vector<char> informativeRef;
			std::vector<char> informativeAlt;
			std::vector<double> informativeLogRR;
			std::vector<std::string> informativePos;

			chrInformativeIndices[chr] = indices;
			chrInformativeRefFreqs[chr] = informativeRefFreqs;
			chrInformativeRef[chr] = informativeRef;
			chrInformativeAlt[chr] = informativeAlt;
			chrInformativeLogRR[chr] = informativeLogRR;
			chrInformativePos[chr] = informativePos;
		}

		// only 1 character refs and alts make it this far, so [0] access is ok
		chrInformativePos[chr].push_back(chr + '\t' + pos + '\t' + str.intToStr(posIndex));
		chrInformativeRef[chr].push_back(ref[0]);
		chrInformativeAlt[chr].push_back(alt[0]);
		chrInformativeRefFreqs[chr].push_back(refFreq);
		chrInformativeLogRR[chr].push_back(logRR);
		chrInformativeIndices[chr].push_back(posIndex);	// not using chrPosIndex
		bafInformativeIndices.push_back(posIndex);

//		std::cout << "position " << posIndex << "\t" << trimmedLine << "\n";
		posIndex += 1;

	}
	std::cout << "done loading " + filename + "\n";
	infile.close();

	return posIndex;
}

unsigned long InputProcessor::readVcfInformativeGenotypesChr(
								std::vector< std::vector<char> > &phasedAlleles,
								const std::string &filename,
								const std::string &vcfSampleName,
								std::map<std::string, std::vector<std::string> >& vcfStrFields,
								std::map<std::string, std::vector<unsigned int> >& vcfIntFields,
								std::map<std::string, std::vector<double> >& vcfDoubleFields,
								std::vector<std::string> &orderedChrs,
								unsigned int &vcfMinDepth,
								std::vector<unsigned long> &vcfInformativeIndices,
								std::map<std::string, std::vector<unsigned long> > &chrVcfInformativeIndices,
								std::map<std::string, std::vector<double> > &chrInformativeRefFrequencies,
								std::map<std::string, std::vector<char> > &chrInformativeRef,
								std::map<std::string, std::vector<char> > &chrInformativeAlt,
								std::map<std::string, std::vector<unsigned int> > &chrInformativeDepth,
								std::map<std::string, std::vector<std::string> > &chrInformativePos) {

	// VCF columns (standard)
	std::string headerPrefix = "#CHROM";
	std::map<std::string, unsigned int> headerIndices;
	unsigned int CHR			= 0;
	unsigned int POS			= 1;
	unsigned int REF			= 3;
	unsigned int ALT			= 4;
	unsigned int FORMAT_FLAGS	= 8;
	unsigned int FORMAT_VALUES	= 9;	// assuming 1 sample

	std::string ALLELE_DEPTH	= "AD";
	std::string GENOTYPE		= "GT";
	std::string DEPTH			= "DP";
//	std::string REF_DEPTH		= "REFD";
//	std::string ALT_DEPTH		= "ALTD";
//	std::string REF_FREQ		= "REF_FREQ";

	char delim = '\t';
	char formatDelim = ':';
	StringUtil str;
	std::ifstream infile(filename.c_str(), std::ifstream::in);

	std::string line;
	unsigned long posIndex = 0;
	std::string previousChr = "NA";

//	std::stringstream buffer;
//	buffer << infile.rdbuf();
//	while(std::getline(buffer, line)) {
	while(std::getline(infile, line)) {
//		std::map<std::string, std::string> vcfEntry;
		std::string trimmedLine = str.trim(line);
		if (trimmedLine[0] == '#' || trimmedLine == "") {
			if (trimmedLine.substr(0,headerPrefix.size()) == headerPrefix) {
				std::vector<std::string> headerTokens;
				str.split(line, delim, headerTokens);
				unsigned int i = 0;
				BOOST_FOREACH(std::string col, headerTokens) {
					headerIndices[col] = i++;
				}
			}
			continue;	// comment line
		}
		std::vector<std::string> lineTokens;
		str.split(line, delim, lineTokens);
		std::vector<std::string> formatFlags;
		std::vector<std::string> formatValues;

		// split format flags and put in formatFlags vector
		str.split(lineTokens[FORMAT_FLAGS], formatDelim, formatFlags);
		// formatValues - order corresponds with formatFlags

		// by default, the first sample is used
		if (vcfSampleName != "") {
			FORMAT_VALUES = headerIndices[vcfSampleName];
		}

		str.split(lineTokens[FORMAT_VALUES], formatDelim, formatValues);
		// put in a map for easy access
		std::map<std::string, std::string> flagValues;
		int i = 0;
		BOOST_FOREACH(std::string flag, formatFlags) {
			flagValues[flag] = formatValues[i++];
		}

		// [1] COLLECT VCF INFO
		std::string chr = lineTokens[CHR];
		std::string pos = lineTokens[POS];
		std::string ref = lineTokens[REF];
		std::string alt = lineTokens[ALT];
		std::string gt = flagValues[GENOTYPE];

		unsigned int depth = 0;
		std::map<std::string, std::string>::iterator it = flagValues.find(DEPTH);
		if (it != flagValues.end()) {
			depth = atoi((it->second).c_str());	// element found
		}

		double refFrequency = 0.0;
		unsigned int refDepth = 0;
		unsigned int altDepth = 0;
		if (depth != 0) {
			// check for allelic depth - there are cases where depth != 0 and allelic depth is not provided
			it = flagValues.find(ALLELE_DEPTH);
			if (it != flagValues.end()) {
				std::string ad = flagValues.find(ALLELE_DEPTH)->second;
				std::vector<std::string> alleleDepths;
				str.split(ad, ',', alleleDepths);
				altDepth = atoi(alleleDepths[1].c_str());
				if (alleleDepths[0].compare(".") == 0) {
					refDepth = depth - altDepth;
				} else {
					refDepth = atoi(alleleDepths[0].c_str());
				}

				refFrequency = ((double) refDepth)/((double)(refDepth + altDepth));
			}
		}

		// populate vcfFields for later use
		vcfStrFields[VCF_CHR].push_back(chr);
		vcfIntFields[VCF_POS].push_back(atoi(pos.c_str()));
		vcfStrFields[VCF_REF].push_back(ref);
		vcfStrFields[VCF_ALT].push_back(alt);
		vcfIntFields[VCF_DP].push_back(depth);
		vcfDoubleFields[VCF_REF_FREQ].push_back(refFrequency);
		vcfStrFields[VCF_GT].push_back(gt);

		// [2] FILTER UNINFORMATIVE SITES
		if (depth < vcfMinDepth) {  // this will capture depths of 0 (or missing depths)
			posIndex += 1;
			continue;		// unreliable site because lack of coverage
		}

		// DETERMINE HET STATUS FROM PHASED ALLELES (which comes from normal sample in paired mode)
		char allele1 = phasedAlleles[0][posIndex];
		char allele2 = phasedAlleles[1][posIndex];

		// imputed sites are noninformative
		if (allele1 == allele2 || gt == "./." || islower(allele1) || islower(allele2) ) {

		// imputed sites are informative
		//if (std::toupper(allele1) == std::toupper(allele2) || gt == "./.") { // || islower(allele1) || islower(allele2) ) {
			posIndex += 1;
			continue;			// non-informative site
		}

		if (ref.length() > 1 || alt.length() > 1) {
				posIndex += 1;
				continue;	// unreliable site because ref or alt is not a simple SNP (may not be dealt with correctly downstream)
		}

		// [3] CREATE CHR SPECIFIC VECTORS
		// usually previousChr will equal chr and this will get bypassed and the search through the orderedChrs set will get bypassed
		if (previousChr.compare(chr) != 0 && (std::find(orderedChrs.begin(), orderedChrs.end(), chr) == orderedChrs.end())) {

			previousChr = chr;
//			std::cout << "adding " << chr << "\n";
			orderedChrs.push_back(chr);

			std::vector<unsigned long> indices;
			std::vector<double> informativeRefFreqs;
			std::vector<char> informativeRef;
			std::vector<char> informativeAlt;
			std::vector<unsigned int> informativeDepth;
			std::vector<std::string> informativePos;

			chrVcfInformativeIndices[chr] = indices;
			chrInformativeRefFrequencies[chr] = informativeRefFreqs;
			chrInformativeRef[chr] = informativeRef;
			chrInformativeAlt[chr] = informativeAlt;
			chrInformativeDepth[chr] = informativeDepth;
			chrInformativePos[chr] = informativePos;
		}

		// only 1 character refs and alts make it this far, so [0] access is ok
		chrInformativePos[chr].push_back(chr + '\t' + pos + '\t' + str.intToStr(posIndex));
		chrInformativeRef[chr].push_back(ref[0]);
		chrInformativeAlt[chr].push_back(alt[0]);
		chrInformativeRefFrequencies[chr].push_back(refFrequency);
		chrInformativeDepth[chr].push_back(depth);
		chrVcfInformativeIndices[chr].push_back(posIndex);	// not using chrPosIndex
		vcfInformativeIndices.push_back(posIndex);

//		std::cout << "position " << posIndex << "\t" << trimmedLine << "\n";
		posIndex += 1;

	}
	infile.close();

	return posIndex;
}

} /* namespace hlsutil */

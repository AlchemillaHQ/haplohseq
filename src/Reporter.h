/*
 * Reporter.h
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#ifndef REPORTER_H_
#define REPORTER_H_

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "Hmm.h"
#include "InputProcessor.h"
#include "StringUtil.h"
#include "VcfUtil.h"

//#include <cmath>
//#include <fstream>
//#include <iomanip>
//#include <iostream>
//#include <iterator>
//#include <list>

namespace haplohseq {

class Reporter {
public:
	Reporter();
	virtual ~Reporter();

	void writeToFile(std::string& filename, std::vector<double>& values, std::string& delimiter, std::string& header);
	void writeToFile(std::string& filename, std::vector<double>& values, std::string& delimiter);

	void writeToFile(std::string& filename, std::vector<unsigned long>& values, std::string& delimiter, std::string& header);
	void writeToFile(std::string& filename, std::vector<unsigned long>& values, std::string& delimiter);

	void writeToFile(std::string& filename, std::vector<std::string>& values, std::string& delimiter, std::string& header);
	void writeToFile(std::string& filename, std::vector<std::string>& values, std::string& delimiter);

	void reportBafBasedResults(	std::vector<std::string>& orderedChrs,
							std::map<std::string, boost::shared_ptr<Hmm> >& chrHmms,
							std::vector<unsigned long>& informativeIndices,
							std::vector< std::vector<char> > phasedAlleles,
							std::map<std::string, std::vector<std::string> >& bafStrFields,
							std::map<std::string, std::vector<unsigned int> >& bafIntFields,
							std::map<std::string, std::vector<double> >& bafDoubleFields,
							std::string& destinationDir,
							std::string& outputPrefix);

	void reportVcfBasedResults(std::vector<std::string>& orderedChrs,
				std::map<std::string, boost::shared_ptr<Hmm> >& chrHmms,
				std::vector<unsigned long>& informativeIndices,
				std::vector< std::vector<char> > phasedAlleles,
				std::map<std::string, std::vector<std::string> >& vcfStrFields,
				std::map<std::string, std::vector<unsigned int> >& vcfIntFields,
				std::map<std::string, std::vector<double> >& vcfDoubleFields,
				std::string& destinationDir,
				std::string& outputPrefix);
};
} /* namespace haplohseq */

#endif /* REPORTER_H_ */

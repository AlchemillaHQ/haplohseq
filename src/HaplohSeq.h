/*
 * HaplohSeq.h
 *
 *  Created on: Jan 2, 2013
 *      Author: fasaan
 */

#ifndef HAPLOHSEQ_H_
#define HAPLOHSEQ_H_

#include <boost/any.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "Hmm.h"
#include "FreqPhase.h"
#include "InputProcessor.h"
#include "Reporter.h"
#include "StringUtil.h"
#include "ThreadPool.h"
#include "VcfUtil.h"

class HaplohSeq {
public:
	HaplohSeq();
	virtual ~HaplohSeq();

	static const std::string 	INIT, EST_TRANS, EST_ABERRANT_EMIT, EST_NORMAL_EMIT, EVENT_PREVALENCE, EVENT_MB, GENOME_MB,
								INITIAL_PARAM_NORMAL, INITIAL_PARAM_EVENT, NUM_STATES, ITERATIONS, LAMBDA_0, LAMBDA_1,
								MEDIAN, NORMAL_STATE, RANDOMSEED, REF, ALT, B_ALLELE, A_ALLELE, HAPLOTYPE1, HAPLOTYPE2;

	void runBafHaplohseq(	std::string& obsType,
							unsigned int& numStates,
							double& epsilon,
							unsigned int& iterations,
							double& eventPrevalence,
							unsigned int& eventMb,
							unsigned int& genomeMb,
							double& initialParamNormal,
							double& initialparamEvent,
							double& endParamEvent,
							unsigned int& numParamEventStarts,
							std::string& bafFilename,
							unsigned int& numThreads,
							std::string& phasedAllelesFilename,
							std::string& destinationDir,
							std::map<std::string, bool>& estimateParam,
							std::string& outputPrefix,
							const unsigned int& randomSeed);

	void runVcfHaplohseq(	std::string& obsType,
							unsigned int& numStates,
							double& epsilon,
							unsigned int& iterations,
							double& eventPrevalence,
							unsigned int& eventMb,
							unsigned int& genomeMb,
							double& initialParamNormal,
							double& initialparamEvent,
							double& endParamEvent,
							unsigned int& numParamEventStarts,
							std::string& vcfFilename,
							std::string& vcfSampleName,
							unsigned int& vcfDepthThreshold,
							unsigned int& numThreads,
							std::string& phasedAllelesFilename,
							std::string& destinationDir,
							std::map<std::string, bool>& estimateParam,
							std::string& outputPrefix,
							const unsigned int& randomSeed);

	void runBafHaplohseqChr(std::string &chr,
							std::map<std::string, boost::shared_ptr<Hmm> > &chrHmms,
							haplohseq::FreqPhase &fp,
							const std::vector<double> &bafFreqs,
							const std::vector<double> &logRRs,
							std::map<std::string, boost::any> &primitives,
							boost::shared_ptr<std::map<std::string, std::vector<char> > > &allelesPtr);

	void runVcfHaplohseqChr(std::string &chr,
							std::map<std::string, boost::shared_ptr<Hmm> > &chrHmms,
							haplohseq::FreqPhase &fp,
							const std::vector<double> &refFreqs,
							const std::vector<unsigned int> &coverages,
							std::map<std::string, boost::any> &primitives,
							boost::shared_ptr<std::map<std::string, std::vector<char> > > &allelesPtr);

};

#endif /* HAPLOHSEQ_H_ */

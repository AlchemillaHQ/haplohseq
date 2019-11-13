/*
 * HaplohSeq.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */


#include "HaplohSeq.h"

// NOTE:
// The reason why there are so many constants is that boost multi-threading
// has a maximum of 9 parameters for binding methods.  Using the constants
// allows variables to be stored in maps.  Using the maps reduces the number
// of parameters being bound because a map of many variables only counts as
// a single parameter.

const std::string HaplohSeq::INIT = "init";
const std::string HaplohSeq::EST_TRANS = "estTrans";
const std::string HaplohSeq::EST_ABERRANT_EMIT = "estAberrantEmit";
const std::string HaplohSeq::EST_NORMAL_EMIT = "estNormalEmit";
const std::string HaplohSeq::EVENT_PREVALENCE = "eventPrevalence";
const std::string HaplohSeq::EVENT_MB = "eventMb";
const std::string HaplohSeq::GENOME_MB = "genomeMb";
const std::string HaplohSeq::INITIAL_PARAM_NORMAL = "intialParamNormal";
const std::string HaplohSeq::INITIAL_PARAM_EVENT = "initialParamEvent";
const std::string HaplohSeq::ITERATIONS = "iterations";
const std::string HaplohSeq::LAMBDA_0 = "lambda_0";
const std::string HaplohSeq::LAMBDA_1 = "lambda_1";
const std::string HaplohSeq::NUM_STATES = "numStates";
const std::string HaplohSeq::MEDIAN = "median";
const std::string HaplohSeq::NORMAL_STATE = "S0";
const std::string HaplohSeq::RANDOMSEED = "randomSeed";
const std::string HaplohSeq::REF = "ref";
const std::string HaplohSeq::ALT = "alt";
const std::string HaplohSeq::B_ALLELE = "B";
const std::string HaplohSeq::A_ALLELE = "A";
const std::string HaplohSeq::HAPLOTYPE1 = "haplotype1";
const std::string HaplohSeq::HAPLOTYPE2 = "haplotype2";

HaplohSeq::HaplohSeq() {}

HaplohSeq::~HaplohSeq() {}

// A helper function to simplify the main part.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
	return os;
}

void HaplohSeq::runBafHaplohseq(	std::string& obsType,
									unsigned int& numStates,
									double& epsilon,
									unsigned int& iterations,
									double& eventPrevalence,
									unsigned int& eventMb,
									unsigned int& genomeMb,
									double& initialParamNormal,
									double& initialParamEvent,
									double& endParamEvent,
									unsigned int& numParamEventStarts,
									std::string& bafFilename,
									unsigned int& numThreads,
									std::string& phasedAllelesFilename,
									std::string& destinationDir,
									std::map<std::string, bool>& estimateParams,
									std::string& outputPrefix,
									const unsigned int& randomSeed) {
	// load baf
	std::cerr << "loading " + bafFilename + "\n";
	std::map<std::string, std::vector<unsigned long> > chrInformativeIndices;
	std::map<std::string, std::vector<double> > chrInformativeRefFreqs;
	std::map<std::string, std::vector<char> > chrInformativeRef;
	std::map<std::string, std::vector<char> > chrInformativeAlt;
	std::map<std::string, std::vector<double> > chrInformativeLogRR;
	std::map<std::string, std::vector<unsigned int> > chrInformativeSwitch;
	std::map<std::string, std::vector<std::string> > chrInformativePos;

	std::vector<unsigned long> bafInformativeIndices;
	std::vector<std::string> orderedChrs;
	hlsutil::InputProcessor inputProc;
	std::map<std::string, std::vector<std::string> > bafStrFields;
	std::map<std::string, std::vector<unsigned int> > bafIntFields;
	std::map<std::string, std::vector<double> > bafDoubleFields;

	std::vector<std::string> bafChr;
	std::vector<unsigned int> bafPos;
	std::vector<std::string> bafRef;
	std::vector<std::string> bafAlt;
	std::vector<unsigned int> bafLogRR;
	std::vector<double> bafRefAlleleFreq;

	bafStrFields[inputProc.BAF_CHR] = bafChr;
	bafIntFields[inputProc.BAF_POS] = bafPos;
	bafStrFields[inputProc.BAF_REF] = bafRef;
	bafStrFields[inputProc.BAF_ALT] = bafAlt;
	bafDoubleFields[inputProc.BAF_RAF] = bafRefAlleleFreq;

	// load phased genotypes
	std::cerr << "loading " + phasedAllelesFilename + "\n";
	std::vector< std::vector<char> > phasedAlleles;
	std::map<std::string, std::vector< std::vector<char> > > chrPhasedAlleles;

	// read in chrPhasedAlleles
	//unsigned long numPhasedGenotypes = inputProc.readAlleles(phasedAllelesFilename, phasedAlleles, ' ', '\t', 0);
	unsigned long numPhasedGenotypes = inputProc.readAlleles(phasedAllelesFilename, phasedAlleles, '\0', ' ', 2);
	//inputProc.readAlleles(phasedAllelesFilename, phasedAlleles, '\0', ' ', 2);
	unsigned long numBafGenotypes = inputProc.readBafInformativeGenotypesChr(
											phasedAlleles,
											bafFilename,
											bafStrFields,
											bafIntFields,
											bafDoubleFields,
											orderedChrs,
											bafInformativeIndices,
											chrInformativeIndices,
											chrInformativeRefFreqs,
											chrInformativeRef,
											chrInformativeAlt,
											chrInformativeLogRR,
											chrInformativePos);

	// assemble informative phased alleles
	BOOST_FOREACH(std::string chr, orderedChrs) {
		std::cerr << "Num informative indices " << chr << ": " << chrInformativeIndices[chr].size() << std::endl;
		std::vector< std::vector<char> > informativePhasedAlleles;

		// assemble haplotypes of informative indices for each chromosome
		std::vector<char> haplotype0;
		std::vector<char> haplotype1;
		for (unsigned int i = 0; i < chrInformativeIndices[chr].size(); i++) {
			haplotype0.push_back(phasedAlleles[0][chrInformativeIndices[chr][i]]);
			haplotype1.push_back(phasedAlleles[1][chrInformativeIndices[chr][i]]);
		}
		informativePhasedAlleles.push_back(haplotype0);
		informativePhasedAlleles.push_back(haplotype1);
		chrPhasedAlleles[chr] = informativePhasedAlleles;
	}

	std::cerr << "Num baf genotypes in: " << numBafGenotypes << std::endl;
	std::cerr << "Num phased genotypes: " << numPhasedGenotypes << std::endl;
	if (numPhasedGenotypes != numBafGenotypes) {
		std::cerr << "ERROR: the number of phased genotypes does not equal the number of genotypes in the MDF file." << std::endl;
		exit(EXIT_FAILURE);
	}

	// get median ref frequency for informative sites
	haplohseq::FreqPhase fp(inputProc);
	std::vector<double> refFreqs;
	BOOST_FOREACH(std::string chr, orderedChrs) {
		refFreqs.insert(refFreqs.end(), chrInformativeRefFreqs[chr].begin(), chrInformativeRefFreqs[chr].end());
	}
	const double &median = fp.medianValue(refFreqs);
	std::cout << "Median ref allele frequency: " << std::setprecision(10) << median << std::endl;

	// Convert eventLength, genomeMb and eventPrevalence to lambda_0 and lambda_1
	unsigned int eventLengthMarkers = (unsigned int) (refFreqs.size()) * ((double) eventMb/ (double) genomeMb);
	std::cout << "Event length markers " << eventLengthMarkers << "\n";
	double lambda_0 = 1.0/( ((1.0-eventPrevalence)/eventPrevalence) * (double) eventLengthMarkers);
	double lambda_1 = 1.0/((double) eventLengthMarkers);

	// NOW THE MULTI-THREADED PART
	// for each chromosome, perform frequency-based phasing, enumerate switches and run HMM

	// map or primitives to allow for future flexibilitiy
	std::map<std::string, boost::any> primitives;
	primitives[MEDIAN] = median;
	primitives[RANDOMSEED] = randomSeed;
	primitives[EVENT_MB] = eventMb;
	primitives[EVENT_PREVALENCE] = eventPrevalence;
	primitives[GENOME_MB] = genomeMb;
	primitives[INITIAL_PARAM_NORMAL] = initialParamNormal;
	primitives[INITIAL_PARAM_EVENT] = initialParamEvent;
	primitives[NUM_STATES] = numStates;
	primitives[ITERATIONS] = iterations;
	primitives[EST_TRANS] = estimateParams[EST_TRANS];
	primitives[EST_ABERRANT_EMIT] =estimateParams[EST_ABERRANT_EMIT];
	primitives[EST_NORMAL_EMIT] = estimateParams[EST_NORMAL_EMIT];
	primitives[LAMBDA_0] = lambda_0;
	primitives[LAMBDA_1] = lambda_1;

	// keep track of the highest likelihood HMMs#################################################################################
	std::map<std::string, boost::shared_ptr<Hmm> > mleChrHmms;
	std::map<std::string, double> chrMLE;
	std::map<std::string, std::vector<double> > chrLikelihoodProfiles;
	std::vector<double> eventParamStartPoints;

	// run haplohseq on chromosomes at different start points
	hlsutil::StringUtil str;
	double interval = 0;
	if (numParamEventStarts > 1) {
		interval = (endParamEvent - initialParamEvent) / (numParamEventStarts - 1);
	}
	for (unsigned int i = 0; i < numParamEventStarts; i++) {

		// Create thread pool
		boost::asio::io_service ioService;
	//	boost::asio::io_service::work work(ioService);
		boost::shared_ptr<boost::asio::io_service::work> work(new boost::asio::io_service::work(ioService));
		boost::thread_group threadPool;

		// Spawn worker threads
		for (std::size_t t = 0; t < numThreads; t++) {
			threadPool.create_thread(boost::bind((unsigned long int (boost::asio::io_service::*)())&boost::asio::io_service::run, &ioService));
		}

		std::map<std::string, boost::shared_ptr<Hmm> > chrHmms;
		BOOST_FOREACH(std::string chr, orderedChrs) {

			boost::shared_ptr<std::map<std::string, std::vector<char> > > chrInformativeAllelesPtr(new std::map<std::string, std::vector<char> >);
			(*chrInformativeAllelesPtr)[REF] = chrInformativeRef[chr];
			(*chrInformativeAllelesPtr)[ALT] = chrInformativeAlt[chr];
			(*chrInformativeAllelesPtr)[HAPLOTYPE1] = chrPhasedAlleles[chr][0];
			(*chrInformativeAllelesPtr)[HAPLOTYPE2] = chrPhasedAlleles[chr][1];

			boost::shared_ptr<Hmm> hmmPtr(new Hmm(numStates, obsType));
			hmmPtr->setEpsilon(epsilon);
			chrHmms[chr] = hmmPtr;

			std::cout << "Posting job to thread pool for chromosome " << chr << "\n";
			// IMPORTANT: do not pass shared_ptr by reference because the shared_ptr expects a new pointer to the object
			ioService.post(boost::bind(&HaplohSeq::runBafHaplohseqChr, this,
													chr,
													boost::ref(chrHmms),
													boost::ref(fp),
													boost::ref(chrInformativeRefFreqs[chr]),
													boost::ref(chrInformativeLogRR[chr]),
													boost::ref(primitives),
													chrInformativeAllelesPtr
			));
		}

		work.reset();
		threadPool.join_all();	// wait for all threads to finish
		std::cout << "Joined all thread jobs.\n";

		std::string runOutputPrefix = outputPrefix + "_" + str.doubleToStr(boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]));
//		this->report(orderedChrs, chrHmms, vcfInformativeIndices, phasedAlleles, vcfStrFields, vcfIntFields, vcfDoubleFields, destinationDir, runOutputPrefix);
		BOOST_FOREACH(std::string chr, orderedChrs) {
			Hmm &chrHmm = *chrHmms[chr];
			double chrHmmLikelihood = exp(chrHmm.calcLogProbObs());
			if (!chrMLE.count(chr)) {
				chrMLE[chr] = chrHmmLikelihood;
				mleChrHmms[chr] = chrHmms[chr];
			} else {
				if (chrHmmLikelihood > chrMLE[chr]) {
					chrMLE[chr] = chrHmmLikelihood;
					mleChrHmms[chr] = chrHmms[chr];
				}
			}

			chrLikelihoodProfiles[chr].push_back(chrHmmLikelihood);
		}
		eventParamStartPoints.push_back(boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]));
		primitives[INITIAL_PARAM_EVENT] = boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]) + interval;
	}

	haplohseq::Reporter reporter;
	std::cout << "Generating reports...\n";
	reporter.reportBafBasedResults(orderedChrs, mleChrHmms, bafInformativeIndices, phasedAlleles, bafStrFields, bafIntFields, bafDoubleFields, destinationDir, outputPrefix);
	std::cout << "Reports generated.\n";
}

void HaplohSeq::runVcfHaplohseq(	std::string& obsType,
									unsigned int& numStates,
									double& epsilon,
									unsigned int& iterations,
									double& eventPrevalence,
									unsigned int& eventMb,
									unsigned int& genomeMb,
									double& initialParamNormal,
									double& initialParamEvent,
									double& endParamEvent,
									unsigned int& numParamEventStarts,
									std::string& vcfFilename,
									std::string& vcfSampleName,
									unsigned int& vcfDepthThreshold,
									unsigned int& numThreads,
									std::string& phasedAllelesFilename,
									std::string& destinationDir,
									std::map<std::string, bool>& estimateParams,
									std::string& outputPrefix,
									const unsigned int& randomSeed) {
	// load vcfs
	std::map<std::string, std::vector<unsigned long> > chrVcfInformativeIndices;
	std::map<std::string, std::vector<double> > chrInformativeRefFreqs;
	std::map<std::string, std::vector<char> > chrInformativeRef;
	std::map<std::string, std::vector<char> > chrInformativeAlt;
	std::map<std::string, std::vector<unsigned int> > chrInformativeDepth;
	std::map<std::string, std::vector<std::string> > chrInformativePos;

	std::vector<unsigned long> vcfInformativeIndices;
	std::vector<std::string> orderedChrs;
	hlsutil::InputProcessor inputProc;
	std::map<std::string, std::vector<std::string> > vcfStrFields;
	std::map<std::string, std::vector<unsigned int> > vcfIntFields;
	std::map<std::string, std::vector<double> > vcfDoubleFields;

	std::vector<std::string> vcfChr;
	std::vector<unsigned int> vcfPos;
	std::vector<std::string> vcfRef;
	std::vector<std::string> vcfAlt;
	std::vector<unsigned int> vcfRefDepth;
	std::vector<unsigned int> vcfAltDepth;
	std::vector<double> vcfRefFreq;

	vcfStrFields[inputProc.VCF_CHR] = vcfChr;
	vcfIntFields[inputProc.VCF_POS] = vcfPos;
	vcfStrFields[inputProc.VCF_REF] = vcfRef;
	vcfStrFields[inputProc.VCF_ALT] = vcfAlt;
//	vcfIntFields[inputProc.VCF_REF_DP] = vcfRefDepth;
//	vcfIntFields[inputProc.VCF_ALT_DP] = vcfAltDepth;
	vcfDoubleFields[inputProc.VCF_REF_FREQ] = vcfRefFreq;

	// load phased genotypes
	std::vector< std::vector<char> > phasedAlleles;
	std::map<std::string, std::vector< std::vector<char> > > chrPhasedAlleles;

	// read in chrPhasedAlleles
	std::cout << "Loading phased alleles " + phasedAllelesFilename + ".\n";
	unsigned long numPhasedGenotypes = inputProc.readAlleles(phasedAllelesFilename, phasedAlleles, '\0', ' ', 2);
	std::cout << "Loading VCF " + vcfFilename + ".\n";
	unsigned long numVcfGenotypes = inputProc.readVcfInformativeGenotypesChr(
											phasedAlleles,
											vcfFilename,
											vcfSampleName,
											vcfStrFields,
											vcfIntFields,
											vcfDoubleFields,
											orderedChrs,
											vcfDepthThreshold,
											vcfInformativeIndices,
											chrVcfInformativeIndices,
											chrInformativeRefFreqs,
											chrInformativeRef,
											chrInformativeAlt,
											chrInformativeDepth,
											chrInformativePos);
	std::cout << "Done loading input files.\n";

	// assemble informative phased alleles
	BOOST_FOREACH(std::string chr, orderedChrs) {
//		std::cout << "Num informative indices for chromosome " << chr << ": " << chrVcfInformativeIndices[chr].size() << std::endl;
		std::vector< std::vector<char> > informativePhasedAlleles;

		// assemble haplotypes of informative indices for each chromosome
		std::vector<char> haplotype0;
		std::vector<char> haplotype1;
		for (unsigned int i = 0; i < chrVcfInformativeIndices[chr].size(); i++) {
			haplotype0.push_back(phasedAlleles[0][chrVcfInformativeIndices[chr][i]]);
			haplotype1.push_back(phasedAlleles[1][chrVcfInformativeIndices[chr][i]]);
		}
		informativePhasedAlleles.push_back(haplotype0);
		informativePhasedAlleles.push_back(haplotype1);
		chrPhasedAlleles[chr] = informativePhasedAlleles;
	}

	std::cout << "Num VCF genotypes: " << numVcfGenotypes << std::endl;
	std::cout << "Num phased genotypes: " << numPhasedGenotypes << std::endl;
	if (numPhasedGenotypes != numVcfGenotypes) {
		std::cerr << "ERROR: the number of phased genotypes does not equal the number of genotypes in the VCF file." << std::endl;
		exit(EXIT_FAILURE);
	}

	// get median ref frequency for informative sites
	haplohseq::FreqPhase fp(inputProc);
	std::vector<double> refFreqs;
	BOOST_FOREACH(std::string chr, orderedChrs) {
		refFreqs.insert(refFreqs.end(), chrInformativeRefFreqs[chr].begin(), chrInformativeRefFreqs[chr].end());
	}
	const double &median = fp.medianValue(refFreqs);
	std::cout << "Median RAF: " << std::setprecision(10) << median << std::endl;

	// Convert eventLength, genomeMb and eventPrevalence to lambda_0 and lambda_1
	unsigned int eventLengthMarkers = (unsigned int) (refFreqs.size()) * ((double) eventMb/ (double) genomeMb);
//	std::cout << "Event length markers " << eventLengthMarkers << "\n";
	double lambda_0 = 1.0/( ((1.0-eventPrevalence)/eventPrevalence) * (double) eventLengthMarkers);
	double lambda_1 = 1.0/((double) eventLengthMarkers);

	// NOW THE MULTI-THREADED PART
	// for each chromosome, perform frequency-based phasing, enumerate switches and run HMM

	// map or primitives to allow for future flexibilitiy
	std::map<std::string, boost::any> primitives;
	primitives[MEDIAN] = median;
	primitives[RANDOMSEED] = randomSeed;
	primitives[EVENT_MB] = eventMb;
	primitives[EVENT_PREVALENCE] = eventPrevalence;
	primitives[GENOME_MB] = genomeMb;
	primitives[INITIAL_PARAM_NORMAL] = initialParamNormal;
	primitives[INITIAL_PARAM_EVENT] = initialParamEvent;
	primitives[NUM_STATES] = numStates;
	primitives[ITERATIONS] = iterations;
	primitives[EST_TRANS] = estimateParams[EST_TRANS];
	primitives[EST_ABERRANT_EMIT] =estimateParams[EST_ABERRANT_EMIT];
	primitives[EST_NORMAL_EMIT] = estimateParams[EST_NORMAL_EMIT];
	primitives[LAMBDA_0] = lambda_0;
	primitives[LAMBDA_1] = lambda_1;

	// keep track of the highest likelihood HMMs#################################################################################
	std::map<std::string, boost::shared_ptr<Hmm> > mleChrHmms;
	std::map<std::string, double> chrMLE;
	std::map<std::string, std::vector<double> > chrLikelihoodProfiles;
	std::vector<double> eventParamStartPoints;

	// run haplohseq on chromosomes at different start points
	hlsutil::StringUtil str;
	double interval = 0;
	if (numParamEventStarts > 1) {
		interval = (endParamEvent - initialParamEvent) / (numParamEventStarts - 1);
	}
	for (unsigned int i = 0; i < numParamEventStarts; i++) {

		// Create thread pool
		boost::asio::io_service ioService;
	//	boost::asio::io_service::work work(ioService);
		boost::shared_ptr<boost::asio::io_service::work> work(new boost::asio::io_service::work(ioService));
		boost::thread_group threadPool;

		// Spawn worker threads
		for (std::size_t t = 0; t < numThreads; t++) {
			threadPool.create_thread(boost::bind((unsigned long int (boost::asio::io_service::*)())&boost::asio::io_service::run, &ioService));
		}

		std::map<std::string, boost::shared_ptr<Hmm> > chrHmms;
		BOOST_FOREACH(std::string chr, orderedChrs) {

			boost::shared_ptr<std::map<std::string, std::vector<char> > > chrInformativeAllelesPtr(new std::map<std::string, std::vector<char> >);
			(*chrInformativeAllelesPtr)[REF] = chrInformativeRef[chr];
			(*chrInformativeAllelesPtr)[ALT] = chrInformativeAlt[chr];
			(*chrInformativeAllelesPtr)[HAPLOTYPE1] = chrPhasedAlleles[chr][0];
			(*chrInformativeAllelesPtr)[HAPLOTYPE2] = chrPhasedAlleles[chr][1];

			boost::shared_ptr<Hmm> hmmPtr(new Hmm(numStates, obsType));
			hmmPtr->setEpsilon(epsilon);
			chrHmms[chr] = hmmPtr;

//			std::cout << "Posting job to thread pool for chromosome " << chr << "\n";
			// IMPORTANT: do not pass shared_ptr by reference because the shared_ptr expects a new pointer to the object
			ioService.post(boost::bind(&HaplohSeq::runVcfHaplohseqChr, this,
													chr,
													boost::ref(chrHmms),
													boost::ref(fp),
													boost::ref(chrInformativeRefFreqs[chr]),
													boost::ref(chrInformativeDepth[chr]),
													boost::ref(primitives),
													chrInformativeAllelesPtr
			));
		}

		work.reset();
		threadPool.join_all();	// wait for all threads to finish
//		std::cout << "Joined all thread jobs.\n";

		std::string runOutputPrefix = outputPrefix + "_" + str.doubleToStr(boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]));
//		this->report(orderedChrs, chrHmms, vcfInformativeIndices, phasedAlleles, vcfStrFields, vcfIntFields, vcfDoubleFields, destinationDir, runOutputPrefix);
		BOOST_FOREACH(std::string chr, orderedChrs) {
			Hmm &chrHmm = *chrHmms[chr];
			double chrHmmLikelihood = exp(chrHmm.calcLogProbObs());
			if (!chrMLE.count(chr)) {
				chrMLE[chr] = chrHmmLikelihood;
				mleChrHmms[chr] = chrHmms[chr];
			} else {
				if (chrHmmLikelihood > chrMLE[chr]) {
					chrMLE[chr] = chrHmmLikelihood;
					mleChrHmms[chr] = chrHmms[chr];
				}
			}

			chrLikelihoodProfiles[chr].push_back(chrHmmLikelihood);
		}
		eventParamStartPoints.push_back(boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]));
		primitives[INITIAL_PARAM_EVENT] = boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]) + interval;
	}

	haplohseq::Reporter reporter;
	reporter.reportVcfBasedResults(orderedChrs, mleChrHmms, vcfInformativeIndices, phasedAlleles, vcfStrFields, vcfIntFields, vcfDoubleFields, destinationDir, outputPrefix);

	std::cout << "\n**** hapLOHseq output reports generated\n";
	std::cout << destinationDir + outputPrefix + ".posterior.dat\n";
	std::cout <<  destinationDir + outputPrefix + ".param.dat\n";
	std::cout <<  destinationDir + outputPrefix + ".likelihoods.dat\n";
	std::cout <<  destinationDir + outputPrefix + ".events.dat\n";
}

void HaplohSeq::runBafHaplohseqChr(
									std::string &chr,
									std::map<std::string, boost::shared_ptr<Hmm> > &chrHmms,
									haplohseq::FreqPhase &fp,
									const std::vector<double> &refFreqs,
									const std::vector<double> &logRRs,
									std::map<std::string, boost::any> &primitives,
									boost::shared_ptr<std::map<std::string, std::vector<char> > > &allelesPtr) {

	double initialParamNormal = boost::any_cast<double>(primitives[INITIAL_PARAM_NORMAL]);
	double initialParamEvent = boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]);
	double lambda_0 = boost::any_cast<double>(primitives[LAMBDA_0]);
	double lambda_1 = boost::any_cast<double>(primitives[LAMBDA_1]);
	double median = boost::any_cast<double>(primitives[MEDIAN]);
	unsigned int randomSeed = boost::any_cast<unsigned int>(primitives[RANDOMSEED]);
	unsigned int iterations = boost::any_cast<unsigned int>(primitives[ITERATIONS]);
	std::cout << "informative indexes size " << refFreqs.size() << "\n";
	std::cout << "median RAF " << median << "\n";
	std::cout << "randomSeed " << randomSeed << "\n";

	bool estimateTrans = boost::any_cast<bool>(primitives[EST_TRANS]);
	bool estimateAberrantEmit = boost::any_cast<bool>(primitives[EST_ABERRANT_EMIT]);
	bool estimateNormalEmit = boost::any_cast<bool>(primitives[EST_NORMAL_EMIT]);
	std::cout << "estimate trans " << estimateTrans << "\n";
	std::cout << "estimate aberrant emit " << estimateAberrantEmit << "\n";
	std::cout << "estimate normal emit " << estimateNormalEmit << "\n";

	std::map<std::string, std::vector<char> > alleles = *allelesPtr;
	std::vector<char> refs = alleles[REF];
	std::vector<char> alts = alleles[ALT];
	std::vector<char> haplotype1 = alleles[HAPLOTYPE1];
	std::vector<char> haplotype2 = alleles[HAPLOTYPE2];

	std::vector< std::vector<char> > haplotypes;
	haplotypes.push_back(haplotype1);
	haplotypes.push_back(haplotype2);

	Hmm &hmm = *chrHmms[chr];
	std::vector<std::vector<char> > freqPhased;

	if (hmm.getObsType() == Hmm::PHASE_CONCORDANCE) {
		std::cout << "Frequency-based phasing for chromosome " << chr << "\n";
		freqPhased = fp.phase(refFreqs,
								median,
								refs,
								alts,
								haplotypes,
								randomSeed);
		// enumerate switches using only informative sites
		// for VCFs hets have already been filtered for previously (could be optimized)
		std::cout << "Enumerating switches for chromosome " << chr << "\n";
		std::vector<int> switches = fp.enumerateHetSwitchesHaploh(haplotypes, freqPhased);
		double phaseConcordance = fp.meanValue(switches);

		std::cout << "Num phase concordance observations chromosome " << chr << ": "  << switches.size() << std::endl;
		std::cout << "Phase concordance chromosome " << chr << ": " << phaseConcordance << std::endl;

		// initialize HMM
		// using shared_ptr here so that we can keep the hmm in scope after leaving this block
		// if we didn't use shared_ptr, the hmm would be destroyed at the close of this block

		// just for ease of understanding, aliasing these events
		double alpha_0 = initialParamNormal;
		double alpha_1 = initialParamEvent;

		hmm.loadProbs(lambda_0, lambda_1, alpha_0, alpha_1, NORMAL_STATE);

		// add observed events
		std::vector<std::string> observedSequence;
		std::vector<bool> dummyRefAllele0ForNow;
		for (unsigned int i = 0; i < switches.size(); i++) {
			std::string switchStr = boost::lexical_cast<std::string>(switches[i]);
			observedSequence.push_back(switchStr);
			dummyRefAllele0ForNow.push_back(false);
		}
		// NOTE: the baumWelch here and addObservation should not take in coverages and refAllele vectors
		// for the PHASE_CONCORDANCE method
		if (estimateTrans || estimateAberrantEmit || estimateNormalEmit) {
			hmm.baumWelch(observedSequence, logRRs, dummyRefAllele0ForNow, iterations, estimateTrans, estimateAberrantEmit, estimateNormalEmit, NORMAL_STATE);
		} else {
			// add observations to hmm (because baum welch was not run and obs were not added to hmm)
			for (unsigned int i = 0; i < switches.size(); i++) {
				std::string switchStr = boost::lexical_cast<std::string>(switches[i]);
				hmm.addObservation(switchStr);
			}
		}
	}

//	if (hmm.getObsType() == Hmm::READ_COUNTS) {
//		throw an error because BAFs don't have reads

	std::cout << "Calculating posterior probabilities for chromosome " << chr << std::endl;
	hmm.forwardBackward(true);

	// run viterbi
//	std::vector<Transition*> path;
//	double jointProb = hmm.viterbi(path);
//	std::cout << "viterbi path probability " << chr << ": " << exp(jointProb - hmm.calcLogProbObs()) << std::endl;

//	hmm.saveProbs();
}

// This method is limited to 8 params because boost::bind can only handle 9 params (one of which is "this" instance).
void HaplohSeq::runVcfHaplohseqChr(
									std::string &chr,
									std::map<std::string, boost::shared_ptr<Hmm> > &chrHmms,
									haplohseq::FreqPhase &fp,
									const std::vector<double> &refFreqs,
									const std::vector<unsigned int> &coverages,
									std::map<std::string, boost::any> &primitives,
									boost::shared_ptr<std::map<std::string, std::vector<char> > > &allelesPtr) {

	double initialParamNormal = boost::any_cast<double>(primitives[INITIAL_PARAM_NORMAL]);
	double initialParamEvent = boost::any_cast<double>(primitives[INITIAL_PARAM_EVENT]);
	double lambda_0 = boost::any_cast<double>(primitives[LAMBDA_0]);
	double lambda_1 = boost::any_cast<double>(primitives[LAMBDA_1]);
	double median = boost::any_cast<double>(primitives[MEDIAN]);
	unsigned int randomSeed = boost::any_cast<unsigned int>(primitives[RANDOMSEED]);
	unsigned int iterations = boost::any_cast<unsigned int>(primitives[ITERATIONS]);
//	std::cout << "Informative indexes size " << refFreqs.size() << "\n";
//	std::cout << "Median RAF " << median << "\n";
//	std::cout << "Random seed " << randomSeed << "\n";

	bool estimateTrans = boost::any_cast<bool>(primitives[EST_TRANS]);
	bool estimateAberrantEmit = boost::any_cast<bool>(primitives[EST_ABERRANT_EMIT]);
	bool estimateNormalEmit = boost::any_cast<bool>(primitives[EST_NORMAL_EMIT]);
//	std::cout << "Estimate HMM trans " << estimateTrans << "\n";
//	std::cout << "Estimate HMM aberrant emit " << estimateAberrantEmit << "\n";
//	std::cout << "Estimate HMM normal emit " << estimateNormalEmit << "\n";

	std::map<std::string, std::vector<char> > alleles = *allelesPtr;
	std::vector<char> refs = alleles[REF];
	std::vector<char> alts = alleles[ALT];
	std::vector<char> haplotype1 = alleles[HAPLOTYPE1];
	std::vector<char> haplotype2 = alleles[HAPLOTYPE2];

	std::vector< std::vector<char> > haplotypes;
	haplotypes.push_back(haplotype1);
	haplotypes.push_back(haplotype2);

	Hmm &hmm = *chrHmms[chr];
	std::vector<std::vector<char> > freqPhased;

	std::cout << "\n**** Processing chromosome " << chr << ".\n";
	if (hmm.getObsType() == Hmm::PHASE_CONCORDANCE) {
		std::cout << "Frequency-based phasing for chromosome " << chr << ".\n";
		freqPhased = fp.phase(	refFreqs,
								median,
								refs,
								alts,
								haplotypes,
								randomSeed);
		// enumerate switches using only informative sites
		// for VCFs hets have already been filtered for previously (could be optimized)
		std::cout << "Enumerating switches for chromosome " << chr << ".\n";
		std::vector<int> switches = fp.enumerateHetSwitchesHaploh(haplotypes, freqPhased);
		double phaseConcordance = fp.meanValue(switches);

		std::cout << "Num phase concordance observations for chromosome " << chr << ": "  << switches.size() << "." << std::endl;
		std::cout << "Phase concordance for chromosome " << chr << ": " << phaseConcordance << "." << std::endl;

		// initialize HMM
		// using shared_ptr here so that we can keep the hmm in scope after leaving this block
		// if we didn't use shared_ptr, the hmm would be destroyed at the close of this block

		// just for ease of understanding, aliasing these events
		double alpha_0 = initialParamNormal;
		double alpha_1 = initialParamEvent;

		hmm.loadProbs(lambda_0, lambda_1, alpha_0, alpha_1, NORMAL_STATE);

		// add observed events
		std::vector<std::string> observedSequence;
		std::vector<bool> dummyRefAllele0ForNow;
		for (unsigned int i = 0; i < switches.size(); i++) {
			std::string switchStr = boost::lexical_cast<std::string>(switches[i]);
			observedSequence.push_back(switchStr);
			dummyRefAllele0ForNow.push_back(false);
		}
		// NOTE: the baumWelch here and addObservation should not take in coverages and refAllele vectors
		// for the PHASE_CONCORDANCE method
		if (estimateTrans || estimateAberrantEmit || estimateNormalEmit) {
			hmm.baumWelch(observedSequence, coverages, dummyRefAllele0ForNow, iterations, estimateTrans, estimateAberrantEmit, estimateNormalEmit, NORMAL_STATE);
		} else {
			// add observations to hmm (because baum welch was not run and obs were not added to hmm)
			for (unsigned int i = 0; i < switches.size(); i++) {
				std::string switchStr = boost::lexical_cast<std::string>(switches[i]);
				hmm.addObservation(switchStr);
			}
		}
	}

	if (hmm.getObsType() == Hmm::RAF_DEVIATION) {

		// TODO: create observations for all informative sites
		// (1) coverages
		// (2) ref depth = (unsigned int) (coverage * refFreq)
		// (3) is hap0 = R or A? if R then ref = true else ref = false
		// (4) then need to edit getEmitProb() to perform Christina/Paul logic

		// initialize HMM
		// using shared_ptr here so that we can keep the hmm in scope after leaving this block
		// if we didn't use shared_ptr, the hmm would be destroyed at the close of this block

		// just for ease of understanding, aliasing these events
		double alpha_0 = initialParamNormal;
		double alpha_1 = initialParamEvent;

		hmm.loadProbs(lambda_0, lambda_1, alpha_0, alpha_1, NORMAL_STATE);

		// add observed events
		std::vector<unsigned int> refCoverages;
		std::vector<std::string> refCoveragesStr;
		std::vector<bool> refSequence;
		for (unsigned int i = 0; i < refFreqs.size(); i++) {
			double refFreq = refFreqs[i];
			unsigned int cov = coverages[i];
			unsigned int refReadCount = (unsigned int) (cov * refFreq);
			refCoverages.push_back(refReadCount);
			refCoveragesStr.push_back(boost::lexical_cast<std::string>(refReadCount));
			if (haplotype1[i] == refs[i]) {
				refSequence.push_back(true);
			} else {
				refSequence.push_back(false);
			}
		}
//		std::cout << chr << "\t" << refCoverages << "\n";
//		std::cout << chr << "\t" << refSequence << "\n";
//		std::cout << chr << "\t" << coverages << "\n";

		std::string observation = "0";
		if (estimateTrans || estimateAberrantEmit || estimateNormalEmit) {
			// TODO: not implemented yet
			hmm.baumWelch(refCoveragesStr, coverages, refSequence, iterations, estimateTrans, estimateAberrantEmit, estimateNormalEmit, NORMAL_STATE);
		} else {
			// add observations to hmm (because baum welch was not run and obs were not added to hmm)
			for (unsigned int i = 1; i < refFreqs.size(); i++) {
				hmm.addObservation(observation, coverages[i-1], coverages[i], refCoverages[i-1], refCoverages[i], refSequence[i-1], refSequence[i]);
//				std::cout << i << "\t" << observation << "\t" <<  coverages[i-1] << "\t" <<  coverages[i] << "\t" <<  refCoverages[i-1] << "\t" <<  refCoverages[i] << "\t" <<  refSequence[i-1] << "\t" <<  refSequence[i] << std::endl;
			}
		}
	}

	std::cout << "Calculating event posterior probabilities for chromosome " << chr << "." << std::endl;
	hmm.forwardBackward(true);
	std::cout << "Successfully analyzed chromosome " << chr << ".\n";
	// run viterbi
//	std::vector<Transition*> path;
//	double jointProb = hmm.viterbi(path);
//	std::cout << "viterbi path probability " << chr << ": " << exp(jointProb - hmm.calcLogProbObs()) << std::endl;

//	hmm.saveProbs();
}



int main(int argc, char* argv[]) {
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");

	std::string configFilename = "",
				initProbFilename = "",
				transProbFilename = "",
				emitProbFilename = "",
				obsSeqFilename = "",
				phasedAllelesFilename = "",
				alleleFreqFilename = "",
				vcfFilename = "",
				vcfSampleName = "",
				markerDensity = "",
				destinationDir = "",
				obsType = "",
				outputPrefix = "";

	bool 		estimateTrans = false,
				estimateAberrantEmit = false,
				estimateNormalEmit = false;

	unsigned int 	randomSeed = 0,
					emIterations = 0,
					vcfDepthThreshold = 0,
					genomeMb = 0,
					eventMb = 0,
					numParamEventStarts = 0,
					numStates = 0,
					numThreads = 0;

	double		epsilon = 0.0,
				initialParamNormal = 0.0,
				initialParamEvent = 0.0,
				endParamEvent = 0.0,
				eventPrevalence = 0.0;

	try {
		// Declare a group of options that will be allowed only on command line
		po::options_description generic("Command-line parameters");
		generic.add_options()
	            		("version,v", "print version string")
	            		("help,h", "produce help message")
	            		("phased", po::value<std::string>(&phasedAllelesFilename)->default_value(""), "file of phased alleles")
	            		("freq,f", po::value<std::string>(&alleleFreqFilename)->default_value(""), "file of allele frequencies - only for microarray data (1-to-1 correspondence with markers in --phased)")
	            		("vcf,v", po::value<std::string>(&vcfFilename)->default_value(""), "vcf file of genotype calls for sample - inferred from sequencing data")
	            		("vcf_min_depth", po::value<unsigned int>(&vcfDepthThreshold)->default_value(10), "min depth for genotyped sites to be included in HMM")
	            		("vcf_sample_name", po::value<std::string>(&vcfSampleName)->default_value(""), "sample name to run haplohseq on - needed if using a multi-sample VCF")
	            		("obs,o", po::value<std::string>(&obsSeqFilename)->default_value(""), "file of observed sequences used to generate hidden most likely hidden sequences")
	            		("dest,d", po::value<std::string>(&destinationDir)->default_value(""), "destination directory for output files")
	            		("est_trans", po::value<bool>(&estimateTrans)->zero_tokens(), "estimate HMM event state transition params")
	            		("est_aberrant_emissions", po::value<bool>(&estimateAberrantEmit)->zero_tokens(), "estimate HMM aberrant event state emission params")
	            		("est_normal_emissions", po::value<bool>(&estimateNormalEmit)->zero_tokens(), "estimate HMM normal event state emission params")
	            		("marker_density,m", po::value<std::string>(&markerDensity)->default_value("uniform"), "strategy for modeling distances between genotype markers (uniform, genomic_pos)")
	            		("prefix,p", po::value<std::string>(&outputPrefix)->default_value(""), "prefix for output files")
	            		("config,c", po::value<std::string>(&configFilename)->default_value(""), "name of HMM configuration file")
						("threads,t", po::value<unsigned int>(&numThreads)->default_value(1), "number of threads to use for haplohseq")
						("genome_mb", po::value<unsigned int>(&genomeMb)->default_value(3156), "genome size in megabases (defaults to 3156MB - for a whole human genome (hg19) size estimate) use 40MB for an exome.")
						("event_prevalence", po::value<double>(&eventPrevalence)->default_value(0.1), "proportion of the genome expected to be aberrant (defaults to 0.2, or 20% - haplohseq is robust to incorrect estimates of event_prevalence)")
						("event_mb", po::value<unsigned int>(&eventMb)->default_value(20), "expected event size in megabases (defaults to 20MB - haplohseq is robust to incorrect estimates of event_length)")
						("initial_param_normal", po::value<double>(&initialParamNormal)->default_value(0.5), "probability of emitting a 1 in the normal state alpha_0 (defaults to 0.5)")
						("initial_param_event", po::value<double>(&initialParamEvent)->default_value(0.7), "value for the first iteration probability of emitting a 1 in the event state alpha_i (defaults to 0.6)")
						("end_param_event", po::value<double>(&endParamEvent)->default_value(0.7), "value for the last iteration probability of emitting a 1 in the event state alpha_i (defaults to 0.6)")
						("num_param_event_starts", po::value<unsigned int>(&numParamEventStarts)->default_value(1), "number of starts to use for the param_event starting at initial_param_event to end_param_event (defaults to 1 with initial_param_event == end_param_event or to 2 if initial_param_event != end_param_event)")
						("em_iterations", po::value<unsigned int>(&emIterations)->default_value(25), "number of EM iterations to fit HMM parameters (defaults to 10)")
						("phase_error_estimate", po::value<double>(&epsilon)->default_value(0.07), "estimate for phasing error (defaults to 0.07)")
						("obs_type", po::value<std::string>(&obsType)->default_value("phase_concordance"), "observed data for the HMM: phase_concordance or read_counts")
						("num_states", po::value<unsigned int>(&numStates)->default_value(2), "number of states for HMM (defaults to 2).");


		po::options_description cmdline_options;
		cmdline_options.add(generic);

		po::options_description visible("Allowed options");
		visible.add(generic);

		po::positional_options_description p;
//		p.add("input-file", -1);

		po::variables_map vm;
		store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
		notify(vm);

		if (vm.count("help")) {
			std::cout << "haplohseq (version 0.1.2) provides probabilities of allelic imbalance events in sequencing data." << std::endl << std::endl;
			std::cout << generic << std::endl;
			exit(EXIT_SUCCESS);
		}

		if (vm.count("version")) {
			std::cout << "haplohseq, version 0.1.2\n";
			exit(EXIT_SUCCESS);
		}

		// check if destination directory exists
		if (!boost::filesystem::exists(destinationDir)) {
			std::cerr<< "ERROR: Destination directory for output files (" << destinationDir << ") does not exist" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (alleleFreqFilename != "" && vcfFilename != "") {
			std::cerr << "ERROR: Only one sample file is allowed (either a --freq or a --vcf file).";
			exit(EXIT_FAILURE);
		}

		if (alleleFreqFilename == "" && vcfFilename == "") {
			std::cerr << "ERROR: At least one sample file must be provided (either a --freq or a --vcf file).";
			exit(EXIT_FAILURE);
		}

		// check if input files exist
		if (vcfFilename != "" && !boost::filesystem::exists(vcfFilename)) {
			std::cerr << "ERROR: VCF input file (" << vcfFilename << ") does not exist" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (alleleFreqFilename != "" && !boost::filesystem::exists(alleleFreqFilename)) {
			std::cerr << "ERROR: Allele frequency input file (" << alleleFreqFilename << ") does not exist" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!boost::filesystem::exists(phasedAllelesFilename)) {
			std::cerr << "ERROR: Phased haplotype input file (" << phasedAllelesFilename << ") does not exist" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (vm.count("input-file")) {
			std::cerr << "Input files are: "
					<< vm["input-file"].as< std::vector<std::string> >() << std::endl;
		}

		if (obsType != Hmm::PHASE_CONCORDANCE && obsType != Hmm::RAF_DEVIATION) {
			std::cerr << "ERROR: Valid options for obs_type are: phase_concordance or read_counts" << std::endl;
			exit(EXIT_FAILURE);
		}

		// adjust numParamEventStarts if initialParamEvent != endParamEvent
		if (numParamEventStarts == 1 && initialParamEvent != endParamEvent) {
			numParamEventStarts = 2;
		}

	}
	catch(std::exception& e) {
		std::cerr << e.what() << "\n";
		exit(EXIT_FAILURE);
	}

	std::map<std::string, bool> estimateParams;
	estimateParams[HaplohSeq::INIT] = false;
	estimateParams[HaplohSeq::EST_TRANS] = estimateTrans;
	estimateParams[HaplohSeq::EST_ABERRANT_EMIT] = estimateAberrantEmit;
	estimateParams[HaplohSeq::EST_NORMAL_EMIT] = estimateNormalEmit;

	HaplohSeq hs;
	if (alleleFreqFilename != "") {
		std::cout << "Running BAF version of haplohseq." << std::endl;
		hs.runBafHaplohseq(	obsType,
							numStates,
							epsilon,
							emIterations,
							eventPrevalence,
							eventMb,
							genomeMb,
							initialParamNormal,
							initialParamEvent,
							endParamEvent,
							numParamEventStarts,
							alleleFreqFilename,
							numThreads,
							phasedAllelesFilename,
							destinationDir,
							estimateParams,
							outputPrefix,
							randomSeed);
	}

	if (vcfFilename != "") {
		std::cout << "Running VCF version of haplohseq." << std::endl;
		hs.runVcfHaplohseq(	obsType,
							numStates,
							epsilon,
							emIterations,
							eventPrevalence,
							eventMb,
							genomeMb,
							initialParamNormal,
							initialParamEvent,
							endParamEvent,
							numParamEventStarts,
							vcfFilename,
							vcfSampleName,
							vcfDepthThreshold,
							numThreads,
							phasedAllelesFilename,
							destinationDir,
							estimateParams,
							outputPrefix,
							randomSeed);
	}

	return 0;
}


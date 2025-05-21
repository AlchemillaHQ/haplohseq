/*
 * Reporter.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#include "Reporter.h"

namespace haplohseq {

Reporter::Reporter() {
	// TODO Auto-generated constructor stub

}

Reporter::~Reporter() {
	// TODO Auto-generated destructor stub
}


void Reporter::writeToFile(std::string& filename, std::vector<double>& values, std::string& delimiter, std::string& header) {
	std::ofstream file;
	file.open(filename.c_str());
	if (header != "") {
		file << header << std::endl;
	}
	for (unsigned int i = 0; i < values.size()-1; i++) {
		file << values[i] << delimiter;
	}
	file << values[values.size()-1];
	file.close();
}

void Reporter::writeToFile(std::string& filename, std::vector<double>& values, std::string& delimiter) {
	std::string header = "";
	this->writeToFile(filename, values, delimiter, header);
}

void Reporter::writeToFile(std::string& filename, std::vector<unsigned long>& values, std::string& delimiter, std::string& header) {
	std::ofstream file;
	file.open(filename.c_str());
	if (header != "") {
		file << header << std::endl;
	}
	for (unsigned int i = 0; i < values.size()-1; i++) {
		file << values[i] << delimiter;
	}
	file << values[values.size()-1];
	file.close();
}

void Reporter::writeToFile(std::string& filename, std::vector<unsigned long>& values, std::string& delimiter) {
	std::string header = "";
	this->writeToFile(filename, values, delimiter, header);
}

void Reporter::writeToFile(std::string& filename, std::vector<std::string>& values, std::string& delimiter, std::string& header) {
	std::ofstream file;
	file.open(filename.c_str());
	if (header != "") {
		file << header << std::endl;
	}
	for (unsigned int i = 0; i < values.size()-1; i++) {
		file << values[i] << delimiter;
	}
	file << values[values.size()-1];
	file.close();
}

void Reporter::writeToFile(std::string& filename, std::vector<std::string>& values, std::string& delimiter) {
	std::string header = "";
	this->writeToFile(filename, values, delimiter, header);
}

void Reporter::reportBafBasedResults(	std::vector<std::string>& orderedChrs,
						std::map<std::string, boost::shared_ptr<Hmm> >& chrHmms,
						std::vector<unsigned long>& informativeIndices,
						std::vector< std::vector<char> > phasedAlleles,
						std::map<std::string, std::vector<std::string> >& bafStrFields,
						std::map<std::string, std::vector<unsigned int> >& bafIntFields,
						std::map<std::string, std::vector<double> >& bafDoubleFields,
						std::string& destinationDir,
						std::string& outputPrefix) {

	unsigned int numBafEntries = bafStrFields[hlsutil::InputProcessor::BAF_CHR].size();
	std::cout << "MDF entries size = " << numBafEntries << "\n";
	std::string SPACE = " ";
	std::string ENDL = "\n";
	std::string TAB = "\t";
	std::string nullString = "NA";

	if (destinationDir != " " and destinationDir[-1] != '/') {
		destinationDir += "/";
	}

	hlsutil::StringUtil str;
	hlsutil::VcfUtil vcf;

	// INITIALIZE FILES
	std::string posteriorFilename = destinationDir + outputPrefix + ".posterior.dat";
	std::string hmmParamFilename = destinationDir + outputPrefix + ".param.dat";
	std::string hmmLikelihoodsFilename = destinationDir + outputPrefix + ".likelihoods.dat";
	std::string eventFilename = destinationDir + outputPrefix + ".events.dat";

	std::ofstream posteriorFile(posteriorFilename.c_str(), std::ofstream::out);
	std::ofstream hmmParamFile(hmmParamFilename.c_str(), std::ofstream::out);
	std::ofstream hmmLikelihoodsFile(hmmLikelihoodsFilename.c_str(), std::ofstream::out);
	std::ofstream eventFile(eventFilename.c_str(), std::ofstream::out);

	std::string posteriorHeader = "CHR\tPOS\tREF\tALT\tINDEX\tPHASE_CONCORDANCE_SWITCH\tHAP\tGT\tRAF\tLRR";
	std::string eventHeader = "CHR\tBEGIN\tEND\tEVENT_STATE\tNUM_INFORMATIVE_MARKERS\tMEAN_POSTPROB\tPHASE_CONCORDANCE";

	// GET POSSIBLE STATES (possible states for each chr are the same)
	std::string firstChr = orderedChrs[0];
	Hmm &firstHmm = *chrHmms[firstChr];
	std::vector<unsigned long> states = firstHmm.getStates();

	// COLLECT POSTERIOR PROBABILITIES FOR ALL STATES AT ALL INFORMATIVE MARKERS
	std::vector< std::vector<double> > statePosteriors;
	for (unsigned int i=0; i < states.size(); i++) {
		unsigned long stateId = states[i];
		posteriorHeader += "\t" + (*chrHmms[firstChr]).getStr(stateId);
		std::vector<double> posteriors;

		for (unsigned int i = 0; i < orderedChrs.size(); i++) {
			std::string chr = orderedChrs[i];
			Hmm &chrHmm = *chrHmms[chr];
			std::vector<double> chrPosteriors = chrHmm.extractPosteriors(stateId);
			posteriors.insert(posteriors.end(), chrPosteriors.begin(), chrPosteriors.end());
		}
		statePosteriors.push_back(posteriors);
	}

	// COLLECT SWITCHES
	std::vector<std::string> observations;
	for (unsigned int i = 0; i < orderedChrs.size(); i++) {
		std::string chr = orderedChrs[i];
		Hmm &chrHmm = *chrHmms[chr];
		std::vector<std::string> chrObservations = chrHmm.getObservations();
		observations.insert(observations.end(), chrObservations.begin(), chrObservations.end());
	}

	// TRACK PUTATIVE EVENTS
	unsigned int eventState = 0;
	bool inEvent = false;							// in possible event if in stretch where posterior >= eventPosteriorThreshold
	bool validEvent = false;						// event is valid if any posterior >= validEventPosteriorThreshold
	unsigned int eventBegin = 0;
	unsigned int eventEnd = 0;
	unsigned int eventMarkerCount = 0;
	double validEventPosteriorThreshold = 0.9;
	double eventPosteriorThreshold = 0.5;
	std::vector<double> eventPosteriors;
	std::vector<int> eventSwitches;

	// OUTPUT POSTERIOR AND EVENT DAT FILE (TODO: refactor this into its own method)
	posteriorFile << posteriorHeader << std::endl;
	eventFile << eventHeader << std::endl;
	unsigned long posIndex = 0;
	unsigned long informativeIndex = 0;
	bool chrChange = false;
	std::string previousChr = "_initChr";
	for (unsigned int bafIndex = 0; bafIndex < numBafEntries; bafIndex++) {
		std::string gt = bafStrFields[hlsutil::InputProcessor::BAF_GT][bafIndex];
		std::string chr = bafStrFields[hlsutil::InputProcessor::BAF_CHR][bafIndex];
		unsigned int pos = bafIntFields[hlsutil::InputProcessor::BAF_POS][bafIndex];
		std::string ref = bafStrFields[hlsutil::InputProcessor::BAF_REF][bafIndex];
		std::string alt = bafStrFields[hlsutil::InputProcessor::BAF_ALT][bafIndex];
		double logRR = bafDoubleFields[hlsutil::InputProcessor::BAF_LOGRR][bafIndex];
		double refFreq = bafDoubleFields[hlsutil::InputProcessor::BAF_RAF][bafIndex];

		if (posIndex == informativeIndices[informativeIndex]) {
			int haplotype = 0;
			// Each marker is a SNP so each ref string is only 1 character long.
			if (phasedAlleles[0][posIndex] == ref[0]) {
				haplotype = 1;
			} else {
				if (phasedAlleles[1][posIndex] == ref[0]) {
					haplotype = 2;
				}
			}

			// There is 1 less posterior probability than het sites for phase_concordance.
			// In the HMM, the last
			// posterior is stored as 0 across states.  For a more accurate output,
			// we'll output a null string for this last state.
			if (chr != previousChr) {
				chrChange = true;
			} else {
				chrChange = false;
			}

			int phaseConcordanceSwitch = atoi(observations[informativeIndex].c_str());

			std::string gtNum = vcf.abToNumStr(gt, nullString);

			posteriorFile 	<< chr << "\t"
							<< pos << "\t"
							<< ref << "\t"
							<< alt << "\t"
							<< posIndex << "\t"
							<< phaseConcordanceSwitch << "\t"
							<< haplotype << "\t"
							<< gt << "\t"
							<< refFreq << "\t"
							<< logRR;

			for (unsigned int stateIndex = 0; stateIndex < statePosteriors.size(); stateIndex++) {
				if (chrChange) {
					posteriorFile << "\t" << nullString;

					// UPDATE EVENT DETAILS
					// if chrChange then check if was inEvent
					// if inEvent, output event and reset values
					if (inEvent && eventState == stateIndex && validEvent) {
						double posteriorSum = 0;
						for (double ep : eventPosteriors) {
							posteriorSum += ep;
						}
						double eventPosterior = posteriorSum / (double) eventMarkerCount;

						double switchSum = 0;
						for (int ss : eventSwitches) {
							switchSum += ss;
						}
						double eventPhaseConcordance = switchSum / (double) eventMarkerCount;

						eventFile 	<< previousChr << "\t"
									<< eventBegin << "\t"
									<< eventEnd << "\t"
									<< "S" << eventState << "\t"
									<< eventMarkerCount << "\t"
									<< eventPosterior << "\t"
									<< eventPhaseConcordance << "\n";

						// reset values for chr changes
						inEvent = false;
						validEvent = false;
						eventState = 0;
						eventMarkerCount = 0;
						eventPosteriors.clear();
						eventSwitches.clear();
					}

				}
				else {
					double posterior = statePosteriors[stateIndex][informativeIndex];
					posteriorFile << "\t" << posterior;

					if (stateIndex == 0 || (inEvent && eventState != stateIndex)) {
						continue;
					}

					// enter event (really event state shouldn't have to equal 0....if it doesn't, switching event states)
					if (!inEvent && posterior >= eventPosteriorThreshold) {
						inEvent = true;

						eventBegin = pos;
						eventEnd = pos;
						eventMarkerCount += 1;
						eventPosteriors.push_back(posterior);
						eventSwitches.push_back(phaseConcordanceSwitch);
						eventState = stateIndex;
						if (posterior >= validEventPosteriorThreshold) {
							validEvent = true;
						}
					} else {
						// extend event
						if (inEvent && eventState == stateIndex && posterior >= eventPosteriorThreshold) {
							eventEnd = pos;
							eventMarkerCount += 1;
							eventPosteriors.push_back(posterior);
							eventSwitches.push_back(phaseConcordanceSwitch);
							if (!validEvent) {
								if (posterior >= validEventPosteriorThreshold) {
									validEvent = true;
								}
							}
//							std::cout << "event " << stateIndex << "\t" << chr << "\t" << eventBegin << "\t" << eventEnd << "\t" << eventMarkerCount << "\n";
						} else {
							// leave event
							if (inEvent && eventState == stateIndex && posterior < eventPosteriorThreshold) {

								if (validEvent) {
									double posteriorSum = 0;
									for (double ep : eventPosteriors) {
										posteriorSum += ep;
									}
									double eventPosterior = posteriorSum / (double) eventMarkerCount;

									double switchSum = 0;
									for (int ss : eventSwitches) {
										switchSum += ss;
									}
									double eventPhaseConcordance = switchSum / (double) eventMarkerCount;

									eventFile 	<< chr << "\t"
												<< eventBegin << "\t"
												<< eventEnd << "\t"
												<< "S" << eventState << "\t"
												<< eventMarkerCount << "\t"
												<< eventPosterior << "\t"
												<< eventPhaseConcordance << "\n";
								}

								inEvent = false;
								validEvent = false;
								eventState = 0;
								eventMarkerCount = 0;
								eventPosteriors.clear();
								eventSwitches.clear();
							}
						}
					}
				}
			}
			posteriorFile << std::endl;
			informativeIndex += 1;

			if (chrChange) {
				previousChr = chr;
			}
		}
		posIndex += 1;
	}
	posteriorFile.close();

	// EXPORT CHROMOSOME BASED REPORTS
	std::cout << "Reporting data for chromosomes.\n";
	for (unsigned int i = 0; i < orderedChrs.size(); i++) {
		std::string chr = orderedChrs[i];
		std::cout << "Reporting data for chromosome " << chr << "\n";
		Hmm &chrHmm = *chrHmms[chr];
		chrHmm.calcLogProbObs();
		hmmLikelihoodsFile << "> " << chr << "\t(" << chrHmm.getLogProbObs() << ")\n";
		std::vector<double> likelihoodProfile = chrHmm.getLikelihoodProfile();
		for (unsigned int iteration = 0; iteration < likelihoodProfile.size(); iteration++) {
			hmmLikelihoodsFile << likelihoodProfile[iteration] << "\t";
		}
		hmmLikelihoodsFile << "\n";
		// allow writing of some type of header for each chr and matrix type
		hmmParamFile << "> " << chr << "\t(from\tto\tprob)\n";
		chrHmm.writeTrans(hmmParamFile);
		chrHmm.writeEmit(hmmParamFile);
	}
	hmmParamFile.close();
}

void Reporter::reportVcfBasedResults(	std::vector<std::string>& orderedChrs,
						std::map<std::string, boost::shared_ptr<Hmm> >& chrHmms,
						std::vector<unsigned long>& informativeIndices,
						std::vector< std::vector<char> > phasedAlleles,
						std::map<std::string, std::vector<std::string> >& vcfStrFields,
						std::map<std::string, std::vector<unsigned int> >& vcfIntFields,
						std::map<std::string, std::vector<double> >& vcfDoubleFields,
						std::string& destinationDir,
						std::string& outputPrefix) {

	unsigned int numVcfEntries = vcfStrFields[hlsutil::InputProcessor::VCF_CHR].size();
//	std::cout << "VCF entries size = " << numVcfEntries << "\n";
	std::string SPACE = " ";
	std::string ENDL = "\n";
	std::string TAB = "\t";
	std::string nullString = "NA";

	if (destinationDir != " " and destinationDir[-1] != '/') {
		destinationDir += "/";
	}

	hlsutil::StringUtil str;
	hlsutil::VcfUtil vcf;

	// INITIALIZE FILES
	std::string posteriorFilename = destinationDir + outputPrefix + ".posterior.dat";
	std::string hmmParamFilename = destinationDir + outputPrefix + ".param.dat";
	std::string hmmLikelihoodsFilename = destinationDir + outputPrefix + ".likelihoods.dat";
	std::string eventFilename = destinationDir + outputPrefix + ".events.dat";

	std::ofstream posteriorFile(posteriorFilename.c_str(), std::ofstream::out);
	std::ofstream hmmParamFile(hmmParamFilename.c_str(), std::ofstream::out);
	std::ofstream hmmLikelihoodsFile(hmmLikelihoodsFilename.c_str(), std::ofstream::out);
	std::ofstream eventFile(eventFilename.c_str(), std::ofstream::out);

	std::string posteriorHeader = "CHR\tPOS\tREF\tALT\tINDEX\tPHASE_CONCORDANCE_SWITCH\tHAP\tGT\tRAF\tDP";
	std::string eventHeader = "CHR\tBEGIN\tEND\tEVENT_STATE\tNUM_INFORMATIVE_MARKERS\tMEAN_POSTPROB\tPHASE_CONCORDANCE";

	// GET POSSIBLE STATES (possible states for each chr are the same)
	std::string firstChr = orderedChrs[0];
	Hmm &firstHmm = *chrHmms[firstChr];
	std::vector<unsigned long> states = firstHmm.getStates();

	// COLLECT POSTERIOR PROBABILITIES FOR ALL STATES AT ALL INFORMATIVE MARKERS
	std::vector< std::vector<double> > statePosteriors;
	for (unsigned int i=0; i < states.size(); i++) {
		unsigned long stateId = states[i];
		posteriorHeader += "\t" + (*chrHmms[firstChr]).getStr(stateId);
		std::vector<double> posteriors;

		for (unsigned int i = 0; i < orderedChrs.size(); i++) {
			std::string chr = orderedChrs[i];
			Hmm &chrHmm = *chrHmms[chr];
			std::vector<double> chrPosteriors = chrHmm.extractPosteriors(stateId);
			posteriors.insert(posteriors.end(), chrPosteriors.begin(), chrPosteriors.end());
		}
		statePosteriors.push_back(posteriors);
	}

	// COLLECT SWITCHES
	std::vector<std::string> observations;
	for (unsigned int i = 0; i < orderedChrs.size(); i++) {
		std::string chr = orderedChrs[i];
		Hmm &chrHmm = *chrHmms[chr];
		std::vector<std::string> chrObservations = chrHmm.getObservations();
		observations.insert(observations.end(), chrObservations.begin(), chrObservations.end());
	}

	// TRACK PUTATIVE EVENTS
	unsigned int eventState = 0;
	bool inEvent = false;							// in possible event if in stretch where posterior >= eventPosteriorThreshold
	bool validEvent = false;						// event is valid if any posterior >= validEventPosteriorThreshold
	unsigned int eventBegin = 0;
	unsigned int eventEnd = 0;
	unsigned int eventMarkerCount = 0;
	double validEventPosteriorThreshold = 0.9;
	double eventPosteriorThreshold = 0.5;
	std::vector<double> eventPosteriors;
	std::vector<int> eventSwitches;

	// OUTPUT POSTERIOR AND EVENT DAT FILE (TODO: refactor this into its own method)
	posteriorFile << posteriorHeader << std::endl;
	eventFile << eventHeader << std::endl;
	unsigned long posIndex = 0;
	unsigned long informativeIndex = 0;
	bool chrChange = false;
	std::string previousChr = "_initChr";
	for (unsigned int vcfIndex = 0; vcfIndex < numVcfEntries; vcfIndex++) {

		std::string gt = vcfStrFields[hlsutil::InputProcessor::VCF_GT][vcfIndex];
		std::string chr = vcfStrFields[hlsutil::InputProcessor::VCF_CHR][vcfIndex];
		unsigned int pos = vcfIntFields[hlsutil::InputProcessor::VCF_POS][vcfIndex];
		std::string ref = vcfStrFields[hlsutil::InputProcessor::VCF_REF][vcfIndex];
		std::string alt = vcfStrFields[hlsutil::InputProcessor::VCF_ALT][vcfIndex];
		unsigned int depth = vcfIntFields[hlsutil::InputProcessor::VCF_DP][vcfIndex];
		double refFreq = vcfDoubleFields[hlsutil::InputProcessor::VCF_REF_FREQ][vcfIndex];

		std::string gtNum = vcf.gtToNumStr(gt, nullString);

		if (posIndex == informativeIndices[informativeIndex]) {
			int haplotype = 0;
			// Each marker is a SNP so each ref string is only 1 character long.
			if (phasedAlleles[0][posIndex] == ref[0]) { //bAllele[0]) {
				haplotype = 1;
			} else {
				if (phasedAlleles[1][posIndex] == ref[0]) { //bAllele[0]) {
					haplotype = 2;
				}
			}

			// There is 1 less posterior probability than het sites for phase_concordance.
			// In the HMM, the last
			// posterior is stored as 0 across states.  For a more accurate output,
			// we'll output a null string for this last state.
			if (chr != previousChr) {
				chrChange = true;
			} else {
				chrChange = false;
			}

			int phaseConcordanceSwitch = atoi(observations[informativeIndex].c_str());

			std::string gtNum = vcf.abToNumStr(gt, nullString);

			posteriorFile 	<< chr << "\t"
								<< pos << "\t"
								<< ref << "\t"
								<< alt << "\t"
								<< posIndex << "\t"
								<< phaseConcordanceSwitch << "\t"
								<< haplotype << "\t"
								<< gt << "\t"
								<< refFreq << "\t"
								<< depth;

			for (unsigned int stateIndex = 0; stateIndex < statePosteriors.size(); stateIndex++) {
				if (chrChange) {
					posteriorFile << "\t" << nullString;

					// UPDATE EVENT DETAILS
					// if chrChange then check if was inEvent
					// if inEvent, output event and reset values
					if (inEvent && eventState == stateIndex && validEvent) {
						double posteriorSum = 0;
						for (double ep : eventPosteriors) {
							posteriorSum += ep;
						}
						double eventPosterior = posteriorSum / (double) eventMarkerCount;

						double switchSum = 0;
						for (int ss : eventSwitches) {
							switchSum += ss;
						}
						double eventPhaseConcordance = switchSum / (double) eventMarkerCount;

						eventFile 	<< previousChr << "\t"
									<< eventBegin << "\t"
									<< eventEnd << "\t"
									<< "S" << eventState << "\t"
									<< eventMarkerCount << "\t"
									<< eventPosterior << "\t"
									<< eventPhaseConcordance << "\n";

						// reset values for chr changes
						inEvent = false;
						validEvent = false;
						eventState = 0;
						eventMarkerCount = 0;
						eventPosteriors.clear();
						eventSwitches.clear();
					}

				}
				else {
					double posterior = statePosteriors[stateIndex][informativeIndex];
					posteriorFile << "\t" << posterior;

					if (stateIndex == 0 || (inEvent && eventState != stateIndex)) {
						continue;
					}

					// enter event (really event state shouldn't have to equal 0....if it doesn't, switching event states)
					if (!inEvent && posterior >= eventPosteriorThreshold) {
						inEvent = true;

						eventBegin = pos;
						eventEnd = pos;
						eventMarkerCount += 1;
						eventPosteriors.push_back(posterior);
						eventSwitches.push_back(phaseConcordanceSwitch);
						eventState = stateIndex;
						if (posterior >= validEventPosteriorThreshold) {
							validEvent = true;
						}
					} else {
						// extend event unless we are at the last informative index (then go to else)
						if (inEvent && eventState == stateIndex && posterior >= eventPosteriorThreshold && informativeIndex != informativeIndices.size() - 2) {
							eventEnd = pos;
							eventMarkerCount += 1;
							eventPosteriors.push_back(posterior);
							eventSwitches.push_back(phaseConcordanceSwitch);
							if (!validEvent) {
								if (posterior >= validEventPosteriorThreshold) {
									validEvent = true;
								}
							}
//							std::cout << "event " << stateIndex << "\t" << chr << "\t" << eventBegin << "\t" << eventEnd << "\t" << eventMarkerCount << "\n";
						} else {
							// leave event if in event and posterior drops below minimum posterior threshold for an event OR if we are at the last posterior in our dataset and in an event
							if (inEvent && eventState == stateIndex && (posterior < eventPosteriorThreshold || informativeIndex == informativeIndices.size() - 2)) {

								if (validEvent) {
									double posteriorSum = 0;
									for (double ep : eventPosteriors) {
										posteriorSum += ep;
									}
									double eventPosterior = posteriorSum / (double) eventMarkerCount;

									double switchSum = 0;
									for (int ss : eventSwitches) {
										switchSum += ss;
									}
									double eventPhaseConcordance = switchSum / (double) eventMarkerCount;

									eventFile 	<< chr << "\t"
												<< eventBegin << "\t"
												<< eventEnd << "\t"
												<< "S" << eventState << "\t"
												<< eventMarkerCount << "\t"
												<< eventPosterior << "\t"
												<< eventPhaseConcordance << "\n";
								}

								inEvent = false;
								validEvent = false;
								eventState = 0;
								eventMarkerCount = 0;
								eventPosteriors.clear();
								eventSwitches.clear();
							}
						}
					}
				}
			}
			posteriorFile << std::endl;
			informativeIndex += 1;

			if (chrChange) {
				previousChr = chr;
			}
		}
		posIndex += 1;
	}
	posteriorFile.close();

	// EXPORT CHROMOSOME BASED REPORTS
	for (unsigned int i = 0; i < orderedChrs.size(); i++) {
		std::string chr = orderedChrs[i];
		Hmm &chrHmm = *chrHmms[chr];
		chrHmm.calcLogProbObs();
		hmmLikelihoodsFile << "> " << chr << "\t(" << chrHmm.getLogProbObs() << ")\n";
		std::vector<double> likelihoodProfile = chrHmm.getLikelihoodProfile();
		for (unsigned int iteration = 0; iteration < likelihoodProfile.size(); iteration++) {
			hmmLikelihoodsFile << likelihoodProfile[iteration] << "\t";
		}
		hmmLikelihoodsFile << "\n";
		// allow writing of some type of header for each chr and matrix type
		hmmParamFile << "> " << chr << "\t(from\tto\tprob)\n";
		chrHmm.writeTrans(hmmParamFile);
		chrHmm.writeEmit(hmmParamFile);
	}
	hmmParamFile.close();
}

//void Reporter::reportVcfBasedResults(	std::vector<std::string>& orderedChrs,
//						std::map<std::string, boost::shared_ptr<Hmm> >& chrHmms,
//						std::vector<unsigned long>& informativeIndices,
//						std::vector< std::vector<char> > phasedAlleles,
//						std::map<std::string, std::vector<std::string> >& vcfStrFields,
//						std::map<std::string, std::vector<unsigned int> >& vcfIntFields,
//						std::map<std::string, std::vector<double> >& vcfDoubleFields,
//						std::string& destinationDir,
//						std::string& outputPrefix) {
//
//	unsigned int numVcfEntries = vcfStrFields[hlsutil::InputProcessor::VCF_CHR].size();
//	std::cout << "vcf entries size = " << numVcfEntries << "\n";
//	std::string SPACE = " ";
//	std::string ENDL = "\n";
//	std::string TAB = "\t";
//	std::string nullString = "NA";
//
//	if (destinationDir != " " and destinationDir[-1] != '/') {
//		destinationDir += "/";
//	}
//
//	hlsutil::StringUtil str;
//	hlsutil::VcfUtil vcf;
//
////	std::string datFilename = destinationDir + outputPrefix + ".dat";
//	std::string posteriorDatFilename = destinationDir + outputPrefix + ".posterior.dat";
//	std::string hmmParamFilename = destinationDir + outputPrefix + ".param.dat";
////	std::string hmmLikelihoodsFilename = destinationDir + outputPrefix + ".likelihoods.dat";
//	//std::string likelihoodsFilename = destinationDir + outputPrefix + ".likelihoods.dat";
//
////	std::ifstream infile(vcfFilename.c_str(), std::ifstream::in);
//
////	std::ofstream datFile(datFilename.c_str(), std::ofstream::out);
//	std::ofstream posteriorDatFile(posteriorDatFilename.c_str(), std::ofstream::out);
////	std::ofstream hmmLikelihoodsFile(hmmLikelihoodsFilename.c_str(), std::ofstream::out);
//
//	std::string datHeader = "chr\tpos\tref\talt\tgt\tcov\trefFreq\tindex";
////	datFile << datHeader << std::endl;
//
//	// get posteriors for each state
//	std::string posteriorDatHeader = datHeader + "\tref_hap";
//
//	// possible states for each chr are the same
//	std::string firstChr = orderedChrs[0];
//	Hmm &firstHmm = *chrHmms[firstChr];
//	std::vector<unsigned long> states = firstHmm.getStates();
////	std::cout << "states size " << firstHmm.getStates().size()  << "\n";
//
//	std::vector< std::vector<double> > statePosteriors;
//	for (unsigned int i=0; i < states.size(); i++) {
//		unsigned long stateId = states[i];
////		std::cout << "About to access state " << firstChr << "\t" << stateId << "\n";
//		posteriorDatHeader += "\t" + (*chrHmms[firstChr]).getStr(stateId);
//		std::vector<double> posteriors;
//
//		for (unsigned int i = 0; i < orderedChrs.size(); i++) {
//			std::string chr = orderedChrs[i];
//			Hmm &chrHmm = *chrHmms[chr];
////			std::cout << "accessing " << chr << "\n";
//			std::vector<double> chrPosteriors = chrHmm.extractPosteriors(stateId);
//			posteriors.insert(posteriors.end(), chrPosteriors.begin(), chrPosteriors.end());
//		}
////		std::cerr << "Num posteriors " << stateId << ": " << posteriors.size() << std::endl;
//		statePosteriors.push_back(posteriors);
//	}
//	posteriorDatFile << posteriorDatHeader << std::endl;
//
//	unsigned long posIndex = 0;
//	unsigned long informativeIndex = 0;
//	bool chrChange = false;
//	std::string previousChr = "_initChr";
//	for (unsigned int vcfIndex = 0; vcfIndex < numVcfEntries; vcfIndex++) {
//
//		std::string gt = vcfStrFields[hlsutil::InputProcessor::VCF_GT][vcfIndex];
//		std::string chr = vcfStrFields[hlsutil::InputProcessor::VCF_CHR][vcfIndex];
//		unsigned int pos = vcfIntFields[hlsutil::InputProcessor::VCF_POS][vcfIndex];
//		std::string ref = vcfStrFields[hlsutil::InputProcessor::VCF_REF][vcfIndex];
//		std::string alt = vcfStrFields[hlsutil::InputProcessor::VCF_ALT][vcfIndex];
//		unsigned int depth = vcfIntFields[hlsutil::InputProcessor::VCF_DP][vcfIndex];
//		double refFreq = vcfDoubleFields[hlsutil::InputProcessor::VCF_REF_FREQ][vcfIndex];
//
//		std::string gtNum = vcf.gtToNumStr(gt, nullString);
//
//
//		if (posIndex == informativeIndices[informativeIndex]) {
//			int haplotype = 0;
//			// Each marker is a SNP so each ref string is only 1 character long.
//			if (phasedAlleles[0][posIndex] == ref[0]) {
//				haplotype = 1;
//			} else {
//				if (phasedAlleles[1][posIndex] == ref[0]) {
//					haplotype = 2;
//				}
//			}
//			posteriorDatFile 	<< chr << "\t"
//								<< pos << "\t"
//								<< ref << "\t"
//								<< alt << "\t"
//								<< gtNum << "\t"
//								<< depth << "\t"
//								<< refFreq << "\t"
//								<< posIndex << "\t"
//								<< haplotype;
//
//			// There is 1 less posterior probability than het sites for phase_concordance.
//			// In the HMM, the last
//			// posterior is stored as 0 across states.  For a more accurate output,
//			// we'll output a null string for this last state.
//			if (chr != previousChr) {
//				chrChange = true;
//				previousChr = chr;
//			} else {
//				chrChange = false;
//			}
//
//			for (unsigned int stateIndex = 0; stateIndex < statePosteriors.size(); stateIndex++) {
//				if (chrChange == true) {
//					posteriorDatFile << "\t" << nullString;
//				}
//				else {
//					double posterior = statePosteriors[stateIndex][informativeIndex];
//					posteriorDatFile << "\t" << posterior;
//				}
//			}
//
//			posteriorDatFile << std::endl;
//
//			informativeIndex += 1;
//		}
//
//		posIndex += 1;
//
//	}
//
//	posteriorDatFile.close();
////	datFile.close();
//
////	std::ofstream hmmParamFile(hmmParamFilename.c_str(), std::ofstream::out);
////	for (unsigned int i = 0; i < orderedChrs.size(); i++) {
////		std::string chr = orderedChrs[i];
////		Hmm &chrHmm = *chrHmms[chr];
////		chrHmm.calcLogProbObs();
////		hmmLikelihoodsFile << chr << "\t" << chrHmm.getLogProbObs() << "\n";
////		// allow writing of some type of header for each chr and matrix type
////		chrHmm.writeTrans(hmmParamFile);
////		chrHmm.writeEmit(hmmParamFile);
////	}
//}

} /* namespace haplohseq */

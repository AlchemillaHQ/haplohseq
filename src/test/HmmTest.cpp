/*
 * HmmTest.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#define BOOST_TEST_MODULE HmmTest
#include <boost/test/included/unit_test.hpp>
#include <string>
#include "DataStructures.h"
#include "Hmm.h"

BOOST_AUTO_TEST_SUITE(HmmTest)

BOOST_AUTO_TEST_CASE(loadProbsAndObservationsTest) {

	// tolerance level of differences between
	// floating point numbers when testing for equality
	double tolerance = 0.00001;

	// initialize toy HMM
	std::string initState = "_initial";
	std::string state1 = "S1";
	std::string state2 = "S2";
	std::string emit1 = "A";
	std::string emit2 = "C";
	std::string emit3 = "G";
	std::string emit4 = "T";

	std::string obsType = "phase_concordance";
	Hmm hmm(2, obsType);
	hmm.setInitState(initState);
//	std::cout << hmm.getInitState() << std::endl;

	// load transition probabilities
	hmm.getTransProbs().load(initState, state1, 0.5, hmm.getStrIdMap());
	hmm.getTransProbs().load(initState, state2, 0.5, hmm.getStrIdMap());
	hmm.getTransProbs().load(state1, state1, 0.8, hmm.getStrIdMap());
	hmm.getTransProbs().load(state1, state2, 0.2, hmm.getStrIdMap());
	hmm.getTransProbs().load(state2, state1, 0.2, hmm.getStrIdMap());
	hmm.getTransProbs().load(state2, state2, 0.8, hmm.getStrIdMap());

	// check transition probs
	BOOST_CHECK_CLOSE_FRACTION(0.5, exp(hmm.getTransProbs().get(hmm.getId(state1), hmm.getId(initState))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.5, exp(hmm.getTransProbs().get(hmm.getId(state2), hmm.getId(initState))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.8, exp(hmm.getTransProbs().get(hmm.getId(state1), hmm.getId(state1))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.2, exp(hmm.getTransProbs().get(hmm.getId(state2), hmm.getId(state1))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.2, exp(hmm.getTransProbs().get(hmm.getId(state1), hmm.getId(state2))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.8, exp(hmm.getTransProbs().get(hmm.getId(state2), hmm.getId(state2))), tolerance);

	// load emission probabilities
	hmm.getEmitProbs().load(state1, emit1, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit2, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit3, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit4, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit1, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit2, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit3, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit4, 0.4, hmm.getStrIdMap());

	// check emission probs
	BOOST_CHECK_CLOSE_FRACTION(0.4, exp(hmm.getEmitProbs().get(hmm.getId(emit1), hmm.getId(state1))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.1, exp(hmm.getEmitProbs().get(hmm.getId(emit2), hmm.getId(state1))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.4, exp(hmm.getEmitProbs().get(hmm.getId(emit3), hmm.getId(state1))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.1, exp(hmm.getEmitProbs().get(hmm.getId(emit4), hmm.getId(state1))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.1, exp(hmm.getEmitProbs().get(hmm.getId(emit1), hmm.getId(state2))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.4, exp(hmm.getEmitProbs().get(hmm.getId(emit2), hmm.getId(state2))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.1, exp(hmm.getEmitProbs().get(hmm.getId(emit3), hmm.getId(state2))), tolerance);
	BOOST_CHECK_CLOSE_FRACTION(0.4, exp(hmm.getEmitProbs().get(hmm.getId(emit4), hmm.getId(state2))), tolerance);

//	hmm.saveProbs(); // prints transition and emission probabilities

	// add some observations - in this case each observation is a character
	std::string observations = "CGTCAG";
	for (int i = 0; i < observations.length(); i++) {
		hmm.addObservation(boost::lexical_cast<std::string>(observations[i]));
		BOOST_REQUIRE_EQUAL(observations[i], hmm.getObservation(i+1)[0]);
	}
}

BOOST_AUTO_TEST_CASE(viterbiTest) {

	// tolerance level of differences between
	// floating point numbers when testing for equality
	double tolerance = 0.00001;

	// initialize toy HMM
	std::string initState = "_initial";
	std::string state1 = "S1";
	std::string state2 = "S2";
	std::string emit1 = "A";
	std::string emit2 = "C";
	std::string emit3 = "G";
	std::string emit4 = "T";

	std::string obsType = "phase_concordance";
	Hmm hmm(2, obsType);
	hmm.setInitState(initState);
//	std::cout << hmm.getInitState() << std::endl;

	// load transition probabilities
	hmm.getTransProbs().load(initState, state1, 0.5, hmm.getStrIdMap());
	hmm.getTransProbs().load(initState, state2, 0.5, hmm.getStrIdMap());
	hmm.getTransProbs().load(state1, state1, 0.8, hmm.getStrIdMap());
	hmm.getTransProbs().load(state1, state2, 0.2, hmm.getStrIdMap());
	hmm.getTransProbs().load(state2, state1, 0.2, hmm.getStrIdMap());
	hmm.getTransProbs().load(state2, state2, 0.8, hmm.getStrIdMap());

	// load emission probabilities
	hmm.getEmitProbs().load(state1, emit1, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit2, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit3, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit4, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit1, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit2, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit3, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit4, 0.4, hmm.getStrIdMap());

	//	hmm.saveProbs(); // prints transition and emission probabilities

	// add some observations - in this case each observation is a character
	std::string observations = "CGTCAG";
	for (int i = 0; i < observations.length(); i++) {
		hmm.addObservation(boost::lexical_cast<std::string>(observations[i]));
	}

	// TEST VITERBI
	std::string predictedPath;
	std::vector<Transition*> path;
	double jointProb = hmm.viterbi(path);
	double viterbiProb = exp(jointProb - hmm.calcLogProbObs());
	for (unsigned int i = 0; i < path.size(); i++) {
		Transition* trans = path[i];
		if (trans==0) continue;
		predictedPath.append(boost::lexical_cast<std::string>(trans->to->getState()));
	}

	BOOST_REQUIRE_EQUAL(predictedPath, "222211");
	BOOST_CHECK_CLOSE_FRACTION(viterbiProb, 0.269585, tolerance);
}

BOOST_AUTO_TEST_CASE(forwardBackwardTest) {

	// tolerance level of differences between
	// floating point numbers when testing for equality
	double tolerance = 0.00001;

	// initialize toy HMM
	std::string initState = "_initial";
	std::string state1 = "S1";
	std::string state2 = "S2";
	std::string emit1 = "A";
	std::string emit2 = "C";
	std::string emit3 = "G";
	std::string emit4 = "T";

	std::string obsType = "phase_concordance";
	Hmm hmm(2, obsType);
	hmm.setInitState(initState);
//	std::cout << hmm.getInitState() << std::endl;

	// load transition probabilities
	hmm.getTransProbs().load(initState, state1, 0.5, hmm.getStrIdMap());
	hmm.getTransProbs().load(initState, state2, 0.5, hmm.getStrIdMap());
	hmm.getTransProbs().load(state1, state1, 0.8, hmm.getStrIdMap());
	hmm.getTransProbs().load(state1, state2, 0.2, hmm.getStrIdMap());
	hmm.getTransProbs().load(state2, state1, 0.2, hmm.getStrIdMap());
	hmm.getTransProbs().load(state2, state2, 0.8, hmm.getStrIdMap());

	// load emission probabilities
	hmm.getEmitProbs().load(state1, emit1, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit2, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit3, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state1, emit4, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit1, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit2, 0.4, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit3, 0.1, hmm.getStrIdMap());
	hmm.getEmitProbs().load(state2, emit4, 0.4, hmm.getStrIdMap());

	//	hmm.saveProbs(); // prints transition and emission probabilities

	// add some observations - in this case each observation is a character
	std::string observations = "CGTCAG";
	for (int i = 0; i < observations.length(); i++) {
		hmm.addObservation(boost::lexical_cast<std::string>(observations[i]));
	}
//	std::cout << hmm.getObservation(3L);

	// TEST FORWARD-BACKWARD
	hmm.forwardBackward();

//	hmm.forward();
//	hmm.backward();
	// expected values for alphas, betas and posteriors to test for

	std::map<unsigned int, std::map<std::string, double> > expectedAlphas;
	expectedAlphas[1]["S1"] = 0.2;
	expectedAlphas[1]["S2"] = 0.8;
	expectedAlphas[2]["S1"] = 0.653061;
	expectedAlphas[2]["S2"] = 0.346939;
	expectedAlphas[3]["S1"] = 0.266055;
	expectedAlphas[3]["S2"] = 0.733945;
	expectedAlphas[4]["S1"] = 0.123116;
	expectedAlphas[4]["S2"] = 0.876884;
	expectedAlphas[5]["S1"] = 0.601379;
	expectedAlphas[5]["S2"] = 0.398621;
	expectedAlphas[6]["S1"] = 0.836281;
	expectedAlphas[6]["S2"] = 0.163719;

	std::map<unsigned int, std::map<std::string, double> > expectedBetas;
	expectedBetas[1]["S1"] = 5.12546;
	expectedBetas[1]["S2"] = 3.71863;
	expectedBetas[2]["S1"] = 3.4965;
	expectedBetas[2]["S2"] = 8.1242;
	expectedBetas[3]["S1"] = 3.8297;
	expectedBetas[3]["S2"] = 4.73673;
	expectedBetas[4]["S1"] = 7.84657;
	expectedBetas[4]["S2"] = 2.80234;
	expectedBetas[5]["S1"] = 6.95804;
	expectedBetas[5]["S2"] = 3.27437;
	expectedBetas[6]["S1"] = 3.72789;
	expectedBetas[6]["S2"] = 3.72789;

	std::map<unsigned int, double > expectedScaleFactors;
	expectedScaleFactors[1] = 4;
	expectedScaleFactors[2] = 5.10204;
	expectedScaleFactors[3] = 4.49541;
	expectedScaleFactors[4] = 3.42337;
	expectedScaleFactors[5] = 5.48966;
	expectedScaleFactors[6] = 3.72789;

	std::map<unsigned int, std::map<std::string, double> > expectedPosteriors;
	expectedPosteriors[1]["S1"] = 0.256273;
	expectedPosteriors[1]["S2"] = 0.743727;
	expectedPosteriors[2]["S1"] = 0.447552;
	expectedPosteriors[2]["S2"] = 0.552448;
	expectedPosteriors[3]["S1"] = 0.226656;
	expectedPosteriors[3]["S2"] = 0.773344;
	expectedPosteriors[4]["S1"] = 0.282188;
	expectedPosteriors[4]["S2"] = 0.717812;
	expectedPosteriors[5]["S1"] = 0.762238;
	expectedPosteriors[5]["S2"] = 0.237762;
	expectedPosteriors[6]["S1"] = 0.836281;
	expectedPosteriors[6]["S2"] = 0.163719;

	for (unsigned int i = 1; i < hmm.getLatentStatesZ().size(); i++) {
		LatentStateZ* zi = hmm.getLatentStatesZ()[i];
		double scaleFactor = zi->getScaleFactor();
		BOOST_CHECK_CLOSE_FRACTION(scaleFactor, expectedScaleFactors[i], tolerance);
		for (LatentStateZ::iterator it = zi->begin(); it != zi->end(); it++) {
			PossibleState *ps = *it;
			std::string obs = hmm.getStr((*it)->getInTrans()[0]->obs);
			std::string state = hmm.getStr(ps->getState());
			double alpha = exp(ps->getLogAlpha());
			double beta = exp(ps->getLogBeta());
			double posterior =	exp(ps->getLogPosterior());

			BOOST_CHECK_CLOSE_FRACTION(alpha, expectedAlphas[i][state], tolerance);
			BOOST_CHECK_CLOSE_FRACTION(beta, expectedBetas[i][state], tolerance);
			BOOST_CHECK_CLOSE_FRACTION(posterior, expectedPosteriors[i][state], tolerance);

			std::cout << i << "\t" 	<< obs << "\t"
									<< state << "\t"
									<< alpha << "\t"
									<< beta << "\t"
									<< scaleFactor << "\t"
									<< posterior << "\t"
									<< std::endl;
		}
	}

}

//BOOST_AUTO_TEST_CASE(estimateProbsTest) {
//
//	// tolerance level of differences between
//	// floating point numbers when testing for equality
//	double tolerance = 0.00001;
//
//	// initialize toy HMM
//	std::string initState = "_initial";
//	std::string state1 = "S1";
//	std::string state2 = "S2";
//	std::string emit1 = "A";
//	std::string emit2 = "C";
//	std::string emit3 = "G";
//	std::string emit4 = "T";
//
//	Hmm hmm;
//	hmm.setInitState(initState);
//	std::cout << hmm.getInitState() << std::endl;
//
//	// load transition probabilities
//	hmm.getTransProbs().load(initState, state1, 0.5, hmm.getStrIdMap());
//	hmm.getTransProbs().load(initState, state2, 0.5, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state1, state1, 0.8, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state1, state2, 0.2, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state2, state1, 0.2, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state2, state2, 0.8, hmm.getStrIdMap());
//
//	// load emission probabilities
//	hmm.getEmitProbs().load(state1, emit1, 0.4, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state1, emit2, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state1, emit3, 0.4, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state1, emit4, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit1, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit2, 0.4, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit3, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit4, 0.4, hmm.getStrIdMap());
//
//	hmm.saveProbs(); // prints transition and emission probabilities
//
//	// add some observations - in this case each observation is a character
//	std::string observations = "CGTCAG";
//	std::vector< std::vector<std::string>* > observedSequences;
//	std::vector<std::string> observedSequence;
//	for (unsigned int i = 0; i < observations.size(); i++) {
//		std::string switchStr = boost::lexical_cast<std::string>(observations[i]);
//		observedSequence.push_back(switchStr);
//	}
//	observedSequences.push_back(&observedSequence);
//
//	// TEST PARAMETER ESTIMATION
//	for (unsigned int i=0; i<observedSequences.size(); i++) {
//		std::vector<std::string>& seq = *observedSequences[i];
//		for (unsigned int j=0; j<seq.size(); j++) {
//			hmm.addObservation(seq[j]);
//		}
//	}
//
//	hmm.estimateProbs();
//	hmm.saveProbs();
//
//}

//BOOST_AUTO_TEST_CASE(baumWelchTest) {
//
//	// tolerance level of differences between
//	// floating point numbers when testing for equality
//	double tolerance = 0.00001;
//
//	// initialize toy HMM
//	std::string initState = "_initial";
//	std::string state1 = "S1";
//	std::string state2 = "S2";
//	std::string emit1 = "A";
//	std::string emit2 = "C";
//	std::string emit3 = "G";
//	std::string emit4 = "T";
//
//	Hmm hmm;
//	hmm.setInitState(initState);
//	std::cout << hmm.getInitState() << std::endl;
//
//	// load transition probabilities
//	hmm.getTransProbs().load(initState, state1, 0.5, hmm.getStrIdMap());
//	hmm.getTransProbs().load(initState, state2, 0.5, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state1, state1, 0.8, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state1, state2, 0.2, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state2, state1, 0.2, hmm.getStrIdMap());
//	hmm.getTransProbs().load(state2, state2, 0.8, hmm.getStrIdMap());
//
//	// load emission probabilities
//	hmm.getEmitProbs().load(state1, emit1, 0.4, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state1, emit2, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state1, emit3, 0.4, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state1, emit4, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit1, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit2, 0.4, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit3, 0.1, hmm.getStrIdMap());
//	hmm.getEmitProbs().load(state2, emit4, 0.4, hmm.getStrIdMap());
//
//	hmm.saveProbs(); // prints transition and emission probabilities
//
//	// add some observations - in this case each observation is a character
//	std::string observations = "CGTCAG";
//	std::vector< std::vector<std::string>* > observedSequences;
//	std::vector<std::string> observedSequence;
//	for (unsigned int i = 0; i < observations.size(); i++) {
//		std::string switchStr = boost::lexical_cast<std::string>(observations[i]);
//		observedSequence.push_back(switchStr);
//	}
//
//	// TEST BAUM-WELCH
//	int iterations = 2;
//	observedSequences.push_back(&observedSequence);
//	hmm.baumWelch(observedSequences, iterations);
//
//}

BOOST_AUTO_TEST_SUITE_END()


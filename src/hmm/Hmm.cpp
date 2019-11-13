#include "Hmm.h"

const std::string Hmm::PHASE_CONCORDANCE	=	"phase_concordance";
const std::string Hmm::RAF_DEVIATION		=	"raf_deviation";
const std::string Hmm::BAF_DEVIATION		=	"baf_deviation";

Hmm::Hmm(unsigned int numStates, std::string& obsType) {
	this->minLogProb = log(0.000001);
	this->initState = -1;
	this->logProbObs = 0;
	this->backwardHasBeenRun = false;
	this->forwardHasBeenRun = false;
	this->forwardBackwardHasBeenRun = false;
	this->numStates = numStates;
	this->obsType = obsType;
	this->epsilon = 1.0;
}

double Hmm::getTransProb(Transition* trans) {
	return this->transProbs.get(trans->to->getState(), trans->from->getState());
}

double Hmm::getEmitProb(Transition* trans) {
	if (this->obsType == PHASE_CONCORDANCE) {
		double emissionProb = this->emitProbs.get(trans->obs, trans->to->getState());
		return emissionProb;	// this is alpha_l for the PHASE_CONCORDANCE method
	}

	if (this->obsType == RAF_DEVIATION) {
		// otherwise, obsType == READ_COUNTS
		// P(s|Lm,hg,epsilon) =
		double epsilon = this->getEpsilon();
		unsigned int n1 = trans->n1;
		unsigned int n2 = trans->n2;
		unsigned int r1 = trans->r1;
		unsigned int r2 = trans->r2;
		bool ref1 = trans->ref1;
		bool ref2 = trans->ref2;

		double emissionProb = 0.0;

		double p_norm = 0.5;	// for now.  at least make this configurable

		std::string pParam = "0";
		double stateParam = exp(this->emitProbs.get(strIdMap.getId(pParam), trans->to->getState()));
		double d = stateParam - p_norm; // abs seems to convert d to an int
		double p_rl = p_norm + d;
		double p_al = p_norm - d;

	//	std::cout << strIdMap.getStr(trans->to->getState()) << "\t" << stateParam << "\t" << exp(stateParam) << "\n";
		if (n1 > 100) {
			int factor = n1 / 100;
			n1 = n1 / factor;
			r1 = r1 / factor;
		}
		if (n2 > 100) {
			int factor = n2 / 100;
			n2 = n2 / factor;
			r2 = r2 / factor;
		}

	//	std::cout << n1 << "\t" << r1 << "\t" << n2 << "\t" << r2 << "\n";
		double c1 = boost::math::binomial_coefficient<double>(n1, r1);
		double c2 = boost::math::binomial_coefficient<double>(n2, r2);

	//	return 0.1;

	//// LOG VERSION INCORRECT
	//	// hg = RR/AA
	//	if (ref1 == ref2) {
	//		emissionProb =
	//				(log(0.5) + (r1*log(p_rl) + (n1-r1)*log(1-p_rl) + r2*log(p_rl) + (n2-r2)*log(1-p_rl)) *
	//						    (r1*log(p_al) + (n1-r1)*log(1-p_al) + r2*log(p_al) + (n2-r2)*log(1-p_al)) +
	//					        log(1-epsilon)) *
	//				(log(0.5) + (r1*log(p_rl) + (n1-r1)*log(1-p_rl) + r2*log(p_al) + (n2-r2)*log(1-p_al)) *
	//					        (r1*log(p_al) + (n1-r1)*log(1-p_al) + r2*log(p_rl) + (n2-r2)*log(1-p_rl)) +
	//					        log(epsilon));
	//	}
	//
	//	// hg = RA/AR
	//	if (ref1 != ref2) {
	//		emissionProb =
	//				(log(0.5) + (r1*log(p_rl) + (n1-r1)*log(1-p_rl) + r2*log(p_rl) + (n2-r2)*log(1-p_rl)) *
	//						    (r1*log(p_al) + (n1-r1)*log(1-p_al) + r2*log(p_al) + (n2-r2)*log(1-p_al)) +
	//					        log(epsilon)) *
	//				(log(0.5) + (r1*log(p_rl) + (n1-r1)*log(1-p_rl) + r2*log(p_al) + (n2-r2)*log(1-p_al)) *
	//					        (r1*log(p_al) + (n1-r1)*log(1-p_al) + r2*log(p_rl) + (n2-r2)*log(1-p_rl)) +
	//					        log(1-epsilon));
	//	}
	//
	//	double testA = (r1*log(p_rl) + (n1-r1)*log(1-p_rl) + r2*log(p_rl) + (n2-r2)*log(1-p_rl));
	//	double testB = (r1*log(p_al) + (n1-r1)*log(1-p_al) + r2*log(p_al) + (n2-r2)*log(1-p_al));
	//	std::cout << testA << std::endl;
	//	std::cout << testB << std::endl;
	//	std::cout << log(0.5) << "\t" << log(epsilon) << std::endl;
	//	std::cout << (log(0.5) + testA * testB + log(epsilon)) << std::endl;
	//	std::cout << testA * testB * testA * testB << std::endl;
	//
	//	std::cout << strIdMap.getStr(trans->to->getState()) << "\t" << stateParam << "\t" << d << "\t" << p_rl << "\t" << p_al << "\t" << n1 << "\t" << n2 << "\t" << r1 << "\t" << r2 << "\t" << emissionProb << std::endl;
	//	std::cout << log(c1) + log(c2) + emissionProb << std::endl;
	//	return emissionProb;
	//	return log(c1) + log(c2) + emissionProb;

	//  NON-LOG VERSION
	//	n1 = 100;
	//	n2 = 100;
	//	r1 = 50;
	//	r2 = 50;
		// hg = RR/AA
		if (ref1 == ref2) {
			emissionProb =
					0.5 * ((pow(p_rl,r1) * pow(1 - p_rl, n1-r1)) * (pow(p_rl,r2) * pow(1 - p_rl, n2-r2)) +
						  (pow(p_al,r1) * pow(1 - p_al, n1-r1)) * (pow(p_al,r2) * pow(1 - p_al, n2-r2))) *
						  (1 - epsilon) +
					0.5 * ((pow(p_rl,r1) * pow(1 - p_rl, n1-r1)) * (pow(p_al,r2) * pow(1 - p_al, n2-r2)) +
						  (pow(p_al,r1) * pow(1 - p_al, n1-r1)) * (pow(p_rl,r2) * pow(1 - p_rl, n2-r2))) *
						  epsilon;
		}

		// hg = RA/AR
		if (ref1 != ref2) {
			emissionProb =
					0.5 * ((pow(p_rl,r1) * pow(1 - p_rl, n1-r1)) * (pow(p_rl,r2) * pow(1 - p_rl, n2-r2)) +
						  (pow(p_al,r1) * pow(1 - p_al, n1-r1)) * (pow(p_al,r2) * pow(1 - p_al, n2-r2))) *
						  epsilon +
					0.5 * ((pow(p_rl,r1) * pow(1 - p_rl, n1-r1)) * (pow(p_al,r2) * pow(1 - p_al, n2-r2)) +
						  (pow(p_al,r1) * pow(1 - p_al, n1-r1)) * (pow(p_rl,r2) * pow(1 - p_rl, n2-r2))) *
						  (1- epsilon);
		}
	//	std::cout << emissionProb << std::endl;
	//	std::cout << p_rl << "\t" << p_al << "\t" << d << "\t" << n1 << "\t" << n2 << "\t" << r1 << "\t" << r2 << "\t" << c1 << "\t" << c2 << "\t" << log(emissionProb) << std::endl;
		return log(c1) + log(c2) + log(emissionProb);
	}

	else {
		// BAF_DEVIATION: to implement
		return 0.0;
	}

}

void Hmm::forward(bool force) {
	// don't run if forward() has been run before and force is false
	if (!force && this->hasForwardBeenRun()) {
		return;
	}

	// compute forward probabilities at time 0
	LatentStateZ* t0 = this->latentStatesZ[0];
	PossibleState* init = (*t0)[0];
	init->setLogAlpha(0);
	t0->setScaleFactor(1);

	// compute forward probabilities at time t using the alpha values for time t-1
	for (unsigned int t = 1; t < this->latentStatesZ.size(); t++) {
		LatentStateZ* ts = this->latentStatesZ[t];

		// loop through possible states at time t
		// scaleFactor = 1/sum(probs)
		double scaleFactorDenom = 0;
		for (LatentStateZ::iterator it = ts->begin(); it != ts->end(); it++) {
			std::vector<Transition*>& ins = (*it)->getInTrans();
			std::vector<double> logProbs(ins.size());
			for (unsigned int i = 0; i < ins.size(); i++) {
				Transition* trans = ins[i];
				double logProb = trans->from->getLogAlpha() +
						getTransProb(trans) +
						getEmitProb(trans);
				logProbs[i] = logProb;
//				std::cout << "forward: " << trans->from->getLogAlpha() << "\t" << getTransProb(trans) << "\t" << getEmitProb(trans) << "\t" << logProb << "\n";
				scaleFactorDenom += exp(logProb);
			}

			(*it)->setLogAlpha(sumLogProb(logProbs));
		}

//		std::cout << t << " *************** saving scale factors\n";
		double scaleFactor = 1/scaleFactorDenom;
		ts->setScaleFactor(scaleFactor);	// needed for backward algorithm

		// scale logAlphas
//		std::cout << t <<  "*************** scaling log alphas\n";
		for (LatentStateZ::iterator it = ts->begin(); it != ts->end(); it++) {
			double scaledLogAlpha = (*it)->getLogAlpha() + log(ts->getScaleFactor());
			(*it)->setLogAlpha(scaledLogAlpha);
		}
	}

//	std::cout << "************** done forward\n";
	this->forwardHasBeenRun = true;
}

void Hmm::backward(bool force) {
	// don't run if backward() has been run before and force is false
	if (!force && this->hasBackwardBeenRun()) {
		return;
	}

	// run forward if it hasn't been run because this is dependent
	// on scaling factors generated by the forward process
	if (force || !this->hasForwardBeenRun()) {
		this->forward(true);
	}

	int T = this->latentStatesZ.size()-1;
	if (T<1) // no observation
		return;
	for (int t = T; t>=0; t--) {
		LatentStateZ* zt = this->latentStatesZ[t];

		double scaleFactor = zt->getScaleFactor();
		for (LatentStateZ::iterator it = zt->begin(); it!=zt->end(); it++) {
			PossibleState* node = *it;
			if (t==T)
				node->setLogBeta(0 + log(scaleFactor));	// 0 here just for documentation (log(0) = 1)
			else {
				std::vector<Transition*>& outs = node->getOutTrans();
				std::vector<double> logProbs(outs.size());
				for (unsigned int i = 0; i<outs.size(); i++) {
					Transition* trans = outs[i];
					double logProb = trans->to->getLogBeta()+getTransProb(trans)+getEmitProb(trans);
					logProbs[i] =logProb;
				}
				node->setLogBeta(sumLogProb(logProbs) + log(scaleFactor));
			}
		}
	}

	this->backwardHasBeenRun = true;
}

void Hmm::forwardBackward(bool force) {

	// run forward and backward if they haven't been run before
	// or force run them
	if (force) {
//		std::cout << "******************RUNNING FORWARD\n";
		this->forward(true);
//		std::cout << "******************RAN FORWARD\n";
		this->backward(true);
//		std::cout << "******************RAN BACKWARD\n";
	} else {
		if (this->hasForwardBackwardBeenRun()) {
			return;
		}
		if (!this->hasForwardBeenRun()) {
			this->forward();
		}
		if (!this->hasBackwardBeenRun()) {
			this->backward();
		}
	}

	// compute forward probabilities at time t using the alpha values for time t-1
	double totalProb;
	for (unsigned int i = 1; i < this->latentStatesZ.size(); i++) {
		LatentStateZ* zi = this->latentStatesZ[i];
		totalProb = 0;

		// loop through possible states
		for (LatentStateZ::iterator it = zi->begin(); it != zi->end(); it++) {
			if (this->obsType == PHASE_CONCORDANCE) {
				// see Rabiner (eq 27, 96, 97)
				// in eq 97 substitute t+1 with t (results in extra scale factor)
				// probObs is implicitly in the scaling coefficients of alpha and beta
				// but there is one extra scale factor at zi that needs to be removed
				double logPosterior =	(*it)->getLogAlpha() +
										(*it)->getLogBeta() -
										log(zi->getScaleFactor());

				(*it)->setLogPosterior(logPosterior);
//				std::cout << exp((*it)->getLogAlpha()) << " " << exp((*it)->getLogBeta()) << " " << exp(this->getLogProbObs()) << " " << exp(logPosterior) << std::endl;
			} else {
				if (this->obsType == RAF_DEVIATION) {
					// for READ_COUNTS method, we are normalizing to avoid calculating
					// the binomial coefficient for efficiency
					double logPosterior = (*it)->getLogAlpha() + (*it)->getLogBeta();
					totalProb += exp(logPosterior);
				}
			}
		}

		if (this->obsType == RAF_DEVIATION) {
			for (LatentStateZ::iterator it = zi->begin(); it != zi->end(); it++) {
				double posterior =	exp((*it)->getLogAlpha() +
										(*it)->getLogBeta()) / totalProb;
				(*it)->setLogPosterior(log(posterior));

			}
		}
	}

	this->forwardBackwardHasBeenRun = true;
}

double Hmm::viterbi(std::vector<Transition*>& path) {
	// set nodes at time 0 according to initial probabilities.
	LatentStateZ* ts = this->latentStatesZ[0];
	PossibleState* init = (*ts)[0];
	init->setLogAlpha(0);

	// find the best path up to path t;
	for (unsigned int t = 1; t < this->latentStatesZ.size(); t++) {
		ts = this->latentStatesZ[t];
		for (LatentStateZ::iterator it = ts->begin(); it!=ts->end(); it++) {
			PossibleState* node = *it;
			std::vector<Transition*>& ins = node->getInTrans();
			double maxProb = log(0.0);
			Transition* bestTrans = 0;
			for (unsigned int i = 0; i<ins.size(); i++) {
				Transition* trans = ins[i];
				double logProb = trans->from->getLogAlpha()+getTransProb(trans)+getEmitProb(trans);
				if (bestTrans==0 || maxProb<logProb) {
					bestTrans = trans;
					maxProb = logProb;
				}
			}
			node->setLogAlpha(maxProb); // store the highest probability in logAlpha
			node->setPsi(bestTrans); // store the best transition in psi
		}
	}

	// Find the best node at time T. It will be the last node in the best path
	ts = this->latentStatesZ[this->latentStatesZ.size()-1];
	PossibleState* best = 0;
	for (LatentStateZ::iterator it = ts->begin(); it!=ts->end(); it++) {
		PossibleState* node = *it;
		if (best==0 || best->getLogAlpha()<node->getLogAlpha())
			best = node;
	}

	// retrieve the nodes in the best path
	for (PossibleState* nd = best; nd;) {
		if (nd->getPsi()) {
			path.push_back(nd->getPsi());
			nd = nd->getPsi()->from;
		}
		else
			nd = 0;
	}

	// reverse the path
	for (int i = 0, j=path.size()-1; i<j; i++, j--) {
		Transition* tmp = path[i];
		path[i] = path[j];
		path[j] = tmp;
	}

	// Done
	return best->getLogAlpha();
}

// Compute P(e_1:T) = sum_s P(e_1:T, x_T=s) = sum_s alpha_s(T);
double Hmm::calcLogProbObs() {
	if (this->latentStatesZ.size()<1)
		return 1; // no observations

	// see Rabiner (eq 103)
	// logProbObs = -sum(logScalingFactors)
	forward();
//	backward();
	double logProbObs = 0;
	for (unsigned int i = 1; i < this->latentStatesZ.size(); i++) {
		LatentStateZ* zi = this->latentStatesZ[i];
		double scaleFactor = zi->getScaleFactor();
		logProbObs += log(scaleFactor);
	}

	logProbObs *= -1;
	this->setLogProbObs(logProbObs);
	return logProbObs;
}

std::vector<double> Hmm::extractPosteriors(unsigned long state) {
	std::vector<double> posteriors;

	// had t=1 to ignore the initial state but easier to send it along and let the caller deal handle the one less posterior probability.
	for (unsigned int t = 0; t < this->latentStatesZ.size(); t++) {
		LatentStateZ* ts = this->latentStatesZ[t];
		// LatentStateZ::iterator it = ts->begin();
		bool posteriorExists = false;
		for (LatentStateZ::iterator it = ts->begin(); it != ts->end(); it++) {
			if ((*it)->getState() == state) {
				posteriors.push_back(exp((*it)->getLogPosterior()));
				posteriorExists = true;
				break;
			}
		}
		if (!posteriorExists) {
			posteriors.push_back(0);
		}
	}

	return posteriors;
}

void Hmm::estimateProbs(bool estimateTrans, bool estimateAberrantEmit, bool estimateNormalEmit, std::string normalState) {

	// (1, 2, 3) Calculate forward/backward probabilities
	//forwardBackward(true);	// this overwrites previous alphas and betas

	// (4) Compute diGammas and gammas
	double denom, numer = 0;

	// NOTE: t = 1 is skipping initial jump into chain (TODO: 1 or 0?) !!!!!!!!!WAS t=0 and size-1 !!!!!!
	for (unsigned int t = 0; t < this->latentStatesZ.size()-1; t++) {
		denom = 0;

		LatentStateZ* ls = this->latentStatesZ[t];
		LatentStateZ::iterator i = ls->begin();
		for (i = ls->begin(); i!=ls->end(); i++) {
			PossibleState* node_i = *i;
			std::vector<Transition*>& outs = node_i->getOutTrans();
			for (unsigned int j = 0; j < outs.size(); j++) {
				Transition* trans = outs[j];
				PossibleState* node_j = trans->to;
				denom += exp(node_i->getLogAlpha() +
							getTransProb(trans) +
							getEmitProb(trans) +
							node_j->getLogBeta());
			}
		}
		for (i = ls->begin(); i!=ls->end(); i++) {
			PossibleState* node_i = *i;
			std::vector<Transition*>& outs = node_i->getOutTrans();

			double gamma = 0;
			for (unsigned int j = 0; j < outs.size(); j++) {
				Transition* trans = outs[j];
				PossibleState* node_j = trans->to;
				double diGamma = 	exp(node_i->getLogAlpha() +
									getTransProb(trans) +
									getEmitProb(trans) +
									node_j->getLogBeta()) /
									denom;
				gamma += diGamma;
//				std::cout << t << "(" << this->getStr(trans->obs) << "): " << this->getStr(trans->from->getState()) << " to " << this->getStr(trans->to->getState()) << ": digamma=" << diGamma << "\t" << this->getStr(trans->from->getState()) << " gamma=" << gamma << std::endl;
				trans->setLogDiGamma(log(diGamma));
				node_i->setLogGamma(log(gamma));
			}
		}
	}

	// (5) Re-estimate probs
	// re-estimate initial probs
	// not going to do this

	// re-estimate alphas and betas
	HmmKeyedMaps& transProbs = this->getTransProbs();
	HmmKeyedMaps& emitProbs = this->getEmitProbs();

	// update transition probabilities
	if (estimateTrans) {
		unsigned long initStateId = this->getInitState();
		for (HmmKeyedMaps::iterator from = transProbs.begin(); from!=transProbs.end(); from++) {
			HmmMap& transitions = *from->second;
			for (HmmMap::iterator to = transitions.begin(); to!=transitions.end(); to++) {
				numer = 0;
				denom = 0;

				unsigned long fromId = from->first;
				unsigned long toId = to->first;

				//std::cerr << "re-estimating prob for " << strIdMap.getStr(fromId) << " to " << strIdMap.getStr(toId) << "(" << exp(to->second) << ")" << std::endl;

				// NOTE: don't use for emission re-estimation
				if (fromId == initStateId) {
					//std::cerr << "skipping " << strIdMap.getStr(fromId) << " to " << strIdMap.getStr(toId) << std::endl;
					continue;
				}

				// NOTE: would be 0 for emission prob re-estimation
				for (unsigned int t = 1; t < this->latentStatesZ.size() - 1; t++) {
					LatentStateZ* ls = this->latentStatesZ[t];

					for (LatentStateZ::iterator i = ls->begin(); i!=ls->end(); i++) {
						PossibleState* node_i = *i;
						std::vector<Transition*>& outs = node_i->getOutTrans();

						for (unsigned int j = 0; j < outs.size(); j++) {
							Transition* trans = outs[j];
							PossibleState* node_j = trans->to;

							// accumulate probabilistic counts for transition from i to j
							if (node_i->getState() == fromId && node_j->getState() == toId) {
								numer += exp(trans->getLogDiGamma());
								denom += exp(node_i->getLogGamma());
	//							std::cerr << t << "\tcalculating probs " << strIdMap.getStr(fromId) << " to " << strIdMap.getStr(toId) << "\t" << exp(trans->getLogDiGamma()) << "\t" << exp(node_i->getLogGamma()) << "\t" << numer << "\t" << denom << std::endl;

								//break;  // optimization
							}
						}  // end j
					}  // end i
				}  // end t
				double updatedLogTransProb = log(numer/denom);
				to->second = updatedLogTransProb;
//				std::cerr << "updated prob for " << strIdMap.getStr(fromId) << " to " << strIdMap.getStr(toId) << " (" << exp(to->second) << ")" << std::endl;

			}  // end to
		}  // end from
	}

	// update emission probabilities
	if (estimateAberrantEmit || estimateNormalEmit) {
		for (HmmKeyedMaps::iterator to = emitProbs.begin(); to!=emitProbs.end(); to++) {
			HmmMap& emissions = *to->second;
			for (HmmMap::iterator emit = emissions.begin(); emit!=emissions.end(); emit++) {
				numer = 0;
				denom = 0;

				unsigned long toId = to->first;
				unsigned long emitId = emit->first;

				//std::cerr << "re-estimating prob for " << strIdMap.getStr(fromId) << " to " << strIdMap.getStr(toId) << "(" << exp(to->second) << ")" << std::endl;

				// NOTE: would be 0 for emission prob re-estimation
				for (unsigned int t = 0; t < this->latentStatesZ.size() - 1; t++) {
					LatentStateZ* ls = this->latentStatesZ[t];
					LatentStateZ::iterator i = ls->begin();
					for (i = ls->begin(); i!=ls->end(); i++) {
						PossibleState* node_i = *i;
						std::vector<Transition*>& outs = node_i->getOutTrans();

						for (unsigned int j = 0; j < outs.size(); j++) {
							Transition* trans = outs[j];
							PossibleState* node_j = trans->to;

							// accumulate probabilistic counts for transition from i to j
							if (node_j->getState() == toId) {
								if (trans->obs == emitId) {
									numer += exp(node_j->getLogGamma());
								}
								denom += exp(node_j->getLogGamma());
							}
						}  // end j
					}  // end i
				}  // end t
				double updatedLogEmitProb = log(numer/denom);
				if ((estimateAberrantEmit && strIdMap.getStr(toId) != normalState) || (estimateNormalEmit && strIdMap.getStr(toId) == normalState)) {
					emit->second = updatedLogEmitProb;
//					std::cerr << "updated prob for " << strIdMap.getStr(toId) << " emits " << strIdMap.getStr(emitId) << "(" << exp(emit->second) << ")" << std::endl;
				}

			}  // end to
		}  // end from
	}
}

double Hmm::getPseudoCounts(PseudoCounts& counts) {
	double probObs = calcLogProbObs(); // this call includes a forward() call.
	backward();
	forwardBackward();

	// Compute the pseudo counts of transitions, emissions, and initializations
	for (unsigned int t = 0; t < this->latentStatesZ.size() - 1; t++) {
		LatentStateZ* ls = this->latentStatesZ[t];
		LatentStateZ::iterator it = ls->begin();

		// add the pseudo counts into counts
		for (it = ls->begin(); it!=ls->end(); it++) {
			PossibleState* node = *it;

			// state "total probabilities"
			double stateCount = node->getLogPosterior();
			counts.getStateCount().add(node->getState(), stateCount);

			// emission "total probabilities" (i.e., compute gamma(i)'s)
			std::vector<Transition*>& ins = node->getInTrans();
			unsigned int k;
			for (k = 0; k < ins.size(); k++) {
				Transition* trans = ins[k];
				PossibleState* from = trans->from;
				double emitCount = from->getLogAlpha()+getTransProb(trans)+getEmitProb(trans)+node->getLogBeta();
				counts.getEmitCount().add(node->getState(), trans->obs, emitCount);
			}

			// transition "total probabilities" (i.e., compute di-gamma(i)'s
			std::vector<Transition*>& outs = node->getOutTrans();
			for (k = 0; k < outs.size(); k++) {
				Transition* trans = outs[k];
				PossibleState* to = trans->to;
				double transCount = node->getLogAlpha()+getTransProb(trans)+getEmitProb(trans)+to->getLogBeta();
				counts.getTransCount().add(node->getState(), to->getState(), transCount);
			}
		}
	}

	return probObs;
}

void Hmm::baumWelch(std::vector<std::string>& observedSequence, const std::vector<double>& logRRSequence, std::vector<bool>& bAlleleSequence, unsigned int iterations, bool estimateTrans, bool estimateAberrantEmit, bool estimateNormalEmit, std::string normalState) {

	if (!estimateTrans && !estimateAberrantEmit && !estimateNormalEmit) {
		return;
	}

//	std::cerr << "Training with Baum-Welch for up to " << iterations << " iterations, using "
//			<< observedSequence.size() << " sequences." << std::endl;
	double prevTotalLogProb = -std::numeric_limits<double>::infinity();

	std::string READ_COUNTS_OBS = "0"; // TODO put in HMM constants
	// initialize hmm with sequences and calculate likelihood of data

	if (this->obsType == PHASE_CONCORDANCE) {
		for (unsigned int i = 0; i < observedSequence.size(); i++) {
			addObservation(observedSequence[i]);
		}
	} else {
//		if (this->obsType == READ_COUNTS) {
//			for (unsigned int i=1; i< observedSequence.size(); i++) {
//				addObservation(READ_COUNTS_OBS, logRRSequence[i-1], logRRSequence[i], atoi(observedSequence[i-1].c_str()), atoi(observedSequence[i].c_str()), bAlleleSequence[i-1], bAlleleSequence[i]);
//			}
//		}
	}
	forwardBackward();

	// now re-estimate parameters and calculate likelihood
//	std::cerr << "***** Initial probs: " << std::endl;
//	this->saveProbs("");
	for (unsigned int k = 0; k < iterations; k++) {
		double totalLogProb = 0;
		this->estimateProbs(estimateTrans, estimateAberrantEmit, estimateNormalEmit, normalState);
		totalLogProb = this->calcLogProbObs();
		this->getLikelihoodProfile().push_back(totalLogProb);
//		std::cerr << "Iteration " << k+1 << ' ' << "prevTotalLogProb=" << prevTotalLogProb << std::endl;
//		std::cerr << "Iteration " << k+1 << ' ' << "totalLogProb=" << totalLogProb << std::endl;

		// if improvement is negligible (specified by "delta") stop the algorithm
		// TODO: parameterize "delta"
		double delta = 0.001;
		if (prevTotalLogProb!=-std::numeric_limits<double>::infinity() && prevTotalLogProb + delta >= totalLogProb)
			break;
		else
			prevTotalLogProb = totalLogProb;

//		std::cerr << "***** Updated probs: " << std::endl;
//		this->saveProbs("");

		// forces re-calculation of ci's based on estimated probabilities
		this->forwardBackward(true);
	}
	// not sure if needed
	this->getLikelihoodProfile().push_back(this->calcLogProbObs());

}

void Hmm::baumWelch(std::vector<std::string>& observedSequence, const std::vector<unsigned int>& covSequence, std::vector<bool>& refSequence, unsigned int iterations, bool estimateTrans, bool estimateAberrantEmit, bool estimateNormalEmit, std::string normalState) {

	if (!estimateTrans && !estimateAberrantEmit && !estimateNormalEmit) {
		return;
	}

//	std::cerr << "Training with Baum-Welch for up to " << iterations << " iterations, using "
//			<< observedSequence.size() << " sequences." << std::endl;
	double prevTotalLogProb = -std::numeric_limits<double>::infinity();

	std::string READ_COUNTS_OBS = "0"; // TODO put in HMM constants
	// initialize hmm with sequences and calculate likelihood of data

	if (this->obsType == PHASE_CONCORDANCE) {
		for (unsigned int i = 0; i < observedSequence.size(); i++) {
			addObservation(observedSequence[i]);
		}
	} else {
		if (this->obsType == RAF_DEVIATION) {
			for (unsigned int i=1; i< observedSequence.size(); i++) {
				addObservation(READ_COUNTS_OBS, covSequence[i-1], covSequence[i], atoi(observedSequence[i-1].c_str()), atoi(observedSequence[i].c_str()), refSequence[i-1], refSequence[i]);
			}
		}
	}
	forwardBackward();

	// now re-estimate parameters and calculate likelihood
//	std::cerr << "***** Initial probs: " << std::endl;
//	this->saveProbs("");
	for (unsigned int k = 0; k < iterations; k++) {
		double totalLogProb = 0;
		this->estimateProbs(estimateTrans, estimateAberrantEmit, estimateNormalEmit, normalState);
		totalLogProb = this->calcLogProbObs();
//		std::cerr << "Iteration " << k+1 << ' ' << "prevTotalLogProb=" << prevTotalLogProb << std::endl;
//		std::cerr << "Iteration " << k+1 << ' ' << "totalLogProb=" << totalLogProb << std::endl;

		// if improvement is negligible (specified by "delta") stop the algorithm
		// TODO: parameterize "delta"
		double delta = 0.001;
		if (prevTotalLogProb!=-std::numeric_limits<double>::infinity() && prevTotalLogProb + delta >= totalLogProb)
			break;
		else
			prevTotalLogProb = totalLogProb;

//		std::cerr << "***** Updated probs: " << std::endl;
//		this->saveProbs("");

		// forces re-calculation of ci's based on estimated probabilities
		this->forwardBackward(true);
	}
}

void Hmm::updateProbs(PseudoCounts& counts) {
	this->transProbs.clear();
	this->emitProbs.clear();
	for (HmmKeyedMaps::iterator i = counts.getTransCount().begin(); i!=counts.getTransCount().end(); i++) {
		unsigned long from = i->first;
		double fromCount = counts.getStateCount().get(from);
		HmmMap& cnts = *i->second;
		for (HmmMap::iterator j = cnts.begin(); j!=cnts.end(); j++) {
			//      if (j->second-fromCount>=_minLogProb)
			this->transProbs.add(from, j->first, j->second-fromCount);
		}
	}
	for (HmmKeyedMaps::iterator s = counts.getEmitCount().begin(); s!=counts.getEmitCount().end(); s++) {
		unsigned long state = s->first;
		double stateCount = counts.getStateCount().get(state);
		HmmMap& cnts = *s->second;
		for (HmmMap::iterator o = cnts.begin(); o!=cnts.end(); o++) {
			//      if (o->second-stateCount>_minLogProb)
			this->emitProbs.add(state, o->first, o->second-stateCount);
		}
	}
}

PossibleState::PossibleState(int time, unsigned long state, Hmm* hmm) {
	this->time = time;
	this->state = state;
	this->logAlpha = 1;		// log(0) = 1 (NOTE: previously initialized to 0)
	this->logBeta = 1;
	this->logPosterior = 1;
	this->logGamma = 1;
	this->psi = NULL;
	this->hmm = hmm;
}

void PossibleState::print() {
	std::cerr << hmm->getStr(getState()) << '\t'
			<< "alpha=" << exp(getLogAlpha()) << '\t'
			<< "beta=" << exp(getLogBeta());
	unsigned int i;
	std::cerr << " (in";
	for (i = 0; i < this->inTrans.size(); i++) {
		std::cerr << ' ' << this->hmm->getStr(this->inTrans[i]->from->getState());
	}
	std::cerr << ") ";
	std::cerr << " (out";
	for (i = 0; i < this->outTrans.size(); i++) {
		std::cerr << ' ' << this->hmm->getStr(this->outTrans[i]->to->getState());
	}
	std::cerr << ")";
	std::cerr << std::endl;
}

PossibleState::~PossibleState() {
	for (unsigned int i = 0; i < this->outTrans.size(); i++) {
		delete this->outTrans[i];
	}
}

LatentStateZ::~LatentStateZ() {
	for (LatentStateZ::iterator it = begin(); it!=end(); it++) {
		delete (*it);
	}
}

void Hmm::reset() {
	for (unsigned int t = 0; t < this->latentStatesZ.size(); t++){
		delete this->latentStatesZ[t];
	}
	this->latentStatesZ.clear();
	this->forwardHasBeenRun = false;
	this->backwardHasBeenRun = false;
	this->forwardBackwardHasBeenRun = false;
}

Hmm::~Hmm() {
	reset();
}

// TODO need logic for PHASE_CONCORDANCE and READ_COUNTS
void Hmm::addObservation(unsigned long o, unsigned long n1, unsigned long n2, unsigned long r1, unsigned long r2, bool ref1, bool ref2) {
	std::vector<unsigned long> stateIds;

	// for READ_COUNTS method (and not explicitly specifying matrices), I believe cntxt should be 0
//	if (this->obsType == PHASE_CONCORDANCE) { was commented out
//		ULSet* cntx = this->emitProbs.getCntx(o);

// COMMENTING THIS IF/ELSE LOGIC OUT AND FORCING THE IF FIXED AN ISSUE WHERE WE WERE HAVING A POSTERIOR WITHOUT A VALUE
//		if (cntx==0) {	// if event doesn't have any emission sources, push back all emission states
			for (HmmKeyedMaps::iterator it = this->emitProbs.begin(); it != this->emitProbs.end(); it++) {
				stateIds.push_back(it->first);  // for READ_COUNTS method emitProbs = emitParam
			}
//		}
//		else {			// if event does have emission source, push back those sources
//			for (ULSet::iterator it = cntx->begin(); it!=cntx->end(); it++) {
//				stateIds.push_back(*it);
//			}
//		}
//	}

	// if this is the first observation, create a node for init
	if (this->latentStatesZ.empty()) {
		// create a special state for time 0;
		LatentStateZ* t0 = new LatentStateZ;
		t0->push_back(new PossibleState(0, this->initState, this));
		this->latentStatesZ.push_back(t0);
	}

	// create a time step for the observation to add
	LatentStateZ* ts = new LatentStateZ;
	int time = this->latentStatesZ.size();

	// for every possible emitting state, create a node at current time step
//	std::cout << "time " << time << " has " << stateIds.size() << " states\n";
	for (unsigned int i = 0; i < stateIds.size(); i++) {
		PossibleState* node = new PossibleState(time, stateIds[i], this);
		ts->push_back(node);
		LatentStateZ* prev = this->latentStatesZ[time-1];

		// connect every new node to all previous nodes where trans prob exists and in lattice diagram
		for (LatentStateZ::iterator it = prev->begin(); it!=prev->end();it++) {
			ULSet* possibleSrc = this->transProbs.getCntx(node->getState());
			if (possibleSrc && possibleSrc->find((*it)->getState())!=possibleSrc->end()) {
				new Transition(*it, node, o, n1, n2, r1, r2, ref1, ref2);
			}
		}
	}
	this->latentStatesZ.push_back(ts);
}

std::string Hmm::getObservation(unsigned long time) {
	LatentStateZ* ts = this->latentStatesZ[time];
	PossibleState* ps = (*ts)[0L];
	std::vector<Transition*>& ins = ps->getInTrans();
	return this->getStr(ins[0]->obs);
}

// there is no observation at time 0 because observations are inter-marker events
std::vector<std::string> Hmm::getObservations() {
	std::string nullString = "NA";
	std::vector<std::string> observations;
//	std::cout << "Num latent states " << this->latentStatesZ.size() << "\n";
	observations.push_back(nullString);
	for (unsigned int t = 1; t < this->latentStatesZ.size(); t++) {
		LatentStateZ* ts = this->latentStatesZ[t];
		PossibleState* ps = (*ts)[0L];
		std::vector<Transition*>& ins = ps->getInTrans();
		observations.push_back(this->getStr(ins[0]->obs));
	}

	return observations;
}

//void Hmm::loadProbs(double& eventPrevalence, unsigned int& eventLengthMarkers, double& initialParamNormal, double& initialParamEvent, std::string normalState) {
//	std::string initState = "_init";
//
//	double lambda_0 = 1.0/( ((1.0-eventPrevalence)/eventPrevalence) * (double) eventLengthMarkers);
//	double lambda_1 = 1.0/((double) eventLengthMarkers);
//	double probStayEventState = 1.0 - lambda_1;
//
//	this->initState = this->strIdMap.getId(initState);
//	this->transProbs.load(initState, normalState, 1 - eventPrevalence, this->strIdMap);
//	this->transProbs.load(normalState, normalState, 1 - lambda_0, this->strIdMap);
//
//	for (unsigned int i = 1; i < this->numStates; i++) {
//		std::string state = "S" + boost::lexical_cast<std::string>(i);
//		this->transProbs.load(initState, state, eventPrevalence/(this->numStates - 1), this->strIdMap);
//		this->transProbs.load(normalState, state, lambda_0/(this->numStates - 1), this->strIdMap);
//		this->transProbs.load(state, state, probStayEventState, this->strIdMap);
//		this->transProbs.load(state, normalState, 1-probStayEventState, this->strIdMap);
//	}
//
//	// keep track of what the states are
//	for (HmmKeyedMaps::iterator it = transProbs.begin(); it != transProbs.end(); it++) {
//		unsigned long state = it->first;
//		if(std::find(this->states.begin(), this->states.end(), state) == this->states.end()) {
//		    /* v does not contain x */
//			if (state != this->initState) {
//				this->states.push_back(state);
//			}
//		}
//	}
//
//	std::string concordantEmission = "1";
//	std::string discordantEmission = "0";
//	initialParamNormal = 0.5;
//	initialParamEvent = 0.5;
//	double deviationFromNormal = (1-initialParamNormal)/(this->numStates);
//	this->emitProbs.load(normalState, concordantEmission, initialParamNormal, this->strIdMap);
//	this->emitProbs.load(normalState, discordantEmission, 1-initialParamNormal, this->strIdMap);
//	for (unsigned int i = 1; i < this->numStates; i++) {
//		std::string state = "S" + boost::lexical_cast<std::string>(i);
//		double concordantProb = initialParamEvent + deviationFromNormal * i;
////		// just adding a little state-specific deviation
////		if (concordantProb != 0.99) {
////			concordantProb += 0.01 * i;
////		} else {
////			concordantProb -= 0.01 * i;
////		}
//		double discordantProb = 1 - concordantProb;
//		this->emitProbs.load(state, concordantEmission, concordantProb, this->strIdMap);
//		this->emitProbs.load(state, discordantEmission, discordantProb, this->strIdMap);
//	}
//}

void Hmm::loadProbs(double& lambda_0, double& lambda_1, double& alpha_0, double& alpha_1, std::string normalState) {
	std::string initState = "_init";

	double probStayEventState = 1.0 - lambda_1;

	this->initState = this->strIdMap.getId(initState);
	this->transProbs.load(initState, normalState, 1 - lambda_0, this->strIdMap);
	this->transProbs.load(normalState, normalState, 1 - lambda_0, this->strIdMap);

	for (unsigned int i = 1; i < this->numStates; i++) {
		std::string state = "S" + boost::lexical_cast<std::string>(i);
		this->transProbs.load(initState, state, lambda_0/(this->numStates - 1), this->strIdMap);
		this->transProbs.load(normalState, state, lambda_0/(this->numStates - 1), this->strIdMap);
		this->transProbs.load(state, state, probStayEventState, this->strIdMap);
		this->transProbs.load(state, normalState, 1-probStayEventState, this->strIdMap);
	}

	// keep track of what the states are
	for (HmmKeyedMaps::iterator it = transProbs.begin(); it != transProbs.end(); it++) {
		unsigned long state = it->first;
		if(std::find(this->states.begin(), this->states.end(), state) == this->states.end()) {
		    /* v does not contain x */
			if (state != this->initState) {
				this->states.push_back(state);
			}
		}
	}

	std::string concordantEmission = "1";
	std::string discordantEmission = "0";
	std::string readCountsParam = "0";	// this is just a name for parameter 0.  if we have more parameters, would call them 1...n

//	double deviationFromNormal = (1-alpha_0)/(this->numStates);
	// double stateStagger = 0.01;
	if (this->obsType == PHASE_CONCORDANCE) {
		this->emitProbs.load(normalState, concordantEmission, alpha_0, this->strIdMap);
		this->emitProbs.load(normalState, discordantEmission, 1-alpha_0, this->strIdMap);
	} else {
		if (this->obsType == RAF_DEVIATION) {
			double dState = alpha_0;
			this->emitProbs.load(normalState, readCountsParam, dState, this->strIdMap);
		}
	}

	for (unsigned int i = 1; i < this->numStates; i++) {
		std::string state = "S" + boost::lexical_cast<std::string>(i);
		double concordantProb = alpha_1 + (1 - alpha_1) / this->numStates * (i-1);	// staggering probabilities for event states
//		// just adding a little state-specific deviation
//		if (concordantProb != (1 - stateStagger)) {
//			concordantProb += stateStagger * i;
//		} else {
//			concordantProb -= stateStagger * i;
//		}
		double discordantProb = 1 - concordantProb;

		if (this->obsType == PHASE_CONCORDANCE) {
			this->emitProbs.load(state, concordantEmission, concordantProb, this->strIdMap);
			this->emitProbs.load(state, discordantEmission, discordantProb, this->strIdMap);
		} else {
			double dState = alpha_1 + + (1 - alpha_1) / this->numStates * (i-1);  // staggering probabilities for event states
			if (this->obsType == RAF_DEVIATION) {
				this->emitProbs.load(state, readCountsParam, dState, this->strIdMap);
			}
		}

	}
}

void Hmm::loadProbs(std::string initProbFilename, std::string transProbFilename, std::string emitProbFilename) {

	std::ifstream initProb(initProbFilename.c_str());
	std::string initState;
	initProb >> initState;	// initState will contain first line/token of initProbFile
	this->initState = this->strIdMap.getId(initState);	// stores init state in nameIdMap
	this->transProbs.load(initProb, this->strIdMap);   // transProbs contains init and trans

	std::ifstream transProb(transProbFilename.c_str());
	this->transProbs.load(transProb, this->strIdMap);

	// keep track of what the states are
	for (HmmKeyedMaps::iterator it = transProbs.begin(); it != transProbs.end(); it++) {
		unsigned long state = it->first;
		if(std::find(this->states.begin(), this->states.end(), state) == this->states.end()) {
		    /* v does not contain x */
			if (state != this->initState) {
				this->states.push_back(state);
			}
		}
	}

	std::ifstream emitProb(emitProbFilename.c_str());
	this->emitProbs.load(emitProb, this->strIdMap);
}

void Hmm::saveProbs(std::string name) {
	if (name=="") {
		std::cerr << "---------------------------------" << std::endl;
		std::cerr << "transition probabilities:" << std::endl;
		this->transProbs.save(std::cerr, this->strIdMap);
		std::cerr << "---------------------------------" << std::endl;
		std::cerr << "emission probabilities:" << std::endl;
		this->emitProbs.save(std::cerr, this->strIdMap);
		std::cerr << "---------------------------------" << std::endl;
	}
	else {
		std::string s = name+".transitions";
		std::ofstream transProb(s.c_str());
		transProb << this->strIdMap.getStr(this->initState) << std::endl;
		this->transProbs.save(transProb, this->strIdMap);

		s = name+".emissions";
		std::ofstream emitProb(s.c_str());
		this->emitProbs.save(emitProb, this->strIdMap);
	}
}

void Hmm::writeTrans(std::ostream& file) {
	file << this->strIdMap.getStr(this->initState) << std::endl;
	this->transProbs.save(file, this->strIdMap);
}

void Hmm::writeEmit(std::ostream& file) {
	this->emitProbs.save(file, this->strIdMap);
}

void Hmm::print() {
	for (unsigned int i = 0; i < this->latentStatesZ.size(); i++) {
		LatentStateZ* ts = this->latentStatesZ[i];
		std::cerr << "TIME=" << i << std::endl;
		for (unsigned int s = 0; s<ts->size(); s++) {
			(*ts)[s]->print();
		}
	}
}

void Hmm::readSeqs(std::string obsSeqFilename, std::vector<std::vector<unsigned long>*>& sequences) {
	std::string line;
	const std::string delims(" ");
	std::ifstream obs(obsSeqFilename.c_str());
	while (getline(obs, line)) {
		std::vector<unsigned long>* seq = new std::vector<unsigned long>;
		std::string::size_type begIdx, endIdx;
		begIdx = line.find_first_not_of(delims);
		while (begIdx!=std::string::npos) {
			if (line[begIdx]=='#') // the rest of the line are comments
				break;
			endIdx = line.find_first_of(delims, begIdx);
			if (endIdx==std::string::npos) {
				endIdx = line.length();
			}
			std::string word = line.substr(begIdx, endIdx-begIdx);
			seq->push_back(getId(word));
			begIdx = line.find_first_not_of(delims, endIdx);
		}
		if (seq->size()>0)
			sequences.push_back(seq);
		else
			delete seq;
	}
}


void Hmm::genSeqs(std::ostream& ostrm, int seqs) {
	std::vector<unsigned long> seq;
	for (int i = 0; i<seqs; i++) {
		genSeq(seq);
		for (unsigned int k = 0; k<seq.size(); k++) {
			if (k)
				ostrm << ' ';
			ostrm << this->strIdMap.getStr(seq[k]);
		}
		ostrm << std::endl;
		seq.clear();
	}
}


void Hmm::genSeq(std::vector<unsigned long>& seq) {
	unsigned long state = this->initState, next, obs;
	while (true) {
		if (!this->transProbs.rand(state, next) || !this->emitProbs.rand(next, obs))
			break;
		state = next;
		seq.push_back(obs);
	}
}

void PseudoCounts::print(StrIdMap& strIdMap) {
	std::cerr << "TRANSITION"<< std::endl;
	this->transCount.save(std::cerr, strIdMap);
	std::cerr << "*********************" << std::endl;
	std::cerr << "EMISSION"<< std::endl;
	this->emitCount.save(std::cerr, strIdMap);
	std::cerr << "*********************" << std::endl;
	std::cerr << "STATE"<< std::endl;
	this->stateCount.save(std::cerr, strIdMap);
}

Transition::Transition(PossibleState* from, PossibleState* to, unsigned long obs, unsigned int n1, unsigned int n2, unsigned int r1, unsigned int r2, bool ref1, bool ref2) {
	this->from = from;
	this->to = to;
	this->obs = obs;
	if (this->from && this->to) {
		this->from->getOutTrans().push_back(this);
		this->to->getInTrans().push_back(this);
	}
	this->logDiGamma = 1;	// log(0) = 1
	this->n1 = n1;
	this->n2 = n2;
	this->r1 = r1;
	this->r2 = r2;
	this->ref1 = ref1;
	this->ref2 = ref2;
}

Transition::~Transition() {}

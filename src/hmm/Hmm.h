#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>


#include "MathUtil.h"
#include "StrIdMap.h"
#include "DataStructures.h"

#ifndef HMM_H_
#define HMM_H_

class PossibleState;
class Hmm;

/** A transition between two consecutive Hmm nodes. */
class Transition {
public:
	unsigned long obs;						// observed value transitioning to: this is 1 or 0 for the PHASE_CONCORDANCE method
											// it is the reference allele count for the READ_COUNTS method

	// extra observations needed for the READ_COUNTS method
	unsigned int n1,n2,r1,r2;						// this is coverage at a site and ref_coverage: used for the READ_COUNTS method
	bool ref1,ref2;								// true if emitting reference allele, otherwise false

	PossibleState* from;					// source node
	PossibleState* to;						// destination node

	double logDiGamma;						// used for parameter re-estimation
	void setLogDiGamma(double logDiGamma) 	{ this->logDiGamma = logDiGamma;	}
	double getLogDiGamma() const 			{ return this->logDiGamma;}
	Transition(PossibleState* from, PossibleState* to, unsigned long obs, unsigned int n1, unsigned int n2, unsigned int r1, unsigned int r2, bool ref1, bool ref2);
	~Transition();
};

/** A node in an Hmm object. */
class PossibleState {
	int time; 								// The time slot for this node.
	Hmm* hmm; 								// The hmm that this node belongs to
	unsigned long state; 					// the state
	std::vector<Transition*> inTrans; 		// incoming transitions
	std::vector<Transition*> outTrans;		// outgoing transitions
	double logAlpha; 						// alpha_t(s) = P(e_1:t, x_t=s);
	double logBeta;  						// beta_t(s)  = P(e_t+1:T | x_t=s);
	double logPosterior;					// posterior probability that node is in this state given observations
	double logGamma;						// used for parameter re-estimation
	Transition* psi; 						// the last transition of the most probable path that reaches this node
public:
	int getTime() 							{ return this->time; 			}
	unsigned long getState() const 			{ return this->state;			}
	void setLogAlpha(double logAlpha) 		{ this->logAlpha = logAlpha;	}
	double getLogAlpha() const 				{ return this->logAlpha;		}
	void setLogBeta(double logBeta) 		{ this->logBeta = logBeta;		}
	double getLogBeta() const 				{ return this->logBeta;			}
	void setLogPosterior(double logPosterior) 	{ this->logPosterior = logPosterior;	}
	double getLogPosterior() const 			{ return this->logPosterior;			}
	void setLogGamma(double logGamma) 	{ this->logGamma = logGamma;	}
	double getLogGamma() const 			{ return this->logGamma;}

	Transition* getPsi() 					{ return this->psi;				}
	void setPsi(Transition* psi) 			{ this->psi = psi;				}
	std::vector<Transition*>& getInTrans()	{ return this->inTrans;			}
	std::vector<Transition*>& getOutTrans()	{ return this->outTrans;		}

	void print();

	PossibleState(int time, unsigned long state, Hmm* hmm);
	~PossibleState();
};

/** The possible states at a particular time. */
// This kind of notation allows you to use an iterator
// on LatentStateZ to loop through possible states.
class LatentStateZ : public std::vector<PossibleState*> {
	double scaleFactor;
public:
	void setScaleFactor(double scaleFactor) { this->scaleFactor = scaleFactor;	}
	double getScaleFactor() const 			{ return this->scaleFactor;			}

	~LatentStateZ();
};

/** Pseudo Counts */
class PseudoCounts {
	HmmMap stateCount;
	HmmKeyedMaps transCount;
	HmmKeyedMaps emitCount;
public:
	HmmMap& getStateCount() 			{ return this->stateCount;}
	HmmKeyedMaps& getTransCount() 		{ return this->transCount;}
	HmmKeyedMaps& getEmitCount()  		{ return this->emitCount;}
	void print(StrIdMap& strIdMap);
};

/** An Hmm object implements the Hidden Markov Model. */
class Hmm {
	unsigned long initState; 					// the initial state
	std::vector<unsigned long> states;			// states (not including the init state)
	HmmKeyedMaps transProbs;    				// transition probabilities
	HmmKeyedMaps emitProbs;      				// emission probabilities
	StrIdMap strIdMap;        					// mapping between strings and integers
	std::vector<LatentStateZ*> latentStatesZ; 	// the time steps
	std::vector<double> likelihoodProfile;		// likelihood profile (for diagnostics)
	double minLogProb;       					// log probabilities lower than this are set to 0
	double logProbObs;
	bool forwardHasBeenRun;
	bool backwardHasBeenRun;
	bool forwardBackwardHasBeenRun;
	unsigned int numStates;						// number of possible HMM states
	std::string obsType;
	double epsilon;								// phase error (currently treated as a constant)

public:
	// possible observed types for the HMM
	static const std::string PHASE_CONCORDANCE, RAF_DEVIATION, BAF_DEVIATION;

	void setInitState(std::string initState)	{ this->initState = this->strIdMap.getId(initState);}
	void setObsType(std::string obsType)		{ this->obsType = obsType;	}
	std::string& getObsType()					{ return this->obsType;		}
	std::vector<unsigned long>& getStates()		{ return this->states;		}
	unsigned long getInitState()				{ return this->initState;	}
	HmmKeyedMaps& getEmitProbs()				{ return this->emitProbs;	}
	HmmKeyedMaps& getTransProbs()				{ return this->transProbs;	}
	StrIdMap& getStrIdMap()						{ return this->strIdMap;	}
	void setLogProbObs(double logProbObs)		{ this->logProbObs = logProbObs;	}
	double getLogProbObs()						{ return this->logProbObs;	}
	std::vector<LatentStateZ*>& getLatentStatesZ()	{ return this->latentStatesZ;	}
	std::vector<double>& getLikelihoodProfile()	{ return this->likelihoodProfile;	}

	void setEpsilon(double epsilon)		{ this->epsilon = epsilon;	}
	double getEpsilon()						{ return this->epsilon;	}

	double getTransProb(Transition* trans);
	double getEmitProb(Transition* trans);

	bool hasForwardBeenRun()					{ return this->forwardHasBeenRun;	}
	bool hasBackwardBeenRun()					{ return this->backwardHasBeenRun;	}
	bool hasForwardBackwardBeenRun()			{ return this->forwardBackwardHasBeenRun;	}

	void forward(bool force=false);  		// compute the forward probabilities P(e_1:t, X_t=s)
	void backward(bool force=false); 		// compute the backward probabilities P(e_t+1:T | X_t=s)
	void forwardBackward(bool force=false);	// compute forward/backward probabilities

	/** Retrieves posterior probabilities in a vector */
	std::vector<double> extractPosteriors(unsigned long state);

	/** Re-compute the transition and emission probabilities according
      to the pseudo counts. */
	void updateProbs(PseudoCounts& counts);

	/** Accumulate pseudo counts using the BaumWelch algorithm.  The
      return value is the probability of the observations according to
      the current model.  */
	double getPseudoCounts(PseudoCounts& counts);

	void estimateProbs(bool estimateTrans, bool estimateAberrantEmit, bool estimateNormalEmit, std::string normalState);

	/** Add an observation into the Hmm after the current last time
      slot. The states that have non-zero probability of generating
      the observation will be created. The transition between the new
      states and the states in the previous time slots will also be
      created.*/
	void addObservation(unsigned long obs, unsigned long n1, unsigned long n2, unsigned long r1, unsigned long r2, bool ref1, bool ref2);

	/** Same as void addObservation(unsigned long obs) above, except
      with a different form of parameter. */
	// TODO: move implementation to Hmm.cpp
	void addObservation(std::string obs, unsigned long n1, unsigned long n2, unsigned long r1, unsigned long r2, bool ref1, bool ref2) {
		addObservation(atoi(obs.c_str()), n1, n2, r1, r2, ref1, ref2);
	}

	// Dummy values are inserted for the PHASE_CONCORDANCE method
	void addObservation(std::string obs) {
		addObservation(strIdMap.getId(obs), 0, 0, 0, 0, false, false);
	}

	std::string getObservation(unsigned long time);
	std::vector<std::string> getObservations();

	/** Read the transition and emission probability tables from the
      files NAME.trans and NAME.emit, where NAME is the value of the
      variable name.*/
	void loadProbs(std::string initProbFilename, std::string transProbFilename, std::string emitProbFilename);
	void loadProbs(double& lambda_0, double& lambda_1, double& alpha_0, double& alpha_1, std::string normalState);
//	void loadProbs(double& eventPrevalence, unsigned int& eventLengthMarkers, double& initialParamNormal, double& initialParamEvent, std::string normalState);

	/** Save the transition and emission probability tables into the
      files NAME.trans and NAME.emit, where NAME is the value of the
      variable name. If name is "", both tables are printed on the
      standard output. */
	void saveProbs(std::string name="");
	void writeTrans(std::ostream& file);
	void writeEmit(std::ostream& file);

	/** Read the training data from the input stream. Each line in the
      input stream is an observation sequence. */
	void readSeqs(
			std::string obsSeqFilename, std::vector<std::vector<unsigned long>*>& sequences
	);


	/** Find the state sequence (a path) that has the maximum
      probability given the sequence of observations:
         max_{x_1:T} P(x_1:T | e_1:T);
      The return value is the logarithm of the joint probability 
      of the state sequence and the observation sequence:
        log P(x_1:T, e_1:T)
	 */
	double viterbi(std::vector<Transition*>& path);

	/** return the logarithm of the observation sequence: log P(e_1:T) */
	double calcLogProbObs();

	/** Train the model with the given observation sequences using the
      Baum-Welch algorithm. */
	void baumWelch(std::vector<std::string>& observedSequence, const std::vector<unsigned int>& covSequence, std::vector<bool>& refSequence, unsigned int maxIterations, bool estimateTrans, bool estimateAberrantEmit, bool estimateNormalEmit, std::string normalState);
	void baumWelch(std::vector<std::string>& observedSequence, const std::vector<double>& logRRSequence, std::vector<bool>& bAlleleSequence, unsigned int iterations, bool estimateTrans, bool estimateAberrantEmit, bool estimateNormalEmit, std::string normalState);

	/** Conversion between the integer id and string form of states and
      observations. */
	std::string getStr(unsigned long id) { return this->strIdMap.getStr(id);}
	unsigned long getId(std::string str) { return this->strIdMap.getId(str);}

	/** Clear all time slots to get ready to deal with another
      sequence. */
	void reset();

	/** Print the states at all time slots and the alpha/beta values at
      these states. */
	void print();

	/** Generate seqs observation sequences according to the model. */
	void genSeqs(std::ostream& ostrm, int seqs);

	/** Generate an observation sequence with up to maxlen elements
      according to the model. */
	void genSeq(std::vector<unsigned long>& seq);

	Hmm(unsigned int numStates, std::string& obsType);
	~Hmm();
};

#endif /* HMM_H_ */

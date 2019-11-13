#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <cmath>
#include <vector>

/** The input array contains a set of log probabilities lp1, lp2, lp3
    ... The return value should be the log of the sum of the
    probabilities: log(e^lp1 + e^lp2 + e^lp3 + ...) */
double sumLogProb(std::vector<double>& logprobs);

/** returns log (e^logprob1 + e^logprob2). */
double sumLogProb(double logprob1, double logprob2);

#endif

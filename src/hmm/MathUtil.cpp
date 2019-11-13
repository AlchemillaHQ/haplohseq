#include "MathUtil.h"

/** The input array contains a set of log probabilities lp1, lp2, lp3
    ... The return value should be the log of the sum of the
    probabilities: log(e^lp1 + e^lp2 + e^lp3 + ...) */
double sumLogProb(std::vector<double>& logprobs)
{
  double max = 0;
  unsigned int i;
  for (i = 0; i<logprobs.size(); i++) {
    if (i==0 || logprobs[i]>max)
      max = logprobs[i];
  }
  if (std::isinf(max)) // the largest probability is 0 (log prob= -inf)
    return max;   // return log 0
  double p = 0;
  for (i = 0; i<logprobs.size(); i++) {
    p += exp(logprobs[i]-max);
  }
  return max + log(p);
}

/** returns log (e^logprob1 + e^logprob2). */
double sumLogProb(double logprob1, double logprob2)
{
  if (std::isinf(logprob1) && std::isinf(logprob2))
    return logprob1; // both prob1 and prob2 are 0, return log 0.
  if (logprob1>logprob2)
    return logprob1+log(1+exp(logprob2-logprob1));
  else
    return logprob2+log(1+exp(logprob1-logprob2));
}


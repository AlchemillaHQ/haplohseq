#ifndef TABLES_H
#define TABLES_H

#include <cmath>
#include <iostream>
#include <set>
#include <stdlib.h>

#include "MathUtil.h"
#include "StrIdMap.h"

/** One dimensional table for storing first degree probabilities or
    frequency counts (in log scale).  */
class HmmMap : public std::map<unsigned long, double> {
	double _smoothedZeroCount;
public:
	HmmMap();
	double smoothedZeroCount() { return _smoothedZeroCount;}
	double get(unsigned long event);
	void add(unsigned long event, double count);
	void load(std::istream& file, StrIdMap& strIdMap);
	void save(std::ostream& file, StrIdMap& strIdMap);

	/** Randomly generate the next value according to the probabilities
      in this table. Return true if the next value is filled. */
	bool rand(unsigned long& next);
	~HmmMap() {}
};

typedef std::set<unsigned long> ULSet; // set of unsigned long integers

/** Two-dimensional table for storing second degree probabilities or
    frequency counts (in log scale).  */
class HmmKeyedMaps : public std::map<unsigned long, HmmMap*> {
	std::map<unsigned long, ULSet*> possibleContexts;
	HmmMap _backoff;
public:
	/** Retrieve the value associated with event and context. */
	double get(unsigned long event, unsigned long context);

	/** Add the count to the  */
	void add(unsigned long event, unsigned long context, double count);

	/** Return the set of contexts such that (event, context) has an
      associated values. */
	ULSet* getCntx(unsigned long event);

	/** fill the table with data in the file. */
	void load(std::istream& file, StrIdMap& strIdMap);
	void load(std::string& name1, std::string& name2, const double value, StrIdMap& strIdMap);

	/** save the contents in the table in the file. */
	void save(std::ostream& file, StrIdMap& strIdMap);

	/** clear the contents of this table */
	void clear();

	/** Randomly generate the next value according to the probabilities
      in this table. Return true if the next value is filled. */
	bool rand(unsigned long curr, unsigned long& next);

	~HmmKeyedMaps();
};

#endif

#include "DataStructures.h"

const double SMOOTHEDZEROCOUNT=-40;

void HmmKeyedMaps::clear() {
	for (HmmKeyedMaps::iterator i = begin(); i!=end(); i++) {
		delete i->second;
	}
	std::map<unsigned long, HmmMap*>::clear();
	for (std::map<unsigned long, ULSet*>::iterator j = possibleContexts.begin(); j!=possibleContexts.end(); j++) {
		delete j->second;
	}
	possibleContexts.clear();
	_backoff.clear();
}

HmmKeyedMaps::~HmmKeyedMaps() {
	this->clear();
	// Note: the possibleContexts object will get deleted automatically because it is not a pointer,
	// but its contents are pointers, so what they point to needs to be deleted.  This is accomplished
	// with the clear() function call.
	for (HmmKeyedMaps::iterator it = begin(); it!=end(); it++) {
		delete it->second;
	}
}

ULSet* HmmKeyedMaps::getCntx(unsigned long event) {

	std::map<unsigned long, ULSet*>::iterator g = possibleContexts.find(event);
	if (g==possibleContexts.end()) {// this means event was not found
//		std::cout << "returning 0 for cntx for event " << event << "\n";
		return 0;
	}
	else {
//		std::cout << "returning a value for cntx of size " << (g->second)->size() << " for event " << event << "\n";
		return g->second;	// a copy (?) of the second object (ULSet) in pair
	}
}

double HmmKeyedMaps::get(unsigned long event, unsigned long context) {
	HmmKeyedMaps::iterator f = find(context);
	if (f==end())
		return _backoff.get(event);
	else
		return f->second->get(event);
}

void HmmKeyedMaps::add(unsigned long context, unsigned long event, double value) {
//	std::cout << "CONTEXT " << context << ":\tadding " << event << " " << value << "\n";
	HmmKeyedMaps::iterator f = this->find(context);
	HmmMap* entry = NULL;
	if (f==end()) {		// if context not in matrix, add a new HmmMap for it
		entry = new HmmMap;
//		std::cout << "ADDING CONTEXT TO HmmKeyedMaps\n";
		(*this)[context] = entry;
//		(*this).insert(std::pair<unsigned long,HmmMap*>(context,entry));
	}
	else				// else it exists and get the ptr to the HmmMap
		entry = f->second;
	entry->add(event, value);	// add event and value to the HmmMap (actually a simple map)
	std::map<unsigned long, ULSet*>::iterator g = possibleContexts.find(event);
	ULSet* possCntx = NULL;
	if (g==possibleContexts.end()) {	// if event not in possible contexts, create possCntx
		possCntx = new ULSet;
		possibleContexts[event] = possCntx;
//		possibleContexts.insert(std::pair<unsigned long,ULSet*>(event,possCntx));
//		std::cout << event << " does not exist, creating and adding context\n";
	}
	else {
		possCntx = g->second;
//		std::cout << event << " does exist, adding context\n";
	}
	possCntx->insert(context);	// add event to possible contexts (which is a set)
}

void HmmKeyedMaps::load(std::istream& file, StrIdMap& strIdMap) {
	std::string name1, name2;
	double value;

	while (file>>name1>>name2>>value) {
		this->load(name1, name2, value, strIdMap);
	}
}

void HmmKeyedMaps::load(std::string& name1, std::string& name2, double value, StrIdMap& strIdMap) {
	add(strIdMap.getId(name1), strIdMap.getId(name2), log(value));
}

void HmmKeyedMaps::save(std::ostream& file, StrIdMap& strIdMap) {
	for (HmmKeyedMaps::iterator i = begin(); i!=end(); i++) {
		HmmMap& vals = *i->second;
		for (HmmMap::iterator j = vals.begin(); j!=vals.end(); j++) {
			file << strIdMap.getStr(i->first) << '\t'
					<< strIdMap.getStr(j->first) << '\t'
					<< exp(j->second) << std::endl;
		}
	}
}

bool HmmKeyedMaps::rand(unsigned long curr, unsigned long& next) {
	HmmKeyedMaps::iterator it = find(curr);
	if (it==end())
		return false;
	return it->second->rand(next);
}


HmmMap::HmmMap() {
	_smoothedZeroCount = SMOOTHEDZEROCOUNT;
}

void HmmMap::load(std::istream& file, StrIdMap& strIdMap) {
	std::string name;
	double c;
	while (file>>name>>c) {
		add(strIdMap.getId(name), log(c));
	}
}

void HmmMap::save(std::ostream& file, StrIdMap& strIdMap) {
	for (HmmMap::iterator i = begin(); i!=end(); i++) {
		file << strIdMap.getStr(i->first) << ' '
				<< exp(i->second) << std::endl;
	}
}

double HmmMap::get(unsigned long event) {
	HmmMap::iterator f = find(event);
	if (f==end())
		return smoothedZeroCount();
	else
		return f->second;
}

void HmmMap::add(unsigned long event, double value) {
	HmmMap::iterator f = find(event);
	if (f==end()) {
//		std::cout << "Adding " << event << " " << value << "\n";
		(*this)[event] += value;
	}
	else {
//		std::cout << "Setting " << event << " " << value << "\n";
		f->second = sumLogProb(f->second, value);

	}
}

bool HmmMap::rand(unsigned long& next) {
	double p = ((double) ::rand())/RAND_MAX;
	double total = 0;
	for (HmmMap::iterator it = begin(); it!=end(); it++) {
		total += exp(it->second);
		if (total >= p) {
			next = it->first;
			return true;
		}
	}
	return false;
}


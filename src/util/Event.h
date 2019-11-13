#ifndef EVENT_H_
#define EVENT_H_

#include <boost/foreach.hpp>

namespace hlsutil {

class Event {
private:
	std::vector<double> posteriors;
	std::vector<int> switches;
public:
	std::string state;
	bool isValid;
	double posterior;
	double phaseConcordance;
	unsigned int numMarkers;

	Event();
	virtual ~Event();
	void setPosteriors(std::vector<double> posteriors);
	void setSwitches(std::vector<int> switches);
};

} /* namespace hlsutil */
#endif /* EVENT_H_ */

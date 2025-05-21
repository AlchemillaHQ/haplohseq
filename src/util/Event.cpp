#include "Event.h"

namespace hlsutil {

/// @brief Constructor.
Event::Event() {
	this->numMarkers = 0;
	this->posterior = 0;
	this->phaseConcordance = 0;
	this->isValid = false;
}

/// @brief Destructor.
Event::~Event() {
}

void Event::setPosteriors(std::vector<double> posteriors) {
	this->posteriors = posteriors;
	double posteriorSum = 0;
	for (double ep : this->posteriors) {
		posteriorSum += ep;
	}
	this->posterior = posteriorSum / (double) this->numMarkers;
}

void Event::setSwitches(std::vector<int> switches) {
	this->switches = switches;
	double switchSum = 0;
	for (int ss : this->switches) {
		switchSum += ss;
	}
	this->phaseConcordance = switchSum / (double) numMarkers;
}

}

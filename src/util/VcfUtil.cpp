/*
 * VcfUtil.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#include "VcfUtil.h"

namespace hlsutil {

VcfUtil::VcfUtil() {}
VcfUtil::~VcfUtil() {}

// return a count of alt
std::string VcfUtil::gtToNumStr(std::string& gt, std::string& nullString) {
	if (gt == "0/0") {
		return "0";
	}
	if (gt == "1/0" || gt == "0/1") {
		return "1";
	}
	if (gt == "1/1") {
		return "2";
	}
	return nullString;
}

// return a count of B
std::string VcfUtil::abToNumStr(std::string& gt, std::string& nullString) {
	if (gt == "A/A") {
		return "0";
	}
	if (gt == "B/A" || gt == "A/B") {
		return "1";
	}
	if (gt == "B/B") {
		return "2";
	}
	return nullString;
}

} /* namespace hlsutil */

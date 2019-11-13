/*
 * StringUtil.cpp
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#include "StringUtil.h"

namespace hlsutil {

StringUtil::StringUtil() {}

StringUtil::~StringUtil() {}

std::string StringUtil::intToStr(int i) {
	std::ostringstream result;
	result << i;
	return result.str();
}

std::string StringUtil::doubleToStr(double d) {
	std::ostringstream result;
	result << d;
	return result.str();
}

std::string StringUtil::charToStr(char c) {
	std::ostringstream result;
	result << c;
	return result.str();
}

std::vector<std::string> &StringUtil::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string>StringUtil::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

// trim from start
std::string &StringUtil::ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
std::string &StringUtil::rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
std::string &StringUtil::trim(std::string &s) {
        return ltrim(rtrim(s));
}

} /* namespace hlsutil */

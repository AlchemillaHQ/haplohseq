/*
 * StringUtil.h
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#ifndef STRINGUTIL_H_
#define STRINGUTIL_H_

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

namespace hlsutil {

class StringUtil {
public:
	StringUtil();
	virtual ~StringUtil();

	// type conversions
	static std::string intToStr(int i);
	static std::string doubleToStr(double d);
	static std::string charToStr(char c);

	// splitting strings
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	std::vector<std::string> split(const std::string &s, char delim);

	// trim from start
	std::string &ltrim(std::string &s);

	// trim from end
	std::string &rtrim(std::string &s);

	// trim from both ends
	std::string &trim(std::string &s);
};

} /* namespace hlsutil */
#endif /* STRINGUTIL_H_ */

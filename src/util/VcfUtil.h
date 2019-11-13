/*
 * VcfUtil.h
 *
 * Author: Anthony San Lucas
 * Email:  sanlucas@gmail.com
 */

#ifndef VCFUTIL_H_
#define VCFUTIL_H_

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

namespace hlsutil {

class VcfUtil {
public:
	VcfUtil();
	virtual ~VcfUtil();

	// type conversions
	static std::string gtToNumStr(std::string& gt, std::string& nullString);
	static std::string abToNumStr(std::string& gt, std::string& nullString);
};

} /* namespace hlsutil */
#endif /* VCFUTIL_H_ */

#ifndef STRIDMAP_H
#define STRIDMAP_H

#include <vector>
#include <string>
#include <map>

class StrIdMap {
	std::map<std::string, unsigned long> _toId;
	std::vector<std::string> _toStr;
public:
	std::string getStr(unsigned long id) {
		return _toStr[id];
	}

	// If the map does not contain the element, the element is added to the end
	// of the _toStr vector and the vector index is stored in the _toId map.
	// Or if the map does contain the element, "second" (i.e., the id) is returned.
	unsigned long getId(std::string str) {
		std::map<std::string, unsigned long>::iterator f = _toId.find(str);
		unsigned long id;
		if (f == _toId.end()) {
			id = _toId.size();
			_toId[str] = id;
			_toStr.push_back(str);
			return id;
		}
		else {
			return f->second;
		}
	}
};

#endif

#ifndef YUTILS_H
#define YUTILS_H

#include <vector>
#include <string>
#include <sstream>

// tokenizer
inline
void Tokenize(const std::string &str, std::vector<std::string> &elems,
    char delimiter = '\t') {
    // http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
    // NOTE: this approach intentionally allows consecutive delimiters
    // NOTE: delimiter in the end of the str is ignored, which works fine for parsing SA tag
    std::stringstream ss(str);
    std::string item;
    while(getline(ss, item, delimiter)) {
        elems.push_back(item);  
    }
}

#endif // YUTILS_H 
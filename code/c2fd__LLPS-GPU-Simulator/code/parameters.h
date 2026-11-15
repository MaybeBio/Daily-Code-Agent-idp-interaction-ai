#pragma once
#include <string>
#include <map>


class Parameters {
    std::string filename;
public:
    std::map<std::string,std::string> m;
    void read_file(const std::string &filename, const std::string &selector);
    inline bool inited() {return filename!="";}

    double operator[](const std::string &s) const;
    std::string as_string(const std::string &s) const;
    int as_int(const std::string &s) const;
    bool as_bool(const std::string &s) const;
};


#include "parameters.h"
#include "tools.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <set>
#include <cctype>


using namespace std;




int Parameters::as_int(const string &s) const {
    return atoi(this->as_string(s).c_str());
}
bool Parameters::as_bool(const string &s) const {
    string v = this->as_string(s);
    if (v.length()==0) return false;
    if (tolower(v)=="false") return false;

    istringstream ss(v);
    double r;
    if (!(ss >> r).fail())
        return bool(r);

    return true;
}
double Parameters::operator[](const string &s) const {
    return atof(this->as_string(s).c_str());
}
string Parameters::as_string(const string &s) const {
    map<string,string>::const_iterator i =  m.find(s);

    if(i == m.end()) {
        cout << "Can't find parameter : "<<filename<<" : "<<s<<endl;
        my_exit(-1);
    }
    return i->second;
}


void Parameters::read_file(const string &filename, const string &selectors) {
    this->filename = filename;
    m.clear();

    // selectors "foo,bar"
    // --> sset = {foo},  sval = {bar:3}
    vector<string> sset_ = split(selectors, ",");
    set<string> sset;
    foreach(i, sset_) {
        string token = trim(*i);
        sset.insert(token);
        if (DEBUG & DEBUG_VERBOSE_PARAMS)
            cout<<" :: "<<token<<endl;
    }


    
    ifstream fin(filename.c_str());
    string s;
    int lineno=0;


    while(!fin.eof()) {
        lineno++;
        getline(fin, s);
        s = s.substr(0, s.find('#')); // remove comments
        s = trim(s);

        size_t i = s.find('=');
        if (i==string::npos) {
            if (s=="")
                continue;
            cout << "Syntax error in parameter file : " << filename
                 << " : " << lineno << endl
                 << s << endl;
            my_exit(-1);
        }

        string selector = "";
        string name  = trim(s.substr(0,i));
        string value = trim(s.substr(i+1,string::npos));

        i = name.find(":");
        if (i!=string::npos) {
            selector = trim(name.substr(0,i));
            name = trim(name.substr(i+1,string::npos));
        }


        if (selector=="" || sset.count(selector)) {
            if (DEBUG & DEBUG_VERBOSE_PARAMS)
                cout<<name<<" = "<<value <<endl;
            m[name] = value;
        }
    }

    cout<<endl;
}


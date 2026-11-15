#include "tools.h"

#include <execinfo.h>
#include <iostream>
#include <iomanip>
#include <cctype>
#include <cuda_runtime.h>

using namespace std;

int DEBUG = (1<<3)-1;
int OUTPUT = (1<<2)-1;





void my_exit(int ret, const string &message) {
    // From 3.0 - 3.2 release note for Linux
    // > It is a known issue that cudaThreadExit() may not be called
    // > implicitly on host thread exit. Due to this, developers are
    // > recommended to explicitly call cudaThreadExit() while the
    // > issue is being resolved.
    cudaThreadExit();

    cerr << flush;
    cout << message << flush;
    exit(ret);
}


void print_backtrace() {
    size_t size = 100;
    void *array[size];

    // get void*'s for all entries on the stack
    size = backtrace(array, size);

    // Write stack frames to file descriptor 2 (stderr)
    // To get function names, need linker flag "-rdynamic".
    backtrace_symbols_fd(array, size, 2);
    cerr << flush;
}



//================================================================


string tolower(const string s) {
    string r = s;
    for (int i=0; i<(int)r.size(); i++)
	r[i] = tolower(r[i]);
    return r;
}


// string + int
string operator+(const string &s, int i) {
    ostringstream os;
    os<<s<<setw(3)<<setfill('0')<<i;
    return os.str();
}


string itoa(int i, int width) {
    ostringstream os;
    os<<setw(width)<<setfill('0')<<i;
    return os.str();
}


// input:  " foo  bar  "
// output: "foo  bar"
string trim(const string &s) {
    size_t i=0;
    while (i<s.size() && isspace(s[i]))  i++;
    size_t j=s.size();
    while (i<j && isspace(s[j-1]))  j--;
    return s.substr(i,j-i);
}


vector<string> split(const string &s, const string &token, bool allow_blank_item) {
    vector<string> ret;
    size_t i=0;
    while(i<s.size()) {
        size_t j = s.find(token, i);
        if (j==string::npos)
            j = s.size();
        if (i<j || allow_blank_item)
            ret.push_back(trim(s.substr(i,j-i)));
        i = j+token.size();
    }
    return ret;
}



double uniform_rand(double a, double b) {
    // Return uniform random var from the interval [a,b) or [b,a),
    // whichever make sense.
    double r = (rand()/(double)RAND_MAX);
    return a + r*(b-a);
}

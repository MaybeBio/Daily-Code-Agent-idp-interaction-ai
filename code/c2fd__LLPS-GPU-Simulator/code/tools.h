#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>



// This needs compiler flag -arch=compute_13 or greater.
#define DOUBLE_PRECISION true	

#define NDIM3 true




enum {
    DEBUG_CHECK           = 1<<0,
    DEBUG_VERBOSE_VERSION = 1<<1,
    DEBUG_VERBOSE_PARAMS  = 1<<2,
    DEBUG_VERBOSE_VELOCITY_ERROR = 1<<3,
    DEBUG_VERBOSE_H5      = 1<<4,
    DEBUG_VERBOSE_STEPS   = 1<<5,
    DEBUG_VERBOSE_MEMORY  = 1<<6,
    DEBUG_LIST_NAN_INF           = 1<<7, // Need to remove +inf from tau_compute_tensor().
    DEBUG_LIST_NAN_INF_VERBOSE   = 1<<8,
    DEBUG_PAD_NAN         = 1<<9
};

//  1 -- io.cu: check()
//  5 -- io.cu: list_nan_inf
// 10 -- tau.inl:  print "tau*"
// 20 -- write_h5: print


// Turn off asserts.  #include "tools.h" needs to be before #include <cassert>
// #define NDEBUG		

enum {
    OUTPUT_PHI   = 1<<0,
    OUTPUT_V     = 1<<1,
    OUTPUT_C     = 1<<2,
    OUTPUT_TAU_N = 1<<3,
    OUTPUT_TAU_P = 1<<4,
    OUTPUT_V_N   = 1<<5, // v_n (also div_v_tau if output_momentum_check)
    OUTPUT_FORCE = 1<<6,
    OUTPUT_P_S   = 1<<7,
    OUTPUT_GN_GC = 1<<8,
    OUTPUT_VISCOUS_STRESS = 1<<9,
    OUTPUT_U              = 1<<10, // has steps
    OUTPUT_DIV            = 1<<11, // has steps, divergence of u, v
    OUTPUT_MOMENTUM_CHECK = 1<<12, // has steps,

    OUTPUT_TAU_DEBUG  = 1<<20,
    OUTPUT_HELMHOLTZ  = 1<<21,  // always output all iteration steps

    OUTPUT_BICG_STEPS = 1<<29,
    // This applies to OUTPUT_* that has steps.
    OUTPUT_ALL_ITERATION_STEPS  = 1<<30 // Output file will be very large.
};


// tools.cpp, inited by biofilm.cu
extern int DEBUG;
extern int OUTPUT;



//================


#if DOUBLE_PRECISION
typedef double myreal;
#else
typedef float myreal;
#endif



#define PI 3.1415926535897932384626433832795028841971693993751



void my_exit(int ret, const std::string &message="");
void print_backtrace();


//================

#define forn(i,n)  for(int i=0; i<int(n); i++)
#define let(a,b) typeof(b) a = (b)
#define foreach(i,c)  for(typeof((c).begin()) i = (c).begin(); i != (c).end(); i++)


inline bool has_element(const std::vector<std::string> &v, const std::string &a) {
    return find(v.begin(), v.end(), a) != v.end();
}



std::string tolower(const std::string s);
std::string operator+(const std::string &s, int i);
inline std::string concat(const std::string &s, int i) {return s+i;}
std::string itoa(int i, int width = 0);
std::string trim(const std::string &s);
std::vector<std::string> split(const std::string &s, const std::string &token, bool allow_blank_item=false);

double uniform_rand(double a=1., double b=0.);


template <typename T>
T next_arg(std::istringstream &ss, T default_value) {
    T ret;
    ss >> ret;
    if (ss.fail())
	ret = default_value;
    return ret;
}

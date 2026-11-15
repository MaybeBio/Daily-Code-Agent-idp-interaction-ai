#include <sys/resource.h>
#include <iostream>
#include <cassert>

#include <csignal>
#include "tools.h"
using namespace std;


void deal_with_rlimit(void) {
    cout<<"RLIM_INFINITY = "<<RLIM_INFINITY<<endl;

    const int N = 4;
    string name[N] = {"RLIMIT_NICE (20-x)", "RLIMIT_FSIZE",
		      "RLIMIT_DATA", "RLIMIT_STACK"};

    int resources[N] = {RLIMIT_NICE, RLIMIT_FSIZE,
			RLIMIT_DATA, RLIMIT_STACK};

    for (int i=0; i<N; i++) {
	rlimit rlim;
	getrlimit(resources[i], &rlim);
	cout<<name[i]<<" = "<<rlim.rlim_cur<<",  "<<rlim.rlim_max<<endl;
    }

    cout<<"\nSetting soft stack limit to 128 MB\n"<<endl;
    rlimit rlim;
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = 1<<27;	// 128 MB
    int result = setrlimit(RLIMIT_STACK, &rlim);
    assert(result==0);
}




void signal_handler(int sig) {
    cerr << "signal_handler(), with signal == " << sig << endl;
    print_backtrace();
    my_exit(-1);
}



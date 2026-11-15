// #include <getopt.h>	// This can't handle whitespace inside an optional argument.
#include <sys/stat.h>
#include <ctime>
#include <iostream>

#include <execinfo.h>
#include <csignal>

#include "state/state.h"
#include "misc.h"
#include "timestamp.h"



using namespace std;


void run(const string &data_dir,
	 const string &para_control,
	 const string &para_global,
	 const string &stage) {
    State sim;
    Parameters GlobalPara;
    Parameters &CP = sim.CP;
    sim.CP.read_file(para_control, stage);
    GlobalPara.read_file(para_global, stage);
    sim.NP.init(sim.CP, GlobalPara);

    gp = GridParams(CP.as_int("Ny"), CP.as_int("Nz"),
                    CP.as_int("dby"), CP.as_int("dbz"),
		    CP["Ly"], CP["Lz"]);
    gp.set_periodicity(true,true);


    globals::data_dir = data_dir;
    globals::write_pad = CP.as_bool("write_pad");
    OUTPUT = CP.as_int("OUTPUT");


    char buffer[1024];
    sprintf(buffer, "cp -fr $PWD %s", data_dir.c_str());
    if(system(buffer)) {
    	cout<<"Fail to execute : "<<buffer<<endl;
    	my_exit(-1);
    }

    if (DEBUG & DEBUG_VERBOSE_PARAMS) {
	sim.NP.display();
	cout<<gp<<endl;
    }

    {
	int device = CP.as_int("device");
	cout << "Using device # "<< device << endl;
        if (device>=0)  {
            if(cudaSuccess != cudaSetDevice(device))
                my_exit(-1);
        } else {
            if(cudaSuccess != cudaGetDevice(&device))
                my_exit(-1);
            cout << "            --> " << device << endl;
        }


	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, device);
	if (DEBUG & DEBUG_VERBOSE_PARAMS) {
	    cout << "    name               : " << prop.name << "\n"
		 << "    totalGlobalMem     : " << prop.totalGlobalMem << "\n"
		 << "    sharedMemPerBlock  : " << prop.sharedMemPerBlock << "\n"
		 << "    regsPerBlock       : " << prop.regsPerBlock << "\n"
		 << "    warpSize           : " << prop.warpSize << "\n"
		 << "    memPitch           : " << prop.memPitch << "\n"
		 << "    maxThreadsPerBlock : " << prop.maxThreadsPerBlock << "\n"
		 << "    totalConstMem      : " << prop.totalConstMem << "\n"
		 << "    major              : " << prop.major << "\n"
		 << "    minor              : " << prop.minor << "\n"
		 << "    clockRate          : " << prop.clockRate << "\n"
		 << "    multiProcessorCount      : " << prop.multiProcessorCount << "\n"
		 << "    kernelExecTimeoutEnabled : " << prop.kernelExecTimeoutEnabled << "\n"
		 << "    computeMode        : " << prop.computeMode << "\n"
		 << "    concurrentKernels  : " << prop.concurrentKernels << "\n"
		 << endl;
	}
    }
    {
	cout << "Test allocation of Vec ... " << flush;
	Vec foo(10, double(0.));
	Vec bar(10, 0.);
	cout << "Pass" << endl;
    }


    sim.init();
    gp.copy_to_device();
    sim.copy_NP_to_device();
    sim.set_initial_condition(); 
    sim.copy_NP_to_device();
    sim.copy_to_device();
    sim.simulate();
}




void print_usage(const string argv0) {
    cout << "Usage: " << argv0
	 << " output_dir [-d DEBUG] [-c para.control] [-g para.global]" << endl;
}


void parse_arguments(int argc, char **argv,
		     string &data_dir,
		     char **para_global,
		     char **para_control,
		     char **stage) {
    int ii=1;
    while (ii<argc) {
	string arg = argv[ii];
	if (data_dir=="" && arg[0]!='-') {
	    data_dir = argv[ii];
            if (data_dir[data_dir.size()-1] != '/')
                data_dir += '/';
	    ii++;
	} else if (arg=="-c" && ii+1<argc) {
	    *para_control = argv[ii+1];
	    ii+=2;
	} else if (arg=="-g" && ii+1<argc) {
	    *para_global = argv[ii+1];
	    ii+=2;
	} else if (arg=="-s" && ii+1<argc) {
	    *stage = argv[ii+1];
	    ii+=2;
	} else if (arg=="-d" && ii+1<argc) {
	    DEBUG = atoi(argv[ii+1]);
	    ii+=2;
	} else {
	    print_built_time();
	    cout<<endl;
	    print_usage(argv[0]);
	    my_exit(-1);
	}
    }
}




int main(int argc, char** argv) {
    string data_dir;
    char *para_global = "para.global";
    char *para_control = "para.control";
    char *stage    = "";

    parse_arguments(argc, argv,
		    data_dir,
		    &para_global,
		    &para_control,
		    &stage);


    if (data_dir=="") {
	print_usage(argv[0]);
	my_exit(-1);
    }
    struct stat st;
    if (stat(data_dir.c_str(), &st) != 0) {
        cout << "Output directory does not exist: " << data_dir <<endl;
	my_exit(-1);
    }

    if (DEBUG & DEBUG_VERBOSE_PARAMS) {
	cout << "\n" << argv[0]
	     << "  " << data_dir
	     << " " << para_control
	     << " " << para_global
	     << " " << stage << endl;
    }
    if (DEBUG & DEBUG_VERBOSE_VERSION) {
	print_built_time();
	cout << endl;
    }

    signal(SIGSEGV, signal_handler);


    //================================================================
    time_t start, end;
    double tdiff;
    int hours, minutes, seconds;

    time(&start);
    // Go, go, go!
    try {
	run(data_dir, para_control, para_global, stage);
    }
    catch (GridException e) {
	save_and_exit();
    }
    time(&end);
    //================================================================

    tdiff = difftime(end, start);
    seconds = (int)(tdiff);

    hours = seconds/3600;
    seconds %= 3600;
    minutes = seconds/60;
    seconds %= 60;
  
    cout << "\n\t\t\t--------------------------------\n";
    cout << "\t\t\t TOTAL TIME ELAPSED: ";
    cout.width(2); cout.fill('0'); cout << hours   << ':';
    cout.width(2); cout.fill('0'); cout << minutes << ':';
    cout.width(2); cout.fill('0'); cout << seconds;
    cout << "\n\t\t\t--------------------------------\n\n\n" << endl;
    return 0;
}

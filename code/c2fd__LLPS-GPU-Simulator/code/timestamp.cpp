#include <iostream>
#include "tools.h"

void print_built_time() {
    if (!(DEBUG & DEBUG_VERBOSE_VERSION))
	return;
    using namespace std;
    cout << "Built : " << __DATE__ << ", " << __TIME__ << endl;
#ifdef FLOW_TEST
    cout << "   with FLOW_TEST" << endl;
#endif
#ifdef TERNARY_MODEL
    cout << "   with TERNARY_MODEL" << endl;
#endif
#if DOUBLE_PRECISION
    cout << "   with DOUBLE_PRECISION" << endl;
#endif
    cout << "   DEBUG == " << DEBUG << endl;
}

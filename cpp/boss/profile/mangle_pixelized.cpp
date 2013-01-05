#include "../mangle2.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <vector>


using namespace std;
using namespace chrono;

struct Data {
	double ra, dec;
	long polyid;
};

int main() {
	// Open the file
	Mangle2::MaskClass mask("../gtests/lowz.ply");
	ifstream ff("../gtests/lowz.random");
	string sbuf;
	double ra, dec, wt;
	long polyid_in {0}, polyid_out {0};

	// Stop as soon as something fails
	// This prevents the code from spewing all over the place.
	vector<Data> dat;
	getline(ff, sbuf);
	while (getline(ff, sbuf)) {
		istringstream(sbuf) >> ra >> dec >> polyid_in;
		dat.push_back({ra, dec, polyid_in});
	}

	wt=0;
	steady_clock::time_point t1 = steady_clock::now();
	for (const Data& x : dat) {
		wt += mask.completeness_radec(x.ra, x.dec, polyid_out);
		if (polyid_out != x.polyid) cout << "Arg!";
	}
	steady_clock::time_point t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << wt << endl;
	cout << "Total evaluation time = "<< time_span.count()*1000.0 << endl;

	ff.close();
}


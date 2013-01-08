#include "src/npgsl/npRandom.h"
#include <iostream>
#include <vector>
#include <chrono>

using namespace std;
using namespace chrono;

int main() {
	double wt;

	{
		wt=0;
		npRandom rng(1234);
		steady_clock::time_point t1 = steady_clock::now();
		for (int ii=0; ii < 100000; ++ii) wt += rng(4)[1];
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		cout << wt << endl;
		cout << "Total evaluation time = "<< time_span.count()*1000.0 << endl;

	}

	{
		wt=0;
		npQuasiRandom qrng(4);
		steady_clock::time_point t1 = steady_clock::now();
		for (int ii=0; ii < 100000; ++ii) wt += qrng()[1];
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		cout << wt << endl;
		cout << "Total evaluation time = "<< time_span.count()*1000.0 << endl;

	}



}

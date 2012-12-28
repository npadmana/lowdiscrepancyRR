#include "ang2d.h"
#include "src/misc/npvecmath.h"
#include <iostream>
#include <sstream>
#include <iterator>
#include <boost/format.hpp>
#include <boost/program_options.hpp>


using namespace std;
using boost::format;
namespace po = boost::program_options;



void Ang2D::vecshift(vector<double> &x, vector<double> &x0) {
	double tmp;
	transform(x0.begin(), x0.end(), x.begin(), x.begin(),
			[&tmp](double a, double b){return modf(a+b,&tmp);});
}

void Ang2D::OutputData::print() {

	double mean, stddev, err;
	int nrand1;
	dvector val1;

	for (const auto& v1 : stats) {
		tie(nrand1, ignore, mean, stddev, err) = v1;
		cout << format("Using %7i terms ... \n") % nrand1;
		cout << format("The mean is %13.10e +/- %13.10e with a scatter of %13.10e, a fractional error of %9.6f percent\n")
				% mean % err % stddev % (stddev/mean * 100);
	}


}

void Ang2D::OutputData::push_back(const int nrand, const dvector& vec) {
	results.push_back(make_tuple(nrand, vec));
}

void Ang2D::OutputData::finalize() {

	// MPI finalization will happen here


	// Compute the statistics
	{
		int nrand1, nsim;
		dvector val1;
		double mean, stddev, err;

		// Compute statistics
		stats.clear();
		for (const auto& v1 : results) {
			tie(nrand1, val1) = v1; // A copy, but this is small
			nsim = val1.size();
			mean = vecmean(val1);
			stddev = vecstddev(val1, mean);
			err = stddev/sqrt(nsim);
			stats.push_back(make_tuple(nrand1, nsim, mean, stddev, err));

		}
	}
}

void Ang2D::OutputData::save(ofstream& out) {
	int nrand1, nsim;
	dvector val1;

	int nelems = results.size();
	out.write((char*)&nelems, sizeof(int));

	for (const auto& v1 : results) {
		tie(nrand1, val1) = v1; // A copy, but this is small
		nsim = val1.size();
		out.write((char*)&nrand1, sizeof(int));
		out.write((char*)&nsim, sizeof(int));
		out.write((char*)&val1[0], sizeof(double)*nsim);
	}
}

double Ang2D::OutputData::mean() {
	return get<2>(stats.back());
}

double Ang2D::OutputData::error() {
	return get<3>(stats.back());
}

template <class T>
vector<T> parsestring(const vector<string>& ss) {
	vector<T> out;

	for (string s1 : ss) {
		istringstream str(s1);
		copy(istream_iterator<T>(str), istream_iterator<T>(), back_inserter(out));
	}

	return out;
}


Ang2D::InputParams::InputParams(string fn) :
		_dump{false},
		_save{false},
		use_prng{false},
		ramin{0},
		ramax{360},
		decmin{-90},
		decmax{90},
		save_schedule{},
		thetabins{}
{
	// Get the mask configuration parameters
	{
		int prng_;
		vector<string> _save_schedule;
		vector<string> _thetabins;
		try {
			po::options_description desc("Allowed options");
			desc.add_options()
	    				("mask1",po::value<string>(&mask1fn), "Mask 1 configuration file")
	    				("mask2",po::value<string>(&mask2fn), "Mask 1 configuration file")
	    				("nrand",po::value<int>(&nrand)->default_value(10000), "Number of points")
	    				("nsim", po::value<int>(&nsim)->default_value(1), "Number of simulations")
	    				("use_prng",po::value<int>(&prng_)->default_value(0), "Use pseudo-random numbers")
	    				("dumpfn", po::value<string>(&dumpfn), "Dump file name")
	    				("savefn", po::value<string>(&savefn),"Save file name")
	    				("save_schedule", po::value< vector<string> >(&_save_schedule),"Save schedule")
	    				("thetabins", po::value< vector<string> >(&_thetabins),"Thetabins")
	    				;

			po::variables_map vm;
			ifstream ifs(fn);
			if (!ifs) {
				cout << "Error opening file " << fn << endl;
				throw runtime_error("Error opening file");
			}
			po::store(po::parse_config_file(ifs, desc), vm);
			po::notify(vm);
			ifs.close();

			if (vm.count("mask1")==0) {
				cout << "ERROR : An acceptance mask is required\n" << endl;
				throw invalid_argument("Missing acceptance mask");
			}

			if (vm.count("mask2")==0) mask2fn = mask1fn;
			if (vm.count("dumpfn")) _dump=true;
			if (vm.count("savefn")) _save=true;
			if (prng_ !=0) use_prng=true;

			if (vm.count("save_schedule")) save_schedule = parsestring<int>(_save_schedule);
			sort(save_schedule.begin(), save_schedule.end());
			if (vm.count("thetabins")) thetabins = parsestring<double>(_thetabins);
			sort(thetabins.begin(), thetabins.end());


		}
		catch (exception &e) {
			cout << e.what() << endl;
			throw;
		}
	}

}

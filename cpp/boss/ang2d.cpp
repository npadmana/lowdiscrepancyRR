#include "ang2d.h"
#include "src/misc/npvecmath.h"
#include <iostream>
#include <sstream>
#include <iterator>
#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include "mpi.h"


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
	MPI_Barrier(MPI_COMM_WORLD); // Barrier for sanity
	{
		dvector val1;
		for (auto& v1 : results) {
			tie(ignore, val1) = v1; // Copy here, don't try to do this in place
			MPI_Allreduce(&val1[0], &(get<1>(v1)[0]), val1.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}


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

double Ang2D::OutputData::error(bool percent) {
	if (percent) {
		return get<3>(stats.back())/mean() * 100.0;
	} else {
		return get<3>(stats.back());
	}
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
		thetabins{},
		verbose{0}
{
	// Get the mask configuration parameters
	{
		int prng_;
		vector<string> _save_schedule;
		vector<string> _thetabins;

		double thetamin, dtheta;
		int nthetabins;

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
	    				("thetamin",po::value<double>(&thetamin),"Theta_min")
	    				("dtheta", po::value<double>(&dtheta), "dtheta")
	    				("nthetabins", po::value<int>(&nthetabins), "number of theta bins")
	    				("verbose", po::value<int>(&verbose)->default_value(0),"Verbosity")
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

			// Thetabins has priority over incremental fills
			if (vm.count("thetabins")) {
				thetabins = parsestring<double>(_thetabins);
				sort(thetabins.begin(), thetabins.end());
			}
			// Incremental fill -- check
			else if ((vm.count("thetamin")>0) &&
					   (vm.count("dtheta")>0) && (vm.count("nthetabins")>0))
			{
				thetabins.resize(nthetabins);
				for (int ii=0; ii<=nthetabins; ++ii) thetabins[ii] = ii*dtheta + thetamin;
			}


		}
		catch (exception &e) {
			cout << e.what() << endl;
			throw;
		}
	}

}

/** Helper function for partition code
 *
 */
tuple<int, int> Ang2D::partition(int n) {
	// Short circuit the case of n=1
	if (n==1) return make_tuple(1,1);

	int sn = static_cast<int>(floor(sqrt(n)+0.01)); // shift slightly to avoid missing perfect squares.
	while ((n % sn)!=0) sn--;

	return make_tuple(sn, n/sn);
}

tuple<double, double> Ang2D::rotatePoint(double ra1, double dec1,
		double costheta, double phi) {

	double cth = sin(dec1*D2R);  // Note that Dec is 90-theta
	double sth = cos(dec1*D2R);  // Note that Dec is 90-theta
	double cphi = cos(ra1*D2R);
	double sphi = sin(ra1*D2R);

	// Define the rotation matrix
	Matrix3d rotmat;
	rotmat << cth * cphi, -sphi, sth*cphi,
			  cth * sphi,  cphi, sth*sphi,
			  -sth      ,     0, cth      ;

	// Define the input and output vectors
	Vector3d x1, x2;
	double sth_ = sqrt(1-costheta*costheta);
	double cphi_ = cos(phi*D2R);
	double sphi_ = sin(phi*D2R);
	x1 << sth_* cphi_, sth_*sphi_, costheta;
	x2.noalias() = rotmat * x1;

	// Now convert back to ra, dec
	double ra2, dec2;
	dec2 = 90 - acos(x2.coeff(2))/D2R; // no range checking
	ra2 = atan2(x2.coeff(1), x2.coeff(0))/D2R; // (-180, 180)
	if (ra2 < 0) ra2 = 360+ra2;

	return make_tuple(ra2,dec2);
}

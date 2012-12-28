#include <iostream>
#include <vector>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv) {
	vector<double> foo;
	{
		try {
			po::options_description desc("Allowed options");
			desc.add_options()
	    				("doubles", po::value< vector<double> >(&foo)->multitoken(), "doubles")
	    				;
			po::variables_map vm;
			po::store(po::parse_command_line(argc, argv, desc), vm);
			po::notify(vm);


		}
		catch (exception &e) {
					cout << e.what() << "\n";
					return 1;
				}
	}

	cout << endl;
	for (auto x : foo) cout << x << " ";
	cout << endl;


}

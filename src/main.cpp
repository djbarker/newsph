#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/mpi/environment.hpp>

#include "core/Simulation.hpp"

using namespace std;
using namespace sim;

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		cerr << "Requires a config file name." << endl;
		return 1;
	}

	cout << std::boolalpha;

	// setup mpi
	boost::mpi::environment(argc,argv);

	Simulation<2> theSim;

	// currently only takes one argument - the config file name
	theSim.loadConfigXML(string(argv[1]));

	theSim.writeOutput(0);

	return 0;
}

#include <iostream>
#include <fstream>
#include <iomanip>

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

	MPI_Init(&argc,&argv);

	Simulation<2> theSim;

	// currently only takes one argument - the config file name
	theSim.loadConfigXML(string(argv[1]));

	theSim.writeOutput();


	MPI_Finalize();

	return 0;
}

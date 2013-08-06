#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/mpi/environment.hpp>

#include "core/Simulation.hpp"
#include "physics/PredictorCorrector.hpp"
#include "physics/Sigma.hpp"
#include "kernels/WendlandQuintic.hpp"

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

	size_t file_number = 0;
	theSim.writeOutput(file_number);

	double tmax = discard_dims(theSim.parameters().tmax);
	for(double t=0;t<tmax; t += discard_dims(theSim.parameters().dt))
	{
		// set values to zero
		theSim.applyFunctions(physics::ResetVals<2>());

		theSim.doSPHSum<kernels::WendlandQuintic,physics::SigmaCalc<2>>(0u,physics::SigmaCalc<2>());

		theSim.writeOutput(file_number);
		++file_number;
	}


	return 0;
}

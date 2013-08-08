#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/mpi/environment.hpp>

#include "core/Simulation.hpp"
#include "physics/PredictorCorrector.hpp"
#include "physics/Sigma.hpp"
#include "kernels/WendlandQuintic.hpp"

#ifndef DIM
#define DIM 2
#endif

using namespace std;
using namespace sim;

int run_main(int argc, char* argv[])
{
	if(argc<2)
	{
		cerr << "Requires a config file name." << endl;
		return 1;
	}

	cout << std::boolalpha;

	// setup mpi
	boost::mpi::environment env(argc,argv,true);

	Simulation<DIM> theSim;

	// currently only takes one argument - the config file name
	theSim.loadConfigXML(string(argv[1]));

	cout << "HERE 0 " << endl;

	size_t file_number = 0;
	theSim.writeOutput(file_number);

	cout << "HERE 1 " << endl;

	double tmax = discard_dims(theSim.parameters().tmax);
	for(double t=0.;t<tmax; t += discard_dims(theSim.parameters().dt))
	{
		// set values to zero
		theSim.applyFunctions(physics::ResetVals<DIM>());

		theSim.placeParticlesIntoLinkedCellGrid(0);

		cout << "HERE 2 " << endl;

		// calculate sigma
		theSim.doSPHSum<kernels::WendlandQuintic>(0u,physics::SigmaCalc<DIM>());

		cout << "HERE 3 " << endl;

		theSim.writeOutput(file_number);
		++file_number;

		if(t>0.) break;
	}


	return 0;
}

int main(int argc, char* argv[])
{
	try
	{
		return run_main(argc,argv);
	}
	catch(std::runtime_error& e)
	{
		cout << "std::runtime_error: " << e.what() << endl;
	}

	return 1;
}

#undef DIM

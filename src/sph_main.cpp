#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/nonblocking.hpp>

#include "core/Simulation.hpp"
#include "physics/PredictorCorrector.hpp"
#include "physics/Sigma.hpp"
#include "kernels/WendlandQuintic.hpp"

#ifndef DIM
#define DIM 2
#endif

using namespace std;
using namespace sim;
using namespace dims;

/*
 * General TODOs which don't fit anywhere specific in the code
 */
// TODO: use extern templates to speed compilation
// TODO: change linked cell grid to boost::multi_array

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
	boost::mpi::communicator comm;

	Simulation<DIM> theSim;

	// currently only takes one argument - the config file name
	theSim.loadConfigXML(string(argv[1]));

	size_t file_number = 0;
	theSim.writeOutput(file_number);
	++file_number;

	double tmax = discard_dims(theSim.parameters().tmax);
	for(double t=0.;t<tmax; t += discard_dims(theSim.parameters().dt))
	{
		if(comm.rank()==0) cout << "t = " << t << endl;

		// set values to zero
		theSim.applyFunctions(physics::ResetVals<DIM>());

		if(comm.rank()==0) cout << "Rest values." << endl;

		theSim.placeParticlesIntoLinkedCellGrid(0);

		if(comm.rank()==0) cout << "Placed into LCG." << endl;

		theSim.exchangeFull();

		if(comm.rank()==0) cout << "Exchanged." << endl;

		// calculate sigma
		theSim.doSPHSum<kernels::WendlandQuintic>(0u,physics::SigmaCalc<DIM>());

		if(comm.rank()==0) cout << "Sigma." << endl;

		if(comm.rank()==0) cout << "Exchange data." << endl;

		// calculate density then pressure
		theSim.applyFunctions(physics::DensityCalc<DIM>(theSim),physics::TaitEquation<DIM>(theSim));

		if(comm.rank()==0) cout << "Density/Tait." << endl;

		// calculate acceleration
		theSim.doSPHSum<kernels::WendlandQuintic>(0u,physics::GradPCalc<DIM>(),physics::ViscCalc<DIM,0>());

		if(comm.rank()==0) cout << "GradP/Visc." << endl;

		theSim.writeOutput(file_number);
		++file_number;

		if(t>0.001) break;
	}

	if(comm.rank()==0) cout << "Finished." << endl;

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

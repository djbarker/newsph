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

		/*
		 * Half-step
		 */

		theSim.applyFunctions(physics::ResetVals<DIM>());	// set values to zero
		theSim.placeParticlesIntoLinkedCellGrid(0);
		theSim.exchangeFull();

		// calculate sigma
		theSim.doSPHSum<kernels::WendlandQuintic>(0u,physics::SigmaCalc<DIM>());
		theSim.exchangeData();

		theSim.writeOutput(file_number);
		file_number++;

		// calculate density then pressure
		theSim.applyFunctions(physics::DensityCalc<DIM>(theSim),physics::TaitEquation<DIM>(theSim));

		// calculate acceleration
		theSim.doSPHSum<kernels::WendlandQuintic>(0u,physics::GradPCalc<DIM>(),physics::ViscCalc<DIM,0>());

		// move particles
		theSim.applyFunctions(physics::PredictorCorrectorUpdater<0>(theSim.parameters().dt));

		/*
		 * Full step
		 */

		theSim.applyFunctions(physics::ResetVals<DIM>());	// set values to zero
		theSim.placeParticlesIntoLinkedCellGrid(1);
		theSim.exchangeFull();

		// calculate sigma
		theSim.doSPHSum<kernels::WendlandQuintic>(1u,physics::SigmaCalc<DIM>());
		theSim.exchangeData();

		// calculate density then pressure
		theSim.applyFunctions(physics::DensityCalc<DIM>(theSim),physics::TaitEquation<DIM>(theSim));

		// calculate acceleration
		theSim.doSPHSum<kernels::WendlandQuintic>(1u,physics::GradPCalc<DIM>(),physics::ViscCalc<DIM,1>());

		// move particles
		theSim.applyFunctions(physics::PredictorCorrectorUpdater<1>(theSim.parameters().dt));

		theSim.writeOutput(file_number);
		++file_number;

		if(t>4*discard_dims(theSim.parameters().dt)) break;
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

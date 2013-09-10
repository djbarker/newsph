#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/nonblocking.hpp>

#include "core/Simulation.hpp"
#include "physics/PredictorCorrector.hpp"
#include "physics/Sigma.hpp"
#include "kernels/WendlandQuintic.hpp"

//if unspecified default to 2 dimensions
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
// TODO: use boost::program_options to get cmd line args

int run_main(int argc, char* argv[], boost::mpi::environment& env)
{
	if(argc<2)
	{
		cerr << "Requires a config file name." << endl;
		return 1;
	}

	cout << std::boolalpha;

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

		comm.barrier();
		theSim.exchangeFull();

		comm.barrier();
		if(comm.rank()==0) cout << "HERE 0" << endl;

		// calculate sigma
		theSim.doSPHSum<kernels::WendlandQuintic>(0,physics::SigmaCalc<DIM>());
		if(comm.rank()==0) cout << "HERE 0.1" << endl;
		theSim.exchangeData();

		if(comm.rank()==0) cout << "HERE 1" << endl;

		// calculate density then pressure
		theSim.applyFunctions(physics::DensityCalc<DIM>(),physics::TaitEquation<DIM>());

		if(comm.rank()==0) cout << "HERE 2" << endl;

		// calculate acceleration
		theSim.doSPHSum<kernels::WendlandQuintic>(0,physics::GradPCalc<DIM>(),physics::ViscCalc<DIM,0>());

		if(comm.rank()==0) cout << "HERE 3" << endl;

		// move particles
		theSim.applyFunctions(physics::PredictorCorrectorUpdater<0,DIM>());

		if(comm.rank()==0) cout << "HERE 4" << endl;

		/*
		 * Full step
		 */

		theSim.exchangeOutOfBounds(1);

		comm.barrier();
		if(comm.rank()==0) cout << "HERE 5" << endl;

		theSim.applyFunctions(physics::ResetVals<DIM>());	// set values to zero
		comm.barrier();
		if(comm.rank()==0) cout << "HERE 5.1" << endl;

		theSim.placeParticlesIntoLinkedCellGrid(1);
		comm.barrier();
		if(comm.rank()==0) cout << "HERE 5.2" << endl;

		theSim.exchangeFull();

		if(comm.rank()==0) cout << "HERE 6" << endl;

		// calculate sigma
		theSim.doSPHSum<kernels::WendlandQuintic>(1,physics::SigmaCalc<DIM>());
		theSim.exchangeData();

		if(comm.rank()==0) cout << "HERE 7" << endl;

		// calculate density then pressure
		theSim.applyFunctions(physics::DensityCalc<DIM>(),physics::TaitEquation<DIM>());

		if(comm.rank()==0) cout << "HERE 8" << endl;

		// calculate acceleration
		theSim.doSPHSum<kernels::WendlandQuintic>(1,physics::GradPCalc<DIM>(),physics::ViscCalc<DIM,1>());

		if(comm.rank()==0) cout << "HERE 9" << endl;

		// move particles
		theSim.applyFunctions(physics::PredictorCorrectorUpdater<1,DIM>());

		if(comm.rank()==0) cout << "HERE 10" << endl;

		theSim.exchangeOutOfBounds(0);

		theSim.writeOutput(file_number);
		++file_number;

		if(t>4*discard_dims(theSim.parameters().dt)) break;
	}

	if(comm.rank()==0) cout << "Finished." << endl;

	return 0;
}

int main(int argc, char* argv[])
{
	// setup mpi
	boost::mpi::environment env(argc,argv,true);
	boost::mpi::communicator comm;

	// run the program
	try
	{
		return run_main(argc,argv,env);
	}
	catch(std::runtime_error& e)
	{
		cout << "[Proc " << comm.rank() << "] std::runtime_error (or derived): " << e.what() << endl;
		env.abort(1);
	}

	return 1;
}

#undef DIM

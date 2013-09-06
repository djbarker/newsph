#ifndef PREDICTORCORRECTOR_HPP_
#define PREDICTORCORRECTOR_HPP_

#include <dims.hpp>

using namespace dims;

namespace sim
{
namespace physics
{

template<size_t Step, size_t Dim>
struct PredictorCorrectorUpdater{};

template<size_t Dim>
struct PredictorCorrectorUpdater<0,Dim>
{
	template<class PType> void operator() (PType& part, Simulation<Dim>& sim)
	{
		part.pos[1] = part.pos[0] + part.vel[0]*sim.parameters().dt/(2.0_number);
		part.vel[1] = part.vel[0] + part.acc*sim.parameters().dt/(2.0_number);

		// limit velocity to h/dt
		quantity<velocity> tmp = (sim.parameters().h)/sim.parameters().dt;
		if(part.vel[1].magnitude()>tmp)
		{
			part.vel[1] = part.vel[1].unit()*tmp;
		}
	}
};


template<size_t Dim>
struct PredictorCorrectorUpdater<1,Dim>
{
	template<class PType> void operator() (PType& part, Simulation<Dim>& sim)
	{
		part.pos[1] = part.pos[0] + part.vel[1]*sim.parameters().dt/(2.0_number);
		part.vel[1] = part.vel[0] + part.acc*sim.parameters().dt/(2.0_number);

		part.pos[0] = 2.0_number*part.pos[1] - part.pos[0];
		part.vel[0] = 2.0_number*part.vel[1] - part.vel[0];

		// limit velocity to h/dt
		quantity<velocity> tmp = (sim.parameters().h)/sim.parameters().dt;
		if(part.vel[0].magnitude()>tmp)
		{
			part.vel[0] = part.vel[0].unit()*tmp;
		}
	}
};

}
}

#endif /* PREDICTORCORRECTOR_HPP_ */

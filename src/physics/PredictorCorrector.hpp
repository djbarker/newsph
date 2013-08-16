#ifndef PREDICTORCORRECTOR_HPP_
#define PREDICTORCORRECTOR_HPP_

#include <dims.hpp>

using namespace dims;

namespace sim
{
namespace physics
{

template<size_t Step>
struct PredictorCorrectorUpdater
{
	PredictorCorrectorUpdater(quantity<dims::time> dt):dt(dt){};

	template<class PType> void operator() (PType& part);

private:
	quantity<dims::time> dt;
};


template<>
template<class PType>
void PredictorCorrectorUpdater<0>::operator() (PType& part)
{
		part.pos[1] = part.pos[0] + part.vel[0]*dt/(2.0_number);
		part.vel[1] = part.vel[0] + part.acc*dt/(2.0_number);
}

template<>
template<class PType>
void PredictorCorrectorUpdater<1>::operator() (PType& part)
{
		part.pos[1] = part.pos[0] + part.vel[1]*dt/(2.0_number);
		part.vel[1] = part.vel[0] + part.acc*dt/(2.0_number);

		part.pos[0] = 2.0_number*part.pos[1] - part.pos[0];
		part.vel[0] = 2.0_number*part.vel[1] - part.vel[0];
}

}
}

#endif /* PREDICTORCORRECTOR_HPP_ */
